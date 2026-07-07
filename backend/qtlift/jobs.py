from __future__ import annotations

import json
import threading
from concurrent.futures import Future, ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from uuid import uuid4

from .pipeline import JobCancelled, run_job


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    temporary.replace(path)


class JobManager:
    def __init__(self, jobs_root: Path):
        self.jobs_root = jobs_root
        self.executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="qtlift-job")
        self._lock = threading.Lock()
        self._events: dict[str, threading.Event] = {}
        self._futures: dict[str, Future] = {}
        self._recover_interrupted()

    def _path(self, job_id: str) -> Path:
        return self.jobs_root / job_id / "summary.json"

    def _recover_interrupted(self) -> None:
        for path in self.jobs_root.glob("*/summary.json"):
            try:
                row = json.loads(path.read_text(encoding="utf-8"))
                if row.get("status") in ("queued", "running", "cancelling"):
                    row.update(status="failed", stage="Interrupted by server restart", error="Server restarted before the job completed.")
                    _write_json(path, row)
            except Exception:
                continue

    def submit(self, payload: dict) -> dict:
        job_id = f"qtlift-{uuid4().hex[:10]}"
        created = datetime.now().isoformat(timespec="seconds")
        payload = {**payload, "job_id": job_id, "_created_at": created}
        row = {"job_id": job_id, "name": payload.get("name") or "Unnamed region", "status": "queued",
               "progress": 0, "stage": "Waiting to start", "created_at": created,
               "source_label": f"{payload['source_ref']} {payload['contig']}:{payload['start']:,}-{payload['end']:,}",
               "final_label": "Pending", "confidence": "Manual check", "warnings": []}
        _write_json(self._path(job_id), row)
        event = threading.Event()
        with self._lock:
            self._events[job_id] = event
            self._futures[job_id] = self.executor.submit(self._run, payload, event)
        return row

    def _run(self, payload: dict, event: threading.Event) -> None:
        job_id = payload["job_id"]
        def progress(percent: int, stage: str) -> None:
            current = self.get(job_id)
            current.update(status="cancelling" if event.is_set() else "running", progress=percent, stage=stage)
            _write_json(self._path(job_id), current)
        try:
            progress(1, "Starting")
            run_job(payload, self.jobs_root, progress, event)
        except JobCancelled as exc:
            current = self.get(job_id); current.update(status="cancelled", stage="Cancelled", error=str(exc))
            _write_json(self._path(job_id), current)
        except Exception as exc:
            current = self.get(job_id); current.update(status="failed", stage="Failed", error=str(exc))
            _write_json(self._path(job_id), current)
        finally:
            with self._lock:
                self._events.pop(job_id, None); self._futures.pop(job_id, None)

    def get(self, job_id: str) -> dict:
        path = self._path(job_id)
        if not path.is_file():
            raise FileNotFoundError(job_id)
        return json.loads(path.read_text(encoding="utf-8"))

    def cancel(self, job_id: str) -> dict:
        with self._lock:
            event, future = self._events.get(job_id), self._futures.get(job_id)
            if not event or not future:
                return self.get(job_id)
            event.set()
            if future.cancel():
                row = self.get(job_id); row.update(status="cancelled", stage="Cancelled before start")
                _write_json(self._path(job_id), row)
            else:
                row = self.get(job_id); row.update(status="cancelling", stage="Cancellation requested")
                _write_json(self._path(job_id), row)
            return row
