"""Deterministic cancellation tests: events control the relevant interleavings."""
import sys
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.jobs import JobManager, _write_json
from qtlift.pipeline import JobCancelled


class JobCancellationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.manager = JobManager(Path(self.directory.name))
        self.addCleanup(self.manager.executor.shutdown, wait=True)

    def submit(self, **overrides):
        return self.manager.submit({"source_ref": "RefB", "contig": "Chr1", "start": 1, "end": 10, **overrides})["job_id"]

    def wait(self, event):
        self.assertTrue(event.wait(5), "worker did not reach the controlled boundary")

    def test_cancel_after_completed_summary_does_not_resurrect_job(self):
        completed, release = threading.Event(), threading.Event()
        self.addCleanup(release.set)

        def run(payload, jobs_root, progress, cancel_event, **kwargs):
            summary = {"job_id": payload["job_id"], "status": "completed", "progress": 100, "stage": "Completed"}
            _write_json(Path(jobs_root) / payload["job_id"] / "summary.json", summary)
            completed.set()
            if not release.wait(5):
                raise RuntimeError("test did not release worker")
            return summary

        with patch("qtlift.jobs.run_job", side_effect=run):
            job_id = self.submit()
            self.wait(completed)
            row = self.manager.cancel(job_id)
            release.set()
            self.manager.executor.shutdown(wait=True)
        self.assertEqual(row["status"], "completed")
        self.assertEqual(self.manager.get(job_id)["status"], "completed")

    def test_cancel_during_report_generation_finishes_cancelled(self):
        writing, release = threading.Event(), threading.Event()
        self.addCleanup(release.set)

        def run(payload, jobs_root, progress, cancel_event, **kwargs):
            progress(95, "Writing reports")
            writing.set()
            if not release.wait(5):
                raise RuntimeError("test did not release worker")
            summary = {"job_id": payload["job_id"], "status": "completed", "progress": 100, "stage": "Completed"}
            if kwargs.get("persist_summary", True):
                _write_json(Path(jobs_root) / payload["job_id"] / "summary.json", summary)
            return summary

        with patch("qtlift.jobs.run_job", side_effect=run):
            job_id = self.submit()
            self.wait(writing)
            self.assertEqual(self.manager.cancel(job_id)["status"], "cancelling")
            release.set()
            self.manager.executor.shutdown(wait=True)
        self.assertEqual(self.manager.get(job_id)["status"], "cancelled")

    def test_cancel_and_worker_terminal_write_cannot_share_temporary_file(self):
        running, cancel_at_replace = threading.Event(), threading.Event()
        worker_update_attempted, cancel_published = threading.Event(), threading.Event()
        original_lock, original_replace = self.manager._lock, Path.replace

        class ObservedLock:
            # Signal before acquisition, so a correctly serialized worker need not
            # acquire the lock while the test holds cancel() at its atomic rename.
            def __enter__(self):
                if threading.current_thread().name.startswith("qtlift-job") and cancel_at_replace.is_set():
                    worker_update_attempted.set()
                original_lock.acquire()

            def __exit__(self, *args):
                original_lock.release()

        self.manager._lock = ObservedLock()

        def run(payload, jobs_root, progress, cancel_event, **kwargs):
            running.set()
            if not cancel_at_replace.wait(5):
                raise RuntimeError("test did not reach cancellation publication")
            raise JobCancelled("Job cancelled by user.")

        def replace(path, target):
            if path.name == "summary.json.tmp" and threading.current_thread() is threading.main_thread():
                cancel_at_replace.set()
                self.wait(worker_update_attempted)
                result = original_replace(path, target)
                cancel_published.set()
                return result
            if cancel_at_replace.is_set():
                worker_update_attempted.set()
                self.wait(cancel_published)
            return original_replace(path, target)

        with patch("qtlift.jobs.run_job", side_effect=run):
            job_id = self.submit()
            self.wait(running)
            future = self.manager._futures[job_id]
            with patch.object(Path, "replace", replace):
                self.manager.cancel(job_id)
                future.result(timeout=5)
        self.assertEqual(self.manager.get(job_id)["status"], "cancelled")

    def test_queued_cancellation_cleans_up_and_never_starts(self):
        running, release = threading.Event(), threading.Event()
        self.addCleanup(release.set)
        started = []

        def run(payload, jobs_root, progress, cancel_event, **kwargs):
            started.append(payload["job_id"])
            running.set()
            if not release.wait(5):
                raise RuntimeError("test did not release worker")
            return {"job_id": payload["job_id"], "status": "completed", "progress": 100}

        with patch("qtlift.jobs.run_job", side_effect=run):
            first_id = self.submit()
            self.wait(running)
            queued_id = self.submit()
            row = self.manager.cancel(queued_id)
            self.assertEqual(row["status"], "cancelled")
            self.assertEqual(row["stage"], "Cancelled before start")
            self.assertEqual(self.manager.cancel(queued_id), row)
            self.assertNotIn(queued_id, self.manager._events)
            self.assertNotIn(queued_id, self.manager._futures)
            release.set()
            self.manager.executor.shutdown(wait=True)
        self.assertEqual(started, [first_id])

    def test_failed_job_is_terminal_for_cancellation(self):
        with patch("qtlift.jobs.run_job", side_effect=ValueError("synthetic failure")):
            job_id = self.submit()
            self.manager.executor.shutdown(wait=True)
        row = self.manager.get(job_id)
        self.assertEqual(row["status"], "failed")
        self.assertEqual(row["error"], "synthetic failure")
        self.assertEqual(self.manager.cancel(job_id), row)

    def test_real_pipeline_keeps_running_summary_until_reports_finish(self):
        from qtlift.reporting import write_outputs
        from scripts.create_sample_data import main

        main()
        reports_written, release = threading.Event(), threading.Event()
        self.addCleanup(release.set)

        def write(*args, **kwargs):
            files = write_outputs(*args, **kwargs)
            reports_written.set()
            if not release.wait(5):
                raise RuntimeError("test did not release report writer")
            return files

        with patch("qtlift.pipeline.write_outputs", side_effect=write):
            job_id = self.submit(genome_root=str(ROOT / "sample_data" / "genomes"),
                                 target_ref="RefA", start=100, end=850, mapping_backend="exact")
            self.wait(reports_written)
            self.assertEqual(self.manager.get(job_id)["status"], "running")
            self.assertEqual(self.manager.get(job_id)["progress"], 95)
            self.assertEqual(self.manager.cancel(job_id)["status"], "cancelling")
            release.set()
            self.manager.executor.shutdown(wait=True)
        row = self.manager.get(job_id)
        self.assertEqual(row["status"], "cancelled")
        self.assertNotIn("final", row)
        self.assertNotIn("files", row)
        self.assertFalse(list(Path(self.directory.name).rglob("*.tmp")))

    def test_cancel_running_job_ends_cancelled(self):
        running = threading.Event()

        def run(payload, jobs_root, progress, cancel_event, **kwargs):
            running.set()
            if not cancel_event.wait(5):
                raise RuntimeError("test did not request cancellation")
            raise JobCancelled("Job cancelled by user.")

        with patch("qtlift.jobs.run_job", side_effect=run):
            job_id = self.submit()
            self.wait(running)
            self.assertEqual(self.manager.cancel(job_id)["status"], "cancelling")
            self.manager.executor.shutdown(wait=True)
        self.assertEqual(self.manager.get(job_id)["status"], "cancelled")
        self.assertNotIn(job_id, self.manager._events)
        self.assertNotIn(job_id, self.manager._futures)


if __name__ == "__main__":
    unittest.main()
