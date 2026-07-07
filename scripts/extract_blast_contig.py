from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main(database: str, contig: str, output: str) -> None:
    process = subprocess.Popen(["/usr/bin/blastdbcmd", "-db", database, "-entry", "all", "-outfmt", "%f"],
                               stdout=subprocess.PIPE, text=True, encoding="utf-8")
    found = False
    destination = Path(output)
    try:
        with destination.open("w", encoding="ascii", newline="\n") as handle:
            assert process.stdout is not None
            for line in process.stdout:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    if found and name != contig:
                        break
                    found = name == contig
                if found:
                    handle.write(line)
    finally:
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
    if not found or destination.stat().st_size == 0:
        destination.unlink(missing_ok=True)
        raise SystemExit(f"Contig not found in BLAST database: {contig}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
