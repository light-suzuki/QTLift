from __future__ import annotations

import gzip
import sys
from pathlib import Path


def main(source_name: str, output_name: str, chunk_size: int) -> None:
    source, output = Path(source_name), Path(output_name)
    opener = gzip.open if source.suffix.lower() == ".gz" else open
    sequence = bytearray()
    name, offset = None, 0
    with opener(source, "rt", encoding="utf-8") as src, output.open("wb") as dst:
        def flush(final: bool = False) -> None:
            nonlocal offset
            while len(sequence) >= chunk_size or (final and sequence):
                size = min(chunk_size, len(sequence))
                dst.write(f">QTLIFT|{name}|{offset}\n".encode())
                chunk = bytes(sequence[:size])
                dst.write(chunk + b"\n")
                del sequence[:size]
                offset += size
        for line in src:
            if line.startswith(">"):
                if name is not None:
                    flush(True)
                name, offset = line[1:].split()[0], 0
            else:
                sequence.extend(line.strip().upper().encode("ascii"))
                flush()
        if name is not None:
            flush(True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
