#!/usr/bin/env python3
"""Write the deterministic FASTQs used by .test_dotmatch."""

from __future__ import annotations

import csv
import gzip
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TEST = ROOT / ".test_dotmatch"
LIBRARY = TEST / "resources" / "guides.csv"
READS = TEST / "reads"

PROFILES = {
    "HT_1": [40, 36, 32, 28, 24, 20, 8, 8, 6, 6, 4, 4],
    "HT_2": [38, 34, 30, 26, 22, 18, 9, 7, 7, 5, 5, 3],
    "noHT_1": [6, 6, 4, 4, 8, 8, 40, 36, 32, 28, 24, 20],
    "noHT_2": [5, 7, 5, 3, 9, 7, 38, 34, 30, 26, 22, 18],
}
UNMATCHED = 5
FILLER = "G" * 30
QUAL = "I" * 50


def main() -> None:
    guides = []
    with LIBRARY.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            guides.append((row["sgRNA"], row["sequence"].upper()))
    if len(guides) != 12:
        raise SystemExit(f"expected 12 guides in {LIBRARY}, found {len(guides)}")

    READS.mkdir(parents=True, exist_ok=True)
    rng = random.Random(16)
    for sample, counts in PROFILES.items():
        records = []
        n = 0
        for (name, seq), count in zip(guides, counts):
            if len(seq) != 20:
                raise SystemExit(f"{name} is not 20 bases")
            for _ in range(count):
                n += 1
                records.append(f"@{sample}.{n}\n{seq}{FILLER}\n+\n{QUAL}\n")
        for _ in range(UNMATCHED):
            n += 1
            junk = "".join(rng.choice("ACGT") for _ in range(20))
            records.append(f"@{sample}.{n}\n{junk}{FILLER}\n+\n{QUAL}\n")
        rng.shuffle(records)
        with gzip.open(READS / f"{sample}.fastq.gz", "wt", encoding="utf-8") as handle:
            handle.writelines(records)
        print(f"{sample}: {n} reads")


if __name__ == "__main__":
    main()
