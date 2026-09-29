"""
Simulate synthetic FASTQ reads for the .test_bagel2 fixture, from the
library built by 01_build_library.py.

Four samples: T0_1/T0_2 (plasmid/reference, control) and T18_1/T18_2
(post-selection, test), matching BAGEL2's own bundled HAP1-TKOv3 example
naming convention. Essential genes (real CEGv2 calls) are given a strong,
per-guide-noisy dropout in T18; non-essential genes (real NEGv1 calls) stay
roughly flat. See ../README.md for the full rationale.

Usage:
    python 02_simulate_reads.py bagel
"""

import csv
import gzip
import random
import sys

SEED = 20260929
BASELINE_MEAN = 400  # ~mean reads/guide in the plasmid/T0 samples
SAMPLES = ["T0_1", "T0_2", "T18_1", "T18_2"]
READ_LEN = 50
QUAL = "I" * READ_LEN  # uniform high quality (Phred+33 'I' = Q40)
# sgRNA scaffold sequence immediately following the 20nt protospacer in
# common lentiviral CRISPR vectors (e.g. lentiCRISPRv2/lentiGuide-Puro),
# used to pad synthetic reads out to a realistic 50bp read length.
SCAFFOLD = "GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"


def main(bagel_repo_dir):
    random.seed(SEED)

    ceg = set(
        l.split("\t")[0]
        for l in open(f"{bagel_repo_dir}/CEGv2.txt").read().splitlines()[1:]
        if l.strip()
    )

    guides = list(csv.DictReader(open("../resources/tkov3_subset.csv")))
    print(f"{len(guides)} guides loaded")

    counts = {s: {} for s in SAMPLES}
    for g in guides:
        gene, sg = g["Gene"], g["sgRNA"]
        is_essential = gene in ceg

        base_factor = random.lognormvariate(0, 0.35)
        base_mean = BASELINE_MEAN * base_factor

        if is_essential:
            gene_log2fc = random.gauss(-4.2, 0.6)
        else:
            gene_log2fc = random.gauss(0.0, 0.3)
        guide_log2fc = gene_log2fc + random.gauss(0, 0.4)

        for s in SAMPLES:
            mean = (
                base_mean
                if s.startswith("T0")
                else max(base_mean * (2**guide_log2fc), 0.5)
            )
            rep_mean = mean * random.lognormvariate(0, 0.2)
            counts[s][sg] = max(0, int(random.gauss(rep_mean, rep_mean**0.5)))

    for s in SAMPLES:
        total = sum(counts[s].values())
        print(f"{s}: total reads = {total}, mean/guide = {total / len(guides):.1f}")

    for s in SAMPLES:
        reads = []
        for g in guides:
            sg, seq = g["sgRNA"], g["sequence"]
            full_seq = (seq + SCAFFOLD)[:READ_LEN]
            reads.extend([full_seq] * counts[s][sg])
        random.shuffle(reads)

        path = f"../reads/{s}.fastq.gz"
        with gzip.open(path, "wt") as fh:
            for i, seq in enumerate(reads):
                fh.write(f"@synthetic_read_{i}\n{seq}\n+\n{QUAL}\n")
        print(f"wrote {path}: {len(reads)} reads")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <path to cloned hart-lab/bagel repo>")
    main(sys.argv[1])
