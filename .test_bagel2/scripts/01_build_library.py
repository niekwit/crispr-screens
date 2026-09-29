"""
Build resources/tkov3_subset.csv for the .test_bagel2 fixture.

Selects a small, real subset of the TKOv3 library (real guide sequences and
genomic coordinates, real essential/non-essential gene calls) from data
bundled in the hart-lab/bagel repository, and writes ONE combined CSV that
serves as both:
  - the sgRNA library for the workflow itself (lib_info.library_file; sgRNA,
    Gene, sequence columns, indices set via config csv.*_column)
  - CRISPRcleanR's annotation library (crisprcleaner.R reads its annotation
    from the SAME file as lib_info.library_file, not a separate one, via its
    literally-named seq/GENES/CODE/CHRM/STARTpos/ENDpos/STRAND columns)

See ../README.md for why this fixture exists and how the genes were chosen.

Usage:
    git clone https://github.com/hart-lab/bagel.git
    python 01_build_library.py bagel
"""

import sys
from collections import defaultdict

SEED = 20260929
N_ESSENTIAL = 50
N_NON_ESSENTIAL = 50
# Gene-dense chromosomes to concentrate the selection on, so each ends up
# with enough guides for CRISPRcleanR's per-chromosome CBS smoothing (very
# sparse chromosomes crash it - see crisprcleaner.R's non-real-chr pooling
# comment for the same underlying issue with unplaced contigs).
DENSE_CHROMS = {"1", "11", "17", "19"}


def main(bagel_repo_dir):
    import random

    random.seed(SEED)

    ceg = set(
        l.split("\t")[0]
        for l in open(f"{bagel_repo_dir}/CEGv2.txt").read().splitlines()[1:]
        if l.strip()
    )
    neg = set(
        l.split("\t")[0]
        for l in open(f"{bagel_repo_dir}/NEGv1.txt").read().splitlines()[1:]
        if l.strip()
    )

    # Real TKOv3 CRISPRcleanR annotation (real guide seqs + coordinates).
    # Columns are [unnamed row-id, CODE, GENES, EXONE, CHRM, STRAND,
    # STARTpos, ENDpos] where the unnamed first column is the actual unique
    # per-guide reagent ID ("{GENE}_{20nt sequence}") - CODE/GENES here are
    # NOT unique per guide (both just hold the gene symbol), so the guide
    # sequence is recovered by stripping the "{gene}_" prefix off the
    # reagent ID instead.
    lib_path = (
        f"{bagel_repo_dir}/pipeline-script-example/"
        "TKOv3_library_forCRISPRcleanR_REAGENT_ID.txt"
    )
    rows = []
    with open(lib_path) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 8:
                continue
            reagent_id, _code, gene, _exon, chrm, strand, start, end = parts
            if not reagent_id.startswith(gene + "_"):
                continue
            seq = reagent_id[len(gene) + 1 :]
            if len(seq) < 15 or not set(seq) <= set("ACGT"):
                continue
            rows.append(
                {
                    "id": reagent_id,
                    "gene": gene,
                    "chrm": chrm,
                    "strand": strand,
                    "start": int(start),
                    "end": int(end),
                    "seq": seq,
                }
            )

    gene_guides = defaultdict(list)
    gene_chrom = {}
    for r in rows:
        gene_guides[r["gene"]].append(r)
        gene_chrom[r["gene"]] = r["chrm"]

    ess_avail = [g for g in ceg if g in gene_guides]
    non_avail = [g for g in neg if g in gene_guides]

    ess_pref = sorted(g for g in ess_avail if gene_chrom[g] in DENSE_CHROMS)
    non_pref = sorted(g for g in non_avail if gene_chrom[g] in DENSE_CHROMS)
    random.shuffle(ess_pref)
    random.shuffle(non_pref)

    ess_genes = sorted(ess_pref[:N_ESSENTIAL])
    non_genes = sorted(non_pref[:N_NON_ESSENTIAL])
    print(f"selected {len(ess_genes)} essential / {len(non_genes)} non-essential genes")

    guides = []
    for g in ess_genes + non_genes:
        guides.extend(gene_guides[g])
    print("total guides:", len(guides))

    with open("../resources/tkov3_subset.csv", "w") as fh:
        fh.write("sgRNA,Gene,sequence,seq,GENES,CODE,CHRM,STARTpos,ENDpos,STRAND\n")
        for r in guides:
            fh.write(
                f"{r['id']},{r['gene']},{r['seq']},{r['seq']},{r['gene']},{r['id']},"
                f"{r['chrm']},{r['start']},{r['end']},{r['strand']}\n"
            )
    print("wrote ../resources/tkov3_subset.csv")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <path to cloned hart-lab/bagel repo>")
    main(sys.argv[1])
