# .test_bagel2

A dedicated, small test fixture sized for BAGEL2, which `.test_mageck_test`
and `.test_mageck_mle` cannot exercise.

## Why this fixture exists

`.test_mageck_test`'s library (`resources/bassik.csv`) is genome-wide:
~211,700 guides across ~20,470 genes, with only 250,000 reads per sample —
about 1 read/guide on average. That is far too sparse for BAGEL2's Bayes
Factor regression, which needs real per-guide variance structure to fit,
and routinely fails there with a regression error (see the
`retries: 3  # Regression sometimes fails` comment already on the
`bagel2bf` rule). It's also blocked by an independent bug: `crisprcleanr.
library_name: TKOv3` in that fixture doesn't match any of CRISPRcleanR's
seven *bundled* library names (`AVANA_Library`, `Brunello_Library`,
`GeCKO_Library_v2`, `KY_Library_v1.0`, `KY_Library_v1.1`,
`MiniLibCas9_Library`, `Whitehead_Library` — see `crisprcleaner.R`), so it
falls through to reading annotation columns straight out of
`bassik.csv`, which doesn't have them.

Getting real genome-wide depth working would need tens of millions of
reads per sample, defeating the point of a fast CI fixture. Instead, this
fixture uses a small (~100 gene) subset with proper depth, so BAGEL2 (and
MAGeCK/DrugZ/STRING-db, also enabled here for full coverage on a fast
dataset) can run end to end without errors.

## How the data was built

Source: the real [hart-lab/bagel](https://github.com/hart-lab/bagel)
repository, which bundles:

- `CEGv2.txt` / `NEGv1.txt` — BAGEL2's own reference essential /
  non-essential gene lists.
- `pipeline-script-example/TKOv3_library_forCRISPRcleanR_REAGENT_ID.txt` —
  a real, full TKOv3 CRISPRcleanR annotation file (real guide sequences and
  real hg38 genomic coordinates).

`scripts/01_build_library.py`:

1. Picks 50 real essential genes (from `CEGv2.txt`) and 50 real
   non-essential genes (from `NEGv1.txt`) that have TKOv3 guide coverage,
   preferring genes on a handful of gene-dense chromosomes (1, 11, 17, 19).
   This keeps enough guides per chromosome for CRISPRcleanR's
   per-chromosome CBS smoothing, which crashes on very sparse chromosomes
   (the same underlying issue `crisprcleaner.R`'s non-standard-chromosome
   pooling already works around for unplaced contigs).
2. Writes `resources/tkov3_subset.csv`: real guide sequences and real
   genomic coordinates for those ~100 genes (~387 guides). This one file
   serves double duty as both the workflow's sgRNA library
   (`lib_info.library_file`; `sgRNA,Gene,sequence` columns) and
   CRISPRcleanR's annotation (`crisprcleaner.R` reads its annotation from
   that *same* file, not a separate one, via its own literally-named
   `seq,GENES,CODE,CHRM,STARTpos,ENDpos,STRAND` columns — `config.yml` sets
   `crisprcleanr.library_name` to a name that isn't one of the seven
   bundled ones, so it reads these columns directly).

`scripts/02_simulate_reads.py` then simulates four synthetic FASTQ samples
(`reads/T0_1.fastq.gz`, `T0_2.fastq.gz`, `T18_1.fastq.gz`, `T18_2.fastq.gz`
— named after BAGEL2's own bundled HAP1-TKOv3 example screen's T0/plasmid
vs T18/post-selection convention):

- **T0** (plasmid/reference, 2 replicates): each guide gets a baseline mean
  of ~400 reads, with a per-guide multiplicative lognormal factor so
  guides aren't perfectly uniform (real libraries never are).
- **T18** (post-selection, 2 replicates): essential-gene guides are given a
  real dropout signal — a per-*gene* target log2FC drawn from
  `Normal(-4.2, 0.6)`, plus per-*guide* noise (`Normal(0, 0.4)`, since not
  every guide targeting an essential gene is equally potent) — while
  non-essential genes stay roughly flat (`Normal(0, 0.3)` gene-level, same
  per-guide noise). This is deliberate: the fixture should catch a real
  regression in the biological result, not just verify the pipeline
  doesn't crash.
- Actual per-replicate counts are then sampled with additional
  replicate-level lognormal noise and a Gaussian draw around that mean, so
  there's real, non-degenerate variance for BAGEL2's regression to fit.
- Each read is the guide's real 20nt sequence followed by the standard
  sgRNA scaffold sequence immediately downstream of the protospacer in
  common lentiviral CRISPR vectors (e.g. lentiCRISPRv2/lentiGuide-Puro),
  padded to a 50bp read, with a uniform high-quality string (no simulated
  sequencing errors).

Both scripts use a fixed random seed (`20260929`), so re-running them
against a fresh clone of `hart-lab/bagel` reproduces this fixture exactly.

## Regenerating

```console
$ git clone https://github.com/hart-lab/bagel.git
$ cd scripts
$ python 01_build_library.py ../../bagel
$ python 02_simulate_reads.py ../../bagel
```
