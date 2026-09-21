# Snakemake workflow: `crispr-screens`

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10286661.svg)](https://doi.org/10.5281/zenodo.10286661)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/crispr-screens/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/crispr-screens/actions/workflows/main.yml)


<p align="center">
  <img src="docs/_static/logo2_small.png" width="200" alt="GPSW Logo" />
</p>


A Snakemake workflow for the analysis of CRISPR screens.

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

Instructions of how to use `crispr-screens` can be found here:

https://crispr-screens.readthedocs.io/en/latest/

## Annotating sgRNAs with genomic coordinates

CRISPRcleanR, which is always run before BAGEL2 (and optionally before MAGeCK and DrugZ), needs the genomic coordinates of every sgRNA. Some libraries do not provide them. The standalone script `annotate_sgrna_coordinates.py` (not part of the Snakemake workflow) finds them by searching the sgRNA sequences in a genome FASTA file and writes a CRISPRcleanR library file.

It needs Python with `pandas` and `numpy` (both are in the `stats` conda environment of this workflow) and three input files:

- a CSV file with the sgRNAs (any layout, you tell the script which columns to use)
- a genome FASTA file (use the primary assembly, `.gz` is fine)
- a GTF file with the gene annotation of the same assembly (`.gz` is fine)

### Usage

```shell
python annotate_sgrna_coordinates.py \
    --library resources/library.csv \
    --name-column "sgRNA" --gene-column "Gene" --sequence-column "sequence" \
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --gtf Homo_sapiens.GRCh38.112.gtf \
    --genome-fallback \
    --output resources/library_crisprcleanr.csv
```

Columns can be given by header name or by 0-based column number. sgRNA names must be unique (if the library has no names, the sequence column can be used for both). Run `python annotate_sgrna_coordinates.py --help` for all options:

| Option | Default | Description |
| --- | --- | --- |
| `--scope` | `locus` | `locus`: only search around annotated genes (fast, recommended); `genome`: search the whole genome (slow, no GTF needed) |
| `--genome-fallback` | off | Also search the whole genome for sgRNAs that were not found in any annotated locus. Needed to place sgRNAs that do not target genes, e.g. safe-targeting/intergenic controls |
| `--no-annotation-fallback` | off | Only accept a position in a locus of the sgRNA's own gene |
| `--feature` | `exon` | GTF feature to search around; use `gene` to search complete gene bodies |
| `--flank` | 30 | Bases added on both sides of each feature (should be at least the sgRNA length) |
| `--gene-attributes` | `gene_name,gene_id` | GTF attributes that the gene names of the library are matched against |
| `--report`, `--log` | next to output | Per-sgRNA mapping report and log file |

### Output

The output file has the columns `CODE,GENES,seq,CHRM,STARTpos,ENDpos,STRAND`, with all sgRNAs in the same order as in the input file:

- Coordinates are 1-based and inclusive and are those of the sgRNA sequence (without PAM).
- `STRAND` is the genomic strand: `+` if the sgRNA sequence is identical to the forward strand of the FASTA file, `-` if its reverse complement is. It is not relative to the gene. CRISPRcleanR does not use the strand.
- sgRNAs that cannot be placed are kept, with empty coordinates. `crisprcleaner.R` treats these as control sgRNAs (gene `CONTROL_GENE`, positioned on a made-up chromosome, so they are never corrected).
- `<output>_mapping_report.csv` lists for each sgRNA how it was placed (`mapped_gene_locus`, `mapped_other_locus`, `mapped_genome`) or why not (`gene_not_in_annotation`, `no_gene`, `sequence_not_found`, `invalid_sequence`) and the number of perfect matches (`n_hits`). If an sgRNA matches at several places, the first position is used, so consider `n_hits > 1` sgRNAs with care.
- Check the report for sgRNAs of real genes that could not be placed, as `crisprcleaner.R` will also treat these as controls.

To use the output in the workflow, set `lib_info: library_file` to the output file (`name_column: 0`, `gene_column: 1`, `sequence_column: 2`) and `stats: crisprcleanr: library_name` to any name that is not one of the CRISPRcleanR libraries.

### Why search around genes instead of the whole genome?

Only searching the exons of annotated genes (plus flanks) covers 177 Mb of the human genome instead of 3,100 Mb, which is roughly 10x faster, but it also gives better results: an sgRNA that matches at several places in the genome (paralogs, pseudogenes, repeats) is assigned to a locus of its own gene. In the test below, the whole genome search chose the intended copy for only 89.7% of the in-gene sgRNAs, against 99.7% for the locus search. sgRNAs of genes whose name is absent from the GTF (e.g. renamed genes) are still found, as the loci of all annotated genes are searched, and `--genome-fallback` finds the rest.

### Tests

All tests used the GRCh38 primary assembly and the Ensembl 112 GTF, on a single thread.

| Test | sgRNAs | Result |
| --- | --- | --- |
| Small synthetic genome (both strands, soft-masked bases, N run, exon boundary spanning sgRNAs, palindromic and duplicate sequences, sgRNAs of 18-19 nt, `chr1`/`1` and `chrM`/`MT` naming mismatches, gene given as name/ID/lower case, renamed genes, controls) | 20 | Correct coordinates for all sgRNAs that can be found in every search mode (`locus`, `genome`, `--genome-fallback`, `--no-annotation-fallback`, `--feature gene`), no wrong coordinates. Same output with gzipped input. |
| Simulated library with known coordinates (sequences extracted with `samtools faidx`): 59,948 in-gene sgRNAs (both strands, 7% spanning an exon boundary, 18-20 nt), 300 with renamed genes, 200 deep intronic, 500 controls | 60,948 | `locus` scope: all in-gene and renamed-gene sgRNAs placed; 99.68% of the in-gene sgRNAs at the exact position; 247 of the 253 placements that differ from the simulated position are sgRNAs with several perfect matches. 116 of 200 deep intronic sgRNAs and all controls unplaced, as expected. 2 min 21 s, 3.7 GB. |
| Same library, `--scope genome` | 60,948 | All 60,448 targeting sgRNAs placed (also the deep intronic ones), but only 89.7% of the in-gene sgRNAs at the intended position (9,190 in-gene sgRNAs have several matches, against 332 in the locus search). 23 min 3 s, 9.1 GB. |
| Brunello library from CRISPRcleanR (real library with GRCh38 coordinates and 2016 gene symbols), `--genome-fallback` | 76,379 | 100% placed (71,470 at their own gene, 4,851 through another gene's locus because of renamed genes, 58 by genome search). 99.52% agree with the CRISPRcleanR coordinates, once the constant offset is taken into account that comes from Brunello using the Cas9 cut site (start position 17 bases lower for `+` sgRNAs and 3 bases lower for `-` sgRNAs than in Brunello). The strand agrees for 71,497 of 71,498 sgRNAs after converting Brunello's `sense`/`antisense` (gene-relative) to genomic strand with the gene strand from the GTF. 5 min 47 s. |
| TKOv3 (annotated in hg19, so only the chromosome and strand can be compared) | 71,090 | 70,901 placed (65,928 at their own gene, 4,973 through another gene's locus); 189 unplaced (142 of these are the LacZ, luciferase and EGFP controls). 99.94% on the same chromosome and 99.42% on the same strand as the hg19 annotation. |
| Bassik library from `.test_mageck_test` (variable sgRNA length of 17-25 nt) | 211,696 | 83.5% at their own gene, 10.0% through another gene's locus, 0.2% not found. The 6.3% with a gene name that is not in the GTF are mostly the `safe` (6,697) and `none` (5,644) controls. 4 min 28 s, 3.8 GB. The `crisprcleanr` rule of this workflow ran to completion on the output (with `min_reads: 0`). |
| Jacquere library (no known coordinates), `--genome-fallback` | 60,550 | 99.8% placed (55,566 at their own gene, 3,279 through another gene's locus, 1,605 by genome search, including all 900 `ONE_SITE_INTERGENIC` controls). The 100 `NO_SITE` (non-targeting) controls are unplaced, as expected. 8 min 55 s, 3.7 GB. |

Note that the `crisprcleanr` rule fails on the small test data set of this repository with the default `min_reads: 10`, because only 86 sgRNAs pass the read filter, no matter which annotation is used.
