# DotMatch vs HISAT2 concordance on the public test data

Counts were generated from `.test_mageck_test` FASTQs after the same
`cutadapt -q 20 -l 20` trim the workflow uses. Bassik is mixed-length
(17–25 nt), so both methods were given the 75,680 length-20 guides only.
DotMatch was `dotmatch==0.2.2`, Hamming `k=0`, `ambiguity_policy=radius`,
window `0:20`. HISAT2 was 2.2.1.

## Workflow HISAT2 (`count.sh` flags)

`hisat2 --no-hd -N 0` with default spliced alignment and soft-clipping:

| sample | guides | HISAT2 assigned | DotMatch assigned | identical guides | Pearson | Spearman |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| HT_1 | 75680 | 63829 | 29833 | 53052 | 0.796 | 0.755 |
| HT_2 | 75680 | 64504 | 30190 | 53037 | 0.797 | 0.756 |
| noHT_1 | 75680 | 63346 | 29958 | 56941 | 0.976 | 0.783 |
| noHT_2 | 75680 | 64292 | 30494 | 57024 | 0.994 | 0.785 |

The extra HISAT2 assignments are not reverse-complement hits. After trimming,
HT_1 has 250,000 reads: 29,833 exact forward 20-mers, 0 exact reverse
complements, and 220,167 sequences that are not a library 20-mer. DotMatch
assigned those 29,833 exact windows and nothing else.

## HISAT2 restricted to the same contract

`hisat2 --no-hd --no-spliced-alignment --no-softclip --norc -N 0`:

| sample | guides | HISAT2 assigned | DotMatch assigned | identical guides | Pearson | Spearman |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| HT_1 | 75680 | 29833 | 29833 | 75680 | 1.000 | 1.000 |
| HT_2 | 75680 | 30190 | 30190 | 75680 | 1.000 | 1.000 |
| noHT_1 | 75680 | 29958 | 29958 | 75680 | 1.000 | 1.000 |
| noHT_2 | 75680 | 30496 | 30494 | 75678 | 1.000 | 1.000 |

On the 12-guide `.test_dotmatch` fixture the two methods also agreed exactly
(Pearson 1.0, 216/216 assigned reads).

Raw tables: `dotmatch-hisat2-concordance.csv` and
`dotmatch-hisat2-concordance-strict.csv`. Reproduce with
`workflow/scripts/compare_count_methods.py`; add `--strict` for the
restricted HISAT2 table.
