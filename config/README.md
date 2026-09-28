# CONFIGURATION

## config.yml

Use config.yml to provide information about your experiment. See `workflow/schemas/config.schema.yaml` for the full, authoritative list of keys.

### lib_info

Under `lib_info`, set `library_file` to the path of the sgRNA library CSV file (in the `resources` folder) and `species` to the species the library targets (e.g. `human`).

### cutadapt_args

Extra arguments passed to `cutadapt` when trimming the raw reads, as a single string (e.g. `"-q 20 -l 20"`). `--minimum-length 10` is always added by the workflow and does not need to be included. See the [cutadapt manual](https://cutadapt.readthedocs.io/en/stable/) for all options.

### csv

A CSV file (in the `resources` folder) must be provided that contains the sgRNA names, sequences, and gene names, each in its own column. Under `name_column`, `sequence_column`, and `gene_column` (0-based), set the column numbers of these columns. If no fasta file is available for building the Bowtie index, one is generated from this CSV file.

### bowtie_args

Extra arguments passed to `bowtie` (v1) when aligning the trimmed reads to the sgRNA index. The read file (`-q`), index (`-x`) and number of threads (`-p`) are set by the workflow and should not be included here.

The default (`-v 1 -m 1`) allows one mismatch (`-v`) and discards reads that align to more than one sgRNA (`-m 1`). Note that Bowtie aligns end-to-end, so reads are expected to be trimmed to the length of the sgRNAs (see `cutadapt_args`). Reads that are longer than the sgRNA they originate from (e.g. in libraries with variable sgRNA lengths) will not align unless they are trimmed accordingly. Use `-v 0` to only allow perfect matches. See the [Bowtie manual](https://bowtie-bio.sourceforge.net/manual.shtml) for all options.

### stats

Each statistical tool (`bagel2`, `mageck`, `drugz`) has its own `run: True`/`False` switch to enable or disable it.

- `crisprcleanr`: settings for the CRISPRcleanR normalisation step. It always runs ahead of BAGEL2, and optionally ahead of MAGeCK/DrugZ (see `apply_crisprcleanr` under those sections). `library_name` is either the name of one of the libraries bundled with CRISPRcleanR (see the comments in `config.yml` for the full list) or any other name, in which case a CRISPRcleanR-formatted library file must be provided (see the main `README.md`, e.g. via `annotate_sgrna_coordinates.py`, for how to generate one). `min_reads` sets the minimum read count in the control sample for an sgRNA to be kept.
- `bagel2`: `custom_gene_lists.essential_genes`/`non_essential_genes` can point to custom essential/non-essential gene list files (`none` uses BAGEL2's own default lists). `extra_args.bf`/`pr` pass extra arguments to the BAGEL2 `bf` and `pr` subcommands respectively.
- `mageck`: `command` is `test` (pairwise, needs `config/stats.csv`) or `mle` (needs one or more design matrices, see `mle.design_matrix`; each file must be placed in the `config` directory). `extra_mageck_arguments` passes extra arguments to the MAGeCK `test`/`mle` command. `mageck_control_genes` is `all` or a path to a file with control gene names, one per line, used to build the null distribution/normalisation instead of the whole library. `apply_CNV_correction` and `cell_line` enable copy-number correction of MAGeCK results.
- `drugz`: `extra` passes extra arguments to the `drugz` command.
- `pathway_analysis`/`string_db`: run pathway (g:Profiler) and/or STRING-db enrichment analysis on the MAGeCK results. `data` selects `enriched`, `depleted`, or `both` gene sets, `fdr` sets the significance threshold, and `top_genes` (if not 0) overrides `fdr` and takes the top N genes instead.

### stats.csv

Pairwise comparisons for MAGeCK (`command: test`), DrugZ, and BAGEL2 are defined in `config/stats.csv`, with `test` and `control` columns naming samples exactly as they appear in the `reads` directory (without the `.fastq.gz`/`.cram` extension). Multiple replicate samples can be combined in one comparison by separating their names with a semicolon. An optional `bagel2_only` column (`y`/`n` per row) can be added to run some comparisons only through BAGEL2 (and CRISPRcleanR) and the rest only through MAGeCK/DrugZ; if the column is absent, every row is available to the enabled tools. See the main documentation for details and examples.

### ngs_tracker

Optional integration with [NGS Tracker](https://github.com/niekwit/ngs-tracker) to register the workflow run and attach output files. Set `enabled: false` to skip this entirely. When enabled, `base_url` and `project_id` identify the NGS Tracker instance/project, and `files` lists the paths (glob patterns allowed) and types of files to attach after a successful run.

### resources

Computational resources (threads, runtime, memory) are set per rule in the workflow itself, not in `config.yml`. To override them, use a Snakemake profile (`--set-resources`/`--set-threads`) or edit the relevant rule.
