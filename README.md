# Snakemake workflow: `crispr-screens`

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10286662.svg)](https://doi.org/10.5281/zenodo.10286662)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.12.0-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/crispr-screens/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/crispr-screens/actions/workflows/main.yml)
[![CodeFactor](https://www.codefactor.io/repository/github/niekwit/crispr-screens/badge)](https://www.codefactor.io/repository/github/niekwit/crispr-screens)


A Snakemake workflow for the analysis of CRISPR screens.

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).


## Installation of required software

### Conda/Mamba and Apptainer

For reproducibility, `crispr-screen` uses (containerized) Conda environments.

Please follow the instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for a detailed guide to install Conda/Mamba, and optionally, but highly recommended, [Apptainer](https://apptainer.org/docs/user/main/quick_start.html).

### Snakemake

To install Snakemake create the following environment with `mamba`:

```shell
$ mamba create -n snakemake snakemake
```

Activate the environment as follows:

```shell
$ mamba activate snakemake
```

If you want to deploy Snakemake on an HPC system using slurm also run:

```shell
$ pip install snakemake-executor-plugin-slurm
```

### Workflow code

The easiest way to obtain the workflow code is to use [snakefetch](https://pypi.org/project/snakefetch/):

```shell
$ pip install snakefetch
$ snakefetch --outdir /path/to/analysis --repo-version v0.6.0 --url https://github.com/niekwit/crispr-screens
Downloading archive file for version v0.6.0 from https://github.com/niekwit/crispr-screens...
Extracting config and workflow directories from tar.gz file to /path/to/analysis...
Done!
```

This will copy the config and workflow directories to the path set with the `--outdir` flag.


## Preparing raw sequencing data

In the directory containing config/workflow create a directory called reads:

```shell
$ cd /path/to/analysis
$ mdkir -p reads 
```

All fastq files can be copied here. Only single end data is processed.

> [!IMPORTANT]  
> All files must have the extension *fastq.gz*


## sgRNA library sequences

In the directory containing config/workflow create a directory called resources:

```shell
$ mdkir -p resources
```
Copy the fasta file with sgRNA names as a header to this directory. 

If you only have a text file containing the sgRNA names and sequences in different columns, this can be used instead (csv format).

> [!IMPORTANT]  
> All sgRNA names must be in the following format: GENE_sgGENE_sgRNAnumber (e.g. B2M_sgB2M_1).

If no fasta file has been generated yet, but a csv file is available that contains sgRNA sequences and gene names in separate column, then this can be copied into the resources directory instead. The workflow will generate the fasta file from this automatically. The columns that contains this information must be included in config.yaml under "csv".

## sgRNA library meta data and analysis settings

Analysis settings are stored in config/config.yml:

```yaml
lib_info:
    sg_length: 20
    vector: "N" # Vector sequence to be removed from reads (N for none)
    left_trim : 0 # Trim n bases from 5' end of reads
    species: human
csv: # If no fasta is available, provide a csv file with sgRNA sequences
  name_column: 3 # Column number with gene names
  sequence_column: 2 # Column number with sgRNA sequences
mismatch: 0 # Mismatches allowed during alignment
stats: 
  skip: none # Skip mageck, bagel2, both, or none
  mageck:
    extra_mageck_arguments: "" 
    mageck_control_genes: all # All or file with control genes
    fdr: 0.25 # FDR threshold for downstream mageck analysis
    apply_CNV_correction: False # Apply CNV correction to mageck results
    cell_line: K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE # Cell line for CNV correction
  pathway_analysis: 
      run: True # Perform pathway analysis on mageck results
      data: both # enriched, depleted, or both
      dbs: ["GO_Molecular_Function_2023","GO_Biological_Process_2023","Reactome_2022"]
      top_genes: 50 # Number of top genes to consider for pathway analysis (overrides fdr, use 0 to disable)
      terms: 10 # Number of terms to plot
resources:
  trim:
    cpu: 4
    time: 60
  fastqc:
    cpu: 4
    time: 60
  count:
    cpu: 8
    time: 120
  stats:
    cpu: 2
    time: 60

```

If your sgRNA lengths are uniform, set that length (*l*) with `sg_length` and set `vector` to "N". This will only keep the the first *l* nucleotides of the reads during trimming.

When sgRNA lengths vary within the library, then provide the vector sequence (i.e. the DNA sequence (spacer/scaffold) immediately downstream of the last nucleotide of the sgRNA sequence). This will trim off that sequence from the reads during trimming.

With the `left_trim` option, one can set the number of bases the be trimmed at the 5' end of the read. This is done before any other trimming events.


## Sample comparisons for MAGeCK and/or BAGEL2

In config/stats.csv the pairwise comparisons between samples for MAGeCK and/or BAGEL2 can be defined:

test | control | bagel2
--- | --- | ---
day18_1 | plasmid | n
day18_2 | plasmid | n
day18_1;day18_2 | plasmid | y

Test and control samples should be included in the test and control columns, respectively. Replicate conditions can be put in the sample column, separated by a semi-colon (;). 

If you want to disable BAGEL2 analysis of some comparisons, add "n" to the bagel2 column.


## Configuration of Snakemake

Running Snakemake can entail quite a few command line flags. To make this easier these can be set in a global profile that is defined in a user-specific configuration directory in order to simplify this process.

For example, a profile `config.yaml` can be stored at /home/user/.config/snakemake/profile:
```yaml
cores: 40
latency-wait: 20
use-conda: True
use-apptainer: True
keep-going: False
rerun-incomplete: True
printshellcmds: True
cache: True
show-failed-logs: True
```

When running on a slurm-based HPC, the following lines can be included in `config.yaml`:
```yaml
executor: slurm
jobs: 100
apptainer-args: "--bind '/parent_dir/of/analysis'" # if analysis in not in /home/$USER
default-resources:
        slurm_partition: icelake
        slurm_account: <ACCOUNT>
```

Some system have limited space allocated to `/tmp`, which can be problematic when using Apptainer. Add the following line to `~/.bashrc` to set a different temporary directory location:

```shell
export APPTAINER_TMPDIR=~/rds/hpc-work/apptainer_tmp
```

## Dry-run of the analysis

Before running the actual analyis with your own data, a dry-run can be performed:

```shell
$ cd path/to/analysis/directory
$ snakemake -np
```

Snakemake will create the DAG of jobs and print the shell command, but it will not execute anything.

## Visualization of workflow

To visualize the workflow run (this command excludes the target rule):
```shell
$ mkdir -p images
$ snakemake --forceall --rulegraph | grep -v '\-> 0\|0\[label = \"all\"' | dot -Tpng > images/rule_graph.png
```

## Running the analysis

Once you know that the test and/or dry run has worked, the actual analysis can be initiated as follows:
```shell
$ snakemake --profile /home/user/.config/snakemake/profile --directory .test/
```

> [!IMPORTANT]  
> Always make sure to use the absolute path (i.e. /home/user/.config/...) rather than the relative path (~/.config/...) when providing the path for the profile file.

## Report of the results

When the analysis has finished succesfully, an HTML report can be created as follows:

```shell
$ snakemake --report report.html
```

This report will contain run time information for the Snakemake rules, as well as figures generated by the workflow, and the code used to create these.
