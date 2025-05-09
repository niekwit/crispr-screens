Preparing analysis data directory
=================================

To run `crispr-screens` on your own data, first prepare a data directory (where the workflow code is located) as follows:

.. code-block:: console

    $ cd my_experiment
    $ mkdir reads resources

Copy the fastq files to the `reads` directory and the resources to the `resources` directory. 

.. note::
    
    All fastq files must have the extension *.fastq.gz* or *.cram*


Copy the the fasta file with sgRNA sequences to the `resources` directory. If no fasta file is available, provide a csv file with sgRNA sequences and gene names in separate columns. Set the column numbers in `config/config.yml` (see below). 

The final directory structure should look like this:

.. code-block:: console

    .
    ├── config
    │   ├── config.yml
    │   └── stats.csv
    ├── reads
    │   ├── HT_1.fastq.gz
    │   ├── HT_2.fastq.gz
    │   ├── noHT_1.fastq.gz
    │   └── noHT_2.fastq.gz
    ├── resources
    │   └── bassik.csv
    └── workflow
        ├── envs
        │   └── stats.yaml
        ├── report
        │   ├── alignment-rates.rst
        │   ├── bagel2_plots.rst
        │   ├── drugz.rst
        │   ├── gini-index.rst
        │   ├── lfc_neg.rst
        │   ├── lfc_pos.rst
        │   ├── mageck.rst
        │   ├── missed-rgrnas.rst
        │   ├── multiqc.rst
        │   ├── pathway_analysis.rst
        │   ├── plot-coverage.rst
        │   ├── sample-correlation.rst
        │   ├── sgrank.rst
        │   └── workflow.rst
        ├── rules
        │   ├── bagel2.smk
        │   ├── count.smk
        │   ├── drugz.smk
        │   ├── mageck.smk
        │   ├── qc.smk
        │   └── trim.smk
        ├── schemas
        │   ├── config.schema.yaml
        │   └── stats.schema.yaml
        ├── scripts
        │   ├── aggregate_counts.py
        │   ├── bagel2bf.py
        │   ├── bagel2pr.py
        │   ├── cnv_cell_lines.txt
        │   ├── count.sh
        │   ├── crisprcleaner.R
        │   ├── csv_to_fasta.py
        │   ├── general_functions.smk
        │   ├── gprofiler.R
        │   ├── mageck.py
        │   ├── plot_alignment_rate.R
        │   ├── plot_bf.R
        │   ├── plot_coverage.R
        │   ├── plot_drugz_results.R
        │   ├── plot_gini_index.R
        │   ├── plot_lfc.R
        │   ├── plot_missed_sgrnas.R
        │   ├── plot_pr.R
        │   └── plot_sgrank.R
        └── Snakefile

    9 directories, 50 files


Experiment meta data
====================

Experiment meta data is described in `config/config.yml`:

.. code-block:: yaml

    lib_info:
    library_file: resources/bassik.csv # Path to library file with sgRNA sequences and gene names

    cutadapt:
        g: "" # 5' adapter sequence to trim
        a: "" # 3' adapter sequence to trim
        u: 0 # trim u bases (before a/g trimming)
        l: 20 # shorten reads to l bases
        extra: "" # Extra arguments for cutadapt

    species: human

    csv: 
        # 0-based column numbers
        name_column: 0 # Column number with sgRNA names 
        gene_column: 1 # Column number with gene names
        sequence_column: 2 # Column number with sgRNA sequences

    mismatch: 0 # Mismatches allowed during alignment

    stats: 
    crisprcleanr:
        # For BAGEL2, crisprcleanr is always run as it allows for combining replicates better
        # It is optional for MAGeCK and DrugZ
        # Path to library file for crisprcleanr with sgRNA annotations
        # With column names: GENE,seq,CODE,CHRM,STARTpos,ENDpos,EXONE(optional),STRAND
        # If library name is one of the following, the library info is loaded from the crisprcleanr library database:
        # AVANA_Library (https://doi.org/10.1038/ng.3984)
        # Brunello_Library (https://doi.org/10.1038/nbt.3437)
        # GeCKO_Library_v2 (https://doi.org/10.1038/nmeth.3047)
        # KY_Library_v1.0 (https://doi.org/10.1016/j.celrep.2016.09.079)
        # KY_Library_v1.1 (https://doi.org/10.1016/j.celrep.2016.09.079)
        # MiniLibCas9_Library (https://doi.org/10.1186/s13059-021-02268-4)
        # Whitehead_Library (https://doi.org/10.1126/science.aac7041)
        library_name: TKOv3
        min_reads: 30 # Keep sgRNAs with at least this many reads in control sample

    bagel2:
        run: True # Perform bagel2 analysis
        custom_gene_lists: 
        # Paths to custom gene lists for bagel2 analysis
        # Use "none" to use BAGEL2 default gene lists
            essential_genes: none
            non_essential_genes: none
        extra_args: # Extra arguments for bagel2 subcommands
        bf: ""
        pr: ""

    mageck:
        run: True # Perform mageck analysis
        command: test # test or mle
        mle:
            design_matrix: ["config/matrix.txt"] # Design matrix for mageck mle
        # It is recommended to disable crisprcleanr when using non-genome-wide sgRNA libraries
        apply_crisprcleanr: False 
        extra_mageck_arguments: "" 
        mageck_control_genes: all # All or file with control genes
        apply_CNV_correction: False # Apply CNV correction to mageck results
        cell_line: K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE # Cell line for CNV correction
  
  drugz:
    run: True # Perform drugZ analysis
    # It is recommended to disable crisprcleanr when using non-genome-wide sgRNA libraries
    apply_crisprcleanr: False
    extra: "" # Extra arguments for drugZ

  pathway_analysis: 
    run: False # Perform pathway analysis on mageck results
    data: both # enriched, depleted, or both
    fdr: 0.25 # FDR threshold for significant genes
    top_genes: 50 # Number of top genes to consider for pathway analysis (overrides fdr, use 0 to disable)


When CRISPRcleanR is not applied, the library csv file will have to contain columns for just the gene names and sgRNA sequences and names, e.g.:

.. csv-table:: Example library csv file
    :header: "sgRNA", "Gene", "sequence"

    ENSG00000121410_A1BG_PROT_195964.1,A1BG,GCTGACGGGTGACACCCA
    ENSG00000121410_A1BG_PROT_195965.2,A1BG,GACTTCCAGCTGTTCAAGAA
    ENSG00000121410_A1BG_PROT_195966.3,A1BG,GCAGGTGAGTCAAGGTGCAC
    ENSG00000121410_A1BG_PROT_195967.4,A1BG,GCCGCTCGGGCTTGTCCAC

However, when CRISPRcleanR is applied, the library csv file will have to be constructed as follows (exon column is optional):

.. csv-table:: Example library csv file for use with CRISPRcleanR
    :header: CODE, GENES, seq, CHRM, STARTpos, ENDpos, EXONE, STRAND

    "chr19\:58864777\-58864796\_A1BG\_\+",A1BG,CAAGAGAAAGACCACGAGCA,chr19,58864777,58864796,ex1,"\+"
    "chr19\:58864319\-58864338\_A1BG\_\+",A1BG,GCTCAGCTGGGTCCATCCTG,chr19,58864319,58864338,ex3,"\+"
    "chr19\:58863885\-58863904\_A1BG\_\+",A1BG,ACTGGCGCCATCGAGAGCCA,chr19,58863885,58863904,ex4,"\+"
    "chr19\:58862759\-58862778\_A1BG\_\-",A1BG,GTCGAGCTGATTCTGAGCGA,chr19,58862759,58862778,ex5,"\-"

In both cases, the column numbers (0-based) for the gene names, sgRNA sequences, and sgRNA names must be set in `config/config.yml` (under the csv section).

.. note::
    
    Some libraries are available in the CRISPRcleanR library database. If the library name is mentioned in the `config/config.yml` file, the library info is loaded from the database, if that name is set as the library name.


BAGEL2 analysis
===============

`BAGEL2 <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00809-3>`_ can be used to perform gene essentiality analysis. It requires a file containing the pair-wise comparisons of the samples. The file should be in the `config` directory and have the following format:

.. csv-table:: Example stats.csv file
    :header: "test", "control"

    HT_1,noHT_1
    HT_2,noHT_2
    HT_1;HT_2,noHT_1;noHT_2

.. note::

    Replicate samples can be used by separating the sample names with a semicolon (see above). The sample names must match the fastq file names in the `reads` directory.


`CRISPRcleanR <https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4989-y>`_ is used to create a count table as input for BAGEL2.

BAGEL2 can be customized as follows:

.. code-block:: yaml

    bagel2:
        run: True # Perform bagel2 analysis
        custom_gene_lists: 
        # Paths to custom gene lists for bagel2 analysis
        # Use "none" to use BAGEL2 default gene lists
            essential_genes: none
            non_essential_genes: none
        extra_args: # Extra arguments for bagel2 subcommands
            bf: ""
            pr: ""

Custom essential and non-essential gene lists that BAGEL2 requires for Bayes Factor calculation can be provided in the `resources` directory. The files should be in the following format:

.. csv-table:: Example essential gene list

    A1BG
    A1CF
    A2M
    A3GALT2
    A4GALT
    A4GNT
    A4GNTL1
    A4GNTL2

If no custom gene lists are provided, the default gene lists from BAGEL2 will be used. This can be done by setting the `essential_genes` and `non_essential_genes` parameters to `none`.

Extra arguments for the BAGEL2 `bf` and `pr` subcommands can be provided in the `extra_args` section.


MAGeCK analysis
===============

`MAGeCK <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0554-4>`_ can be used to perform multiple types of analyses: RRA and MLE.


RRA
---

To perform pair-wise analyses, MAGeCK RRA can be run with the `test` command. When the `test` command is used, a file (config/stats.csv) with these pairwise comparisons of the samples must be provided. The file format is described in the BAGEL2 section.


MLE
---

For more complex designs, the `mle` command can be used. A design matrix must be provided in the `config` directory. An example of a design matrix (tab-delimited) is shown below. :

.. csv-table:: Example design matrix
    :header: "Samples", "baseline", "noHT_vs_HT"

    HT_1,1,0
    HT_2,1,0
    noHT_1,1,1
    noHT_2,1,1


The MAGeCK analysis can be customized as follows:

.. code-block:: yaml

    mageck:
        run: True # Perform mageck analysis
        command: test # test or mle
        mle:
            design_matrix: ["config/matrix.txt"] # Design matrix for mageck mle
        # It is recommended to disable crisprcleanr when using non-genome-wide sgRNA libraries
        apply_crisprcleanr: False 
        extra_mageck_arguments: "" 
        mageck_control_genes: all # All or file with control genes
        apply_CNV_correction: False # Apply CNV correction to mageck results
        cell_line: K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE # Cell line for CNV correction

CRISPRcleanR can be used to create a normalised count table as input for MAGeCK. This will disable normalisation by MAGeCK.

A list of control genes can be provided in the `resources` directory. The file should be in the following format:

.. csv-table:: Example control gene list

    A1BG
    A1CF
    A2M
    A3GALT2
    A4GALT
    A4GNT
    A4GNTL1
    A4GNTL2

These genes will then be used for normalisation and for generating the null distribution of RRA.

Extra arguments for the MAGeCK `test` and `mle` commands can be provided in the `extra_mageck_arguments` section. 


drugZ analysis
==============

`drugZ <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0665-3>`_ can be used to perform chemogenetic analysis of CRISPR screens. It requires a file containing the pair-wise comparisons of the samples. The file format is described in the BAGEL2 section.

The drugZ analysis can be customized as follows:

.. code-block:: yaml

    drugz:
        run: True # Perform drugZ analysis
        # It is recommended to disable crisprcleanr when using non-genome-wide sgRNA libraries
        apply_crisprcleanr: False
        extra: "" # Extra arguments for drugZ

CRISPRcleanR can be used to create a normalised count table as input for MAGeCK.

Extra arguments for the drugZ command can be provided in the `extra` section.


Setup global Snakemake profile
==============================

To setup a profile for custom command line arguments, create a new profile (`config.yaml`) in `$HOME/.config/snakemake/standard/`:

.. code-block:: yaml

    cores: 40
    latency-wait: 10
    use-conda: True # Recommended
    rerun-incomplete: True
    printshellcmds: True
    show-failed-logs: True
    use-apptainer: True # Recommended
    apptainer-args: "--bind '/parent_dir/of/analysis'" # If analysis is not in /home/$USER
    # For execution on a SLURM cluster add:
    executor: slurm
    jobs: 100 # Maximum number of jobs to run in parallel
    local-cores: 2 # Limit core usage for local rules
    default-resources:
        slurm_partition: icelake
        slurm_account: <ACCOUNT>


With the above settings, `Snakemake` will download a `Docker` image from `Docker Hub` and convert it to an `Apptainer` image. In this image all required `Conda` environments are pre-installed. The `--bind` option is required when the analysis directory is not in `/home/$USER`. If the analysis directory is in `/home/$USER`, this option can be omitted.


Dry-run on experimental data
============================

To test if the workflow is working with your data and settings, run a dry-run:

.. code-block:: console

    $ snakemake -np


Execution of workflow
=====================

To execute the workflow:

.. code-block:: console

    $ snakemake --profile $HOME/.config/snakemake/standard/


Report of analysis
==================

When the workflow has finished, a report of the results can be generated (HTML format):

.. code-block:: console

    $ snakemake --report report.html

