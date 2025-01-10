Preparing analysis data directory
=================================

To run `crispr-screens` on your own data, first prepare a data directory (where the workflow code is located) as follows:

.. code-block:: console

    $ cd my_experiment
    $ mkdir reads resources

Copy the fastq files to the `reads` directory and the resources to the `resources` directory. 

.. note::
    
    All fastq files must have the extension *.fastq.gz*


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

    10 directories, 53 files


Experiment meta data
====================

Experiment meta data is described in `config/config.yml`:

.. code-block:: yaml

    lib_info:
        library_file: resources/bassik.csv
        sg_length: 17
        vector: "N" # Vector sequence to be removed from reads (N for none)
        left_trim : 1 # Trim n bases from 5' end of reads
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

    bagel2:
        run: False # Perform bagel2 analysis
        custom_gene_lists: # Paths to custom gene lists for bagel2 analysis
        essential_genes: none
        non_essential_genes: none
    
    mageck:
        run: True # Perform mageck analysis
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


Setup global Snakemake profile
==============================

To setup a profile for custom command line arguments, create a new profile (`config.yaml`) in `$HOME/.config/snakemake/standard/`:

.. code-block:: yaml

    cores: 40
    latency-wait: 10
    use-conda: True
    rerun-incomplete: True
    printshellcmds: True
    cache: False
    show-failed-logs: True
    use-apptainer: True


Dry-run on experimental data
============================

To test if the workflow is defined properly and to estimate the amount of required computational resources, run a dry-run:

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

