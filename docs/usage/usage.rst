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

.. note::
    
    The sgRNA sequence names must follow the following pattern: GENE_sgGENE_number (checked with regex: [A-Za-z0-9\\-]+_sg[A-Za-z0-9\\-]+_[0-9]+).

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
    │   └── bassik_lib.fasta
    └── workflow
        ├── envs
        │   └── stats.yaml
        ├── report
        │   ├── alignment-rates.rst
        │   ├── bagel2_plots.rst
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
        │   ├── count.smk
        │   ├── qc.smk
        │   ├── stats.smk
        │   └── trim.smk
        ├── schemas
        │   ├── config.schema.yaml
        │   └── stats.schema.yaml
        ├── scripts
        │   ├── aggregate_counts.py
        │   ├── bagel2bf.py
        │   ├── bagel2fc.py
        │   ├── bagel2pr.py
        │   ├── convert_count_table.py
        │   ├── count.sh
        │   ├── general_functions.smk
        │   ├── lfc_plots.R
        │   ├── mageck_plots.R
        │   ├── mageck.py
        │   ├── missed_sgrnas.R
        │   ├── normalise_count_table.py
        │   ├── pathway_analysis.R
        │   ├── plot_alignment_rate.R
        │   ├── plot_bf.R
        │   ├── plot_coverage.R
        │   ├── plot_gini_index.R
        │   ├── plot_pr.R
        │   └── sgrank_plot.R
        └── Snakefile

9 directories, 47 files


Experiment meta data
====================

Experiment meta data is described in `config/config.yml`:

.. code-block:: yaml

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
    
        bagel2:
            run: True # Perform bagel2 analysis
            custom_gene_lists: # Paths to custom gene lists for bagel2 analysis
            essential_genes: none
            non_essential_genes: none
        
        mageck:
            run: True # Perform mageck analysis
            extra_mageck_arguments: "" 
            mageck_control_genes: all # All or file with control genes
            fdr: 0.25 # FDR threshold for downstream mageck analysis
            apply_CNV_correction: False # Apply CNV correction to mageck results
            cell_line: K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE # Cell line for CNV correction
        
        drugz:
            run: False # Perform drugZ analysis
            extra: "" # Extra arguments for drugZ

        pathway_analysis: 
            run: True # Perform pathway analysis on mageck results
            data: both # enriched, depleted, or both
            dbs: ["GO_Molecular_Function_2023","GO_Biological_Process_2023","Reactome_2022"]
            top_genes: 50 # Number of top genes to consider for pathway analysis (overrides fdr, use 0 to disable)
            terms: 10 # Number of terms to plot



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

