Running test data
=================

To run `crispr-screens` on test data, first clone the repository:

.. code-block:: console

    $ git clone https://github.com/niekwit/crispr-screens.git
    $ cd crispr-screens
    $ mamba activate snakemake
    $ snakemake --use-conda --use-apptainer --show-failed-logs --directory .test/
    $ snakemake --report .test/report.html --directory .test/


Preparing analysis data directory
=================================

To run `crispr-screens` on your own data, first prepare a data directory as follows:

.. code-block:: console

    $ mkdir my_experiment
    $ cd my_experiment
    $ mkdir reads resources

Copy the fastq files to the `reads` directory and the resources to the `resources` directory. 

.. note::
    
    All fastq files must have the extension *fastq.gz*


Obtaining workflow code
=======================

To obtain the workflow code, run snakefetch from my_experiment/:

.. code-block:: console

    $ snakefetch --outdir . --repo-version v0.6.1 --url https://github.com/niekwit/crispr-screens
    Downloading archive file for version v0.6.1 from https://github.com/niekwit/crispr-screens...
    Extracting config and workflow directories from tar.gz file to /path/to/analysis...
    Done!


Experiment meta data
====================

TO DO

Dry-run on experimental data
============================

If the test data runs successfully, you can run `crispr-screens` on your own data. To do this, first create a new directory and copy the test data:

.. code-block:: console

    $ mkdir my_experiment
    $ snakemake -np


