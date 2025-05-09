Conda
===========

For reproducibility, `crispr-screens` uses (containerized) Conda environments.

To install the `(Miniforge) <https://github.com/conda-forge/miniforge>`_  run on a Unix-like platform:

.. code-block:: console

    $ curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    $ bash Miniforge3-$(uname)-$(uname -m).sh


Snakemake
=========

To install Snakemake create the following environment with `conda`:

.. code-block:: console

    $ conda create -n snakemake snakemake=8.25.5

Activate the environment as follows:

.. code-block:: console

    $ conda activate snakemake

If you want to deploy Snakemake on an HPC system using slurm also run:

.. code-block:: console

    $ pip install snakemake-executor-plugin-slurm


.. note::
   Depending on your HPC/cloud storage you may also need to install a storage plugin. For more information on this see: https://snakemake.github.io/snakemake-plugin-catalog/index.html.


Apptainer
=========

For (optional) containerization `crispr-screens` uses `Apptainer <https://apptainer.org>`_. If this is applied, Snakemake will pull a `Docker image <https://hub.docker.com/repository/docker/niekwit/crispr-screens/general>`_ and convert it to an Apptainer image. This image contains all pre-made Conda environments.


Installation on a local machine
-------------------------------

The easiest, and recommended, way to install `Apptainer` on a local Machine is via `Conda`:

.. code-block:: console

    $ conda install -c conda-forge apptainer



Installation on HPC systems
---------------------------

As the installation of Apptainer requires administrator privileges, please contact your HPC administrator to install Apptainer if it is not available yet.


Snakefetch
=======================================

The easiest way to obtain the workflow code is to use `(snakefetch) <https://pypi.org/project/snakefetch/>`_:

.. code-block:: console

    $ pip install snakefetch


Obtaining workflow code
=======================

To obtain the workflow code, run snakefetch from the analysis directory:

.. code-block:: console

   $ mkdir my_experiment
   $ cd my_experiment
   $ snakefetch --outdir . --repo-version v0.8.2 --url https://github.com/niekwit/crispr-screens
   Downloading archive file for version v0.8.2 from https://github.com/niekwit/crispr-screens...
   Extracting config and workflow directories from tar.gz file to /path/to/analysis...
   Done!

Alternatively, you can install the latest development version from source:

.. code-block:: console

    
    $ git clone https://github.com/niekwit/crispr-screens.git
    $ cd crispr-screens
    $ cp -r workflow /path/to/analysis

To download a specific version, go to the `release page <https://github.com/niekwit/crispr-screens/releases>`_ on GitHub and download the source code there.