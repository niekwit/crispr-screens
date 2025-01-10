Conda
===========

For reproducibility, `crispr-screens` uses (containerized) Conda environments.

To install the `(Miniforge) <https://github.com/conda-forge/miniforge>`_  run on Unix-like platform:

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

On a local machine Apptainer can be `installed <https://github.com/apptainer/apptainer/blob/release-1.3/INSTALL.md>`_ as follows:

First, install all required dependencies:

On Debian-based systems, including Ubuntu:

.. code-block:: console

    $ sudo apt-get update
    $ sudo apt-get install -y build-essential libseccomp-dev pkg-config uidmap squashfs-tools fakeroot  cryptsetup tzdata dh-apparmor curl wget git

On CentOS/RHEL:

.. code-block:: console

    $ sudo yum groupinstall -y 'Development Tools'
    $ sudo yum install -y epel-release
    $ sudo yum install -y libseccomp-devel squashfs-tools fakeroot cryptsetup wget git

Finally, install Apptainer with conda:

.. code-block:: console

    $ conda activate snakemake
    $ conda install conda-forge::apptainer

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
   $ snakefetch --outdir . --repo-version v0.8.0 --url https://github.com/niekwit/crispr-screens
   Downloading archive file for version v0.8.0 from https://github.com/niekwit/crispr-screens...
   Extracting config and workflow directories from tar.gz file to /path/to/analysis...
   Done!

Alternatively, you can install the latest development version from source:

.. code-block:: console

    
    $ git clone https://github.com/niekwit/crispr-screens.git
    $ cd crispr-screens
    $ cp -r workflow /path/to/analysis

To download a specific version, go to the `release page <https://github.com/niekwit/crispr-screens/releases>`_ on GitHub and download the source code there.