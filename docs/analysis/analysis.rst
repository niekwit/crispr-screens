There are several options for analysing the sequencing data:

Snakemake workflow for CRISPR screens
=====================================

Our lab has developped a Snakemake_ pipeline for the analysis of CRISPR screens. The pipeline is available at https://github.com/niekwit/crispr-screens, along with its documentatio.

MAGeCK
======

MAGeCK_ (Li et al 2014) is a software package for identifying genes in CRISPR/Cas9 screens, but also has a built-in pipeline for the complete analysis (i.e. starting from raw fastq files).

The pipeline can be run as follows (in this case with the demo data included in the package):

.. code-block:: console

    $ mageck count -l library.txt -n demo --sample-label L1,CTRL  --fastq test1.fastq test2.fastq 
    $ mageck test -k demo.count.txt -t L1 -c CTRL -n demo


References
==========

#. Li W, Xu H, Xiao T, Cong L, Love MI, Zhang F, Irizarry RA, Liu JS, Brown M, Liu XS. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biol. 2014 Dec 19;15(12):554. doi: 10.1186/s13059-014-0554-4.






.. _Snakemake: https://snakemake.readthedocs.io
.. _MAGeCK: https://sourceforge.net/p/mageck/wiki/Home/

