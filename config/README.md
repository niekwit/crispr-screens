# CONFIGURATION

## config.yml

Use config.yml to provide information about your experiment:

### lib_info

Under `lib_info` information about the sgRNA library must be provided. If your library has a fixed sgRNA length, then provide this with *sg_length*. This way the reads in the raw data will be trimmed from the 3' end until the read length is equal to *sg_length*

Some libraries have variable sgRNA length. In this case provide a vector sequence under *vector*, so that this sequence is removed, instead of trimming to a fixed length.

Finally, with some sequencing strategies, the first base sequenced can be the same for all sgRNA sequences. In this case the base calling can be of very poor quality. By setting `left_trim` to 1, one can remove the first 5' base for all reads in order to improve the quality of the read. This step will be performed before 3' end trimming or removal of the vector sequence.

### csv

Inside the resources folder, a *fasta* file should be placed that contains unique sgRNA names and sequences, which will be used to build an index for alignment using Bowtie.

This Snakemake workflow can be run without a fasta file, as long as a CSV file (also in the resources folder) is provided that contains the unique sgRNA sequences and corresponding gene names in separate columns. Under `name_column` the column number of the gene names, and under `sequence_column` the column number of the sequence column have to be set.

### bowtie_args

Extra arguments passed to `bowtie` (v1) when aligning the trimmed reads to the sgRNA index. The read file (`-q`), index (`-x`) and number of threads (`-p`) are set by the workflow and should not be included here.

The default (`-v 1 -m 1`) allows one mismatch (`-v`) and discards reads that align to more than one sgRNA (`-m 1`). Note that Bowtie aligns end-to-end, so reads are expected to be trimmed to the length of the sgRNAs (see `cutadapt_args`). Reads that are longer than the sgRNA they originate from (e.g. in libraries with variable sgRNA lengths) will not align unless they are trimmed accordingly. Use `-v 0` to only allow perfect matches. See the [Bowtie manual](https://bowtie-bio.sourceforge.net/manual.shtml) for all options.

### stats

With `skip` one can skip statistical analyses with MAGeCK, and/or BAGEL2.

Any extra argument to be parsed to MAGeCK can be defined with `extra_mageck_arguments`.

Normally MAGeCK builds the statistical model using all sgRNAs in the library. However, sometimes this is not optimal when many sgRNAs change, for example when using smaller libraries. In this case, the user can provide a file, with each gene name on a new line, that contain genes that should be used to build the model instead. 

### resources

Under resources the computational requirements can be set. The CPU count will be used locally and on any HPC/cloud platform, while the time will only be relevant for the latter.

