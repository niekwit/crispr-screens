#!/usr/bin/env bash


zcat ${snakemake_input[fq]} | hisat2 --no-hd -p ${snakemake[threads]} -t -N ${snakemake_params[mm]} -x ${snakemake_params[idx]} - 2> ${snakemake_log[0]} | sed '/XS:/d' | cut -f3 | sort | uniq -c | sed 's/^ *//' | sed '1d' > ${snakemake_output[0]}
