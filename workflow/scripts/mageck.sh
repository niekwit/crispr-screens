#!/usr/bin/env bash

mageck test --normcounts-to-file -k ${snakemake_input[0]} -t $(echo ${snakemake_wildcards[mcomparison]} | awk -F '_vs_' '{print $1}' | sed 's/-/,/') -c $(echo ${snakemake_wildcards[mcomparison]} | awk -F '_vs_' '{print $2}' | sed 's/-/,/') -n results/mageck/${snakemake_wildcards[mcomparison]}/${snakemake_wildcards[mcomparison]} ${snakemake_params[control]} ${snakemake_params[extra]} 2> ${snakemake_log[0]}
