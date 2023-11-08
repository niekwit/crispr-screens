rule hisat2_index: 
    output:
        files=multiext(idx_prefix, 
                ".1.ht2",
                ".2.ht2",
                ".3.ht2",
                ".4.ht2",
                ".5.ht2",
                ".6.ht2",
                ".7.ht2",
                ".8.ht2"),
    params:
        fa=fasta,
        idx=idx_prefix,
    threads: config["resources"]["count"]["cpu"]
    resources:
        runtime=config["resources"]["count"]["time"]
    log:
        "logs/hisat2/index.log"
    conda:
        "../envs/count.yaml"
    shell:
        "hisat2-build -p {threads} {params.fa} {params.idx} 2> {log}"


rule count:
    input: 
        fq="results/trimmed/{sample}.fastq.gz",
        idx_files=multiext(idx_prefix, 
                ".1.ht2",
                ".2.ht2",
                ".3.ht2",
                ".4.ht2",
                ".5.ht2",
                ".6.ht2",
                ".7.ht2",
                ".8.ht2"),
    output:
        "results/count/{sample}.guidecounts.txt"
    params:
        mm=config["mismatch"],
        idx=idx_prefix,
    threads: config["resources"]["count"]["cpu"]
    resources:
        runtime=config["resources"]["count"]["time"]
    log:
        "logs/count/{sample}.log"
    conda:
        "../envs/count.yaml"
    shell:
        "zcat {input.fq} | hisat2 --no-hd -p {threads} -t -N {params.mm} -x {params.idx} - 2> {log} | "
        "sed '/XS:/d' | cut -f3 | sort | uniq -c | sed 's/^ *//' | sed '1d' > {output}"


rule aggregated_counts:
    input:
        files=expand("results/count/{sample}.guidecounts.txt", sample=SAMPLES)
    output:
        "results/count/counts-aggregated.tsv"
    params:
        fa=fasta,
    conda:
        "../envs/count.yaml"
    log:
        "logs/count/aggregate_counts.log"
    script:
        "../scripts/join.py"


rule normalise_count_table:
        input:
            counts="results/count/counts-aggregated.tsv"
        output:
            norm_counts=temp("results/count/counts-aggregated_normalised.csv")
        run:
            df = pd.read_table(input.counts)
            column_range = range(2,len(df.columns))
            for i in column_range:
                column_sum = df.iloc[:,i].sum()
                df.iloc[:,i] = df.iloc[:,i] / column_sum * 1E8
                df.iloc[:,i] = df.iloc[:,i].astype(int)
            df.to_csv(output.norm_counts,
                    index = False,
                    header = True)


