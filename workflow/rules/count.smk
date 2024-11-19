rule hisat2_index:
    input:
        fasta = fasta,
    output:
        multiext(
            "resources/index/index",
            ".1.ht2",
            ".2.ht2",
            ".3.ht2",
            ".4.ht2",
            ".5.ht2",
            ".6.ht2",
            ".7.ht2",
            ".8.ht2",
        ),
    params:
        extra="",
        prefix = lambda wildcard, output: output[0].replace(".1.ht2", ""),
    log:
        "logs/hisat2/index.log"
    threads: 4
    resources:
        runtime=30
    wrapper:
        "v5.1.0/bio/hisat2/index"


rule count:
    input: 
        fq="results/trimmed/{sample}.fastq.gz",
        idx=multiext(
            "resources/index/index",
            ".1.ht2",
            ".2.ht2",
            ".3.ht2",
            ".4.ht2",
            ".5.ht2",
            ".6.ht2",
            ".7.ht2",
            ".8.ht2",
        ),
    output:
        "results/count/{sample}.guidecounts.txt"
    params:
        mm=config["mismatch"],
        idx=lambda wildcard, input: input.idx[0].replace(".1.ht2", ""),
    threads: 6
    resources:
        runtime=45
    log:
        "logs/count/{sample}.log"
    conda:
        "../envs/count.yaml"
    script:
        "../scripts/count.sh"


rule aggregate_counts:
    input:
        files=expand("results/count/{sample}.guidecounts.txt", sample=SAMPLES)
    output:
        "results/count/counts-aggregated.tsv"
    params:
        fa=fasta,
    threads: 1
    resources:
        runtime=10
    conda:
        "../envs/count.yaml"
    log:
        "logs/count/aggregate_counts.log"
    script:
        "../scripts/aggregate_counts.py"


rule normalise_count_table:
    input:
        counts="results/count/counts-aggregated.tsv"
    output:
        norm_counts=temp("results/count/counts-aggregated_normalised.csv")
    conda:
        "../envs/count.yaml"
    threads: 1
    resources:
        runtime=5
    log:
        "logs/count/normalise_counts.log"
    script:
        "../scripts/normalise_count_table.py"
