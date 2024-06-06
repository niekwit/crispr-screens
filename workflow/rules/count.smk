rule hisat2_index:
    input:
        fasta = fasta,
    output:
        directory(f"resources/index/"),
    params:
        extra="",
        prefix = "resources/index/index",
    log:
        "logs/hisat2/index.log"
    threads: config["resources"]["count"]["cpu"]
    resources:
        runtime=config["resources"]["count"]["time"]
    wrapper:
        f"{wrapper_version}/bio/hisat2/index"


rule count:
    input: 
        fq="results/trimmed/{sample}.fastq.gz",
        idx=f"resources/index/",
    output:
        "results/count/{sample}.guidecounts.txt"
    params:
        mm=config["mismatch"],
        idx="resources/index/index",
    threads: config["resources"]["count"]["cpu"]
    resources:
        runtime=config["resources"]["count"]["time"]
    log:
        "logs/count/{sample}.log"
    conda:
        "../envs/count.yaml"
    shell:
        "zcat {input.fq} | "
        "hisat2 --no-hd -p {threads} -t -N {params.mm} -x {params.idx} - 2> {log} | "
        "sed '/XS:/d' | "
        "cut -f3 | "
        "sort | "
        "uniq -c | "
        "sed 's/^ *//' | "
        "sed '1d' > {output}"


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
        conda:
            "../envs/count.yaml"
        log:
            "logs/count/normalise_counts.log"
        script:
            "../scripts/normalise_count_table.py"
