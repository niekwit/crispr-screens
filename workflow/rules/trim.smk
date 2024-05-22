rule cutadapt:
    input:
        "reads/{sample}.fastq.gz",
    output:
        fastq=temp("results/trimmed/{sample}.fastq.gz"),
        qc="results/trimmed/{sample}.qc.txt",
    params:
        extra=f"-q 20 {cut_adapt_arg(config)}",
    log:
        "logs/cutadapt/{sample}.log",
    threads: config["resources"]["trim"]["cpu"]  # set desired number of threads here
    resources:
        runtime=config["resources"]["trim"]["time"]
    wrapper:
        f"{wrapper_version}/bio/cutadapt/se"