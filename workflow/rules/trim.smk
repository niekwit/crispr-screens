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
    threads: 4  # set desired number of threads here
    resources:
        runtime=25,
    wrapper:
        "v5.2.1/bio/cutadapt/se"
