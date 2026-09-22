rule cutadapt:
    input:
        "reads/{sample}.fastq.gz",
    output:
        fastq="results/trimmed/{sample}.fastq.gz",
        qc="results/trimmed/{sample}.qc.txt",
    log:
        "logs/cutadapt/{sample}.log",
    threads: 4  # set desired number of threads here
    resources:
        runtime=25,
    params:
        # Empty reads (after trimming) crash Bowtie; user arguments come last so they take precedence
        extra=f"--minimum-length 1 {config['cutadapt_args']}",
    wrapper:
        "v5.2.1/bio/cutadapt/se"
