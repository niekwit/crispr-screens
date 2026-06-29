rule fastqc:
    input:
        "results/trimmed/{sample}.fastq.gz",
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip",
    log:
        "logs/fastqc/{sample}.log",
    threads: 4
    resources:
        runtime=15,
        mem_mb=2048,
    params:
        extra="--quiet",
    wrapper:
        "v5.2.1/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_fastqc.zip", sample=SAMPLES),
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="MultiQC",
        ),
    log:
        "logs/multiqc/multiqc.log",
    threads: 2
    resources:
        runtime=30,
        mem_mb=2048,
    params:
        extra="",  # Optional: extra parameters for multiqc
    wrapper:
        "v5.2.1/bio/multiqc"


rule plot_alignment_rate:
    input:
        expand("logs/count/{sample}.log", sample=SAMPLES),
    output:
        report(
            "results/qc/alignment-rates.pdf",
            caption="../report/alignment-rates.rst",
            category="Alignment rates",
        ),
        csv="results/qc/alignment-rates.csv",
    log:
        "logs/plot-alignment-rate.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    script:
        "../scripts/plot_alignment_rate.R"


rule plot_coverage:
    input:
        "results/count/counts-aggregated.tsv",
    output:
        report(
            "results/qc/sequence-coverage.pdf",
            caption="../report/plot-coverage.rst",
            category="Sequence coverage",
        ),
    log:
        "logs/plot-coverage.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    params:
        fasta=fasta,
    script:
        "../scripts/plot_coverage.R"


rule plot_gini_index:
    input:
        "results/count/counts-aggregated.tsv",
    output:
        report(
            "results/qc/gini-index.pdf",
            caption="../report/gini-index.rst",
            category="Gini index",
        ),
    log:
        "logs/gini-index.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    params:
        yaml="workflow/envs/plot_settings.yaml",
    script:
        "../scripts/plot_gini_index.R"


rule plot_missed_sgrnas:
    input:
        "results/count/counts-aggregated.tsv",
    output:
        report(
            "results/qc/missed-rgrnas.pdf",
            caption="../report/missed-rgrnas.rst",
            category="Missed sgRNAs",
        ),
    log:
        "logs/missed-rgrnas.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    params:
        yaml="workflow/envs/plot_settings.yaml",
    script:
        "../scripts/plot_missed_sgrnas.R"
