rule link_fastq_for_fastqc:
    # FastQC embeds the name of the file it was given inside its report, and
    # MultiQC reads that embedded name (not the output file name) to tell
    # samples apart. Raw and trimmed reads share the same file name
    # ({sample}.fastq.gz in different directories), so without this rename
    # step MultiQC would treat pre- and post-trimming FastQC reports for the
    # same sample as one sample and silently drop one of them.
    input:
        fastqc_input,
    output:
        temp("results/qc/fastqc_input/{sample}_{stage}.fastq.gz"),
    log:
        "logs/fastqc/link_{sample}_{stage}.log",
    wildcard_constraints:
        stage="raw|trimmed",
    threads: 1
    resources:
        runtime=2,
    shell:
        "ln -sr {input} {output} 2> {log}"


rule fastqc:
    input:
        "results/qc/fastqc_input/{sample}_{stage}.fastq.gz",
    output:
        html="results/qc/fastqc/{sample}_{stage}.html",
        zip="results/qc/fastqc/{sample}_{stage}_fastqc.zip",
    log:
        "logs/fastqc/{sample}_{stage}.log",
    wildcard_constraints:
        stage="raw|trimmed",
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
        expand(
            "results/qc/fastqc/{sample}_{stage}_fastqc.zip",
            sample=SAMPLES,
            stage=["raw", "trimmed"],
        ),
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
        "v9.18.0/bio/multiqc"


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


rule plot_count_distribution:
    input:
        "results/count/counts-aggregated.tsv",
    output:
        report(
            "results/qc/count-distribution.pdf",
            caption="../report/count-distribution.rst",
            category="Count distribution",
        ),
        csv="results/qc/count-distribution.csv",
    log:
        "logs/count-distribution.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    script:
        "../scripts/plot_count_distribution.R"


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
