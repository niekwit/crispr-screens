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
        # MultiQC's current template defaults to dark mode; force light mode
        # since this report also gets embedded (as-is, no runtime theme
        # override) into qc-report.html, which is otherwise light-themed.
        extra='--cl-config "template_dark_mode: false"',
    wrapper:
        "v9.18.0/bio/multiqc"


rule qc_report:
    input:
        logs=expand("logs/count/{sample}.log", sample=SAMPLES),
        counts="results/count/counts-aggregated.tsv",
        multiqc="results/qc/multiqc.html",
        crisprcleanr_roc=expand(
            "results/crisprcleanr/roc_{comparison}.csv",
            comparison=crisprcleanr_comparisons(),
        ),
        crisprcleanr_pr=expand(
            "results/crisprcleanr/pr_{comparison}.csv",
            comparison=crisprcleanr_comparisons(),
        ),
        crisprcleanr_drnk_genes=expand(
            "results/crisprcleanr/depletion_rank_{comparison}_genes.csv",
            comparison=crisprcleanr_comparisons(),
        ),
        crisprcleanr_drnk_signatures=expand(
            "results/crisprcleanr/depletion_rank_{comparison}_signatures.csv",
            comparison=crisprcleanr_comparisons(),
        ),
    output:
        html=report(
            "results/qc-report.html",
            caption="../report/qc-report.rst",
            category="QC report",
        ),
        csv="results/qc/alignment-rates.csv",
        count_csv="results/qc/count-distribution.csv",
    log:
        "logs/qc-report.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=10,
    params:
        fasta=fasta,
        crisprcleanr_comparisons=crisprcleanr_comparisons(),
        sorted_samples=sorted_samples(),
    script:
        "../scripts/qc_report.py"
