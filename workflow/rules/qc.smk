rule fastqc:
    input:
        "results/trimmed/{sample}.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip",
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}.log"
    threads: 4
    resources:
        runtime=15,
        mem_mb = 2048,
    wrapper:
        f"{wrapper_version}/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_fastqc.zip", sample=SAMPLES)
    output:
        report("results/qc/multiqc.html", caption="../report/multiqc.rst", category="MultiQC"),
    params:
        extra="",  # Optional: extra parameters for multiqc
    threads: 2
    resources:
        runtime=30,
        mem_mb = 2048,
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        f"{wrapper_version}/bio/multiqc"


rule plot_alignment_rate:
    input:
        expand("logs/count/{sample}.log", sample=SAMPLES)
    output:
        report("results/qc/alignment-rates.pdf", caption="../report/alignment-rates.rst", category="Alignment rates")
    log:
        "logs/plot-alignment-rate.log"
    threads: 1
    resources:
        runtime=5
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot_alignment_rate.R"


rule plot_coverage:
    input:
        "results/count/counts-aggregated.tsv"
    output:
        report("results/qc/sequence-coverage.pdf", caption="../report/plot-coverage.rst", category="Sequence coverage")
    params:
        fasta=fasta
    threads: 1
    resources:
        runtime=5
    log:
        "logs/plot-coverage.log"
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot_coverage.R"


rule plot_gini_index:
    input:
        "results/count/counts-aggregated.tsv"
    output:
        report("results/qc/gini-index.pdf", caption="../report/gini-index.rst", category="Gini index")
    params:
        yaml="workflow/envs/plot_settings.yaml"
    threads: 1
    resources:
        runtime=5
    log:
        "logs/gini-index.log"
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot_gini_index.R"


rule plot_sample_correlation:
    input:
        "results/count/counts-aggregated_normalised.csv"
    output:
        report("results/qc/sample-correlation.pdf", caption="../report/sample-correlation.rst", category="Sample correlation")
    threads: 1
    resources:
        runtime=5
    log:
        "logs/sample-correlation.log"
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/sample_correlation.R"


rule plot_missed_sgrnas:
    input:
        "results/count/counts-aggregated.tsv"
    output:
        report("results/qc/missed-rgrnas.pdf", caption="../report/missed-rgrnas.rst", category="Missed sgRNAs")
    params:
        yaml="workflow/envs/plot_settings.yaml"
    threads: 1
    resources:
        runtime=5
    log:
        "logs/missed-rgrnas.log"
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/missed_sgrnas.R"