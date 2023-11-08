rule fastqc:
    input:
        "results/trimmed/{sample}.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    wrapper:
        "v2.0.0/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_fastqc.zip", sample=SAMPLES)
    output:
        "results/qc/multiqc.html",
    params:
        extra="",  # Optional: extra parameters for multiqc
    threads: config["resources"]["fastqc"]["cpu"]
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        "v2.9.0/bio/multiqc"


rule plot_alignment_rate:
    input:
        expand("logs/count/{sample}.log", sample=SAMPLES)
    output:
        report("results/qc/alignment-rates.pdf", caption="workflow/report/alignment-rates.rst", category="Alignment rates")
    log:
        "logs/plot-alignment-rate.log"
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot_alignment_rate.R"


rule plot_coverage:
    input:
        "results/count/counts-aggregated.tsv"
    output:
        report("results/qc/sequence-coverage.pdf", caption="workflow/report/plot-coverage.rst", category="Sequence coverage")
    params:
        fasta=fasta
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
        report("results/qc/gini-index.pdf", caption="workflow/report/gini-index.rst", category="Gini index")
    params:
        yaml="workflow/envs/plot_settings.yaml"
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
        report("results/qc/sample-correlation.pdf", caption="workflow/report/sample-correlation.rst", category="Sample correlation")
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
        report("results/qc/missed-rgrnas.pdf", caption="workflow/report/missed-rgrnas.rst", category="Missed sgRNAs")
    params:
        yaml="workflow/envs/plot_settings.yaml"
    log:
        "logs/missed-rgrnas.log"
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/missed_sgrnas.R"