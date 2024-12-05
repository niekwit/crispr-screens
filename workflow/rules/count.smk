rule hisat2_index:
    input:
        fasta=fasta,
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
        "v5.2.1/bio/hisat2/index"


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
        "../envs/stats.yaml"
    script:
        "../scripts/count.sh"


rule aggregate_counts:
    input:
        files=expand("results/count/{sample}.guidecounts.txt", sample=SAMPLES),
        fasta=fasta,
    output:
        "results/count/counts-aggregated.tsv"
    threads: 1
    resources:
        runtime=10
    conda:
        "../envs/stats.yaml"
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
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5
    log:
        "logs/count/normalise_counts.log"
    script:
        "../scripts/normalise_count_table.py"


if config["stats"]["crisprcleanr"]["run"] or config["stats"]["bagel2"]["run"]:
    rule crisprcleanr:
        input:
            counts="results/count/counts-aggregated.tsv",
        output:
            # Input for MAGeCK (if run)
            corr_counts="results/count/crisprcleanr/corrected_counts_{bcomparison}.tsv",
            # Input for BAGEL2
            corr_lfc="results/count/crisprcleanr/corrected_lfc_{bcomparison}.foldchange",
            # QC plots
            roc="results/plots/crisprcleanr/roc_{bcomparison}.pdf",
            pr="results/plots/crisprcleanr/pr_{bcomparison}.pdf",
            drnk="results/plots/crisprcleanr/depletion_rank_{bcomparison}.pdf",
        params:
            lib_name=config["stats"]["crisprcleanr"]["library_name"],
            lib=config["stats"]["crisprcleanr"]["library_file"],
            control=lambda wc, output: wc.comparison.split("_vs_")[1].replace("-", ","),
            test=lambda wc, output: wc.comparison.split("_vs_")[0].replace("-", ","),
        conda:
            "../envs/stats.yaml"
        threads: 2
        resources:
            runtime=30
        log:
            "logs/count/crisprcleanr_{comparison}.log"
        script:
            "../scripts/crisprcleanr.R"