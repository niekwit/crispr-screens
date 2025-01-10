rule create_fasta:
    input:
        csv=csv,
    output:
        fasta=fasta
    conda:
        "../envs/stats.yaml"
    log:
        "logs/create_fasta.log"
    script:
        "../scripts/csv_to_fasta.py"


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
        csv=csv,
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

if config["stats"]["bagel2"]["run"] or config["stats"]["mageck"]["apply_crisprcleanr"] or config["stats"]["drugz"]["apply_crisprcleanr"]:
    rule crisprcleanr:
        input:
            counts="results/count/counts-aggregated.tsv",
            fasta=fasta,
        output:
            # Input for MAGeCK/DrugZ (if run)
            corr_counts="results/count/crisprcleanr/corrected_counts_{comparison}.tsv",
            # Input for BAGEL2
            corr_lfc="results/count/crisprcleanr/corrected_lfc_{comparison}.foldchange",
            # QC plots
            roc="results/plots/crisprcleanr/roc_{comparison}.pdf",
            pr="results/plots/crisprcleanr/pr_{comparison}.pdf",
            drnk="results/plots/crisprcleanr/depletion_rank_{comparison}.pdf",
        params:
            lib_name=config["stats"]["crisprcleanr"]["library_name"],
            lib=config["stats"]["crisprcleanr"]["library_file"],
            control=lambda wc, output: wc.comparison.split("_vs_")[1].replace("-", ","),
            test=lambda wc, output: wc.comparison.split("_vs_")[0].replace("-", ","),
            ceg=config["stats"]["bagel2"]["custom_gene_lists"]["essential_genes"],
            cneg=config["stats"]["bagel2"]["custom_gene_lists"]["non_essential_genes"]
        conda:
            "../envs/stats.yaml"
        threads: 2
        resources:
            runtime=30
        log:
            "logs/crisprcleanr/{comparison}.log"
        script:
            "../scripts/crisprcleaner.R"