rule create_fasta:
    input:
        csv=csv,
    output:
        fasta=fasta,
    log:
        "logs/create_fasta.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/csv_to_fasta.py"


if cram():

    rule cram2fastq:
        input:
            "reads/{sample}.cram",
        output:
            "reads/{sample}.fastq.gz",
        log:
            "logs/samtools-fastq/{sample}.interleaved.log",
        conda:
            "../envs/stats.yaml"
        threads: 1
        params:
            " ",
        shell:
            "samtools fastq {input} 2> {log} | gzip -c > {output}"


rule bowtie_index:
    input:
        fasta=fasta,
    output:
        multiext(
            "resources/index/index",
            ".1.ebwt",
            ".2.ebwt",
            ".3.ebwt",
            ".4.ebwt",
            ".rev.1.ebwt",
            ".rev.2.ebwt",
        ),
    log:
        "logs/bowtie/index.log",
    conda:
        "../envs/stats.yaml"
    threads: 4
    resources:
        runtime=30,
    params:
        prefix=lambda wildcard, output: output[0].replace(".1.ebwt", ""),
    shell:
        "bowtie-build --threads {threads} {input.fasta} {params.prefix} > {log} 2>&1"


rule count:
    input:
        fq="results/trimmed/{sample}.fastq.gz",
        idx=multiext(
            "resources/index/index",
            ".1.ebwt",
            ".2.ebwt",
            ".3.ebwt",
            ".4.ebwt",
            ".rev.1.ebwt",
            ".rev.2.ebwt",
        ),
    output:
        "results/count/{sample}.guidecounts.txt",
    log:
        "logs/count/{sample}.log",
    conda:
        "../envs/stats.yaml"
    threads: 6
    resources:
        runtime=45,
    params:
        extra=config["bowtie_args"],
        idx=lambda wildcard, input: input.idx[0].replace(".1.ebwt", ""),
    shell:
        "zcat {input.fq} | "
        "bowtie -q {params.extra} -p {threads} -x {params.idx} - 2> {log} | "
        "cut -f3 | sort | uniq -c | sed 's/^ *//' > {output}"


rule aggregate_counts:
    input:
        files=expand("results/count/{sample}.guidecounts.txt", sample=SAMPLES),
        csv=csv,
    output:
        "results/count/counts-aggregated.tsv",
    log:
        "logs/count/aggregate_counts.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=10,
    script:
        "../scripts/aggregate_counts.py"


if (
    config["stats"]["bagel2"]["run"]
    or config["stats"]["mageck"]["apply_crisprcleanr"]
    or config["stats"]["drugz"]["apply_crisprcleanr"]
):

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
        log:
            "logs/crisprcleanr/{comparison}.log",
        conda:
            "../envs/stats.yaml"
        threads: 2
        resources:
            runtime=30,
        params:
            lib_name=config["stats"]["crisprcleanr"]["library_name"],
            lib=csv,
            control=lambda wc, output: wc.comparison.split("_vs_")[1].replace("-", ","),
            test=lambda wc, output: wc.comparison.split("_vs_")[0].replace("-", ","),
            ceg=config["stats"]["bagel2"]["custom_gene_lists"]["essential_genes"],
            cneg=config["stats"]["bagel2"]["custom_gene_lists"]["non_essential_genes"],
            min_reads=config["stats"]["crisprcleanr"]["min_reads"],
        script:
            "../scripts/crisprcleaner.R"
