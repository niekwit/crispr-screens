if config["stats"]["mageck"]["apply_CNV_correction"]:
    logger.info("Applying CNV correction to MAGeCK analysis")
    check_sum = "33abf0c446f3e5116f35bda5483de28c3e504b1740457e1658b3fc2d22fbfa58"

    rule get_CNV_data:
        output:
            ensure("resources/cnv_data.txt.xz", sha256=check_sum),
        params:
            url="https://github.com/niekwit/crispr-screens/raw/main/.resources/CCLE_copynumber_byGene_2013-12-03.txt.xz",
        threads: 1
        resources:
            runtime=10,
        log:
            "logs/get_cnv_data.log",
        conda:
            "../envs/stats.yaml"
        shell:
            "wget -O {output} {params.url}  2> {log}"

    rule unpack_CNV_data:
        input:
            "resources/cnv_data.txt.xz",
        output:
            "resources/cnv_data.txt",
        threads: 1
        resources:
            runtime=10,
        log:
            "logs/unpack_cnv_data.log",
        conda:
            "../envs/stats.yaml"
        shell:
            "xz -dv {input} 2> {log}"


rule mageck_test:
    input:
        unpack(mageck_input),
    output:
        rnw="results/mageck/{comparison}/{cnv}/{comparison}_summary.Rnw",
        gs=report(
            "results/mageck/{comparison}/{cnv}/{comparison}.gene_summary.txt",
            caption="../report/mageck.rst",
            category="MAGeCK",
        ),
        ss="results/mageck/{comparison}/{cnv}/{comparison}.sgrna_summary.txt",
    params:
        command=config["stats"]["mageck"]["command"],
        control=mageck_control(),
        dir_name=lambda wc, output: os.path.dirname(output["gs"]),
        extra=extra_mageck_args(),
    threads: 2
    resources:
        runtime=20,
    conda:
        "../envs/stats.yaml"
    log:
        "logs/mageck/test/{comparison}_{cnv}.log",
    script:
        "../scripts/mageck.py"


rule create_mageck_mle_count_table:
    input:
        counts="results/count/counts-aggregated.tsv",
        matrix="config/{matrix}.txt",
    output:
        "results/count/counts-aggregated_{matrix}.tsv",
    threads: 1
    resources:
        runtime=5,
    conda:
        "../envs/stats.yaml"
    log:
        "logs/count/create_count_table_{matrix}.log",
    script:
        "../scripts/create_mageck_mle_count_table.py"


rule mageck_mle:
    input:
        unpack(mageck_input),
    output:
        gs="results/mageck/mle/{cnv}/{matrix}.gene_summary.txt",
        ss="results/mageck/mle/{cnv}/{matrix}.sgrna_summary.txt",
    params:
        command=config["stats"]["mageck"]["command"],
        control=mageck_control(),
        dir_name=lambda wc, output: os.path.dirname(output["gs"]),
        extra=extra_mageck_args(),
    threads: 24
    resources:
        runtime=45,
    conda:
        "../envs/stats.yaml"
    log:
        "logs/mageck/mle_{matrix}_{cnv}.log",
    script:
        "../scripts/mageck.py"


rule lfc_plots:
    input:
        "results/mageck/{comparison}/{cnv}/{comparison}.gene_summary.txt",
    output:
        pos=report(
            "results/plots/mageck/{comparison}/{cnv}/{comparison}.lfc_pos.pdf",
            caption="../report/lfc_pos.rst",
            category="MAGeCK plots",
            subcategory="{comparison}",
            labels={"Comparison": "{comparison}", "Figure": "lfc plot enriched genes"},
        ),
        neg=report(
            "results/plots/mageck/{comparison}/{cnv}/{comparison}.lfc_neg.pdf",
            caption="../report/lfc_neg.rst",
            category="MAGeCK plots",
            subcategory="{comparison}",
            labels={"Comparison": "{comparison}", "Figure": "lfc plot depleted genes"},
        ),
    threads: 1
    resources:
        runtime=5,
    conda:
        "../envs/stats.yaml"
    log:
        "logs/mageck_plots/lfc_{comparison}_{cnv}.log",
    script:
        "../scripts/plot_lfc.R"


rule sg_rank_plot:
    input:
        sg="results/mageck/{comparison}/{cnv}/{comparison}.sgrna_summary.txt",
        gene="results/mageck/{comparison}/{cnv}/{comparison}.gene_summary.txt",
    output:
        report(
            "results/plots/mageck/{comparison}/{cnv}/{comparison}.sgrank.pdf",
            caption="../report/sgrank.rst",
            category="MAGeCK plots",
            subcategory="{comparison}",
            labels={"Comparison": "{comparison}", "Figure": "sgrank plot"},
        ),
    threads: 1
    resources:
        runtime=5,
    conda:
        "../envs/stats.yaml"
    log:
        "logs/mageck_plots/sgrank_{comparison}_{cnv}.log",
    script:
        "../scripts/plot_sgrank.R"


rule gprofiler_mageck:
    input:
        txt="results/mageck/{comparison}/{cnv}/{comparison}.gene_summary.txt",
    output:
        csv="results/mageck/gprofiler/{comparison}/{cnv}/{pathway_data}.csv",
        pdf=report(
            "results/plots/mageck/gprofiler/{comparison}/{cnv}/{pathway_data}.pdf",
            caption="../report/pathway_analysis.rst",
            category="gprofiler plots",
            subcategory="{comparison}",
            labels={"Comparison": "{comparison}", "Figure": "pathway analysis"},
        ),
    params:
        fdr=config["stats"]["pathway_analysis"]["fdr"],
        top_genes=config["stats"]["pathway_analysis"]["top_genes"],
        data="mageck",
    threads: 1
    resources:
        runtime=10,
    conda:
        "../envs/stats.yaml"
    log:
        "logs/gprofiler/mageck/{comparison}_{cnv}_{pathway_data}.log",
    script:
        "../scripts/gprofiler.R"
