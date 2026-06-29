if config["stats"]["mageck"]["apply_CNV_correction"]:
    logger.info("Applying CNV correction to MAGeCK analysis")
    check_sum = "33abf0c446f3e5116f35bda5483de28c3e504b1740457e1658b3fc2d22fbfa58"

    rule get_CNV_data:
        output:
            ensure("resources/cnv_data.txt.xz", sha256=check_sum),
        log:
            "logs/get_cnv_data.log",
        conda:
            "../envs/stats.yaml"
        threads: 1
        resources:
            runtime=10,
        params:
            url="https://github.com/niekwit/crispr-screens/raw/main/.resources/CCLE_copynumber_byGene_2013-12-03.txt.xz",
        shell:
            "wget -O {output} {params.url}  2> {log}"

    rule unpack_CNV_data:
        input:
            "resources/cnv_data.txt.xz",
        output:
            "resources/cnv_data.txt",
        log:
            "logs/unpack_cnv_data.log",
        conda:
            "../envs/stats.yaml"
        threads: 1
        resources:
            runtime=10,
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
    log:
        "logs/mageck/test/{comparison}_{cnv}.log",
    conda:
        "../envs/stats.yaml"
    threads: 2
    resources:
        runtime=20,
    params:
        command=config["stats"]["mageck"]["command"],
        control=mageck_control(),
        dir_name=lambda wc, output: os.path.dirname(output["gs"]),
        extra=extra_mageck_args(),
    script:
        "../scripts/mageck.py"


rule create_mageck_mle_count_table:
    input:
        counts="results/count/counts-aggregated.tsv",
        matrix="config/{matrix}.txt",
    output:
        "results/count/counts-aggregated_{matrix}.tsv",
    log:
        "logs/count/create_count_table_{matrix}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    script:
        "../scripts/create_mageck_mle_count_table.py"


rule mageck_mle:
    input:
        unpack(mageck_input),
    output:
        gs="results/mageck/mle/{cnv}/{matrix}.gene_summary.txt",
        ss="results/mageck/mle/{cnv}/{matrix}.sgrna_summary.txt",
    log:
        "logs/mageck/mle_{matrix}_{cnv}.log",
    conda:
        "../envs/stats.yaml"
    threads: 24
    resources:
        runtime=45,
    params:
        command=config["stats"]["mageck"]["command"],
        control=mageck_control(),
        dir_name=lambda wc, output: os.path.dirname(output["gs"]),
        extra=extra_mageck_args(),
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
    log:
        "logs/mageck_plots/lfc_{comparison}_{cnv}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
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
    log:
        "logs/mageck_plots/sgrank_{comparison}_{cnv}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
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
    log:
        "logs/gprofiler/mageck/{comparison}_{cnv}_{pathway_data}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=10,
    params:
        fdr=config["stats"]["pathway_analysis"]["fdr"],
        top_genes=config["stats"]["pathway_analysis"]["top_genes"],
        data="mageck",
    script:
        "../scripts/gprofiler.R"


rule string_db:
    input:
        txt="results/mageck/{comparison}/{cnv}/{comparison}.gene_summary.txt",
    output:
        svg="results/mageck/stringdb/{cnv}/{comparison}/{pathway_data}/pathway_analysis.svg",
        csv="results/mageck/stringdb/{cnv}/{comparison}/{pathway_data}/pathway_analysis.csv",
    log:
        "logs/stringdb/mageck/{cnv}/{comparison}_{pathway_data}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=10,
    script:
        "../scripts/string_db.py"
