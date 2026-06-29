rule install_bagel2:
    output:
        directory("resources/bagel2"),
    log:
        "logs/bagel2/install.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    shell:
        "git clone https://github.com/hart-lab/bagel.git {output} 2> {log}"


rule bagel2bf:
    input:
        b2dir="resources/bagel2/",
        fc="results/count/crisprcleanr/corrected_lfc_{comparison}.foldchange",
    output:
        bf="results/bagel2/{comparison}/{comparison}.bf",
    log:
        command="logs/bagel2/bf/{comparison}_command.log",
        stdout="logs/bagel2/bf/{comparison}_stdout.log",
        stderr="logs/bagel2/bf/{comparison}_stderr.log",
    retries: 3  # Regression sometimes fails
    conda:
        "../envs/stats.yaml"
    threads: 2
    resources:
        runtime=15,
    params:
        species=config["lib_info"]["species"],
        ceg=config["stats"]["bagel2"]["custom_gene_lists"]["essential_genes"],
        cneg=config["stats"]["bagel2"]["custom_gene_lists"]["non_essential_genes"],
        extra=config["stats"]["bagel2"]["extra_args"]["bf"],
    script:
        "../scripts/bagel2bf.py"


rule bagel2pr:
    input:
        b2dir="resources/bagel2/",
        bf="results/bagel2/{comparison}/{comparison}.bf",
    output:
        "results/bagel2/{comparison}/{comparison}.pr",
    log:
        command="logs/bagel2/pr/{comparison}_command.log",
        stdout="logs/bagel2/pr/{comparison}_stdout.log",
        stderr="logs/bagel2/pr/{comparison}_stderr.log",
    conda:
        "../envs/stats.yaml"
    threads: 2
    resources:
        runtime=15,
    params:
        species=config["lib_info"]["species"],
        ceg=config["stats"]["bagel2"]["custom_gene_lists"]["essential_genes"],
        cneg=config["stats"]["bagel2"]["custom_gene_lists"]["non_essential_genes"],
        extra=config["stats"]["bagel2"]["extra_args"]["pr"],
    script:
        "../scripts/bagel2pr.py"


rule plot_bf:
    input:
        "results/bagel2/{comparison}/{comparison}.bf",
    output:
        report(
            "results/plots/bagel2/{comparison}/{comparison}.bf.pdf",
            caption="../report/bagel2_plots.rst",
            category="BAGEL2 plots",
            subcategory="{comparison}",
            labels={"Comparison": "{comparison}", "Figure": "BF plot"},
        ),
    log:
        "logs/bagel2/plot/bf_{comparison}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    script:
        "../scripts/plot_bf.R"


rule plot_pr:
    input:
        "results/bagel2/{comparison}/{comparison}.pr",
    output:
        report(
            "results/plots/bagel2/{comparison}/{comparison}.pr.pdf",
            caption="../report/bagel2_plots.rst",
            category="BAGEL2 plots",
            subcategory="{comparison}",
            labels={"Comparison": "{comparison}", "Figure": "Precision-recall plot"},
        ),
    log:
        "logs/bagel2/plot/pr_{comparison}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    script:
        "../scripts/plot_pr.R"


rule gprofiler_bagel2:
    input:
        txt="results/bagel2/{comparison}/{comparison}.pr",
    output:
        csv="results/bagel2/gprofiler/{comparison}/{pathway_data}.csv",
        pdf=report(
            "results/plots/bagel2/gprofiler/{comparison}/{pathway_data}.pdf",
            caption="../report/pathway_analysis.rst",
            category="gprofiler plots",
            subcategory="{comparison}",
            labels={"Comparison": "{comparison}", "Figure": "pathway analysis"},
        ),
    log:
        "logs/gprofiler/bagel2/{comparison}_{pathway_data}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=10,
    params:
        fdr=config["stats"]["pathway_analysis"]["fdr"],
        top_genes=config["stats"]["pathway_analysis"]["top_genes"],
        data="bagel2",
    script:
        "../scripts/gprofiler.R"
