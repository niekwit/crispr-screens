rule install_drugz:
    output:
        directory("resources/drugz"),
    log:
        "logs/drugz/install.log"
    threads: 1
    resources:
        runtime=5
    conda:
        "../envs/stats.yaml"
    shell:
        "git clone https://github.com/hart-lab/drugz.git {output} 2> {log}"


rule drugz:
    input:
        unpack(drugz_input)
    output:
        report("results/drugz/{comparison}.txt", caption="../report/drugz.rst", category="DrugZ", subcategory="{comparison}", labels={"Comparison":"{comparison}", "Figure":"DrugZ output"})
    params:
        test=lambda wc, output: wc.comparison.split("_vs_")[0].replace("-", ","),
        control=lambda wc, output: wc.comparison.split("_vs_")[1].replace("-", ","),
        extra=config["stats"]["drugz"]["extra"]
    threads: 2
    resources:
        runtime=15
    conda:
        "../envs/stats.yaml"
    log:
        "logs/drugz/{comparison}.log"
    shell:
        "python {input.drugz}/drugz.py "
        "-i {input.counts} "
        "-c {params.control} "
        "-x {params.test} "
        "{params.extra} "
        "-o {output} 2> {log} "


rule plot_drugz_results:
    input:
        txt="results/drugz/{comparison}.txt"
    output:
        pdf=report("results/plots/drugz/dot_plot_{comparison}.pdf", caption="../report/drugz.rst", category="DrugZ plots", subcategory="{comparison}", labels={"Comparison":"{comparison}", "Figure":"DrugZ output"})
    params:
        fdr=config["stats"]["pathway_analysis"]["fdr"], 
    threads: 1
    resources:
        runtime=5
    conda:
        "../envs/stats.yaml"
    log:
        "logs/drugz_plots/{comparison}.log"
    script:
        "../scripts/plot_drugz_results.R"


rule gprofiler_drugz:
    input:
        txt="results/drugz/{comparison}.txt",
    output:
        csv="results/drugz/gprofiler/{comparison}/{pathway_data}.csv",
        pdf=report("results/plots/drugz/gprofiler/{comparison}/{pathway_data}.pdf", caption="../report/pathway_analysis.rst", category="gprofiler plots", subcategory="{comparison}", labels={"Comparison":"{comparison}","Figure": "pathway analysis"})
    params:
        fdr=config["stats"]["pathway_analysis"]["fdr"],
        top_genes=config["stats"]["pathway_analysis"]["top_genes"],
        data="drugz"
    threads: 1
    resources:
        runtime=10
    conda:
        "../envs/stats.yaml"
    log:
        "logs/gprofiler/drugz/{comparison}_{pathway_data}.log"
    script:
        "../scripts/gprofiler.R"

