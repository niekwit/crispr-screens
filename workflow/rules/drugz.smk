rule install_drugz:
    output:
        directory("resources/drugz"),
    log:
        "logs/drugz/install.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    shell:
        "git clone https://github.com/hart-lab/drugz.git {output} 2> {log}"


rule drugz:
    input:
        unpack(drugz_input),
    output:
        drugz=report(
            "results/drugz/{comparison}.txt",
            caption="../report/drugz.rst",
            category="DrugZ",
            subcategory="{comparison}",
            labels={"Comparison": "{comparison}", "Figure": "DrugZ output"},
        ),
        fc="results/drugz/{comparison}.foldchange.txt",
    log:
        "logs/drugz/{comparison}.log",
    conda:
        "../envs/stats.yaml"
    threads: 2
    resources:
        runtime=15,
    params:
        test=lambda wc, output: wc.comparison.split("_vs_")[0].replace("-", ","),
        control=lambda wc, output: wc.comparison.split("_vs_")[1].replace("-", ","),
        extra=config["stats"]["drugz"]["extra"],
    shell:
        "python {input.drugz}/drugz.py "
        "-i {input.counts} "
        "-c {params.control} "
        "-x {params.test} "
        "-f {output.fc}"
        "{params.extra} "
        "-o {output.drugz} 2> {log} "


rule plot_drugz_results:
    input:
        txt="results/drugz/{comparison}.txt",
    output:
        pdf=report(
            "results/plots/drugz/dot_plot_{comparison}.pdf",
            caption="../report/drugz.rst",
            category="DrugZ plots",
            subcategory="{comparison}",
            labels={"Comparison": "{comparison}", "Figure": "DrugZ output"},
        ),
    log:
        "logs/drugz_plots/{comparison}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    params:
        fdr=config["stats"]["string_db"]["fdr"],
    script:
        "../scripts/plot_drugz_results.R"


rule string_db_drugz:
    input:
        txt="results/drugz/{comparison}.txt",
    output:
        svg="results/drugz/stringdb/{comparison}/{pathway_data}/pathway_analysis.svg",
        csv="results/drugz/stringdb/{comparison}/{pathway_data}/pathway_analysis.csv",
    log:
        "logs/stringdb/drugz/{comparison}_{pathway_data}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=10,
    params:
        data="drugz",
    script:
        "../scripts/string_db.py"


rule interactive_drugz_report:
    input:
        txt="results/drugz/{comparison}.txt",
        string_enriched=lambda wc: (
            f"results/drugz/stringdb/{wc.comparison}/enriched/pathway_analysis.csv"
            if config["stats"]["string_db"]["run"]
            and config["stats"]["string_db"]["data"] in ("enriched", "both")
            else []
        ),
        string_depleted=lambda wc: (
            f"results/drugz/stringdb/{wc.comparison}/depleted/pathway_analysis.csv"
            if config["stats"]["string_db"]["run"]
            and config["stats"]["string_db"]["data"] in ("depleted", "both")
            else []
        ),
    output:
        html="results/drugz/interactive/{comparison}.html",
    log:
        "logs/drugz/interactive_{comparison}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    script:
        "../scripts/interactive_drugz_report.py"
