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
    params:
        # Pinned, so that a workflow version always runs the same code
        commit="eb15d34e4dd172965e618d5bb662c053066da799",
    shell:
        "(git clone https://github.com/hart-lab/drugz.git {output} && "
        "git -C {output} checkout -q {params.commit}) 2> {log}"


rule drugz:
    input:
        unpack(drugz_input),
    output:
        drugz="results/drugz/{comparison}.txt",
        fc="results/drugz/{comparison}.foldchange.txt",
    log:
        "logs/drugz/{comparison}.log",
    # Not the files in results/drugz/stringdb/
    wildcard_constraints:
        comparison="[^/]+",
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
        "-f {output.fc} "
        "{params.extra} "
        "-o {output.drugz} 2> {log} "


rule plot_drugz_results:
    input:
        txt="results/drugz/{comparison}.txt",
    output:
        pdf="results/drugz/dot_plot_{comparison}.pdf",
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
        id_mapping="results/drugz/stringdb/{comparison}/{pathway_data}/id_mapping.csv",
        # STRING version these queries used (for the methods)
        version="results/drugz/stringdb/{comparison}/{pathway_data}/string_version.txt",
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
