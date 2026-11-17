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
    params:
        # Pinned, so that a workflow version always runs the same code
        commit="53388adbb4fb0931e5c9dda135502be19e4555f0",
    shell:
        "(git clone https://github.com/hart-lab/bagel.git {output} && "
        "git -C {output} checkout -q {params.commit}) 2> {log}"


rule bagel2bf:
    input:
        b2dir="resources/bagel2/",
        fc="results/crisprcleanr/corrected_lfc_{comparison}.foldchange",
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
        "results/bagel2/{comparison}/{comparison}.bf.pdf",
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
        "results/bagel2/{comparison}/{comparison}.pr.pdf",
    log:
        "logs/bagel2/plot/pr_{comparison}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=5,
    script:
        "../scripts/plot_pr.R"


rule string_db_bagel2:
    input:
        txt="results/bagel2/{comparison}/{comparison}.bf",
    output:
        svg="results/bagel2/stringdb/{comparison}/{pathway_data}/pathway_analysis.svg",
        csv="results/bagel2/stringdb/{comparison}/{pathway_data}/pathway_analysis.csv",
        id_mapping="results/bagel2/stringdb/{comparison}/{pathway_data}/id_mapping.csv",
        # STRING version these queries used (for the methods)
        version="results/bagel2/stringdb/{comparison}/{pathway_data}/string_version.txt",
    log:
        "logs/stringdb/bagel2/{comparison}_{pathway_data}.log",
    conda:
        "../envs/stats.yaml"
    threads: 1
    resources:
        runtime=10,
    params:
        data="bagel2",
    script:
        "../scripts/string_db.py"
