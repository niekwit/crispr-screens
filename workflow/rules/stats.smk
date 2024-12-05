if config["stats"]["mageck"]["run"]:
    if config["stats"]["mageck"]["apply_CNV_correction"]:
        logger.info("Applying CNV correction to MAGeCK analysis")
        check_sum = "33abf0c446f3e5116f35bda5483de28c3e504b1740457e1658b3fc2d22fbfa58"
        rule get_CNV_data:
            output:
                ensure("resources/cnv_data.txt.xz", sha256=check_sum)
            params:
                url="https://github.com/niekwit/crispr-screens/raw/main/.resources/CCLE_copynumber_byGene_2013-12-03.txt.xz"
            threads: 1
            resources:
                runtime=10
            log:
                "logs/get_cnv_data.log"
            conda:
                "../envs/stats.yaml"
            shell:
                "wget -O {output} {params.url}  2> {log}"

        rule unpack_CNV_data:
            input:
                "resources/cnv_data.txt.xz"
            output:
                "resources/cnv_data.txt"
            threads: 1
            resources:
                runtime=10
            log:
                "logs/unpack_cnv_data.log"
            conda:
                "../envs/stats.yaml"
            shell:
                "xz -dv {input} 2> {log}"
    
    rule mageck:
        input: 
            unpack(mageck_input),
        output:
            rnw="results/mageck/{comparison}/{cnv}/{comparison}_summary.Rnw",
            gs=report("results/mageck/{comparison}/{cnv}/{comparison}.gene_summary.txt", caption="../report/mageck.rst", category="MAGeCK"),
            ss="results/mageck/{comparison}/{cnv}/{comparison}.sgrna_summary.txt",
            norm="results/mageck/{comparison}/{cnv}/{comparison}.normalized.txt",
        params:
            control=mageck_control(),
            dir_name=lambda wc, output: os.path.dirname(output["rnw"]),
            extra=extra_mageck_args(),
        threads: 2
        resources:
            runtime=30
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck/{comparison}_{cnv}.log"
        script:
            "../scripts/mageck.py"
            

    rule lfc_plots:
        input:
            "results/mageck/{comparison}/{cnv}/{comparison}.gene_summary.txt",
        output:
            pos=report("results/plots/mageck/{comparison}/{cnv}/{comparison}.lfc_pos.pdf", caption="../report/lfc_pos.rst", category="MAGeCK plots", subcategory="{comparison}", labels={"Comparison":"{comparison}","Figure": "lfc plot enriched genes"}),
            neg=report("results/plots/mageck/{comparison}/{cnv}/{comparison}.lfc_neg.pdf", caption="../report/lfc_neg.rst", category="MAGeCK plots", subcategory="{comparison}", labels={"Comparison":"{comparison}","Figure": "lfc plot depleted genes"}),
        threads: 1
        resources:
            runtime=5
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck_plots/lfc_{comparison}_{cnv}.log"
        script:
            "../scripts/lfc_plots.R"


    rule sg_rank_plot:
        input:
            "results/mageck/{comparison}/{cnv}/{comparison}.sgrna_summary.txt",
        output:
            report("results/plots/mageck/{comparison}/{cnv}/{comparison}.sgrank.pdf", caption="../report/sgrank.rst", category="MAGeCK plots", subcategory="{comparison}", labels={"Comparison":"{comparison}","Figure": "sgrank plot"})
        params:
            fdr=config["stats"]["mageck"]["fdr"],
        threads: 1
        resources:
            runtime=5
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck_plots/sgrank_{comparison}_{cnv}.log"
        script:
            "../scripts/sgrank_plot.R"

    
    rule pathway_analysis:
        input:
            txt="results/mageck/{comparison}/{cnv}/{comparison}.gene_summary.txt",
        output:
            csv=expand("results/mageck/{{comparison}}/{{cnv}}/pathway_analysis/{dbs}_{pathway_data}.csv", dbs=DBS, pathway_data=PATHWAY_DATA),
            plots=report(expand("results/plots/mageck/{{comparison}}/{{cnv}}/pathway_analysis/{dbs}_{pathway_data}.pdf", dbs=DBS, pathway_data=PATHWAY_DATA), caption="../report/pathway_analysis.rst", category="Pathway analysis MAGeCK plots", subcategory="{comparison}", labels={"Comparison":"{comparison}","Figure": "pathway analysis"})
        params:
            dbs=DBS,
            fdr=config["stats"]["mageck"]["fdr"],
            top_genes=config["stats"]["pathway_analysis"]["top_genes"],
            data_type=PATHWAY_DATA,
            terms=config["stats"]["pathway_analysis"]["terms"],
        threads: 1
        resources:
            runtime=5
        conda:
            "../envs/stats.yaml"
        log:
            "logs/pathway_analysis_{comparison}_{cnv}.log"
        script:
            "../scripts/pathway_analysis.R"
    
if config["stats"]["bagel2"]["run"]:
    if COMPARISONS:
        rule install_bagel2:
            output:
                directory("resources/bagel2"),
            log:
                "logs/bagel2/install.log"
            threads: 1
            resources:
                runtime=5
            conda:
                "../envs/stats.yaml"
            shell:
                "git clone https://github.com/hart-lab/bagel.git {output} 2> {log}"
            

        rule bagel2bf:
            input:
                b2dir="resources/bagel2/",
                fc="results/count/crisprcleanr/corrected_lfc_{comparison}.foldchange",
            output:
                bf="results/bagel2/{comparison}/{comparison}.bf"
            params:
                species=config["lib_info"]["species"],
                ceg=config["stats"]["bagel2"]["custom_gene_lists"]["essential_genes"],
                cneg=config["stats"]["bagel2"]["custom_gene_lists"]["non_essential_genes"]
            threads: 2
            resources:
                runtime=15
            conda:
                "../envs/stats.yaml"
            log:
                command="logs/bagel2/bf/{comparison}_command.log",
                stdout="logs/bagel2/bf/{comparison}_stdout.log",
                stderr="logs/bagel2/bf/{comparison}_stderr.log"
            script:
                "../scripts/bagel2bf.py"


        rule bagel2pr:
            input:
                b2dir="resources/bagel2/",
                bf="results/bagel2/{comparison}/{comparison}.bf",
            output:
                "results/bagel2/{comparison}/{comparison}.pr",
            params:
                species=config["lib_info"]["species"],
                ceg=config["stats"]["bagel2"]["custom_gene_lists"]["essential_genes"],
                cneg=config["stats"]["bagel2"]["custom_gene_lists"]["non_essential_genes"],
            threads: 2
            resources:
                runtime=15
            conda:
                "../envs/stats.yaml"
            log:
                command="logs/bagel2/pr/{comparison}_command.log",
                stdout="logs/bagel2/pr/{comparison}_stdout.log",
                stderr="logs/bagel2/pr/{comparison}_stderr.log"
            script:
                "../scripts/bagel2pr.py"


        rule plot_bf:
            input:
                "results/bagel2/{comparison}/{comparison}.bf"
            output:
                report("results/plots/bagel2/{comparison}/{comparison}.bf.pdf", caption="../report/bagel2_plots.rst", category="BAGEL2 plots", subcategory="{comparison}", labels={"Comparison":"{comparison}", "Figure":"BF plot"})
            threads: 1
            resources:
                runtime=5
            conda:
                "../envs/stats.yaml"
            log:
                "logs/bagel2/plot/bf_{comparison}.log"
            script:
                "../scripts/plot_bf.R"


        rule plot_pr:
            input:
                "results/bagel2/{comparison}/{comparison}.pr"
            output:
                report("results/plots/bagel2/{comparison}/{comparison}.pr.pdf", caption="../report/bagel2_plots.rst", category="BAGEL2 plots", subcategory="{comparison}", labels={"Comparison":"{comparison}", "Figure":"Precision-recall plot"})
            threads: 1
            resources:
                runtime=5
            conda:
                "../envs/stats.yaml"
            log:
                "logs/bagel2/plot/pr_{comparison}.log"
            script:
                "../scripts/plot_pr.R"


if config["stats"]["drugz"]["run"]:
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
