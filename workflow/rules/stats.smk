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
            **mageck_input(),
        output:
            rnw="results/mageck/{mcomparison}/{cnv}/{mcomparison}_summary.Rnw",
            gs=report("results/mageck/{mcomparison}/{cnv}/{mcomparison}.gene_summary.txt", caption="../report/mageck.rst", category="MAGeCK"),
            ss="results/mageck/{mcomparison}/{cnv}/{mcomparison}.sgrna_summary.txt",
            norm="results/mageck/{mcomparison}/{cnv}/{mcomparison}.normalized.txt",
        params:
            control=mageck_control(),
            dir_name=lambda wc, output: os.path.dirname(output["rnw"]),
            extra=config["stats"]["mageck"]["extra_mageck_arguments"],
        threads: 2
        resources:
            runtime=30
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck/{mcomparison}_{cnv}.log"
        script:
            "../scripts/mageck.py"
            

    rule lfc_plots:
        input:
            "results/mageck/{mcomparison}/{cnv}/{mcomparison}.gene_summary.txt",
        output:
            pos=report("results/mageck_plots/{mcomparison}/{cnv}/{mcomparison}.lfc_pos.pdf", caption="../report/lfc_pos.rst", category="MAGeCK plots", subcategory="{mcomparison}", labels={"Comparison":"{mcomparison}","Figure": "lfc plot enriched genes"}),
            neg=report("results/mageck_plots/{mcomparison}/{cnv}/{mcomparison}.lfc_neg.pdf", caption="../report/lfc_neg.rst", category="MAGeCK plots", subcategory="{mcomparison}", labels={"Comparison":"{mcomparison}","Figure": "lfc plot depleted genes"}),
        threads: 1
        resources:
            runtime=5
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck_plots/lfc_{mcomparison}_{cnv}.log"
        script:
            "../scripts/lfc_plots.R"


    rule sg_rank_plot:
        input:
            "results/mageck/{mcomparison}/{cnv}/{mcomparison}.sgrna_summary.txt",
        output:
            report("results/mageck_plots/{mcomparison}/{cnv}/{mcomparison}.sgrank.pdf", caption="../report/sgrank.rst", category="MAGeCK plots", subcategory="{mcomparison}", labels={"Comparison":"{mcomparison}","Figure": "sgrank plot"})
        params:
            fdr=config["stats"]["mageck"]["fdr"],
        threads: 1
        resources:
            runtime=5
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck_plots/sgrank_{mcomparison}_{cnv}.log"
        script:
            "../scripts/sgrank_plot.R"

    
    rule pathway_analysis:
        input:
            txt="results/mageck/{mcomparison}/{cnv}/{mcomparison}.gene_summary.txt",
        output:
            csv=expand("results/mageck/{{mcomparison}}/{{cnv}}/pathway_analysis/{dbs}_{pathway_data}.csv", dbs=DBS, pathway_data=PATHWAY_DATA),
            plots=report(expand("results/mageck_plots/{{mcomparison}}/{{cnv}}/pathway_analysis/{dbs}_{pathway_data}.pdf", dbs=DBS, pathway_data=PATHWAY_DATA), caption="../report/pathway_analysis.rst", category="Pathway analysis MAGeCK plots", subcategory="{mcomparison}", labels={"Comparison":"{mcomparison}","Figure": "pathway analysis"})
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
            "logs/pathway_analysis_{mcomparison}_{cnv}.log"
        script:
            "../scripts/pathway_analysis.R"
    
if config["stats"]["bagel2"]["run"]:
    if B_COMPARISONS:
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
            

        rule convert_count_table:
            input:
                "results/count/counts-aggregated.tsv"
            output:
                "results/count/counts-aggregated-bagel2.tsv"
            params:
                fa=fasta,
            threads: 1
            resources:
                runtime=5
            log:
                "logs/bagel2/convert_count_table.log"
            conda:
                "../envs/stats.yaml"
            script:
                "../scripts/convert_count_table.py"


        rule bagel2fc:
            input:
                b2dir="resources/bagel2/",
                ct="results/count/counts-aggregated-bagel2.tsv",
            output:
                fc="results/bagel2/{bcomparison}/{bcomparison}.foldchange"
            threads: 2
            resources:
                runtime=15
            conda:
                "../envs/stats.yaml"
            log:
                command="logs/bagel2/fc/{bcomparison}_command.log",
                stdout="logs/bagel2/fc/{bcomparison}_stdout.log",
                stderr="logs/bagel2/fc/{bcomparison}_stderr.log"
            script:
                "../scripts/bagel2fc.py"


        rule bagel2bf:
            input:
                b2dir="resources/bagel2/",
                fc="results/bagel2/{bcomparison}/{bcomparison}.foldchange",
            output:
                bf="results/bagel2/{bcomparison}/{bcomparison}.bf"
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
                command="logs/bagel2/bf/{bcomparison}_command.log",
                stdout="logs/bagel2/bf/{bcomparison}_stdout.log",
                stderr="logs/bagel2/bf/{bcomparison}_stderr.log"
            script:
                "../scripts/bagel2bf.py"


        rule bagel2pr:
            input:
                b2dir="resources/bagel2/",
                bf="results/bagel2/{bcomparison}/{bcomparison}.bf",
            output:
                "results/bagel2/{bcomparison}/{bcomparison}.pr",
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
                command="logs/bagel2/pr/{bcomparison}_command.log",
                stdout="logs/bagel2/pr/{bcomparison}_stdout.log",
                stderr="logs/bagel2/pr/{bcomparison}_stderr.log"
            script:
                "../scripts/bagel2pr.py"


        rule plot_bf:
            input:
                "results/bagel2/{bcomparison}/{bcomparison}.bf"
            output:
                report("results/bagel2_plots/{bcomparison}/{bcomparison}.bf.pdf", caption="../report/bagel2_plots.rst", category="BAGEL2 plots", subcategory="{bcomparison}", labels={"Comparison":"{bcomparison}", "Figure":"BF plot"})
            threads: 1
            resources:
                runtime=5
            conda:
                "../envs/stats.yaml"
            log:
                "logs/bagel2/plot/bf_{bcomparison}.log"
            script:
                "../scripts/plot_bf.R"


        rule plot_pr:
            input:
                "results/bagel2/{bcomparison}/{bcomparison}.pr"
            output:
                report("results/bagel2_plots/{bcomparison}/{bcomparison}.pr.pdf", caption="../report/bagel2_plots.rst", category="BAGEL2 plots", subcategory="{bcomparison}", labels={"Comparison":"{bcomparison}", "Figure":"Precision-recall plot"})
            threads: 1
            resources:
                runtime=5
            conda:
                "../envs/stats.yaml"
            log:
                "logs/bagel2/plot/pr_{bcomparison}.log"
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
            counts="results/count/counts-aggregated.tsv",
            drugz="resources/drugz",
        output:
            report("results/drugz/{mcomparison}.txt", caption="../report/drugz.rst", category="DrugZ", subcategory="{mcomparison}", labels={"Comparison":"{mcomparison}", "Figure":"DrugZ output"})
        params:
            test=lambda wc, output: wc.mcomparison.split("_vs_")[0].replace("-", ","),
            control=lambda wc, output: wc.mcomparison.split("_vs_")[1].replace("-", ","),
            extra=config["stats"]["drugz"]["extra"]
        threads: 2
        resources:
            runtime=15
        conda:
            "../envs/stats.yaml"
        log:
            "logs/drugz/{mcomparison}.log"
        shell:
            "python {input.drugz}/drugz.py "
            "-i {input.counts} "
            "-c {params.control} "
            "-x {params.test} "
            "{params.extra} "
            "-o {output} 2> {log} "
