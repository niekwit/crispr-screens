if skip_stats != "mageck" and skip_stats !="both":
    rule mageck:
        input: 
            "results/count/counts-aggregated.tsv"
        output:
            "results/mageck/{mcomparison}/{mcomparison}_summary.Rnw",
            report("results/mageck/{mcomparison}/{mcomparison}.gene_summary.txt", caption="../report/mageck.rst", category="MAGeCK"),
            "results/mageck/{mcomparison}/{mcomparison}.sgrna_summary.txt",
            "results/mageck/{mcomparison}/{mcomparison}.normalized.txt"
        params:
            control=mageck_control(config),
            extra=config["stats"]["extra_mageck_arguments"],
        resources:
            runtime=config["resources"]["stats"]["time"]
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck/{mcomparison}.log"
        script:
            "../scripts/mageck.sh"
            

    rule lfc_plots:
        input:
            "results/mageck/{mcomparison}/{mcomparison}.gene_summary.txt",
        output:
            pos=report("results/mageck_plots/{mcomparison}/{mcomparison}.lfc_pos.pdf", caption="../report/lfc_pos.rst", category="MAGeCK plots", subcategory="{mcomparison}", labels={"Comparison":"{mcomparison}","Figure": "lfc plot enriched genes"}),
            neg=report("results/mageck_plots/{mcomparison}/{mcomparison}.lfc_neg.pdf", caption="../report/lfc_neg.rst", category="MAGeCK plots", subcategory="{mcomparison}", labels={"Comparison":"{mcomparison}","Figure": "lfc plot depleted genes"}),
        resources:
            runtime=config["resources"]["stats"]["time"]
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck_plots/lfc_{mcomparison}.log"
        script:
            "../scripts/lfc_plots.R"


    rule sg_rank_plot:
        input:
            "results/mageck/{mcomparison}/{mcomparison}.sgrna_summary.txt",
        output:
            report("results/mageck_plots/{mcomparison}/{mcomparison}.sgrank.pdf", caption="../report/sgrank.rst", category="MAGeCK plots", subcategory="{mcomparison}", labels={"Comparison":"{mcomparison}","Figure": "sgrank plot"})
        params:
            fdr=config["stats"]["fdr"],
        resources:
            runtime=config["resources"]["stats"]["time"]
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck_plots/sgrank_{mcomparison}.log"
        script:
            "../scripts/sgrank_plot.R"

    
if skip_stats != "bagel2" and skip_stats !="both" and B_COMPARISONS != None:
    rule install_bagel2:
        output:
            directory("workflow/scripts/bagel2"),
        log:
            "logs/bagel2/install.log"
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
        resources:
            runtime=config["resources"]["stats"]["time"]
        log:
            "logs/bagel2/convert_count_table.log"
        conda:
            "../envs/stats.yaml"
        script:
            "../scripts/convert_count_table.py"


    rule bagel2fc:
        input:
            b2dir="workflow/scripts/bagel2/",
            ct="results/count/counts-aggregated-bagel2.tsv",
        output:
            fc="results/bagel2/{bcomparison}/{bcomparison}.foldchange"
        resources:
            runtime=config["resources"]["stats"]["time"]
        conda:
            "../envs/stats.yaml"
        log:
            "logs/bagel2/fc/{bcomparison}.log"
        script:
            "../scripts/bagel2fc.py"


    rule bagel2bf:
        input:
            b2dir="workflow/scripts/bagel2/",
            fc="results/bagel2/{bcomparison}/{bcomparison}.foldchange",
        output:
            bf="results/bagel2/{bcomparison}/{bcomparison}.bf"
        params:
            species=config["lib_info"]["species"],
        resources:
            runtime=config["resources"]["stats"]["time"]
        conda:
            "../envs/stats.yaml"
        log:
            "logs/bagel2/bf/{bcomparison}.log"
        script:
            "../scripts/bagel2bf.py"


    rule bagel2pr:
        input:
            b2dir="workflow/scripts/bagel2/",
            bf="results/bagel2/{bcomparison}/{bcomparison}.bf",
        output:
            "results/bagel2/{bcomparison}/{bcomparison}.pr",
        params:
            species=config["lib_info"]["species"]
        resources:
            runtime=config["resources"]["stats"]["time"]
        conda:
            "../envs/stats.yaml"
        log:
            "logs/bagel2/pr/{bcomparison}.log"
        script:
            "../scripts/bagel2pr.py"


    rule plot_bf:
        input:
            "results/bagel2/{bcomparison}/{bcomparison}.bf"
        output:
            report("results/bagel2_plots/{bcomparison}/{bcomparison}.bf.pdf", caption="../report/bagel2_plots.rst", category="BAGEL2 plots", subcategory="{bcomparison}", labels={"Comparison":"{bcomparison}", "Figure":"BF plot"})
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
        conda:
            "../envs/stats.yaml"
        log:
            "logs/bagel2/plot/pr_{bcomparison}.log"
        script:
            "../scripts/plot_pr.R"

            