import os
import glob
import datetime
import pandas as pd
from snakemake.utils import min_version, validate
from snakemake.logging import logger
from snakemake.shell import shell

def targets():
    TARGETS = [
    "results/qc/multiqc.html",
    "results/qc/alignment-rates.pdf",
    "results/qc/sequence-coverage.pdf",
    "results/qc/gini-index.pdf",
    "results/qc/missed-rgrnas.pdf",
]
    if config["stats"]["mageck"]["run"]:
        # Extend targets with MAGeCK files    
        TARGETS.extend([
            expand("results/plots/mageck/{comparison}/{cnv}/{comparison}.lfc_pos.pdf", comparison=COMPARISONS, cnv=CNV),
            expand("results/plots/mageck/{comparison}/{cnv}/{comparison}.lfc_neg.pdf", comparison=COMPARISONS, cnv=CNV),
            expand("results/plots/mageck/{comparison}/{cnv}/{comparison}.sgrank.pdf", comparison=COMPARISONS, cnv=CNV),
        ])
        if config["stats"]["pathway_analysis"]["run"]:
            TARGETS.extend([
                expand("results/mageck/gprofiler/{comparison}/{cnv}/{pathway_data}.csv",pathway_data=PATHWAY_DATA, comparison=COMPARISONS, cnv=CNV),
                expand("results/plots/mageck/gprofiler/{comparison}/{cnv}/{pathway_data}.pdf", pathway_data=PATHWAY_DATA, comparison=COMPARISONS, cnv=CNV),
            ])
    if config["stats"]["bagel2"]["run"]:
        if COMPARISONS:
            # Extend targets with BAGEL2 files 
            TARGETS.extend([
                expand("results/plots/bagel2/{comparison}/{comparison}.bf.pdf", comparison=COMPARISONS),
                expand("results/plots/bagel2/{comparison}/{comparison}.pr.pdf", comparison=COMPARISONS),
            ])
        if config["stats"]["pathway_analysis"]["run"]:
            TARGETS.extend([
                expand("results/bagel2/gprofiler/{comparison}/{pathway_data}.csv", pathway_data=["depleted"], comparison=COMPARISONS),
                expand("results/plots/bagel2/gprofiler/{comparison}/{pathway_data}.pdf", pathway_data=["depleted"], comparison=COMPARISONS),
            ])
    if config["stats"]["drugz"]["run"]:
        # Extend targets with DrugZ files 
        TARGETS.extend([
            expand("results/drugz/{comparison}.txt", comparison=COMPARISONS),
            expand("results/plots/drugz/dot_plot_{comparison}.pdf", comparison=COMPARISONS),
        ])
        if config["stats"]["pathway_analysis"]["run"]:
            TARGETS.extend([
                expand("results/drugz/gprofiler/{comparison}/{pathway_data}.csv", pathway_data=PATHWAY_DATA, comparison=COMPARISONS),
                expand("results/plots/drugz/gprofiler/{comparison}/{pathway_data}.pdf", pathway_data=PATHWAY_DATA, comparison=COMPARISONS),
            ])
    return TARGETS


def fasta():

    csv = glob.glob("resources/*.csv")
    
    if len(csv) == 0:
        raise ValueError("No csv file in resources directory")
    elif len(csv) > 1:
        raise ValueError("More than one csv file in resources directory")
    else:
        csv = csv[0]
        return csv.replace(".csv",".fasta"), csv
       

def cut_adapt_arg(config):
    """
    Generates Cutadapt argument for removing vector sequence
    """
    left_trim = config["lib_info"]["left_trim"]
    vector = config["lib_info"]["vector"]
    if vector.lower() == "n":
        sg_length = config["lib_info"]["sg_length"]
        cut_arg = f"-l {str(sg_length)}"
    else:
        cut_arg = f"-a {vector}"
    
    if left_trim != 0:
        cut_arg += f" -u {str(left_trim)}"

    return cut_arg


def sample_names():
    """
    Get sample names from fastq files and check for invalid characters
    """
    fastq = glob.glob("reads/*.fastq.gz")
    
    # Check if fastq files are present
    assert len(fastq) != 0, "No fastq files (.fastq.gz) found in reads directory"
    
    sample_names = [os.path.basename(x).replace(".fastq.gz","") for x in fastq]
    
    # Check if this matches the samples in stats.csv
    # Check this now, otherwise it will fail later with MAGeCK
    stats_csv = pd.read_csv("config/stats.csv")
    samples = stats_csv["test"].tolist()
    samples.extend(stats_csv["control"].tolist())
    
    not_found = []
    for sample in samples:
        if sample not in sample_names:
            if not ";" in sample:
                not_found.append(sample)
    if not_found:
        raise ValueError(f"Sample(s) {', '.join(not_found)} not found in reads directory")
        
    return sample_names


def comparisons():
    """
    Load comparisons for MAGeCK and BAGEL2
    """
    CSV = pd.read_csv("config/stats.csv")
    COMPARISONS = CSV[["test","control"]].agg('_vs_'.join, axis=1).tolist()
    
    COMPARISONS = [x.replace(";","-") for x in COMPARISONS] # snakemake report does not support ; in filenames

    if "--paired" in config["stats"]["mageck"]["extra_mageck_arguments"]:
        # Remove comparisons with unequal number of test and control samples
        # i.e. the number of - is not zero or an even number
        COMPARISONS = [x for x in COMPARISONS if x.count("-") % 2 == 0]
                    
    return COMPARISONS


def mageck_control():
    """
    Load control genes for MAGeCK
    """
    if config["stats"]["mageck"]["mageck_control_genes"] == "all": # Use all genes as controls
        control = ""
    else: # Use genes from file set in config
        file = config["stats"]["mageck"]["mageck_control_genes"]

        # Check if file exists
        assert os.path.exists(file), f"Control gene file ({file}) does not exist"
        control = f"--control-gene {file}" 
        
    return control


def extra_mageck_args():
    """
    Defines extra arguments for MAGeCK, including disabling normalisation 
    if CRISPRcleanR is used
    """
    # Base args
    args = config["stats"]["mageck"]["extra_mageck_arguments"]
    if config["stats"]["mageck"]["apply_crisprcleanr"]:
        args += " --norm-method none " # Disable normalisation in MAGeCK
    else:
        args += "--normcounts-to-file "
    return args


def mageck_input(wildcards):
    """
    Defines MAGeCK rule input file(s) 
    """
    input_data = {}

    if config["stats"]["mageck"]["apply_CNV_correction"]:
        input_data["cnv"] = "resources/cnv_data.txt"

    if config["stats"]["mageck"]["apply_crisprcleanr"]:
        input_data["counts"] = "results/count/crisprcleanr/corrected_counts_{wildcards.comparison}.tsv".format(wildcards=wildcards)
    else:
        input_data["counts"] = "results/count/counts-aggregated.tsv"

    return input_data


def drugz_input(wildcards):
    """
    Defines DrugZ rule input file(s) 
    
    """
    # Base input
    input_data = {"drugz": "resources/drugz"}
    
    if config["stats"]["drugz"]["apply_crisprcleanr"]:
        input_data["counts"] = "results/count/crisprcleanr/corrected_counts_{wildcards.comparison}.tsv".format(wildcards=wildcards)
    else:
        input_data["counts"] = "results/count/counts-aggregated.tsv"
    
    return input_data


def cnv():
    """
    Returns value for CNV wildcard.
    """
    if config["stats"]["mageck"]["apply_CNV_correction"]:
        return ["CNV-corrected"]
    else:
        return ["not-CNV-corrected"]


def pathway_data():
    """
    Returns value(s) for PATHWAY_DATA wildcard.
    """
    data = config["stats"]["pathway_analysis"]["data"]
    if data == "both":
        return ["enriched", "depleted"]
    else:
        return [data]