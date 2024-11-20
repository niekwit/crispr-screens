import os
import glob
import sys
import re
import pandas as pd

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
            expand("results/mageck_plots/{mcomparison}/{cnv}/{mcomparison}.lfc_pos.pdf", mcomparison=M_COMPARISONS, cnv=CNV),
            expand("results/mageck_plots/{mcomparison}/{cnv}/{mcomparison}.lfc_neg.pdf", mcomparison=M_COMPARISONS, cnv=CNV),
            expand("results/mageck_plots/{mcomparison}/{cnv}/{mcomparison}.sgrank.pdf", mcomparison=M_COMPARISONS, cnv=CNV),
        ])
        if config["stats"]["pathway_analysis"]["run"]:
            TARGETS.extend([
                expand("results/mageck/{mcomparison}/{cnv}/pathway_analysis/{dbs}_{pathway_data}.csv", mcomparison=M_COMPARISONS, cnv=CNV, pathway_data=PATHWAY_DATA, dbs=DBS),
                expand("results/mageck_plots/{mcomparison}/{cnv}/pathway_analysis/{dbs}_{pathway_data}.pdf", mcomparison=M_COMPARISONS, cnv=CNV, dbs=DBS, pathway_data=PATHWAY_DATA),
            ])
    if config["stats"]["bagel2"]["run"]:
        if B_COMPARISONS:
            # Extend targets with BAGEL2 files 
            TARGETS.extend([
                expand("results/bagel2_plots/{bcomparison}/{bcomparison}.bf.pdf", bcomparison=B_COMPARISONS),
                expand("results/bagel2_plots/{bcomparison}/{bcomparison}.pr.pdf", bcomparison=B_COMPARISONS),
            ])
    if config["stats"]["drugz"]["run"]:
        # Extend targets with DrugZ files 
        TARGETS.extend([
            expand("results/drugz/{mcomparison}.txt", mcomparison=M_COMPARISONS),
        ])
    return TARGETS


def fasta():
    fasta = glob.glob("resources/*.*a") # gets both .fa and .fasta files

    # If no fasta file is available, convert csv file to fasta file
    if len(fasta) == 0:
        csv = glob.glob("resources/*.csv")
        
        try:
            if len(csv) == 0:
                print("ERROR: No fasta/csv file in resources directory")
                sys.exit(1)
            elif len(csv) > 1:
                print("ERROR: More than one csv file in resources directory")
                sys.exit(1)
            else:
                csv_to_fasta(csv[0], 
                            config["csv"]["name_column"], 
                            config["csv"]["sequence_column"])
                fasta = glob.glob("resources/*.*a")
        except KeyError:
            print("ERROR: No fasta file in resources directory and no csv file information specified in config.yml")
            sys.exit(1)
    
    # Check that only one fasta file is given in resources directory
    assert len(fasta) == 1, "ERROR: There should be only one fasta file in the resources folder"
    
    # Check if fasta file is correct (same amount of lines starting with and without >)
    with open(fasta[0], "r") as f:
        lines = f.readlines()
        lines_seq = [x for x in lines if not x.startswith(">")]
        lines_name = [x for x in lines if x.startswith(">")]
    assert len(lines_seq) == len(lines_name), "Fasta file is not correct"
    
    # Check if sgRNA names follow correct format (GENE_sgGENE_number)
    bad_sgrna_names = []
    for line in lines_name:
        if not re.match(r">[A-Za-z0-9\-\.]+_sg[A-Za-z0-9\-\.]+_[0-9]+", line):
            bad_sgrna_names.append(line)
    if len(bad_sgrna_names) != 0:
        bad_sgrna_names = "".join(bad_sgrna_names)
        print(f"ERROR: Incorrect sgRNA names in fasta file (format as GENE_sgGENE_number):\n{bad_sgrna_names}")
        sys.exit(1)

    return fasta[0]


def csv_to_fasta(csv, column_gene, column_seq):
    """
    Convert csv file to fasta file. Gene and sequence column numbers must be specified in config.yml.
    """
    df = pd.read_csv(csv)
    fasta = csv.replace(".csv",".fasta")
    
    # Enumerate number of sequences per gene in new column in csv
    df["seq_number"] = df.groupby(df.columns[column_gene - 1]).cumcount() + 1
    seq_number_loc = df.columns.get_loc("seq_number")
        
    # Create sgRNA names
    sgrna = []
    for row in zip(df.iloc[:, column_gene - 1], df.iloc[:, seq_number_loc]):
        sgrna.append(f">{row[0]}_sg{row[0]}_{row[1]}")
    df["sgrna"] = sgrna
    
    # Only keep sgRNA name and sequence columns
    df = df.iloc[:, [seq_number_loc + 1, column_seq - 1]]
    
    # Write fasta file
    df.to_csv(fasta, 
              sep="\n", 
              index=False, 
              header=False)
        

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
        
    return sample_names


def comparisons():
    """
    Load comparisons for MAGeCK and/or BAGEL2
    """
    COMPARISONS = pd.read_csv("config/stats.csv")
    M_COMPARISONS = COMPARISONS[["test","control"]].agg('_vs_'.join, axis=1).tolist()
    
    M_COMPARISONS = [x.replace(";","-") for x in M_COMPARISONS] # snakemake report does not support ; in filenames

    if "--paired" in config["stats"]["mageck"]["extra_mageck_arguments"]:
        # Remove comparisons with unequal number of test and control samples
        # i.e. the number of - is not zero or an even number
        M_COMPARISONS = [x for x in M_COMPARISONS if x.count("-") % 2 == 0]

    # Get comparisons for BAGEL2
    B_COMPARISONS = COMPARISONS[COMPARISONS["bagel2"] == "y"]
    
    if len(B_COMPARISONS) != 0:
        B_COMPARISONS = B_COMPARISONS[["test","control"]].agg('_vs_'.join, axis=1).tolist()
        B_COMPARISONS = [x.replace(";","-") for x in B_COMPARISONS]
    
        # Remove comparisons with pooled control samples (not supported by BAGEL2)
        B_COMPARISONS = [x for x in B_COMPARISONS if not "-" in x.split("_vs_")[1]]
    else:
        B_COMPARISONS = None
                
    return M_COMPARISONS, B_COMPARISONS


def mageck_control():
    """Load control genes for MAGeCK
    """
    if config["stats"]["mageck"]["mageck_control_genes"] == "all": # Use all genes as controls
        control = ""
    else: # Use genes from file set in config
        file = config["stats"]["mageck"]["mageck_control_genes"]

        # Check if file exists
        assert os.path.exists(file), f"Control gene file ({file}) does not exist"
        control = f"--control-gene {file}" 
        
    return control
  

def mageck_input():
    """
    Defines MAGeCK rule input file(s) 
    """
    # Base input
    input_data = {"counts": "results/count/counts-aggregated.tsv"}

    if config["stats"]["mageck"]["apply_CNV_correction"]:
        input_data["cnv"] = "resources/cnv_data.txt"

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