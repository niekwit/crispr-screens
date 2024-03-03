import os
import glob
import sys
import pandas as pd

def targets():
    TARGETS = [
    "results/qc/multiqc.html",
    "results/qc/alignment-rates.pdf",
    "results/qc/sequence-coverage.pdf",
    "results/qc/gini-index.pdf",
    "results/qc/sample-correlation.pdf",
    "results/qc/missed-rgrnas.pdf",
]

    if skip_stats != "mageck" and skip_stats !="both":
        # extend target rule with MAGecK targets     
        TARGETS.extend([
            expand("results/mageck_plots/{mcomparison}/{mcomparison}.lfc_pos.pdf", mcomparison=M_COMPARISONS),
            expand("results/mageck_plots/{mcomparison}/{mcomparison}.lfc_neg.pdf", mcomparison=M_COMPARISONS),
            expand("results/mageck_plots/{mcomparison}/{mcomparison}.sgrank.pdf", mcomparison=M_COMPARISONS),
        ])

    if skip_stats != "bagel2" and skip_stats !="both" and B_COMPARISONS != None:
        # extend target rule with BAGEL2 targets
        TARGETS.extend([
            expand("results/bagel2_plots/{bcomparison}/{bcomparison}.bf.pdf", bcomparison=B_COMPARISONS),
            expand("results/bagel2_plots/{bcomparison}/{bcomparison}.pr.pdf", bcomparison=B_COMPARISONS),
        ])
    return TARGETS


def fasta(config):
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
        

def hisat2_index_path(fasta):
    """
    Generate HISAT2 index path from fasta file
    """
    lib = os.path.basename(fasta).split(".",1)[0]
    index = f"resources/index_{lib}/index_{lib}"
    
    return index


def cut_adapt_arg(config):
    """
    Generates Cutadapt argument for removing vector sequence
    """
    try:
        vector = config["lib_info"]["vector"]

        # Check if vector is DNA sequence or empty
        if vector == "":
            try:
                sg_length = config["lib_info"]["sg_length"]
            except KeyError:
                print("No sg_length set in config.yml")
                sys.exit(1)
        
            assert type(sg_length) == int, "sg_length in config.yml is not an integer"
        
            cut_arg = f"-l {str(sg_length)}"
        elif not set(vector.lower()).issubset(set("atcg")):
            print(f"Vector sequence ({vector}) in config.yml is not a DNA sequence")
            sys.exit(1)
        else:
            cut_arg = f"-a {vector}" 
    except KeyError:
        cut_arg = f"-l {sg_length}"

    return cut_arg


def sample_names():
    """
    Get sample names from fastq files and check for invalid characters
    """
    fastq = glob.glob("reads/*.fastq.gz")
    
    # Check if fastq files are present
    assert len(fastq) != 0, "No fastq files (.fastq.gz) found in reads directory"
    
    sample_names = [os.path.basename(x).replace(".fastq.gz","") for x in fastq]
    
    # Check for invalid characters (.,:;-) in sample names
    invalid = [".",":",";","-"]
    for name in sample_names:
        assert not any(x in name for x in invalid), f"Invalid character(s) (.,:;-) in sample {name}"
    
    return sample_names


def comparisons():
    """
    Load comparisons for MAGeCK and/or BAGEL2
    """
    try:
        COMPARISONS = pd.read_csv("config/stats.csv")
        M_COMPARISONS = (COMPARISONS["test"] + "_vs_" + COMPARISONS["control"]).tolist()
        M_COMPARISONS = [x.replace(";","-") for x in M_COMPARISONS] # snakemake report does not support ; in filenames

        # Get comparisons for BAGEL2
        B_COMPARISONS = COMPARISONS[COMPARISONS["bagel2"] == "y"]
        
        if len(B_COMPARISONS) != 0:
            B_COMPARISONS = (B_COMPARISONS["test"] + "_vs_" + B_COMPARISONS["control"]).tolist()
            B_COMPARISONS = [x.replace(";","-") for x in B_COMPARISONS]
        
            # Remove comparisons with pooled control samples (not supported by BAGEL2)
            B_COMPARISONS = [x for x in B_COMPARISONS if not "-" in x.split("_vs_")[1]]
        else:
            B_COMPARISONS = None
    except FileNotFoundError:
        print("No stats.csv file found")
        sys.exit(1)
            
    return M_COMPARISONS, B_COMPARISONS


def mageck_control(config):
    """Load control genes for MAGeCK
    """
    if config["stats"]["mageck_control_genes"] == "all": # Use all genes as controls
        control = ""
    else: # Use genes from file set in config
        file = config["stats"]["mageck_control_genes"]

        # Check if file exists
        assert os.path.exists(file), f"ERROR: control gene file ({file}) does not exist"
        control = f"--control-gene {file}" 
        
    return control

    
def gene_number(fasta):
    """
    Get number of genes from sgRNA fasta file
    """
    # Read fasta file
    f = open(fasta, "r").readlines()
    
    # Select lines with sgRNA names (start with >)
    f = [x.strip() for x in f if x.startswith(">")]
    
    # Return number of genes
    genes = [x.split("_sg")[0].replace(">","") for x in f]
    
    return len(set(genes))


def lib_name(fasta):
    """
    Get sgRNA library name
    """
    name = os.path.basename(fasta).rsplit(".", 1)[0]
    
    return name
    