import os
import glob
import sys
import pandas as pd

def fasta():
    fasta = glob.glob("resources/*.*a") # gets both .fa and .fasta files

    # Check that only one fasta file is given in resources directory
    assert len(fasta) == 1, "ERROR: There should be only one fasta file in the resources folder"
    
    # Check if fasta file is correct (same amount of lines starting with and without >)
    with open(fasta[0], "r") as f:
        lines = f.readlines()
        
        lines_seq = [x for x in lines if not x.startswith(">")]
        lines_name = [x for x in lines if x.startswith(">")]
    
    assert len(lines_seq) == len(lines_name), "Fasta file is not correct"
    
    return fasta[0]


def hisat2_index_path(fasta):
    """Generate HISAT2 index path from fasta file
    """
    lib = os.path.basename(fasta).split(".",1)[0]
    index = f"resources/index_{lib}/index_{lib}"
    
    return index


def cut_adapt_arg(config):
    """Generates Cutadapt argument for removing vector sequence
    """
    try:
        vector = config["lib_info"]["vector"]

        # check if vector is DNA sequence or empty
        if vector == "":
            try:
                sg_length = config["lib_info"]["sg_length"]
            
            except KeyError:
                print("ERROR: No sg_length in config.yml")
                sys.exit(1)
        
            assert type(sg_length) == int, "ERROR: sg_length in config.yml is not an integer"
        
            cut_arg = f"-l {str(sg_length)}"

        elif not set(vector.lower()).issubset(set("atcg")):
            print(f"ERROR: vector sequence ({vector}) in config.yml is not a DNA sequence")
            
            sys.exit(1)
    
        else:
            cut_arg = f"-a {vector}" 
        
    except KeyError:
    
        cut_arg = f"-l {sg_length}"

    return cut_arg


def sample_names():
    """Get sample names from fastq files and check for invalid characters
    """
    fastq = glob.glob("reads/*.fastq.gz")
    
    # check if fastq files are present
    assert len(fastq) != 0, "ERROR: No fastq files (.fastq.gz) found in reads directory"
    
    sample_names = [os.path.basename(x).replace(".fastq.gz","") for x in fastq]
    
    # check for invalid characters (.,:;-) in sample names
    invalid = [".",":",";","-"]
    for name in sample_names:
        assert not any(x in name for x in invalid), f"ERROR: Invalid character(s) (.,:;-) in sample {name}"
    
    return sample_names


def comparisons():
    """Load comparisons for MAGeCK and/or BAGEL2
    """
    try:
        COMPARISONS = pd.read_csv("config/stats.csv")
        M_COMPARISONS = (COMPARISONS["test"] + "_vs_" + COMPARISONS["control"]).tolist()
        M_COMPARISONS = [x.replace(";","-") for x in M_COMPARISONS] # snakemake report does not support ; in filenames

        # get comparisons for BAGEL2
        B_COMPARISONS = COMPARISONS[COMPARISONS["bagel2"] == "y"]
        
        if len(B_COMPARISONS) != 0:
            B_COMPARISONS = (B_COMPARISONS["test"] + "_vs_" + B_COMPARISONS["control"]).tolist()
            B_COMPARISONS = [x.replace(";","-") for x in B_COMPARISONS]
        
            # remove comparisons with pooled control samples (not supported by BAGEL2)
            B_COMPARISONS = [x for x in B_COMPARISONS if not "-" in x.split("_vs_")[1]]
        else:
            B_COMPARISONS = None

    except FileNotFoundError:
        print("ERROR: No stats.csv file found")
        sys.exit(1)
            
    return M_COMPARISONS, B_COMPARISONS


def mageck_control(config):
    """Load control genes for MAGeCK
    """
    if config["stats"]["mageck_control_genes"] == "all": #use all genes as controls
        control = ""

    else: # use genes from file set in config
        file = config["stats"]["mageck_control_genes"]

        # check if file exists
        assert os.path.exists(file), f"ERROR: control gene file ({file}) does not exist"
        
        control = f"--control-gene {file}" 
        
    return control

    

