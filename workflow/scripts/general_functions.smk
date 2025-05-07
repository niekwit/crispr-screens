import os
import glob
import datetime
import pandas as pd
import numpy as np
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
        # Check csv file for empty columns
        check_csv(csv)

        # Extend targets with MAGeCK files
        if config["stats"]["mageck"]["command"] == "test":
            logger.info("Running MAGeCK test...")
            TARGETS.extend(
                [
                    expand(
                        "results/plots/mageck/{comparison}/{cnv}/{comparison}.lfc_pos.pdf",
                        comparison=COMPARISONS,
                        cnv=CNV,
                    ),
                    expand(
                        "results/plots/mageck/{comparison}/{cnv}/{comparison}.lfc_neg.pdf",
                        comparison=COMPARISONS,
                        cnv=CNV,
                    ),
                    expand(
                        "results/plots/mageck/{comparison}/{cnv}/{comparison}.sgrank.pdf",
                        comparison=COMPARISONS,
                        cnv=CNV,
                    ),
                ]
            )
            if config["stats"]["pathway_analysis"]["run"]:
                TARGETS.extend(
                    [
                        expand(
                            "results/mageck/gprofiler/{comparison}/{cnv}/{pathway_data}.csv",
                            pathway_data=PATHWAY_DATA,
                            comparison=COMPARISONS,
                            cnv=CNV,
                        ),
                        expand(
                            "results/plots/mageck/gprofiler/{comparison}/{cnv}/{pathway_data}.pdf",
                            pathway_data=PATHWAY_DATA,
                            comparison=COMPARISONS,
                            cnv=CNV,
                        ),
                    ]
                )
        else:
            logger.info("Running MAGeCK mle...")
            TARGETS.extend(
                [
                    expand(
                        "results/mageck/mle/{cnv}/results.gene_summary.txt",
                        comparison=COMPARISONS,
                        cnv=CNV,
                    ),
                    expand(
                        "results/mageck/mle/{cnv}/results.sgrna_summary.txt",
                        comparison=COMPARISONS,
                        cnv=CNV,
                    ),
                ]
            )
    if config["stats"]["bagel2"]["run"]:
        if COMPARISONS:
            # Extend targets with BAGEL2 files
            TARGETS.extend(
                [
                    expand(
                        "results/plots/bagel2/{comparison}/{comparison}.bf.pdf",
                        comparison=COMPARISONS,
                    ),
                    expand(
                        "results/plots/bagel2/{comparison}/{comparison}.pr.pdf",
                        comparison=COMPARISONS,
                    ),
                ]
            )
        if config["stats"]["pathway_analysis"]["run"]:
            TARGETS.extend(
                [
                    expand(
                        "results/bagel2/gprofiler/{comparison}/{pathway_data}.csv",
                        pathway_data=["depleted"],
                        comparison=COMPARISONS,
                    ),
                    expand(
                        "results/plots/bagel2/gprofiler/{comparison}/{pathway_data}.pdf",
                        pathway_data=["depleted"],
                        comparison=COMPARISONS,
                    ),
                ]
            )
    if config["stats"]["drugz"]["run"]:
        # Check csv file for empty columns
        check_csv(csv)

        # Extend targets with DrugZ files
        TARGETS.extend(
            [
                expand("results/drugz/{comparison}.txt", comparison=COMPARISONS),
                expand(
                    "results/plots/drugz/dot_plot_{comparison}.pdf",
                    comparison=COMPARISONS,
                ),
            ]
        )
        if config["stats"]["pathway_analysis"]["run"]:
            TARGETS.extend(
                [
                    expand(
                        "results/drugz/gprofiler/{comparison}/{pathway_data}.csv",
                        pathway_data=PATHWAY_DATA,
                        comparison=COMPARISONS,
                    ),
                    expand(
                        "results/plots/drugz/gprofiler/{comparison}/{pathway_data}.pdf",
                        pathway_data=PATHWAY_DATA,
                        comparison=COMPARISONS,
                    ),
                ]
            )
    return TARGETS


def fasta():
    """
    Get fasta and csv file paths
    """
    csv = config["lib_info"]["library_file"]
    if not os.path.exists(csv):
        raise ValueError(f"Library file {csv} does not exist")
    return csv.replace(".csv", ".fasta"), csv


def cut_adapt_arg(config):
    """
    Generates Cutadapt argument for removing vector sequence
    """
    """
    This is the original code, but it kept just as a comment
    left_trim = config["lib_info"]["left_trim"]
    vector = config["lib_info"]["vector"]
    if vector.lower() == "n":
        sg_length = config["lib_info"]["sg_length"]
        cut_arg = f"-l {str(sg_length)}"
    else:
        cut_arg = f"-a {vector}"
    
    if left_trim != 0:
        cut_arg += f" -u {str(left_trim)}"
    """
    cut_arg = ""
    if config["lib_info"]["cutadapt"]["g"]:
        cut_arg = f"-g {config['lib_info']['cutadapt']['g']}".strip()
    if config["lib_info"]["cutadapt"]["a"]:
        cut_arg = f"{cut_arg} -a {config['lib_info']['cutadapt']['a']}".strip()
    if config["lib_info"]["cutadapt"]["u"]:
        cut_arg = f"{cut_arg} -u {config['lib_info']['cutadapt']['u']}".strip()
    if config["lib_info"]["cutadapt"]["l"]:
        cut_arg = f"{cut_arg} -l {config['lib_info']['cutadapt']['l']}".strip()
    if config["lib_info"]["cutadapt"]["extra"]:
        cut_arg = f"{cut_arg} {config['lib_info']['cutadapt']['extra']}".strip()

    return cut_arg


def sample_names():
    """
    Get sample names from fastq files and check for invalid characters
    """
    fastq = glob.glob("reads/*.fastq.gz")
    cram = glob.glob("reads/*.cram")
    if len(fastq) == 0 and len(cram) == 0:
        raise ValueError("No fastq or cram files found in reads directory")

    # Check which format is used
    if len(fastq) > 0 and len(cram) > 0:
        raise ValueError("Both fastq and cram files found in reads directory")
    elif len(fastq) > 0:
        ext = ".fastq.gz"
        files = fastq
    else:
        ext = ".cram"
        files = cram

    sample_names = [os.path.basename(x).replace(ext, "") for x in files]

    if config["stats"]["mageck"]["command"] == "test":
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
            raise ValueError(
                f"Sample(s) {', '.join(not_found)} not found in reads directory"
            )
    else:
        # Check if sample names are in design matrix
        matrix = config["stats"]["mageck"]["mle"]["design_matrix"]
        df = pd.read_csv(matrix, sep="\t")
        
        # Check if sample names are in first column
        sample_names_in_matrix = df.iloc[:, 0].tolist()
        not_found = []
        for sample in sample_names:
            if sample not in sample_names_in_matrix:
                not_found.append(sample)
        if not_found:
            raise ValueError(
                f"Sample(s) {', '.join(not_found)} not found in design matrix"
            )

    return sample_names


def cram():
    """
    Check if raw data is in CRAM format
    """
    # Check if CRAM files are present
    cram = glob.glob("reads/*.cram")

    if len(cram) > 0:
        return True
    else:
        return False


def comparisons():
    """
    Load comparisons for MAGeCK and BAGEL2
    """
    if config["stats"]["mageck"]["command"] == "test":
        # Load comparisons from stats.csv
        # Check if file exists
        assert os.path.exists("config/stats.csv"), "config/stats.csv file does not exist"
    else:
        return None
    CSV = pd.read_csv("config/stats.csv")
    COMPARISONS = CSV[["test", "control"]].agg("_vs_".join, axis=1).tolist()

    COMPARISONS = [
        x.replace(";", "-") for x in COMPARISONS
    ]  # snakemake report does not support ; in filenames

    if "--paired" in config["stats"]["mageck"]["extra_mageck_arguments"]:
        # Remove comparisons with unequal number of test and control samples
        # i.e. the number of - is not zero or an even number
        COMPARISONS = [x for x in COMPARISONS if x.count("-") % 2 == 0]

    return COMPARISONS


def mageck_control():
    """
    Load control genes for MAGeCK
    """
    if (
        config["stats"]["mageck"]["mageck_control_genes"] == "all"
    ):  # Use all genes as controls
        control = ""
    else:  # Use genes from file set in config
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
        args += " --norm-method none "  # Disable MAGeCK normalisation
    else:
        if config["stats"]["mageck"]["command"] == "test":
            args += "--normcounts-to-file "
    return args


def mageck_input(wildcards):
    """
    Defines MAGeCK rule input file(s)
    """
    input_data = {}

    if config["stats"]["mageck"]["apply_CNV_correction"]:
        input_data["cnv"] = "resources/cnv_data.txt"

    if (
        config["stats"]["mageck"]["apply_crisprcleanr"]
        and config["stats"]["mageck"]["command"] == "test"
    ):
        input_data[
            "counts"
        ] = "results/count/crisprcleanr/corrected_counts_{wildcards.comparison}.tsv".format(
            wildcards=wildcards
        )
    else:
        input_data["counts"] = "results/count/counts-aggregated.tsv"

    if config["stats"]["mageck"]["command"] == "mle":
        if config["stats"]["mageck"]["apply_crisprcleanr"]:
            logger.info("Skipping CRISPRcleanR normalisation for MAGeCK mle...")

        matrix = config["stats"]["mageck"]["mle"]["design_matrix"]
        # Check if matrix file exists
        assert os.path.exists(
            matrix
        ), f"MAGeCK mle design matrix file ({matrix}) does not exist"

        # Check if design matrix is valid
        if not design_matrix_valid():
            raise ValueError(f"Design matrix {matrix} is not valid")
        input_data["matrix"] = matrix

    return input_data


def drugz_input(wildcards):
    """
    Defines DrugZ rule input file(s)

    """
    # Base input
    input_data = {"drugz": "resources/drugz"}

    if config["stats"]["drugz"]["apply_crisprcleanr"]:
        input_data[
            "counts"
        ] = "results/count/crisprcleanr/corrected_counts_{wildcards.comparison}.tsv".format(
            wildcards=wildcards
        )
    else:
        input_data["counts"] = "results/count/counts-aggregated.tsv"

    return input_data


def cnv():
    """
    Returns value for CNV wildcard.
    """
    if config["stats"]["mageck"]["apply_CNV_correction"]:
        # First check if CNV data is present for selected cell line
        check_cnv_cell_line()
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


def check_cnv_cell_line():
    """
    Check if CNV data is present for selected cell line
    """
    # Read cnv_cell_lines.txt to list
    with open("workflow/scripts/cnv_cell_lines.txt") as f:
        cell_lines = f.read().splitlines()
    cell_line = config["stats"]["mageck"]["cell_line"]
    if cell_line not in cell_lines:
        raise ValueError(f"CNV data not found for cell line {cell_line}")


def check_csv(csv):
    """
    Check if csv file has any empty columns in columns set in config.
    This might crash MAGeCK or DrugZ.
    """
    # Read csv file
    df = pd.read_csv(csv, low_memory=False)

    # Fetch columns from config to check
    column_indices = [
        config["csv"]["name_column"],
        config["csv"]["gene_column"],
        config["csv"]["sequence_column"],
    ]

    # Check if the indices are within the DataFrame's column range
    if any(index >= df.shape[1] for index in column_indices):
        raise ValueError(
            f"One or more column indices in the config are out of bounds for the CSV file '{csv}' with {df.shape[1]} columns."
        )

    # Check if these columns (accessed by index) contain any empty values
    for index in column_indices:
        column_name = df.columns[
            index
        ]  # Get the actual column name for the error message
        if df.iloc[:, index].isnull().values.any():
            raise ValueError(
                f"Column at index {index} ('{column_name}') in {csv} contains empty values. Please check the file."
            )


def design_matrix_valid():
    """
    Checks if a design matrix is valid by
    detecting perfect multicollinearity.

    """
    df = pd.read_csv(config["stats"]["mageck"]["mle"]["design_matrix"], sep="\t")

    if df.empty:
        logger.error("The design matrix DataFrame is empty.")
        return False

    # Sample names are not needed for the analysis, so drop first column
    df = df.iloc[:, 1:]

    # Ensure all data is numeric
    try:
        numeric_matrix = df.astype(float).to_numpy()
    except ValueError as e:
        logger.error(f"Ensure that the design matrix only contains numerical values.")
        return False

    # Get the number of columns (predictors)
    num_columns = numeric_matrix.shape[1]

    if num_columns == 0:
        print("Error: The design matrix has no columns.")
        return False

    # Calculate the rank of the matrix
    matrix_rank = np.linalg.matrix_rank(numeric_matrix)

    # Check for multicollinearity
    if matrix_rank < num_columns:
        logger.error("Design Matrix appears invalid...")
        logger.error("At least one column is a linear combination of the others.")
        return False
    else:
        logger.info("Design matrix appears valid...")
        return True
