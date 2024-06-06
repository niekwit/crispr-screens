from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Get snakemake variables
config = snakemake.config
control = snakemake.params["control"]
extra = snakemake.params["extra"]
comparison = snakemake.wildcards["mcomparison"]
count_table = snakemake.input["counts"]
dir_name = snakemake.params["dir_name"]

if config["stats"]["mageck"]["apply_CNV_correction"]:
    cnv_file = snakemake.input["cnv"]
    # Get first line of cnv _file
    with open(cnv_file, "r") as f:
        header = f.readline()

    # Check if cell line is in CNV file
    cell_line = config["stats"]["mageck"]["cell_line"]
    if not cell_line in header:
        raise ValueError(f"Cell line {cell_line} not found in CNV file...")
    
    cnv = f"--cnv-norm {cnv_file} --cell-line {cell_line} "
else:
    cnv = ""

# Run MAGeCK command
test_sample = comparison.split("_vs_")[0].replace("-", ",")
control_sample = comparison.split("_vs_")[1].replace("-", ",")

shell(
    "mageck test "
    "--normcounts-to-file "
    "-k {count_table} "
    "-t {test_sample} "
    "-c {control_sample} "
    "-n {dir_name}/{comparison} "
    "{cnv} "
    "{control} "
    "{extra} "
    "{log}"
)