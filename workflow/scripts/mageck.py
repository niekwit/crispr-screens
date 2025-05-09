from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Get snakemake variables
extra = snakemake.params["extra"]
count_table = snakemake.input["counts"]
dir_name = snakemake.params["dir_name"]
control = snakemake.params["control"]

if snakemake.config["stats"]["mageck"]["apply_CNV_correction"]:
    cnv_file = snakemake.input["cnv"]
    # Get first line of cnv _file
    with open(cnv_file, "r") as f:
        header = f.readline()

    # Check if cell line is in CNV file
    cell_line = snakemake.config["stats"]["mageck"]["cell_line"]
    if not cell_line in header:
        raise ValueError(f"Cell line {cell_line} not found in CNV file...")

    cnv = f"--cnv-norm {cnv_file} --cell-line {cell_line} "
else:
    cnv = ""

if snakemake.config["stats"]["mageck"]["command"] == "test":
    # Run MAGeCK test command
    comparison = snakemake.wildcards["comparison"]

    test_sample = comparison.split("_vs_")[0].replace("-", ",")
    control_sample = comparison.split("_vs_")[1].replace("-", ",")

    shell(
        "mageck test "
        "-k {count_table} "
        "-t {test_sample} "
        "-c {control_sample} "
        "-n {dir_name}/{comparison} "
        "{cnv} "
        "{control} "
        "{extra} "
        "{log}"
    )
else:
    # Run MAGeCK mle command
    matrix = snakemake.input["matrix"]
    threads = snakemake.threads
    matrix_name = snakemake.wildcards["matrix"]

    shell(
        "export OMP_NUM_THREADS=1; "  # https://sourceforge.net/p/mageck/wiki/Home/#q-and-a
        "mageck mle "
        "-k {count_table} "
        "-d {matrix} "
        "-n {dir_name}/{matrix_name} "
        "{cnv} "
        "{control} "
        "{extra} "
        "--threads {threads} "
        "{log}"
    )
