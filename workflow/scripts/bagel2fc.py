import os
import pandas as pd
from snakemake.shell import shell

b2dir = snakemake.input["b2dir"]
count_table_bagel2 = snakemake.input["ct"]
fc_table = snakemake.output["fc"]

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# create dictionary to store the sample column numbers in count table
# (samples cannot be refered to by column name)
count_table = pd.read_csv(count_table_bagel2, sep="\t")
column_names = list(count_table.columns)
column_dict = {key: i for i, key in enumerate(column_names)}
column_dict = {key: column_dict[key] - 1 for key in column_dict} #first sample column should have value 1

# get comparison and control sample
comparison = os.path.basename(fc_table.replace(".foldchange",""))
control = comparison.split("_vs_")[1]
control_column = column_dict[control]

# create output directory
shell("mkdir -p bagel2/{comparison}")

# run BAGEL2 fc
shell(
    "python {b2dir}/BAGEL.py fc "
      "-i {count_table_bagel2} "
      "-o results/bagel2/{comparison}/{comparison} "
      "-c {control_column} "
      "{log}"
      )




