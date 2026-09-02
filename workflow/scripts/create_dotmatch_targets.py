import pandas as pd


"""Create the target table consumed by the optional DotMatch counter."""

csv = pd.read_csv(snakemake.input.csv, low_memory=False)
columns = snakemake.config["csv"]
selected = csv.iloc[
    :,
    [
        columns["name_column"],
        columns["sequence_column"],
        columns["gene_column"],
    ],
].copy()
selected.columns = ["target_id", "target_seq", "gene"]

for column in selected.columns:
    if selected[column].isnull().any():
        raise ValueError(f"{column} contains empty values")

selected["target_id"] = selected["target_id"].astype(str)
selected["target_seq"] = selected["target_seq"].astype(str).str.upper()
selected["gene"] = selected["gene"].astype(str)

if selected["target_id"].duplicated().any():
    raise ValueError("target_id values must be unique for DotMatch counting")
duplicate_seqs = selected.loc[selected["target_seq"].duplicated(keep=False)]
if not duplicate_seqs.empty:
    examples = ", ".join(duplicate_seqs["target_id"].head(10).tolist())
    raise ValueError(
        "target_seq values must be unique for DotMatch counting; "
        f"duplicated sequences include: {examples}"
    )
if (~selected["target_seq"].str.fullmatch(r"[ACGT]+", na=False)).any():
    raise ValueError("target_seq values must contain DNA bases A, C, G, or T")
if selected["target_seq"].str.len().nunique() != 1:
    raise ValueError("target_seq values must all have the same length for DotMatch")
target_length = int(selected["target_seq"].str.len().iloc[0])
configured_length = int(snakemake.config["dotmatch"]["target_length"])
if target_length != configured_length:
    raise ValueError(
        f"DotMatch target_length is {configured_length}, but the library sequences are {target_length} bases"
    )

selected.to_csv(snakemake.output.targets, sep="\t", index=False)
