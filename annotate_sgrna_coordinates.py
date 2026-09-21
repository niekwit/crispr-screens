#!/usr/bin/env python3
"""
Find the genomic coordinates of sgRNAs and write a CRISPRcleanR library file.

Some CRISPR libraries do not come with the genomic coordinates of the sgRNAs,
which CRISPRcleanR (and therefore BAGEL2 in this workflow) requires. This
script searches every sgRNA sequence (exact match, both strands) in a genome
FASTA file and writes a CSV file with the columns that CRISPRcleanR expects:

    CODE,GENES,seq,CHRM,STARTpos,ENDpos,STRAND

CODE is the sgRNA name and GENES the gene name, both exactly as in the input
file, and the rows are in the same order as the input file. Coordinates are
1-based and inclusive and are those of the sgRNA sequence itself (without PAM).
STRAND is + if the sgRNA sequence is identical to the genome (forward) strand
and - if its reverse complement is. (The libraries bundled with CRISPRcleanR
use slightly different anchors, e.g. Brunello uses the Cas9 cut site and gene
orientation; this does not matter as CRISPRcleanR only uses the relative order
and distance of sgRNAs.) sgRNAs that cannot be placed (e.g. non-targeting
controls) are kept, with empty coordinates: crisprcleaner.R treats these as
control sgRNAs. Use the primary assembly FASTA, as sgRNAs of multi-copy genes
can also match unplaced scaffolds.

Search strategy
---------------
By default the search is restricted to the loci of the annotated genes
(--scope locus). All features of type --feature (default: exon) in the GTF
file are padded with --flank bases (so that sgRNAs spanning an exon-intron
boundary are found), overlapping features are merged, and each base of the
resulting loci is searched once. For human exons +-30 bases this is 177 Mb,
i.e. 5.7% of the genome. An sgRNA is assigned, in order of preference, to:

  1. a locus of its own gene (mapped_gene_locus)
  2. a locus of any other annotated gene (mapped_other_locus), which rescues
     sgRNAs of genes that have been renamed since the library was designed
  3. only with --genome-fallback: anywhere in the genome (mapped_genome)

--scope genome searches the whole genome instead (no GTF file needed).

If an sgRNA matches perfectly at several positions of its best tier, the first
one (in FASTA order) is used and the number of matches is given as n_hits in
the mapping report, so these can be inspected or removed.

All k-mers of the searched sequence are encoded as integers (2 bits per base)
and looked up in the sorted sgRNA codes with numpy, so no Python loop runs over
bases or sgRNAs. A per-sgRNA mapping report is written next to the output.

Example
-------
    python annotate_sgrna_coordinates.py \\
        --library resources/bassik.csv \\
        --name-column sgRNA --gene-column Gene --sequence-column sequence \\
        --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \\
        --gtf Homo_sapiens.GRCh38.112.gtf.gz \\
        --output resources/bassik_crisprcleanr.csv

Using the output in crispr-screens: set lib_info: library_file to the output
file (its columns are in the order name, gene, sequence, so name_column: 0,
gene_column: 1, sequence_column: 2) and crisprcleanr: library_name to any name
that is not one of the CRISPRcleanR libraries.
"""

import argparse
import gzip
import logging
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

# 2 bits per base have to fit in a signed 64-bit integer
MAX_LENGTH = 31
MIN_LENGTH = 8

# Base to integer lookup table (A=0, C=1, G=2, T=3, everything else 4)
# With this coding the complement of base b is 3 - b
BASE_LUT = np.full(256, 4, dtype=np.uint8)
for _base, _code in zip(b"ACGT", range(4)):
    BASE_LUT[_base] = _code
    BASE_LUT[_base + 32] = _code  # lower case (soft-masked) bases

# Search tiers, lower is preferred
TIER_GENE, TIER_OTHER, TIER_GENOME = 0, 1, 2
TIER_STATUS = {
    TIER_GENE: "mapped_gene_locus",
    TIER_OTHER: "mapped_other_locus",
    TIER_GENOME: "mapped_genome",
}


def setup_logging(log_file):
    """Log to file (same format as the other scripts in this repo) and console"""
    logging.basicConfig(
        format="%(levelname)s:%(asctime)s:%(message)s",
        level=logging.DEBUG,
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
        force=True,
    )


def fail(message):
    """Log an error and raise it, as in the other scripts in this repo"""
    logging.error(message)
    raise ValueError(message)


def resolve_column(df, spec, what):
    """Column can be given as header name or as 0-based column number"""
    if spec in df.columns:
        return spec
    if spec.isdigit() and int(spec) < df.shape[1]:
        return df.columns[int(spec)]
    fail(
        f"Column for {what} ('{spec}') not found in library file. "
        f"Available columns: {', '.join(map(str, df.columns))}"
    )


def normalise_chrom(chrom):
    """Make chromosome names comparable (chr1 = 1, chrMT = chrM = MT = M)"""
    chrom = chrom.astype(str).str.replace(r"^chr", "", case=False, regex=True)
    return chrom.str.upper().replace({"MT": "M"})


def chrom_key(name):
    """Scalar version of normalise_chrom"""
    name = re.sub(r"^chr", "", name, flags=re.IGNORECASE).upper()
    return "M" if name == "MT" else name


def load_guides(args):
    """Read library file and normalise the sgRNA sequences"""
    df = pd.read_csv(args.library, sep=args.sep, dtype=str, low_memory=False)
    logging.info(f"Read {len(df)} sgRNAs from {args.library}")

    name_col = resolve_column(df, args.name_column, "sgRNA name")
    seq_col = resolve_column(df, args.sequence_column, "sgRNA sequence")
    guides = pd.DataFrame({"CODE": df[name_col], "seq": df[seq_col]})
    if args.gene_column is not None:
        gene_col = resolve_column(df, args.gene_column, "gene name")
        guides["GENES"] = df[gene_col]
    else:
        guides["GENES"] = pd.Series(np.nan, index=guides.index, dtype=object)

    if guides["CODE"].isna().any() or guides["CODE"].duplicated().any():
        # CRISPRcleanR uses the sgRNA names as row names
        fail("sgRNA names must be present and unique")

    guides["seq_norm"] = guides["seq"].str.strip().str.upper()
    guides["length"] = guides["seq_norm"].str.len()
    guides["gene_key"] = guides["GENES"].str.strip().str.upper()
    guides["valid"] = (
        guides["seq_norm"].str.fullmatch(r"[ACGT]+").fillna(False)
        & guides["length"].between(MIN_LENGTH, MAX_LENGTH)
    )
    n_invalid = int((~guides["valid"]).sum())
    if n_invalid:
        logging.warning(
            f"{n_invalid} sgRNAs have no valid sequence "
            f"(only ACGT and {MIN_LENGTH}-{MAX_LENGTH} nt are supported)"
        )
    lengths = guides.loc[guides["valid"], "length"].value_counts().sort_index()
    logging.info(
        "sgRNA lengths: " + ", ".join(f"{k} nt: {v}" for k, v in lengths.items())
    )
    return guides


def build_queries(guides, guide_mask):
    """
    Encode sgRNAs (and their reverse complements) as integers, per length.
    Returns {length: (sorted unique codes, table with code/guide_idx/strand)}
    """
    queries = {}
    subset = guides[guide_mask & guides["valid"]]
    for k, group in subset.groupby("length"):
        k = int(k)
        bases = BASE_LUT[
            np.frombuffer("".join(group["seq_norm"]).encode(), dtype=np.uint8)
        ]
        bases = bases.reshape(-1, k).astype(np.int64)
        shifts = 2 * np.arange(k - 1, -1, -1, dtype=np.int64)
        forward = (bases << shifts).sum(axis=1)
        reverse = ((3 - bases[:, ::-1]) << shifts).sum(axis=1)
        idx = group.index.to_numpy()
        table = pd.concat(
            [
                pd.DataFrame({"code": forward, "guide_idx": idx, "strand": "+"}),
                pd.DataFrame({"code": reverse, "guide_idx": idx, "strand": "-"}),
            ],
            ignore_index=True,
        )
        queries[k] = (np.unique(table["code"].to_numpy()), table)
    return queries


def merge_intervals(df, group_cols):
    """
    Merge overlapping/adjacent intervals (start, end) per group. Vectorised:
    a new interval starts when a feature begins after the running maximum end
    of the features before it in the same group.
    """
    df = df.sort_values(group_cols + ["start"], ignore_index=True)
    running_end = df.groupby(group_cols, sort=False)["end"].cummax()
    previous_end = running_end.groupby([df[c] for c in group_cols], sort=False).shift()
    interval = (previous_end.isna() | (df["start"] > previous_end + 1)).cumsum()
    return df.groupby(interval, sort=False).agg(
        **{c: (c, "first") for c in group_cols},
        start=("start", "min"),
        end=("end", "max"),
    )


def build_regions(args, gene_keys):
    """
    Read GTF and return (1) the padded and merged loci over all genes as a dict
    with the normalised chromosome name as key, (2) the loci per gene, and
    (3) the set of annotated genes. If gene_keys is not None, only these genes
    are used.
    """
    start = time.time()
    # The header lines (#) have a different number of fields than the data
    # rows, so they are skipped explicitly instead of using comment="#"
    # (which would also cut attributes that contain a #)
    opener = gzip.open if str(args.gtf).endswith(".gz") else open
    n_header = 0
    with opener(args.gtf, "rb") as handle:
        for line in handle:
            if not line.startswith(b"#"):
                break
            n_header += 1
    gtf = pd.read_csv(
        args.gtf,
        sep="\t",
        header=None,
        skiprows=n_header,
        usecols=[0, 2, 3, 4, 8],
        names=["chrom", "feature", "start", "end", "attributes"],
        dtype={"chrom": str, "feature": str, "attributes": str},
        low_memory=False,
    )
    gtf = gtf[gtf["feature"] == args.feature]
    if gtf.empty:
        fail(f"No '{args.feature}' features found in {args.gtf}")
    gtf["start"] = pd.to_numeric(gtf["start"])
    gtf["end"] = pd.to_numeric(gtf["end"])
    logging.info(f"Read {len(gtf)} '{args.feature}' features from {args.gtf}")

    # Genes can be referred to by name or by ID: make a locus for each key
    frames = []
    for attribute in args.gene_attributes.split(","):
        keys = gtf["attributes"].str.extract(
            rf'{re.escape(attribute)} "([^"]+)"', expand=False
        )
        found = keys.notna()
        if not found.any():
            logging.warning(f"Attribute '{attribute}' not found in {args.gtf}")
            continue
        frames.append(
            pd.DataFrame(
                {
                    "chrom": normalise_chrom(gtf.loc[found, "chrom"]),
                    "start": gtf.loc[found, "start"],
                    "end": gtf.loc[found, "end"],
                    "key": keys[found].str.upper(),
                }
            )
        )
    if not frames:
        fail(f"None of the gene attributes ({args.gene_attributes}) were found")
    regions = pd.concat(frames, ignore_index=True).drop_duplicates()
    annotated_keys = set(regions["key"])
    if gene_keys is not None:
        regions = regions[regions["key"].isin(gene_keys)]

    # Pad the features, then merge them per gene (to check afterwards if an
    # sgRNA was found at a locus of its own gene) and over all genes (the
    # sequence that is searched, so that every base is only searched once)
    regions["start"] = (regions["start"] - args.flank).clip(lower=1)
    regions["end"] = regions["end"] + args.flank
    gene_loci = merge_intervals(regions, ["chrom", "key"])
    loci = merge_intervals(regions[["chrom", "start", "end"]], ["chrom"])
    span = int((loci["end"] - loci["start"] + 1).sum())
    logging.info(
        f"{len(gene_loci)} loci for {gene_loci['key'].nunique()} genes, "
        f"{len(loci)} non-overlapping loci ({span / 1e6:.1f} Mb to search) "
        f"built in {time.time() - start:.1f}s"
    )
    loci = {chrom: df[["start", "end"]] for chrom, df in loci.groupby("chrom")}
    return loci, gene_loci, annotated_keys


def iter_fasta(path, keep):
    """
    Yield (name, base codes) for every FASTA record for which keep(name) is
    True. Only one contig is in memory at a time. Can be gzip compressed.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    name, lines = None, None
    with opener(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if lines:
                    yield name, BASE_LUT[np.frombuffer(b"".join(lines), np.uint8)]
                name = line[1:].split()[0].decode()
                lines = [] if keep(name) else None
            elif lines is not None:
                lines.append(line.rstrip())
    if lines:
        yield name, BASE_LUT[np.frombuffer(b"".join(lines), np.uint8)]


def kmer_codes(bases, k):
    """
    Integer code of every k-mer of bases (int64, 0-3) in O(log k) vector
    operations: codes of length 1, 2, 4, 8... are built by doubling and then
    combined following the binary representation of k
    """
    n = bases.size - k + 1
    pieces = {1: bases}
    size = 1
    while size * 2 <= k:
        previous = pieces[size]
        m = previous.size - size
        pieces[size * 2] = (previous[:m] << (2 * size)) | previous[size : size + m]
        size *= 2
    codes, offset = None, 0
    for size in sorted(pieces, reverse=True):
        if k & size:
            part = pieces[size][offset : offset + n]
            codes = part if codes is None else (codes << (2 * size)) | part
            offset += size
    return codes


def scan_sequence(seq, k, unique_codes, chunk_size):
    """
    Find all windows of k bases in seq (codes 0-3, 4 = N/separator) that equal
    one of unique_codes. Returns 0-based window positions and window codes.
    """
    positions, matched = [], []
    n_windows = seq.size - k + 1
    for start in range(0, max(n_windows, 0), chunk_size):
        stop = min(start + chunk_size, n_windows)
        part = seq[start : stop + k - 1]
        invalid = part > 3
        codes = kmer_codes(np.where(invalid, 0, part).astype(np.int64), k)
        # Sorted lookup of all windows in the sgRNA codes
        idx = np.searchsorted(unique_codes, codes)
        idx[idx == unique_codes.size] = 0
        hit = unique_codes[idx] == codes
        # Windows with N or crossing a separator are not valid
        n_invalid = np.concatenate(([0], np.cumsum(invalid, dtype=np.int32)))
        hit &= (n_invalid[k:] - n_invalid[:-k]) == 0
        found = np.flatnonzero(hit)
        positions.append(found + start)
        matched.append(codes[found])
    if not positions:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    return np.concatenate(positions), np.concatenate(matched)


def extract_loci(seq, loci):
    """
    Concatenate the sequence of all loci of a contig (vectorised gather),
    separated by an N so that no window can span two loci.
    Returns the sequence and the offset and genomic start of each locus in it.
    """
    starts = (loci["start"].to_numpy() - 1).clip(min=0)
    ends = np.minimum(loci["end"].to_numpy(), seq.size)
    keep = starts < ends
    starts, ends = starts[keep], ends[keep]
    lengths = ends - starts
    offsets = np.concatenate(([0], np.cumsum(lengths + 1)[:-1]))
    total = int(lengths.sum() + lengths.size)
    within = np.arange(total) - np.repeat(offsets, lengths + 1)
    gathered = seq[np.minimum(np.repeat(starts, lengths + 1) + within, seq.size - 1)]
    gathered[offsets + lengths] = 4
    table = pd.DataFrame({"offset": offsets, "start": starts + 1})
    return gathered, table


def search_contig(seq, table, queries, chunk_size):
    """Search all sgRNAs in a contig, return hits with genomic coordinates"""
    offsets = table["offset"].to_numpy()
    frames = []
    for k, (unique_codes, query_table) in queries.items():
        positions, codes = scan_sequence(seq, k, unique_codes, chunk_size)
        if positions.size == 0:
            continue
        locus = np.searchsorted(offsets, positions, side="right") - 1
        hits = pd.DataFrame(
            {
                "code": codes,
                "start": table["start"].to_numpy()[locus] + positions - offsets[locus],
            }
        ).merge(query_table, on="code")
        hits["end"] = hits["start"] + k - 1
        frames.append(hits.drop(columns="code"))
    return pd.concat(frames, ignore_index=True) if frames else None


def search_fasta(fasta, queries, regions, chunk_size):
    """
    Search sgRNAs in the FASTA file, either in the given loci (dict with
    chromosome as key) or in the whole genome (regions is None)
    """

    def keep(name):
        return regions is None or chrom_key(name) in regions

    frames = []
    n_searched = 0
    seen = set()
    for order, (name, seq) in enumerate(iter_fasta(fasta, keep)):
        start = time.time()
        if regions is None:
            table = pd.DataFrame({"offset": [0], "start": [1]})
        else:
            seen.add(chrom_key(name))
            seq, table = extract_loci(seq, regions[chrom_key(name)])
        hits = search_contig(seq, table, queries, chunk_size)
        n_searched += seq.size
        if hits is not None:
            hits["contig"] = name
            hits["contig_order"] = order
            frames.append(hits)
        logging.debug(
            f"{name}: searched {seq.size / 1e6:.1f} Mb, "
            f"{0 if hits is None else len(hits)} hits ({time.time() - start:.1f}s)"
        )
    logging.info(f"Searched {n_searched / 1e6:.1f} Mb of sequence")
    if regions is not None and set(regions) - seen:
        logging.warning(
            "GTF chromosomes not found in the FASTA file: "
            + ", ".join(sorted(set(regions) - seen)[:10])
        )
    return pd.concat(frames, ignore_index=True) if frames else None


def in_own_gene_locus(hits, guides, gene_loci):
    """True for hits that overlap a locus of the gene of the sgRNA (vectorised)"""
    contigs = hits["contig"].unique()
    check = pd.DataFrame(
        {
            "row": np.arange(len(hits)),
            "chrom": hits["contig"].map({c: chrom_key(c) for c in contigs}),
            "key": guides["gene_key"].to_numpy()[hits["guide_idx"].to_numpy()],
            "start": hits["start"].to_numpy(),
            "end": hits["end"].to_numpy(),
        }
    ).merge(
        gene_loci.rename(columns={"start": "locus_start", "end": "locus_end"}),
        on=["chrom", "key"],
    )
    inside = check[
        (check["start"] <= check["locus_end"]) & (check["end"] >= check["locus_start"])
    ]
    own = np.zeros(len(hits), dtype=bool)
    own[inside["row"].unique()] = True
    return own


def resolve_hits(hits, guides, gene_loci):
    """
    Pick one location per sgRNA: the best tier (own gene, other gene, genome),
    then the first position in the genome. Without gene_loci (whole genome
    search) all hits are in the genome tier. Returns a table indexed like guides.
    """
    result = pd.DataFrame(
        {"CHRM": pd.NA, "STARTpos": pd.NA, "ENDpos": pd.NA, "STRAND": pd.NA},
        index=guides.index,
    )
    result["tier"] = np.nan
    result["n_hits"] = 0
    if hits is None or hits.empty:
        return result

    if gene_loci is None:
        hits["tier"] = TIER_GENOME
    else:
        own_gene = in_own_gene_locus(hits, guides, gene_loci)
        hits["tier"] = np.where(own_gene, TIER_GENE, TIER_OTHER)
    # For palindromic sgRNAs (both strands match) prefer the + strand
    hits = hits.sort_values(["guide_idx", "tier", "strand"])
    hits = hits.drop_duplicates(["guide_idx", "contig", "start"])

    best_tier = hits.groupby("guide_idx")["tier"].transform("min")
    hits = hits[hits["tier"] == best_tier]
    n_hits = hits.groupby("guide_idx").size()
    first = hits.sort_values(["guide_idx", "contig_order", "start"]).drop_duplicates(
        "guide_idx"
    )
    first = first.set_index("guide_idx")
    result.loc[first.index, "CHRM"] = first["contig"]
    result.loc[first.index, "STARTpos"] = first["start"]
    result.loc[first.index, "ENDpos"] = first["end"]
    result.loc[first.index, "STRAND"] = first["strand"]
    result.loc[first.index, "tier"] = first["tier"]
    result.loc[n_hits.index, "n_hits"] = n_hits
    return result


def merge_results(current, new):
    """Fill in sgRNAs that were not mapped yet"""
    todo = current["tier"].isna() & new["tier"].notna()
    current.loc[todo] = new.loc[todo]
    return current


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-l", "--library", required=True, help="Input sgRNA CSV")
    parser.add_argument(
        "--name-column",
        required=True,
        help="Column with sgRNA names (header name or 0-based number)",
    )
    parser.add_argument(
        "--gene-column",
        help="Column with gene names (header name or 0-based number)",
    )
    parser.add_argument(
        "--sequence-column",
        required=True,
        help="Column with sgRNA sequences (header name or 0-based number)",
    )
    parser.add_argument("--sep", default=",", help="Input CSV separator")
    parser.add_argument(
        "-f", "--fasta", required=True, help="Genome FASTA (can be .gz)"
    )
    parser.add_argument("-g", "--gtf", help="Gene annotation GTF (can be .gz)")
    parser.add_argument("-o", "--output", required=True, help="Output CSV")
    parser.add_argument(
        "--scope",
        choices=["locus", "genome"],
        default="locus",
        help="Search only around annotated genes (default) or in the whole genome",
    )
    parser.add_argument(
        "--feature",
        default="exon",
        help="GTF feature to search around: exon (default) or gene for the "
        "complete gene body",
    )
    parser.add_argument(
        "--gene-attributes",
        default="gene_name,gene_id",
        help="GTF attributes that sgRNA gene names are matched against",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=30,
        help="Bases added on both sides of every feature (should be at least "
        "the sgRNA length)",
    )
    parser.add_argument(
        "--no-annotation-fallback",
        action="store_true",
        help="Only search the loci of the sgRNA's own gene",
    )
    parser.add_argument(
        "--genome-fallback",
        action="store_true",
        help="Search the whole genome for sgRNAs that were not found in any "
        "annotated locus (slow)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=2**24,
        help="Number of windows processed at once (memory use)",
    )
    parser.add_argument("--report", help="Per-sgRNA mapping report CSV")
    parser.add_argument("--log", help="Log file")
    args = parser.parse_args()

    output = Path(args.output)
    report = (
        Path(args.report)
        if args.report
        else output.with_name(output.stem + "_mapping_report.csv")
    )
    setup_logging(args.log if args.log else output.with_suffix(".log"))
    start = time.time()
    logging.info(f"Arguments: {vars(args)}")

    if args.scope == "locus":
        if args.gtf is None or args.gene_column is None:
            fail("--gtf and --gene-column are required for --scope locus")
    for path in (args.library, args.fasta):
        if not Path(path).exists():
            fail(f"File {path} does not exist")

    guides = load_guides(args)
    if not guides["valid"].any():
        fail(
            "No sgRNA has a valid sequence: check --sequence-column "
            f"('{args.sequence_column}')"
        )
    max_length = int(guides.loc[guides["valid"], "length"].max())
    if args.scope == "locus" and args.flank < max_length:
        logging.warning(
            f"--flank ({args.flank}) is smaller than the longest sgRNA "
            f"({max_length}): sgRNAs spanning exon boundaries can be missed"
        )

    annotated_keys, gene_loci = None, None
    if args.scope == "locus":
        gene_keys = None if not args.no_annotation_fallback else set(
            guides["gene_key"].dropna()
        )
        loci, gene_loci, annotated_keys = build_regions(args, gene_keys)
        n_genes = guides["gene_key"].dropna().nunique()
        n_found = len(set(guides["gene_key"].dropna()) & annotated_keys)
        logging.info(f"{n_found} of {n_genes} genes of the library are in the GTF")
        queries = build_queries(guides, pd.Series(True, index=guides.index))
        logging.info(f"Searching {int(guides['valid'].sum())} sgRNAs in loci")
        hits = search_fasta(args.fasta, queries, loci, args.chunk_size)
    else:
        queries = build_queries(guides, pd.Series(True, index=guides.index))
        logging.info(f"Searching {int(guides['valid'].sum())} sgRNAs in genome")
        hits = search_fasta(args.fasta, queries, None, args.chunk_size)

    result = resolve_hits(hits, guides, gene_loci)
    if args.no_annotation_fallback:
        # Only hits in the own gene's loci are allowed
        other = result["tier"] > TIER_GENE
        result.loc[other, ["CHRM", "STARTpos", "ENDpos", "STRAND", "n_hits"]] = [
            pd.NA,
            pd.NA,
            pd.NA,
            pd.NA,
            0,
        ]
        result.loc[other, "tier"] = np.nan

    if args.scope == "locus" and args.genome_fallback:
        todo = result["tier"].isna() & guides["valid"]
        if todo.any():
            logging.info(
                f"Searching the genome for {int(todo.sum())} sgRNAs that were "
                "not found in annotated loci"
            )
            queries = build_queries(guides, todo)
            genome_hits = search_fasta(args.fasta, queries, None, args.chunk_size)
            result = merge_results(result, resolve_hits(genome_hits, guides, None))

    # Reason for sgRNAs without coordinates
    status = result["tier"].map(TIER_STATUS)
    unmapped = result["tier"].isna()
    reason = pd.Series("sequence_not_found", index=guides.index)
    if annotated_keys is not None:
        # Gene names are only used in the locus search
        reason[~guides["gene_key"].isin(annotated_keys)] = "gene_not_in_annotation"
        reason[guides["gene_key"].isna()] = "no_gene"
    reason[~guides["valid"]] = "invalid_sequence"
    status[unmapped] = reason[unmapped]

    # Write CRISPRcleanR library (all sgRNAs, same order as the input file)
    library = pd.DataFrame(
        {
            "CODE": guides["CODE"],
            "GENES": guides["GENES"],
            "seq": guides["seq"],
            "CHRM": result["CHRM"],
            "STARTpos": result["STARTpos"].astype("Int64"),
            "ENDpos": result["ENDpos"].astype("Int64"),
            "STRAND": result["STRAND"],
        }
    )
    library.to_csv(output, index=False)
    pd.DataFrame(
        {
            "CODE": guides["CODE"],
            "GENES": guides["GENES"],
            "status": status,
            "n_hits": result["n_hits"],
        }
    ).to_csv(report, index=False)

    # Summary
    counts = status.value_counts()
    for name, count in counts.items():
        logging.info(f"{name}: {count} sgRNAs ({100 * count / len(guides):.1f}%)")
    n_multiple = int((result["n_hits"] > 1).sum())
    if n_multiple:
        logging.warning(
            f"{n_multiple} sgRNAs have multiple perfect matches in their best "
            "tier: the first position was used (see n_hits in the report)"
        )
    lost = guides.loc[unmapped & guides["gene_key"].notna() & guides["valid"]]
    if len(lost):
        logging.warning(
            f"{len(lost)} sgRNAs with a gene name could not be placed and will "
            "be treated as control sgRNAs by crisprcleaner.R, e.g.: "
            + ", ".join(lost["CODE"].head(5))
        )
    other = int((result["tier"] == TIER_OTHER).sum())
    if other:
        logging.info(
            f"{other} sgRNAs were found at a locus of a different gene than "
            "listed in the library (renamed genes or paralogs)"
        )
    logging.info(f"Wrote {output} and {report} in {time.time() - start:.1f}s")


if __name__ == "__main__":
    main()
