"""Alignment rate per sample, parsed from Bowtie logs.

Also writes alignment-rates.csv, consumed externally by ngs-tracker — its
column names (sample, mapping_rate) must stay stable.
"""

import os
import re

import pandas as pd

from .common import bar_fig, clean_sample_name

MIN_ALIGNMENT_RATE = 60  # samples below this mapping rate (%) are flagged

_ALIGN_RE = re.compile(r"reads with at least one alignment.*\(([\d.]+)%\)")


def _parse_log(path):
    with open(path) as fh:
        for line in fh:
            match = _ALIGN_RE.search(line)
            if match:
                return float(match.group(1))
    raise ValueError(f"Could not find alignment rate in {path}")


def build(log_files, output_csv):
    rows = [
        {
            "sample": clean_sample_name(os.path.basename(f).replace(".log", "")),
            "mapping_rate": _parse_log(f),
        }
        for f in log_files
    ]
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)

    flagged = df["mapping_rate"] < MIN_ALIGNMENT_RATE
    fig = bar_fig(
        df,
        "sample",
        "mapping_rate",
        "Mapping rate (%)",
        "Alignment rate",
        y_range=[0, 100],
        flagged=flagged,
    )
    return dict(
        id="alignment-rate",
        title="Alignment rate",
        description=(
            "Percentage of reads with at least one alignment to the sgRNA library, "
            f"per sample. Samples below {MIN_ALIGNMENT_RATE}% are flagged."
        ),
        fig=fig,
        df=df,
        filename="alignment-rates",
    )
