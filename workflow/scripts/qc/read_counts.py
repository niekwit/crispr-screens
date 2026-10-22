"""Total (post-trimming) read count per sample, parsed from the same Bowtie
logs used for the alignment rate — independent of alignment success, this
flags sequencing-depth imbalances between samples."""

import os
import re

import pandas as pd

from .common import bar_fig, clean_sample_name

_READS_RE = re.compile(r"reads processed:\s*(\d+)")


def _parse_log(path):
    with open(path) as fh:
        for line in fh:
            match = _READS_RE.search(line)
            if match:
                return int(match.group(1))
    raise ValueError(f"Could not find read count in {path}")


def build(log_files):
    rows = [
        {
            "sample": clean_sample_name(os.path.basename(f).replace(".log", "")),
            "total_reads": _parse_log(f),
        }
        for f in log_files
    ]
    df = pd.DataFrame(rows)

    fig = bar_fig(df, "sample", "total_reads", "Total reads", "Total reads")

    return dict(
        id="total-reads",
        title="Total reads",
        description=(
            "Total number of (post-trimming) reads per sample, independent of "
            "how many aligned. Large imbalances between samples can affect "
            "downstream comparisons even when alignment rate and coverage look "
            "fine individually."
        ),
        fig=fig,
        df=df,
        filename="total-reads",
    )
