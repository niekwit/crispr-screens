import csv
import json
import subprocess
import tempfile
from pathlib import Path


"""Run DotMatch 0.2.2 and adapt its per-target table to this workflow's count input."""

COUNT_COLUMNS = (
    "total_count",  # documented long schema
    "count_total",  # unprefixed fallback
)


def total_count_from_row(row):
    """Return the unique-assignment total from a DotMatch count row.

    DotMatch 0.2.2 ``count`` writes ``{sample}_count_total``. The documented
    long schema names the same field ``total_count``. Accept both and fail
    with the observed headers if neither is present.
    """
    for name in COUNT_COLUMNS:
        if name in row:
            return row[name]
    matching = [name for name in row if name.endswith("_count_total")]
    if len(matching) == 1:
        return row[matching[0]]
    if len(matching) > 1:
        raise ValueError(
            "DotMatch count output has multiple *_count_total columns: "
            + ", ".join(matching)
        )
    raise ValueError(
        "DotMatch count output has no total_count or *_count_total column. "
        f"Columns: {list(row)}"
    )


if "snakemake" in globals():
    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="crispr-screens-dotmatch-") as temporary:
        raw_counts = Path(temporary) / "counts.tsv"
        command = [
            "dotmatch",
            "count",
            "--targets",
            str(snakemake.input.targets),
            "--reads",
            str(snakemake.input.fq),
            "--sample-label",
            str(snakemake.wildcards.sample),
            "--target-start",
            str(snakemake.params.target_start),
            "--target-length",
            str(snakemake.params.target_length),
            "--k",
            str(snakemake.params.k),
            "--metric",
            "hamming",
            "--ambiguity-policy",
            str(snakemake.params.ambiguity_policy),
            "--out",
            str(raw_counts),
            "--summary",
            str(snakemake.output.summary),
        ]
        completed = subprocess.run(command, capture_output=True, text=True)

        with log_path.open("w", encoding="utf-8") as log:
            if completed.stdout:
                log.write(completed.stdout)
                if not completed.stdout.endswith("\n"):
                    log.write("\n")
            if completed.stderr:
                log.write(completed.stderr)
                if not completed.stderr.endswith("\n"):
                    log.write("\n")
            if completed.returncode != 0:
                raise subprocess.CalledProcessError(
                    completed.returncode,
                    command,
                    output=completed.stdout,
                    stderr=completed.stderr,
                )

            summary = json.loads(
                Path(snakemake.output.summary).read_text(encoding="utf-8")
            )
            samples = summary.get("samples", [])
            if len(samples) > 1:
                raise ValueError(
                    "DotMatch count summary must contain exactly one sample"
                )
            sample_summary = samples[0] if samples else summary
            total_reads = int(sample_summary["total_reads"])
            assigned_unique = int(sample_summary["assigned_unique"])
            assignment_rate = 100 * assigned_unique / total_reads if total_reads else 0
            log.write(f"DotMatch unique assignment rate: {assignment_rate:.2f}%\n")
            log.write(f"DotMatch total reads: {total_reads}\n")
            log.write(f"DotMatch unique assignments: {assigned_unique}\n")
            log.write(f"DotMatch ambiguous reads: {int(sample_summary['ambiguous'])}\n")
            log.write(f"DotMatch unmatched reads: {int(sample_summary['unmatched'])}\n")
            log.write(f"DotMatch invalid reads: {int(sample_summary['invalid'])}\n")

        with raw_counts.open(encoding="utf-8", newline="") as source, Path(
            snakemake.output.counts
        ).open("w", encoding="utf-8") as target:
            for row in csv.DictReader(source, delimiter="\t"):
                target.write(f"{total_count_from_row(row)} {row['target_id']}\n")
