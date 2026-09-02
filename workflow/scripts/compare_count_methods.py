#!/usr/bin/env python3
"""Compare HISAT2 and DotMatch counts on the same trimmed FASTQs."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def read_guide_table(path: Path, length: int) -> list[tuple[str, str, str]]:
    rows = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            seq = row["sequence"].upper()
            if len(seq) == length and set(seq) <= set("ACGT"):
                rows.append((row["sgRNA"], row["Gene"], seq))
    return rows


def write_fasta(guides: list[tuple[str, str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, _gene, seq in guides:
            handle.write(f">{name}\n{seq}\n")


def write_targets(guides: list[tuple[str, str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["target_id", "target_seq", "gene"])
        for name, gene, seq in guides:
            writer.writerow([name, seq, gene])


def parse_hisat2_counts(path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            count, name = line.strip().split(" ", 1)
            if name != "*":
                counts[name] = int(count)
    return counts


def parse_dotmatch_counts(path: Path) -> dict[str, int]:
    from dotmatch_count import total_count_from_row

    counts: dict[str, int] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts[row["target_id"]] = int(total_count_from_row(row))
    return counts


def pearson(xs: list[float], ys: list[float]) -> float:
    mean_x = statistics.fmean(xs)
    mean_y = statistics.fmean(ys)
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if den_x == 0 or den_y == 0:
        return float("nan")
    return num / (den_x * den_y)


def spearman(xs: list[float], ys: list[float]) -> float:
    def ranks(values: list[float]) -> list[float]:
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            rank = (i + j) / 2 + 1
            for k in range(i, j + 1):
                out[order[k]] = rank
            i = j + 1
        return out

    return pearson(ranks(xs), ranks(ys))


def run(args: argparse.Namespace) -> None:
    work = Path(tempfile.mkdtemp(prefix="crispr-screens-concordance-"))
    guides = read_guide_table(args.library, args.length)
    if not guides:
        raise SystemExit(f"no length-{args.length} guides in {args.library}")
    fasta = work / "guides.fa"
    targets = work / "targets.tsv"
    write_fasta(guides, fasta)
    write_targets(guides, targets)
    index_prefix = work / "index"
    subprocess.run(
        ["hisat2-build", "-q", str(fasta), str(index_prefix)],
        check=True,
    )

    names = [name for name, _gene, _seq in guides]
    rows = []
    for sample, fastq in args.sample:
        trimmed = work / f"{sample}.trimmed.fastq.gz"
        subprocess.run(
            [
                "cutadapt",
                "-q",
                "20",
                "-l",
                str(args.length),
                "-o",
                str(trimmed),
                str(fastq),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        hisat_counts_path = work / f"{sample}.hisat2.txt"
        hisat2_cmd = [
            "hisat2",
            "--no-hd",
            "-p",
            str(args.threads),
            "-t",
            "-N",
            str(args.mismatch),
        ]
        if args.strict:
            hisat2_cmd.extend(
                ["--no-spliced-alignment", "--no-softclip", "--norc"]
            )
        hisat2_cmd.extend(
            [
                "-x",
                str(index_prefix),
                "-U",
                str(trimmed),
            ]
        )
        aligned = subprocess.run(
            hisat2_cmd,
            capture_output=True,
            text=True,
            check=False,
        )
        if aligned.returncode not in (0, 1):
            raise RuntimeError(
                f"HISAT2 failed for {sample}: {aligned.stderr[-1000:]}"
            )
        kept = []
        for line in aligned.stdout.splitlines():
            if "XS:" in line:
                continue
            fields = line.split("\t")
            if len(fields) >= 3 and fields[2] != "*":
                kept.append(fields[2])
        hisat_counts: dict[str, int] = {}
        for name in kept:
            hisat_counts[name] = hisat_counts.get(name, 0) + 1
        hisat_counts_path.write_text(
            "".join(f"{count} {name}\n" for name, count in sorted(hisat_counts.items())),
            encoding="utf-8",
        )

        dotmatch_out = work / f"{sample}.dotmatch.tsv"
        summary = work / f"{sample}.summary.json"
        subprocess.run(
            [
                "dotmatch",
                "count",
                "--targets",
                str(targets),
                "--reads",
                str(trimmed),
                "--sample-label",
                sample,
                "--target-start",
                "0",
                "--target-length",
                str(args.length),
                "--k",
                str(args.mismatch),
                "--metric",
                "hamming",
                "--ambiguity-policy",
                "radius",
                "--out",
                str(dotmatch_out),
                "--summary",
                str(summary),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        dotmatch_counts = parse_dotmatch_counts(dotmatch_out)
        hisat_total = sum(hisat_counts.get(name, 0) for name in names)
        dotmatch_total = sum(dotmatch_counts.get(name, 0) for name in names)
        xs = [float(hisat_counts.get(name, 0)) for name in names]
        ys = [float(dotmatch_counts.get(name, 0)) for name in names]
        equal = sum(1 for x, y in zip(xs, ys) if x == y)
        rows.append(
            {
                "sample": sample,
                "guides": len(names),
                "hisat2_strict_assigned" if args.strict else "hisat2_assigned": hisat_total,
                "dotmatch_assigned": dotmatch_total,
                "identical_guides": equal,
                "pearson": pearson(xs, ys),
                "spearman": spearman(xs, ys),
            }
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    if not args.keep_work:
        shutil.rmtree(work)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--sample", nargs=2, action="append", metavar=("NAME", "FASTQ"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--length", type=int, default=20)
    parser.add_argument("--mismatch", type=int, default=0)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Match DotMatch's contract: --no-spliced-alignment --no-softclip --norc",
    )
    parser.add_argument("--keep-work", action="store_true")
    args = parser.parse_args()
    if not args.sample:
        raise SystemExit("at least one --sample NAME FASTQ pair is required")
    run(args)


if __name__ == "__main__":
    main()
