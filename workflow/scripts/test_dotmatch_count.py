#!/usr/bin/env python3
"""Unit checks for the DotMatch count adapter and schema gate."""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import yaml
from jsonschema import Draft7Validator, ValidationError

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from dotmatch_count import total_count_from_row  # noqa: E402


def test_total_count_columns() -> None:
    assert total_count_from_row({"total_count": "7", "target_id": "g"}) == "7"
    assert total_count_from_row({"count_total": "4", "target_id": "g"}) == "4"
    assert (
        total_count_from_row({"HT_1_count_total": "11", "target_id": "g"}) == "11"
    )
    try:
        total_count_from_row({"target_id": "g", "exact_count": "1"})
    except ValueError as exc:
        assert "total_count" in str(exc)
    else:
        raise AssertionError("expected missing-column error")


def load_schema() -> dict:
    return yaml.safe_load(
        (ROOT / "workflow" / "schemas" / "config.schema.yaml").read_text(
            encoding="utf-8"
        )
    )


def minimal_config() -> dict:
    return {
        "lib_info": {"library_file": "resources/guides.csv", "species": "human"},
        "csv": {"name_column": 0, "sequence_column": 2, "gene_column": 1},
        "mismatch": 0,
        "count_method": "hisat2",
        "stats": {
            "crisprcleanr": {"library_name": "TKOv3", "min_reads": 10},
            "bagel2": {
                "run": False,
                "custom_gene_lists": {
                    "essential_genes": "none",
                    "non_essential_genes": "none",
                },
                "extra_args": {"bf": "", "pr": ""},
            },
            "mageck": {
                "run": False,
                "command": "test",
                "extra_mageck_arguments": "",
                "mageck_control_genes": "all",
                "apply_CNV_correction": False,
                "cell_line": "K562",
            },
            "drugz": {"run": False, "extra": ""},
            "pathway_analysis": {
                "run": False,
                "data": "both",
                "fdr": 0.25,
                "top_genes": 50,
            },
            "string_db": {
                "run": False,
                "data": "both",
                "fdr": 0.25,
                "top_genes": 50,
            },
        },
    }


def test_schema_requires_dotmatch_block() -> None:
    schema = load_schema()
    validator = Draft7Validator(schema)
    config = minimal_config()
    validator.validate(config)

    broken = copy.deepcopy(config)
    broken["count_method"] = "dotmatch"
    try:
        validator.validate(broken)
    except ValidationError as exc:
        assert "dotmatch" in exc.message or "dotmatch" in str(exc)
    else:
        raise AssertionError("expected schema failure without a dotmatch block")

    complete = copy.deepcopy(broken)
    complete["dotmatch"] = {
        "target_start": 0,
        "target_length": 20,
        "ambiguity_policy": "radius",
    }
    validator.validate(complete)

    missing_field = copy.deepcopy(complete)
    del missing_field["dotmatch"]["target_length"]
    try:
        validator.validate(missing_field)
    except ValidationError:
        pass
    else:
        raise AssertionError("expected schema failure without target_length")


if __name__ == "__main__":
    test_total_count_columns()
    test_schema_requires_dotmatch_block()
    print("dotmatch adapter and schema checks passed")
