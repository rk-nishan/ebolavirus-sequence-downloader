#!/usr/bin/env python3
"""
Filter and rank vaccine constructs from all_60_constructs_complete_results.csv.

Default screen:
  - VaxiJen antigenic at or above the row's default threshold
  - AllerTOP non-allergen
  - ToxinPred2 non-toxic
  - ProtParam instability class stable
  - ProtParam GRAVY negative
  - ProtParam half-life: 30 h mammalian reticulocytes, >20 h yeast, >10 h E. coli
  - Soluble by either SoluProt or SOLpro
  - Aliphatic index >= 70
  - Molecular weight between 10 and 110 kDa

The first-stage retained constructs are ranked by highest VaxiJen score first.
The second-stage selection always reads ranked_filtered_constructs.csv and
keeps constructs that pass every metric threshold, capped at the requested top N.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


DEFAULT_SCRIPT_DIR = Path(r"F:\3.Project_Drive\11. Vaccine\000. Version 1.0\12. Vaccine_construct_designs\all_constructs_merged_results")
SCRIPT_DIR = DEFAULT_SCRIPT_DIR if DEFAULT_SCRIPT_DIR.exists() else Path("./construct_results")
DEFAULT_INPUT = SCRIPT_DIR / "all_60_constructs_complete_results.csv"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "filtered_constructs"
DEFAULT_RANKED_INPUT = DEFAULT_OUTPUT_DIR / "ranked_filtered_constructs.csv"

PASS_FIELDS = [
    "pass_antigenic",
    "pass_non_allergen",
    "pass_non_toxin",
    "pass_stable",
    "pass_gravy_negative",
    "pass_half_life",
    "pass_soluble_soluprot_or_solpro",
    "pass_good_aliphatic_index",
    "pass_reasonable_molecular_weight",
]

METRIC_SELECTION_FIELDS = [
    "metric_selection_rank",
    "metric_selection_primary_category",
    "metric_selection_primary_category_rank",
    "metric_selection_primary_metric",
    "metric_selection_primary_metric_value",
    "metric_selection_categories",
]

CATEGORY_HIT_FIELDS = [
    "selection_category",
    "selection_category_description",
    "selection_category_rank",
    "selection_metric",
    "selection_metric_value",
    "selection_direction",
]


def norm(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip().lower())


def to_float(value: object, default: float = math.nan) -> float:
    try:
        text = str(value).strip()
        return float(text) if text else default
    except (TypeError, ValueError):
        return default


def is_finite(value: float) -> bool:
    return not math.isnan(value) and not math.isinf(value)


def best_available_score(row: Dict[str, object], columns: Sequence[str]) -> float:
    values = [to_float(row.get(column)) for column in columns]
    values = [value for value in values if is_finite(value)]
    return max(values) if values else math.nan


def is_antigenic(row: Dict[str, str], antigen_threshold: Optional[float]) -> bool:
    score = to_float(row.get("vaxijen_score"))
    threshold = antigen_threshold
    if threshold is None:
        threshold = to_float(row.get("vaxijen_default_threshold"), 0.4)
    prediction = norm(row.get("vaxijen_prediction"))
    prediction_supports = "antigen" in prediction and "non-antigen" not in prediction
    return score >= threshold and prediction_supports


def is_non_allergen(row: Dict[str, str]) -> bool:
    return re.search(r"\bnon[-\s]?allergen\b", norm(row.get("allertop_prediction"))) is not None


def is_non_toxin(row: Dict[str, str]) -> bool:
    return re.search(r"\bnon[-\s]?toxin\b", norm(row.get("toxinpred2_local_prediction"))) is not None


def has_target_half_life(row: Dict[str, str]) -> bool:
    half_life = norm(row.get("protparam_estimated_half_life"))
    has_mammalian = "30 hours" in half_life and "mammalian reticulocytes" in half_life
    has_yeast = ">20 hours" in half_life and "yeast" in half_life
    has_ecoli = ">10 hours" in half_life and (
        "escherichia coli" in half_life or "e. coli" in half_life
    )
    return has_mammalian and has_yeast and has_ecoli


def solubility_basis(row: Dict[str, str]) -> Tuple[bool, str]:
    source_columns = [
        "soluprot_local_prediction",
        "soluprot_prediction",
        "solpro_local_prediction",
    ]
    favorable = []
    observed = []
    for column in source_columns:
        value = norm(row.get(column))
        if value:
            observed.append(f"{column}={row.get(column)}")
        if value == "soluble":
            favorable.append(f"{column}=soluble")
    return bool(favorable), "; ".join(favorable or observed)


def metric_selection_specs():
    return [
        {
            "key": "vaxijen_ge_0_80",
            "description": "VaxiJen score >= 0.80",
            "metric": "vaxijen_score",
            "direction": "higher",
            "reverse": True,
            "qualifies": lambda row: to_float(row.get("vaxijen_score")) >= 0.80,
            "value": lambda row: to_float(row.get("vaxijen_score")),
        },
        {
            "key": "soluprot_local_prediction_soluble",
            "description": "SoluProt local prediction is soluble",
            "metric": "soluprot_local_prediction",
            "direction": "higher",
            "reverse": True,
            "qualifies": lambda row: norm(row.get("soluprot_local_prediction")) == "soluble",
            "value": lambda row: to_float(row.get("soluprot_local_score"), 0.0),
        },
        {
            "key": "soluprot_local_score_ge_0_70",
            "description": "SoluProt local score >= 0.70",
            "metric": "soluprot_local_score",
            "direction": "higher",
            "reverse": True,
            "qualifies": lambda row: is_finite(to_float(row.get("soluprot_local_score")))
            and to_float(row.get("soluprot_local_score")) >= 0.70,
            "value": lambda row: to_float(row.get("soluprot_local_score")),
        },
        {
            "key": "protparam_instability_index_le_40_0",
            "description": "ProtParam instability index <= 40.0",
            "metric": "protparam_instability_index",
            "direction": "lower",
            "reverse": False,
            "qualifies": lambda row: is_finite(to_float(row.get("protparam_instability_index")))
            and to_float(row.get("protparam_instability_index")) <= 40.0,
            "value": lambda row: to_float(row.get("protparam_instability_index")),
        },
        {
            "key": "protparam_gravy_le_neg_0_10",
            "description": "ProtParam GRAVY <= -0.10",
            "metric": "protparam_gravy",
            "direction": "lower",
            "reverse": False,
            "qualifies": lambda row: is_finite(to_float(row.get("protparam_gravy")))
            and to_float(row.get("protparam_gravy")) <= -0.10,
            "value": lambda row: to_float(row.get("protparam_gravy")),
        },
        {
            "key": "protparam_aliphatic_index_ge_70",
            "description": "ProtParam aliphatic index >= 70",
            "metric": "protparam_aliphatic_index",
            "direction": "higher",
            "reverse": True,
            "qualifies": lambda row: is_finite(to_float(row.get("protparam_aliphatic_index")))
            and to_float(row.get("protparam_aliphatic_index")) >= 70.0,
            "value": lambda row: to_float(row.get("protparam_aliphatic_index")),
        },
    ]


def select_metric_constructs(
    retained: Sequence[Dict[str, object]],
    per_category: int,
    final_top_n: int,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    del per_category
    specs = metric_selection_specs()
    selected: List[Dict[str, object]] = []
    seen: set[str] = set()
    for row in retained:
        record_id = str(row.get("record_id") or row.get("submitted_id") or "")
        if record_id in seen:
            continue
        if not all(spec["qualifies"](row) for spec in specs):
            continue
        selected_row = dict(row)
        selected_row.update(
            {
                "metric_selection_primary_category": "all_metric_thresholds",
                "metric_selection_primary_category_rank": len(selected) + 1,
                "metric_selection_primary_metric": "combined_metric_thresholds",
                "metric_selection_primary_metric_value": "pass",
                "metric_selection_categories": "; ".join(spec["key"] for spec in specs),
            }
        )
        selected.append(selected_row)
        seen.add(record_id)

    if final_top_n > 0:
        selected = selected[:final_top_n]

    for rank, row in enumerate(selected, start=1):
        row["metric_selection_rank"] = rank

    pass_rows: List[Dict[str, object]] = []
    for rank, row in enumerate(selected, start=1):
        hit_row = dict(row)
        hit_row.update(
            {
                "selection_category": "all_metric_thresholds",
                "selection_category_description": "Construct passes every metric threshold",
                "selection_category_rank": rank,
                "selection_metric": "combined_metric_thresholds",
                "selection_metric_value": "pass",
                "selection_direction": "all",
            }
        )
        pass_rows.append(hit_row)

    return selected, pass_rows


def evaluate_row(
    row: Dict[str, str],
    antigen_threshold: Optional[float],
    min_aliphatic_index: float,
    min_mw: float,
    max_mw: float,
) -> Dict[str, object]:
    gravy = to_float(row.get("protparam_gravy"))
    aliphatic_index = to_float(row.get("protparam_aliphatic_index"))
    molecular_weight = to_float(row.get("protparam_molecular_weight"))
    soluble, basis = solubility_basis(row)

    checks = {
        "pass_antigenic": is_antigenic(row, antigen_threshold),
        "pass_non_allergen": is_non_allergen(row),
        "pass_non_toxin": is_non_toxin(row),
        "pass_stable": norm(row.get("protparam_instability_class")) == "stable",
        "pass_gravy_negative": gravy < 0,
        "pass_half_life": has_target_half_life(row),
        "pass_soluble_soluprot_or_solpro": soluble,
        "pass_good_aliphatic_index": aliphatic_index >= min_aliphatic_index,
        "pass_reasonable_molecular_weight": min_mw <= molecular_weight <= max_mw,
    }
    failed = [field.replace("pass_", "") for field, passed in checks.items() if not passed]
    checks["passes_all_filters"] = not failed
    checks["failure_reasons"] = "; ".join(failed)
    checks["solubility_basis"] = basis
    return checks


def read_rows(input_csv: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    with input_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if not rows:
        raise ValueError(f"No data rows found in {input_csv}")
    return rows, fieldnames


def write_csv(path: Path, rows: Sequence[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_csv_with_fallback(
    path: Path, rows: Sequence[Dict[str, object]], fieldnames: Sequence[str]
) -> Path:
    try:
        write_csv(path, rows, fieldnames)
        return path
    except PermissionError:
        fallback = path.with_name(f"{path.stem}_updated{path.suffix}")
        write_csv(fallback, rows, fieldnames)
        return fallback


def make_summary(
    ranked_rows: Sequence[Dict[str, object]],
    metric_rows: Sequence[Dict[str, object]],
    category_hit_rows: Sequence[Dict[str, object]],
    args: argparse.Namespace,
) -> List[Dict[str, object]]:
    total = len(ranked_rows)
    summary = [
        {
            "item": "ranked_input_rows",
            "value": total,
            "notes": str(args.ranked_input),
        },
        {
            "item": "metric_threshold_pass_rows",
            "value": len(category_hit_rows),
            "notes": "Constructs passing all metric thresholds after the top-N cap",
        },
        {
            "item": "unique_metric_selected_rows",
            "value": len(metric_rows),
            "notes": (
                "Unique constructs selected by passing all metric thresholds, "
                f"capped at {args.metric_selection_top_n}"
            ),
        },
        {
            "item": "metric_source_policy",
            "value": "ranked_filtered_constructs.csv",
            "notes": "Metric selection is performed from the ranked filtered CSV, not from the raw all-constructs CSV",
        },
        {
            "item": "aliphatic_index_cutoff",
            "value": args.min_aliphatic_index,
            "notes": "Override with --min-aliphatic-index",
        },
        {
            "item": "molecular_weight_range_da",
            "value": f"{args.min_molecular_weight}-{args.max_molecular_weight}",
            "notes": "Override with --min-molecular-weight / --max-molecular-weight",
        },
    ]
    for field in PASS_FIELDS:
        if not ranked_rows or field not in ranked_rows[0]:
            continue
        passed = sum(1 for row in ranked_rows if row[field] == "TRUE")
        summary.append(
            {
                "item": field,
                "value": passed,
                "notes": f"{passed} of {total} rows passed this independent criterion",
            }
        )
    return summary


def add_rank(rows: List[Dict[str, object]]) -> None:
    for rank, row in enumerate(rows, start=1):
        row["antigenicity_rank"] = rank


def as_audit_value(value: object) -> str:
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return "" if value is None else str(value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter vaccine constructs and rank retained rows by VaxiJen score."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Input CSV path.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for ranked, audit, and summary CSV outputs.",
    )
    parser.add_argument(
        "--ranked-input",
        type=Path,
        default=DEFAULT_RANKED_INPUT,
        help="Ranked filtered constructs CSV used for metric selection.",
    )
    parser.add_argument(
        "--refresh-ranked",
        action="store_true",
        help="Rebuild ranked_filtered_constructs.csv from --input before metric selection.",
    )
    parser.add_argument(
        "--antigen-threshold",
        type=float,
        default=None,
        help="Optional fixed VaxiJen threshold. Default uses each row's vaxijen_default_threshold.",
    )
    parser.add_argument(
        "--min-aliphatic-index",
        type=float,
        default=70.0,
        help="Minimum acceptable ProtParam aliphatic index.",
    )
    parser.add_argument(
        "--min-molecular-weight",
        type=float,
        default=10000.0,
        help="Minimum acceptable molecular weight in Da.",
    )
    parser.add_argument(
        "--max-molecular-weight",
        type=float,
        default=110000.0,
        help="Maximum acceptable molecular weight in Da.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=0,
        help="If set above 0, also write only the top N retained constructs.",
    )
    parser.add_argument(
        "--metric-selection-per-category",
        type=int,
        default=2,
        help="Deprecated; metric selection now requires all threshold criteria.",
    )
    parser.add_argument(
        "--metric-selection-top-n",
        type=int,
        default=12,
        help="Final cap for the de-duplicated unique metric-selected constructs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.ranked_input == DEFAULT_RANKED_INPUT and args.output_dir != DEFAULT_OUTPUT_DIR:
        args.ranked_input = args.output_dir / "ranked_filtered_constructs.csv"

    added_fields = [
        "antigenicity_rank",
        "passes_all_filters",
        *PASS_FIELDS,
        "solubility_basis",
        "failure_reasons",
        "protparam_estimated_half_life",
    ]

    if args.refresh_ranked or not args.ranked_input.exists():
        rows, source_fields = read_rows(args.input)
        audit_rows: List[Dict[str, object]] = []
        for row in rows:
            checks = evaluate_row(
                row,
                args.antigen_threshold,
                args.min_aliphatic_index,
                args.min_molecular_weight,
                args.max_molecular_weight,
            )
            audit_row = dict(row)
            audit_row.update({key: as_audit_value(value) for key, value in checks.items()})
            audit_rows.append(audit_row)

        retained = [row for row in audit_rows if row["passes_all_filters"] == "TRUE"]
        retained.sort(
            key=lambda row: (
                -to_float(row.get("vaxijen_score")),
                str(row.get("record_id") or row.get("submitted_id") or ""),
            )
        )
        add_rank(retained)

        ranked_fields = [field for field in added_fields if field == "antigenicity_rank"] + [
            field for field in source_fields if field not in added_fields
        ] + [field for field in added_fields if field != "antigenicity_rank"]

        write_csv(args.ranked_input, retained, ranked_fields)
        audit_path = args.output_dir / "construct_filter_audit.csv"
        write_csv(audit_path, audit_rows, ranked_fields)
        print(f"Refreshed ranked output: {args.ranked_input}")
        print(f"Audit output: {audit_path}")

    ranked_rows, ranked_fields = read_rows(args.ranked_input)

    metric_selected, category_hits = select_metric_constructs(
        ranked_rows,
        args.metric_selection_per_category,
        args.metric_selection_top_n,
    )

    metric_output_fields = METRIC_SELECTION_FIELDS + [
        field for field in ranked_fields if field not in METRIC_SELECTION_FIELDS
    ]
    category_hit_output_fields = CATEGORY_HIT_FIELDS + [
        field for field in ranked_fields if field not in CATEGORY_HIT_FIELDS
    ]

    top_metric_path = args.output_dir / f"top_{len(category_hits)}_metric_selected_constructs.csv"
    category_hits_path = args.output_dir / "metric_selection_category_hits.csv"
    unique_metric_path = args.output_dir / (
        f"top_{args.metric_selection_top_n}_unique_metric_selected_constructs.csv"
        if args.metric_selection_top_n > 0
        else "unique_metric_selected_constructs.csv"
    )
    audit_path = args.output_dir / "construct_filter_audit.csv"
    summary_path = args.output_dir / "construct_filter_summary.csv"

    actual_top_metric_path = write_csv_with_fallback(
        top_metric_path, category_hits, category_hit_output_fields
    )
    actual_category_hits_path = write_csv_with_fallback(
        category_hits_path, category_hits, category_hit_output_fields
    )
    actual_unique_metric_path = write_csv_with_fallback(
        unique_metric_path, metric_selected, metric_output_fields
    )
    actual_summary_path = write_csv_with_fallback(
        summary_path,
        make_summary(ranked_rows, metric_selected, category_hits, args),
        ["item", "value", "notes"],
    )

    print(f"Ranked input rows: {len(ranked_rows)}")
    print(f"Metric threshold-pass rows: {len(category_hits)}")
    print(f"Unique metric-selected rows: {len(metric_selected)}")
    print(f"Ranked input: {args.ranked_input}")
    print(f"Top metric-selected output: {actual_top_metric_path}")
    print(f"Metric threshold-pass output: {actual_category_hits_path}")
    print(f"Unique metric-selected output: {actual_unique_metric_path}")
    print(f"Summary output: {actual_summary_path}")

    if ranked_rows:
        best = ranked_rows[0]
        print(
            "Top construct: "
            f"{best.get('record_id')} "
            f"(VaxiJen={best.get('vaxijen_score')}, "
            f"MW={best.get('protparam_molecular_weight')} Da, "
            f"AI={best.get('protparam_aliphatic_index')}, "
            f"GRAVY={best.get('protparam_gravy')})"
        )

    if args.top_n > 0:
        top_n_path = args.output_dir / f"top_{args.top_n}_filtered_constructs.csv"
        write_csv(top_n_path, ranked_rows[: args.top_n], ranked_fields)
        print(f"Top {args.top_n} output: {top_n_path}")


if __name__ == "__main__":
    main()
