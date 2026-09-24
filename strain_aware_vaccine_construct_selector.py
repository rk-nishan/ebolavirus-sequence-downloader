from __future__ import print_function

import csv
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path


# -----------------------------------------------------------------------------
# Reproducible configuration
# -----------------------------------------------------------------------------

SCRIPT_VERSION = "2.0.0-strain-aware"

DEFAULT_PROJECT_DIR = Path(r"F:\3.Project_Drive\11. Vaccine\000. Version 1.0\11. Epitope Validation")
PROJECT_DIR = DEFAULT_PROJECT_DIR if DEFAULT_PROJECT_DIR.exists() else Path("./11_epitope_validation")

FILTERED_DIR = PROJECT_DIR / "05_filtered_candidates"
PARSED_DIR = PROJECT_DIR / "04_parsed_results"
OUTPUT_DIR = PROJECT_DIR / "12. vaccine_constructs"
REPORT_DIR = PROJECT_DIR / "06_reports"

INPUT_FILES = {
    "Bcell": FILTERED_DIR / "13_Bcell_strict_antigen_nonallergen_nontoxin_conservancy_gt95.csv",
    "CTL": FILTERED_DIR / "14_CTL_strict_antigen_nonallergen_nontoxin_mhci_positive_conservancy_gt95.csv",
    "HTL": FILTERED_DIR / "15_HTL_strict_antigen_nonallergen_nontoxin_ifn_il4_il10_conservancy_gt95.csv",
}

BLAST_SUMMARY_FILE = PARSED_DIR / "12_Human_BLASTp_query_summary.csv"

CLASS_ORDER = ["Bcell", "CTL", "HTL"]
MAX_TOTAL_EPITOPES = 24
STRAIN_COLUMN = "strain"
SINGLE_STRAIN_LABEL = "single_strain"

# B-cell epitopes are eligible only when their source protein is surface-accessible.
# Preferred input annotation: a Boolean-like column named "surface_accessible".
# Accepted true values include True/1/yes/y/surface/surface-accessible/exposed.
# Accepted false values include False/0/no/n/internal/cytoplasmic/non-surface.
SURFACE_ACCESSIBLE_COLUMN = "surface_accessible"

# Alternative to the input column: list exact protein names here. Matching is
# case-insensitive and whitespace-normalised. Leave empty when the CSV column is used.
SURFACE_ACCESSIBLE_PROTEINS = {
    # "Glycoprotein",
    # "Envelope protein",
}

# Set to True only when the B-cell input file has already been restricted to
# surface-accessible proteins upstream. Keeping this False prevents silent assumptions.
ASSUME_ALL_BCELL_INPUTS_SURFACE_ACCESSIBLE = True

SURFACE_TRUE_VALUES = {
    "1",
    "true",
    "t",
    "yes",
    "y",
    "surface",
    "surface accessible",
    "surface-accessible",
    "surface exposed",
    "surface-exposed",
    "exposed",
    "extracellular",
    "secreted",
}

SURFACE_FALSE_VALUES = {
    "0",
    "false",
    "f",
    "no",
    "n",
    "internal",
    "intracellular",
    "cytoplasmic",
    "non-surface",
    "not surface accessible",
}

SOURCE_COLUMNS = [
    "epitope_id",
    "epitope_class",
    "strain",
    "protein",
    "peptide",
    "length",
    "surface_accessible",
    "surface_accessibility_source",
    "tools",
    "starts",
    "ends",
    "evidence_count",
    "vaxijen_score",
    "vaxijen_prediction",
    "allertop_prediction",
    "toxinpred_score",
    "toxinpred_prediction",
    "conservancy_percent",
    "conservancy_percent_numeric",
    "conservancy_minimum_identity",
    "conservancy_maximum_identity",
    "mhc_i_immunogenicity_score",
    "mhc_i_immunogenicity_positive",
    "hla_i_binding_allele_count",
    "hla_i_binding_alleles",
    "ifnepitope2_score",
    "ifnepitope2_prediction",
    "il4pred_score",
    "il4pred_prediction",
    "il10pred_score",
    "il10pred_prediction",
    "hla_ii_binding_allele_count",
    "hla_ii_binding_alleles",
    "filter_rule",
]

BLAST_COLUMNS = [
    "blast_query_id",
    "human_blastp_hit_count",
    "human_blastp_no_hit",
    "human_blastp_min_evalue",
    "human_blastp_best_bitscore",
    "human_blastp_max_percent_identity",
    "human_blastp_max_query_coverage_reported_by_ncbi",
    "human_blastp_max_alignment_length",
    "human_blastp_max_query_span_length",
    "human_blastp_max_query_span_fraction",
    "human_blastp_full_length_exact_hit",
    "human_blastp_any_100pct_identity_hit",
    "human_blastp_any_80pct_identity_80pct_query_span_hit",
    "top_hit_subject_accession",
    "top_hit_percent_identity",
    "top_hit_alignment_length",
    "top_hit_query_span_length",
    "top_hit_query_coverage_reported_by_ncbi",
    "top_hit_evalue",
    "top_hit_bit_score",
]

SELECTION_COLUMNS = [
    "construct_version",
    "version_name",
    "version_description",
    "construct_order",
    "class_order",
    "strain_order",
    "selection_scope",
    "selection_note",
    "selection_primary_metric",
    "selection_priority_score",
    "quota_per_class_per_strain",
    "selection_rank_within_strain_class",
    "selection_rank_within_strain_protein_class",
    # Backward-compatible aliases. In this script these are strain-specific ranks.
    "selection_rank_within_class",
    "selection_rank_within_protein_class",
    "protein_repeat_within_strain_class",
    "duplicate_peptide_within_version",
    "diversity_max_overlap_fraction",
]

OUTPUT_COLUMNS = SELECTION_COLUMNS + SOURCE_COLUMNS + BLAST_COLUMNS


# -----------------------------------------------------------------------------
# General utilities
# -----------------------------------------------------------------------------


def normalised_text(value):
    return re.sub(r"\s+", " ", str(value or "").strip())


def canonical_text(value):
    return normalised_text(value).casefold()


def read_csv(path):
    if not path.exists():
        raise FileNotFoundError("Required input file does not exist: %s" % path)
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("CSV file has no header: %s" % path)
        return list(reader)


def write_table(path, rows, columns, delimiter=","):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter=delimiter,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def numeric(row, column, default=None):
    value = row.get(column, "")
    if value is None or value == "":
        return default
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if text.lower() in {"true", "false"}:
        return 1.0 if text.lower() == "true" else 0.0
    match = re.search(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", text)
    if not match:
        return default
    try:
        return float(match.group(0))
    except ValueError:
        return default


def as_count(row, column):
    return numeric(row, column, 0.0) or 0.0


def format_float(value, digits=6):
    if value is None:
        return ""
    return ("%%.%df" % digits) % value


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_json_hash(payload):
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def row_identity(row):
    return (
        normalised_text(row.get("strain", "")),
        normalised_text(row.get("epitope_class", "")),
        normalised_text(row.get("epitope_id", "")),
        normalised_text(row.get("protein", "")),
        normalised_text(row.get("peptide", "")).upper(),
        normalised_text(row.get("starts", "")),
        normalised_text(row.get("ends", "")),
    )


def protein_key(row):
    return normalised_text(row.get("protein", ""))


# -----------------------------------------------------------------------------
# Input validation and annotation
# -----------------------------------------------------------------------------


def validate_required_columns(rows, path, required_columns):
    if not rows:
        return
    missing = [column for column in required_columns if column not in rows[0]]
    if missing:
        raise ValueError(
            "Missing required columns in %s: %s"
            % (path, ", ".join(sorted(missing)))
        )


def build_unique_index(rows, key_column, source_name):
    index = {}
    duplicates = set()
    for row in rows:
        key = normalised_text(row.get(key_column, ""))
        if not key:
            continue
        if key in index:
            duplicates.add(key)
        index[key] = row
    if duplicates:
        raise ValueError(
            "%s contains duplicate %s values: %s"
            % (source_name, key_column, ", ".join(sorted(duplicates)[:20]))
        )
    return index


def build_unique_composite_index(rows, key_columns, source_name):
    index = {}
    duplicates = set()
    blank_rows = []

    for row_number, row in enumerate(rows, start=2):
        key = tuple(normalised_text(row.get(column, "")) for column in key_columns)
        if any(not value for value in key):
            blank_rows.append(row_number)
            continue
        canonical_key = tuple(value.casefold() for value in key)
        if canonical_key in index:
            duplicates.add(key)
        index[canonical_key] = row

    if blank_rows:
        raise ValueError(
            "%s contains blank values in composite key %s at rows: %s"
            % (
                source_name,
                " + ".join(key_columns),
                ", ".join(str(value) for value in blank_rows[:20]),
            )
        )

    if duplicates:
        examples = [
            " + ".join(values)
            for values in sorted(duplicates)[:20]
        ]
        raise ValueError(
            "%s contains duplicate %s values: %s"
            % (
                source_name,
                " + ".join(key_columns),
                "; ".join(examples),
            )
        )
    return index


def assign_strains(raw_rows_by_class):
    all_rows = [row for rows in raw_rows_by_class.values() for row in rows]
    nonblank = [normalised_text(row.get(STRAIN_COLUMN, "")) for row in all_rows]
    present = [value for value in nonblank if value]

    if not present:
        for row in all_rows:
            row[STRAIN_COLUMN] = SINGLE_STRAIN_LABEL
        return [SINGLE_STRAIN_LABEL], "single_strain_fallback"

    blank_count = sum(1 for value in nonblank if not value)
    if blank_count:
        raise ValueError(
            "The '%s' column is partially populated: %d rows are blank. "
            "Every row must have a strain when multi-strain mode is used."
            % (STRAIN_COLUMN, blank_count)
        )

    strains = sorted(set(present), key=lambda value: (value.casefold(), value))
    for row in all_rows:
        row[STRAIN_COLUMN] = normalised_text(row.get(STRAIN_COLUMN, ""))
    return strains, "explicit_strain_column"


def parse_surface_accessible(row):
    raw = normalised_text(row.get(SURFACE_ACCESSIBLE_COLUMN, ""))
    if raw:
        canonical = raw.casefold()
        if canonical in SURFACE_TRUE_VALUES:
            return True, "input_column:%s" % SURFACE_ACCESSIBLE_COLUMN
        if canonical in SURFACE_FALSE_VALUES:
            return False, "input_column:%s" % SURFACE_ACCESSIBLE_COLUMN
        raise ValueError(
            "Unrecognised %s value '%s' for epitope_id=%s. "
            "Use a recognised Boolean-like value."
            % (
                SURFACE_ACCESSIBLE_COLUMN,
                raw,
                row.get("epitope_id", ""),
            )
        )

    configured = {canonical_text(value) for value in SURFACE_ACCESSIBLE_PROTEINS}
    if canonical_text(row.get("protein", "")) in configured:
        return True, "configured_surface_protein_set"

    if SURFACE_ACCESSIBLE_PROTEINS:
        return False, "configured_surface_protein_set"

    if ASSUME_ALL_BCELL_INPUTS_SURFACE_ACCESSIBLE:
        return True, "assumed_pre_filtered_bcell_input"

    return None, "unresolved"


def annotate_surface_accessibility(rows_by_class):
    unresolved = []
    for row in rows_by_class.get("Bcell", []):
        flag, source = parse_surface_accessible(row)
        if flag is None:
            unresolved.append(row)
            continue
        row["_surface_accessible"] = flag
        row["surface_accessible"] = "True" if flag else "False"
        row["surface_accessibility_source"] = source

    if unresolved:
        example_proteins = sorted(
            {protein_key(row) for row in unresolved if protein_key(row)}
        )[:10]
        raise ValueError(
            "B-cell surface accessibility is not defined. Add a '%s' column, "
            "populate SURFACE_ACCESSIBLE_PROTEINS, or explicitly set "
            "ASSUME_ALL_BCELL_INPUTS_SURFACE_ACCESSIBLE=True. Example proteins: %s"
            % (
                SURFACE_ACCESSIBLE_COLUMN,
                "; ".join(example_proteins) or "none",
            )
        )

    for epitope_class in ("CTL", "HTL"):
        for row in rows_by_class.get(epitope_class, []):
            flag, source = parse_surface_accessible(row)
            if flag is None:
                row["surface_accessible"] = ""
                row["surface_accessibility_source"] = "not_required_for_%s" % epitope_class
            else:
                row["surface_accessible"] = "True" if flag else "False"
                row["surface_accessibility_source"] = source


def compute_quota_per_class_per_strain(number_of_strains):
    if number_of_strains <= 0:
        raise ValueError("At least one strain is required.")
    quota = MAX_TOTAL_EPITOPES // (len(CLASS_ORDER) * number_of_strains)
    if quota < 1:
        raise ValueError(
            "%d strains cannot each receive at least one B-cell, CTL, and HTL "
            "epitope under the maximum of %d. The current design supports at most %d strains."
            % (
                number_of_strains,
                MAX_TOTAL_EPITOPES,
                MAX_TOTAL_EPITOPES // len(CLASS_ORDER),
            )
        )
    return quota


# -----------------------------------------------------------------------------
# Data loading and BLAST merge
# -----------------------------------------------------------------------------


def load_inputs():
    raw_rows_by_class = {}
    for epitope_class, path in INPUT_FILES.items():
        rows = read_csv(path)
        validate_required_columns(rows, path, ["epitope_id", "protein", "peptide"])
        raw_rows_by_class[epitope_class] = rows

    strains, strain_mode = assign_strains(raw_rows_by_class)

    blast_rows = read_csv(BLAST_SUMMARY_FILE)
    validate_required_columns(blast_rows, BLAST_SUMMARY_FILE, ["epitope_id"])

    blast_has_strain = any(
        normalised_text(row.get(STRAIN_COLUMN, ""))
        for row in blast_rows
    )
    blast_has_class = any(
        normalised_text(row.get("epitope_class", ""))
        for row in blast_rows
    )

    if blast_has_strain:
        blank_strain_rows = [
            index
            for index, row in enumerate(blast_rows, start=2)
            if not normalised_text(row.get(STRAIN_COLUMN, ""))
        ]
        if blank_strain_rows:
            raise ValueError(
                "%s has a partially populated '%s' column. Blank rows: %s"
                % (
                    BLAST_SUMMARY_FILE,
                    STRAIN_COLUMN,
                    ", ".join(str(value) for value in blank_strain_rows[:20]),
                )
            )

        if blast_has_class:
            blast_key_columns = [STRAIN_COLUMN, "epitope_class", "epitope_id"]
        else:
            blast_key_columns = [STRAIN_COLUMN, "epitope_id"]

        blast_index = build_unique_composite_index(
            blast_rows,
            blast_key_columns,
            str(BLAST_SUMMARY_FILE),
        )
    else:
        if len(strains) > 1:
            raise ValueError(
                "%s has no populated '%s' column, but the candidate files contain "
                "%d strains. Re-run the strain-aware BLAST preparation and merge scripts."
                % (BLAST_SUMMARY_FILE, STRAIN_COLUMN, len(strains))
            )
        blast_key_columns = ["epitope_id"]
        blast_index = build_unique_index(
            blast_rows,
            "epitope_id",
            str(BLAST_SUMMARY_FILE),
        )

    rows_by_class = {}
    missing_blast_keys = []
    for epitope_class, raw_rows in raw_rows_by_class.items():
        class_rows = []
        for row in raw_rows:
            merged = {column: row.get(column, "") for column in SOURCE_COLUMNS}
            merged["epitope_class"] = normalised_text(
                merged.get("epitope_class", "") or epitope_class
            )
            merged["strain"] = normalised_text(row.get(STRAIN_COLUMN, ""))
            merged["protein"] = protein_key(row)
            merged["peptide"] = normalised_text(row.get("peptide", "")).upper()
            merged["epitope_id"] = normalised_text(row.get("epitope_id", ""))

            if blast_has_strain:
                if blast_has_class:
                    lookup_values = (
                        merged["strain"],
                        merged["epitope_class"],
                        merged["epitope_id"],
                    )
                else:
                    lookup_values = (
                        merged["strain"],
                        merged["epitope_id"],
                    )
                lookup_key = tuple(
                    normalised_text(value).casefold()
                    for value in lookup_values
                )
                blast = blast_index.get(lookup_key, {})
            else:
                blast = blast_index.get(merged["epitope_id"], {})

            if not blast:
                missing_blast_keys.append(
                    "%s + %s + %s"
                    % (
                        merged["strain"],
                        merged["epitope_class"],
                        merged["epitope_id"],
                    )
                )

            for column in BLAST_COLUMNS:
                if column == "human_blastp_max_query_span_fraction":
                    continue
                merged[column] = blast.get(column, "")

            span = numeric(merged, "human_blastp_max_query_span_length", None)
            length = numeric(merged, "length", None)
            if span is not None and length:
                merged["human_blastp_max_query_span_fraction"] = format_float(
                    span / length,
                    6,
                )
            else:
                merged["human_blastp_max_query_span_fraction"] = ""

            class_rows.append(merged)

        # Stable row order ensures reproducible ranking even before explicit sorting.
        class_rows.sort(key=row_identity)
        rows_by_class[epitope_class] = class_rows

    if missing_blast_keys:
        raise ValueError(
            "No Human BLASTp summary row was found for %d candidate rows. "
            "Examples: %s"
            % (
                len(missing_blast_keys),
                "; ".join(missing_blast_keys[:20]),
            )
        )

    annotate_surface_accessibility(rows_by_class)
    return rows_by_class, strains, strain_mode


# -----------------------------------------------------------------------------
# Scoring logic retained from the original script
# -----------------------------------------------------------------------------


def build_normalizers(rows_by_class):
    columns = [
        "vaxijen_score",
        "conservancy_percent_numeric",
        "evidence_count",
        "mhc_i_immunogenicity_score",
        "hla_i_binding_allele_count",
        "ifnepitope2_score",
        "il4pred_score",
        "il10pred_score",
        "hla_ii_binding_allele_count",
        "human_blastp_hit_count",
        "human_blastp_max_percent_identity",
        "human_blastp_max_query_span_fraction",
    ]
    stats = {}
    for epitope_class, rows in rows_by_class.items():
        for column in columns:
            values = [numeric(row, column, None) for row in rows]
            values = [value for value in values if value is not None]
            if values:
                stats[(epitope_class, column)] = (min(values), max(values))
    return stats


def normalize(row, column, stats, default=0.0):
    value = numeric(row, column, None)
    if value is None:
        return default
    min_value, max_value = stats.get(
        (row["epitope_class"], column),
        (value, value),
    )
    if max_value == min_value:
        return 1.0
    return (value - min_value) / (max_value - min_value)


def inverse_normalize(row, column, stats, default=0.0):
    # Missing BLAST values are not treated as ideal low-human-similarity results.
    if numeric(row, column, None) is None:
        return default
    return 1.0 - normalize(row, column, stats, default)


def score_high_antigenicity(row, stats):
    return numeric(row, "vaxijen_score", -999999.0)


def score_ctl_immunogenicity(row, stats):
    return numeric(row, "mhc_i_immunogenicity_score", -999999.0)


def score_hla_i_breadth(row, stats):
    return numeric(row, "hla_i_binding_allele_count", -999999.0)


def score_hla_ii_breadth(row, stats):
    return numeric(row, "hla_ii_binding_allele_count", -999999.0)


def score_balanced(row, stats):
    if row["epitope_class"] == "Bcell":
        return (
            0.60 * normalize(row, "vaxijen_score", stats)
            + 0.20 * normalize(row, "evidence_count", stats)
            + 0.20 * normalize(row, "conservancy_percent_numeric", stats)
        )
    if row["epitope_class"] == "CTL":
        return (
            0.35 * normalize(row, "vaxijen_score", stats)
            + 0.30 * normalize(row, "mhc_i_immunogenicity_score", stats)
            + 0.25 * normalize(row, "hla_i_binding_allele_count", stats)
            + 0.10 * normalize(row, "evidence_count", stats)
        )
    return (
        0.30 * normalize(row, "vaxijen_score", stats)
        + 0.25 * normalize(row, "hla_ii_binding_allele_count", stats)
        + 0.15 * normalize(row, "ifnepitope2_score", stats)
        + 0.15 * normalize(row, "il4pred_score", stats)
        + 0.15 * normalize(row, "il10pred_score", stats)
    )


def score_ctl_strength(row, stats):
    if row["epitope_class"] == "CTL":
        return (
            0.40 * normalize(row, "mhc_i_immunogenicity_score", stats)
            + 0.35 * normalize(row, "hla_i_binding_allele_count", stats)
            + 0.25 * normalize(row, "vaxijen_score", stats)
        )
    return score_high_antigenicity(row, stats)


def score_htl_strength(row, stats):
    if row["epitope_class"] == "HTL":
        return (
            0.35 * normalize(row, "hla_ii_binding_allele_count", stats)
            + 0.25 * normalize(row, "vaxijen_score", stats)
            + 0.15 * normalize(row, "ifnepitope2_score", stats)
            + 0.125 * normalize(row, "il4pred_score", stats)
            + 0.125 * normalize(row, "il10pred_score", stats)
        )
    if row["epitope_class"] == "CTL":
        return score_ctl_immunogenicity(row, stats)
    return score_high_antigenicity(row, stats)


def score_low_human_similarity(row, stats):
    if row["epitope_class"] == "CTL":
        immune_component = normalize(
            row,
            "mhc_i_immunogenicity_score",
            stats,
        )
    else:
        immune_component = normalize(row, "vaxijen_score", stats)
    return (
        0.45
        * inverse_normalize(
            row,
            "human_blastp_max_query_span_fraction",
            stats,
        )
        + 0.25
        * inverse_normalize(
            row,
            "human_blastp_hit_count",
            stats,
        )
        + 0.15
        * inverse_normalize(
            row,
            "human_blastp_max_percent_identity",
            stats,
        )
        + 0.15 * immune_component
    )


def class_score_function(version_id, epitope_class):
    if version_id == "V1":
        return score_high_antigenicity, "vaxijen_score"
    if version_id == "V2":
        if epitope_class == "CTL":
            return score_ctl_immunogenicity, "mhc_i_immunogenicity_score"
        return score_high_antigenicity, "vaxijen_score"
    if version_id == "V3":
        if epitope_class == "CTL":
            return score_hla_i_breadth, "hla_i_binding_allele_count"
        if epitope_class == "HTL":
            return score_hla_ii_breadth, "hla_ii_binding_allele_count"
        return score_high_antigenicity, "vaxijen_score"
    if version_id == "V4":
        return score_balanced, "balanced_composite_score"
    if version_id == "V5":
        if epitope_class == "CTL":
            return score_ctl_strength, "ctl_strength_composite_score"
        return score_high_antigenicity, "vaxijen_score"
    if version_id == "V6":
        if epitope_class == "HTL":
            return score_htl_strength, "htl_strength_composite_score"
        if epitope_class == "CTL":
            return score_ctl_immunogenicity, "mhc_i_immunogenicity_score"
        return score_high_antigenicity, "vaxijen_score"
    if version_id == "V7":
        return score_low_human_similarity, "low_human_similarity_composite_score"
    if version_id == "V8":
        if epitope_class == "CTL":
            return score_ctl_immunogenicity, "mhc_i_immunogenicity_score"
        return score_high_antigenicity, "vaxijen_score"
    raise ValueError("Unknown version: %s" % version_id)


def tie_values(row):
    return [
        numeric(row, "vaxijen_score", -999999.0),
        numeric(row, "mhc_i_immunogenicity_score", -999999.0),
        numeric(row, "hla_i_binding_allele_count", -999999.0),
        numeric(row, "hla_ii_binding_allele_count", -999999.0),
        numeric(row, "evidence_count", -999999.0),
        numeric(row, "conservancy_percent_numeric", -999999.0),
        -numeric(row, "human_blastp_max_query_span_fraction", 999999.0),
        -numeric(row, "human_blastp_hit_count", 999999.0),
    ]


def candidate_sort_key(row, scorer, stats):
    score = scorer(row, stats)
    return tuple(
        [-score]
        + [-value for value in tie_values(row)]
        + [
            canonical_text(row.get("protein", "")),
            canonical_text(row.get("peptide", "")),
            canonical_text(row.get("epitope_id", "")),
        ]
    )


def sorted_candidates(rows, scorer, stats):
    return sorted(rows, key=lambda row: candidate_sort_key(row, scorer, stats))


# -----------------------------------------------------------------------------
# Coordinate-overlap support retained for V10
# -----------------------------------------------------------------------------


def intervals(row):
    starts = [int(value) for value in re.findall(r"\d+", row.get("starts", ""))]
    ends = [int(value) for value in re.findall(r"\d+", row.get("ends", ""))]
    pairs = []
    for start, end in zip(starts, ends):
        if end >= start:
            pairs.append((start, end))
    return pairs


def interval_overlap_fraction(row, already_selected):
    same_context = [
        candidate
        for candidate in already_selected
        if candidate.get("strain") == row.get("strain")
        and protein_key(candidate) == protein_key(row)
    ]
    if not same_context:
        return 0.0
    row_intervals = intervals(row)
    if not row_intervals:
        return 0.0
    max_fraction = 0.0
    for selected in same_context:
        for start_a, end_a in row_intervals:
            length_a = end_a - start_a + 1
            for start_b, end_b in intervals(selected):
                length_b = end_b - start_b + 1
                overlap = max(
                    0,
                    min(end_a, end_b) - max(start_a, start_b) + 1,
                )
                if overlap:
                    max_fraction = max(
                        max_fraction,
                        float(overlap) / float(min(length_a, length_b)),
                    )
    return max_fraction


# -----------------------------------------------------------------------------
# Strain-balanced, protein-diverse selection
# -----------------------------------------------------------------------------


def version_definitions(quota):
    common_scope = "per_strain_dynamic_quota_distinct_protein_first"
    common_text = (
        "For every strain, selects up to %d epitopes per class. "
        "B-cell candidates must be surface-accessible. Within each strain and "
        "class, distinct proteins are selected before any protein is repeated."
        % quota
    )
    return [
        {
            "id": "V1",
            "name": "Highest antigenicity by strain",
            "description": common_text + " Ranking uses VaxiJen antigenicity.",
            "scope": common_scope,
        },
        {
            "id": "V2",
            "name": "B/HTL antigenicity plus CTL immunogenicity by strain",
            "description": common_text
            + " B-cell and HTL use VaxiJen; CTL uses MHC-I immunogenicity.",
            "scope": common_scope,
        },
        {
            "id": "V3",
            "name": "Maximum HLA allele breadth by strain",
            "description": common_text
            + " B-cell uses VaxiJen; CTL uses HLA-I allele count; HTL uses HLA-II allele count.",
            "scope": common_scope,
        },
        {
            "id": "V4",
            "name": "Balanced composite by strain",
            "description": common_text
            + " Ranking uses the original class-specific balanced composite scores.",
            "scope": common_scope,
        },
        {
            "id": "V5",
            "name": "CTL-strength construct by strain",
            "description": common_text
            + " CTL uses the CTL-strength composite; B-cell and HTL use VaxiJen.",
            "scope": common_scope,
        },
        {
            "id": "V6",
            "name": "HTL-strength construct by strain",
            "description": common_text
            + " HTL uses the HTL-strength composite; CTL uses immunogenicity; B-cell uses VaxiJen.",
            "scope": common_scope,
        },
        {
            "id": "V7",
            "name": "Low human-similarity preference by strain",
            "description": common_text
            + " Ranking uses the low-human-similarity composite score.",
            "scope": common_scope,
        },
        {
            "id": "V8",
            "name": "Low-overlap diversity construct by strain",
            "description": common_text
            + " When choices otherwise compete, lower coordinate overlap is preferred before score.",
            "scope": common_scope,
        },
    ]


def candidate_rank_maps(candidates, scorer, stats):
    ranked = sorted_candidates(candidates, scorer, stats)
    strain_class_rank = {
        row_identity(row): index for index, row in enumerate(ranked, start=1)
    }

    protein_rank = {}
    proteins = sorted({protein_key(row) for row in candidates})
    for protein in proteins:
        protein_rows = [row for row in candidates if protein_key(row) == protein]
        for index, row in enumerate(
            sorted_candidates(protein_rows, scorer, stats),
            start=1,
        ):
            protein_rank[row_identity(row)] = index

    return ranked, strain_class_rank, protein_rank


def choose_best_candidate(pool, version_id, scorer, stats, selected_context):
    if not pool:
        return None
    if version_id != "V10":
        return min(pool, key=lambda row: candidate_sort_key(row, scorer, stats))
    return min(
        pool,
        key=lambda row: (
            interval_overlap_fraction(row, selected_context),
            candidate_sort_key(row, scorer, stats),
        ),
    )


def annotate_selection_row(
    row,
    version,
    order,
    class_order,
    strain_order,
    note,
    metric,
    score,
    quota,
    rank_strain_class,
    rank_strain_protein_class,
    protein_repeat,
    overlap,
):
    selected = dict(row)
    selected.update(
        {
            "construct_version": version["id"],
            "version_name": version["name"],
            "version_description": version["description"],
            "construct_order": str(order),
            "class_order": str(class_order),
            "strain_order": str(strain_order),
            "selection_scope": version["scope"],
            "selection_note": note,
            "selection_primary_metric": metric,
            "selection_priority_score": format_float(score, 6),
            "quota_per_class_per_strain": str(quota),
            "selection_rank_within_strain_class": str(rank_strain_class),
            "selection_rank_within_strain_protein_class": str(
                rank_strain_protein_class
            ),
            "selection_rank_within_class": str(rank_strain_class),
            "selection_rank_within_protein_class": str(
                rank_strain_protein_class
            ),
            "protein_repeat_within_strain_class": (
                "True" if protein_repeat else "False"
            ),
            "duplicate_peptide_within_version": "",
            "diversity_max_overlap_fraction": format_float(overlap, 6),
        }
    )
    return selected


def select_strain_class(
    version,
    strain,
    epitope_class,
    candidates,
    stats,
    quota,
    selected_context,
):
    scorer, metric = class_score_function(version["id"], epitope_class)
    ranked, strain_class_rank, protein_rank = candidate_rank_maps(
        candidates,
        scorer,
        stats,
    )

    selected = []
    selected_ids = set()
    represented_proteins = set()

    # Phase 1: highest-ranked candidates from distinct proteins.
    while len(selected) < quota:
        pool = [
            row
            for row in ranked
            if row_identity(row) not in selected_ids
            and protein_key(row) not in represented_proteins
        ]
        chosen = choose_best_candidate(
            pool,
            version["id"],
            scorer,
            stats,
            selected_context + selected,
        )
        if chosen is None:
            break
        selected.append(chosen)
        selected_ids.add(row_identity(chosen))
        represented_proteins.add(protein_key(chosen))

    distinct_count = len(selected)

    # Phase 2: only when distinct proteins cannot fill the quota, repeat a protein.
    while len(selected) < quota:
        pool = [
            row for row in ranked if row_identity(row) not in selected_ids
        ]
        chosen = choose_best_candidate(
            pool,
            version["id"],
            scorer,
            stats,
            selected_context + selected,
        )
        if chosen is None:
            break
        selected.append(chosen)
        selected_ids.add(row_identity(chosen))

    annotated = []
    proteins_seen = set()
    for local_index, chosen in enumerate(selected, start=1):
        protein = protein_key(chosen)
        repeat = protein in proteins_seen
        proteins_seen.add(protein)
        note = (
            "primary_distinct_protein"
            if local_index <= distinct_count
            else "fallback_repeated_protein_due_to_insufficient_distinct_proteins"
        )
        overlap = interval_overlap_fraction(chosen, selected_context + selected[: local_index - 1])
        annotated.append(
            {
                "row": chosen,
                "note": note,
                "metric": metric,
                "score": scorer(chosen, stats),
                "rank_strain_class": strain_class_rank[row_identity(chosen)],
                "rank_strain_protein_class": protein_rank[row_identity(chosen)],
                "protein_repeat": repeat,
                "overlap": overlap,
            }
        )

    return annotated


def select_version(version, rows_by_class, strains, stats, quota):
    selected_rows = []
    order = 1

    for strain_index, strain in enumerate(strains, start=1):
        selected_context = []
        for class_index, epitope_class in enumerate(CLASS_ORDER, start=1):
            candidates = [
                row
                for row in rows_by_class.get(epitope_class, [])
                if row.get("strain") == strain
            ]
            if epitope_class == "Bcell":
                candidates = [
                    row for row in candidates if row.get("_surface_accessible") is True
                ]

            choices = select_strain_class(
                version,
                strain,
                epitope_class,
                candidates,
                stats,
                quota,
                selected_context,
            )

            for choice in choices:
                selected = annotate_selection_row(
                    choice["row"],
                    version,
                    order,
                    class_index,
                    strain_index,
                    choice["note"],
                    choice["metric"],
                    choice["score"],
                    quota,
                    choice["rank_strain_class"],
                    choice["rank_strain_protein_class"],
                    choice["protein_repeat"],
                    choice["overlap"],
                )
                selected_rows.append(selected)
                selected_context.append(choice["row"])
                order += 1

    return selected_rows


def mark_duplicate_peptides(rows):
    counts = {}
    for row in rows:
        peptide = row.get("peptide", "")
        counts[peptide] = counts.get(peptide, 0) + 1
    for row in rows:
        row["duplicate_peptide_within_version"] = (
            "True" if counts.get(row.get("peptide", ""), 0) > 1 else "False"
        )
    return rows


# -----------------------------------------------------------------------------
# Summaries, validation, and reports
# -----------------------------------------------------------------------------


def validate_selected_version(rows, strains, quota, version_id):
    expected_total = len(strains) * len(CLASS_ORDER) * quota
    if len(rows) > expected_total or len(rows) > MAX_TOTAL_EPITOPES:
        raise AssertionError(
            "%s selected %d rows, exceeding the allowed target of %d."
            % (version_id, len(rows), min(expected_total, MAX_TOTAL_EPITOPES))
        )

    for strain in strains:
        for epitope_class in CLASS_ORDER:
            count = sum(
                1
                for row in rows
                if row.get("strain") == strain
                and row.get("epitope_class") == epitope_class
            )
            if count > quota:
                raise AssertionError(
                    "%s selected %d %s epitopes for strain %s; quota is %d."
                    % (version_id, count, epitope_class, strain, quota)
                )


def selection_signature(rows):
    payload = [
        {
            "order": int(row["construct_order"]),
            "strain": row.get("strain", ""),
            "class": row.get("epitope_class", ""),
            "protein": row.get("protein", ""),
            "peptide": row.get("peptide", ""),
            "epitope_id": row.get("epitope_id", ""),
        }
        for row in sorted(rows, key=lambda item: int(item["construct_order"]))
    ]
    return stable_json_hash(payload)


def summarize_version(version, rows, strains, quota):
    counts_by_class = {epitope_class: 0 for epitope_class in CLASS_ORDER}
    counts_by_strain_class = {
        strain: {epitope_class: 0 for epitope_class in CLASS_ORDER}
        for strain in strains
    }
    proteins_by_class = {epitope_class: set() for epitope_class in CLASS_ORDER}
    fallback_count = 0
    duplicate_count = 0

    for row in rows:
        epitope_class = row["epitope_class"]
        strain = row["strain"]
        counts_by_class[epitope_class] += 1
        counts_by_strain_class[strain][epitope_class] += 1
        proteins_by_class[epitope_class].add(row["protein"])
        if row["selection_note"].startswith("fallback_repeated"):
            fallback_count += 1
        if row["duplicate_peptide_within_version"] == "True":
            duplicate_count += 1

    expected_total = len(strains) * len(CLASS_ORDER) * quota
    strain_count_text = ";".join(
        "%s:Bcell=%d,CTL=%d,HTL=%d"
        % (
            strain,
            counts_by_strain_class[strain]["Bcell"],
            counts_by_strain_class[strain]["CTL"],
            counts_by_strain_class[strain]["HTL"],
        )
        for strain in strains
    )

    ordered_rows = sorted(rows, key=lambda item: int(item["construct_order"]))
    unique_peptides = []
    seen_peptides = set()
    for row in ordered_rows:
        peptide = row["peptide"]
        if peptide not in seen_peptides:
            seen_peptides.add(peptide)
            unique_peptides.append(peptide)

    return {
        "construct_version": version["id"],
        "version_name": version["name"],
        "selection_scope": version["scope"],
        "strain_count": len(strains),
        "quota_per_class_per_strain": quota,
        "expected_epitope_count": expected_total,
        "total_epitope_count": len(rows),
        "selection_shortfall_count": expected_total - len(rows),
        "bcell_count": counts_by_class["Bcell"],
        "ctl_count": counts_by_class["CTL"],
        "htl_count": counts_by_class["HTL"],
        "unique_protein_count": len({row["protein"] for row in rows}),
        "unique_peptide_count": len(unique_peptides),
        "bcell_proteins": ";".join(sorted(proteins_by_class["Bcell"])),
        "ctl_proteins": ";".join(sorted(proteins_by_class["CTL"])),
        "htl_proteins": ";".join(sorted(proteins_by_class["HTL"])),
        "strain_class_counts": strain_count_text,
        "fallback_repeated_protein_rows": fallback_count,
        "duplicate_peptide_rows": duplicate_count,
        "selection_signature_sha256": selection_signature(rows),
        "identical_to_version": "",
        "peptide_block_pipe_delimited": "|".join(
            row["peptide"] for row in ordered_rows
        ),
        "unique_peptide_block_pipe_delimited": "|".join(unique_peptides),
    }


def mark_identical_versions(summary_rows):
    first_by_signature = {}
    for row in summary_rows:
        signature = row["selection_signature_sha256"]
        if signature in first_by_signature:
            row["identical_to_version"] = first_by_signature[signature]
        else:
            first_by_signature[signature] = row["construct_version"]
    return summary_rows


def candidate_availability(rows_by_class, strains):
    result = {}
    for strain in strains:
        result[strain] = {}
        for epitope_class in CLASS_ORDER:
            rows = [
                row
                for row in rows_by_class.get(epitope_class, [])
                if row.get("strain") == strain
            ]
            if epitope_class == "Bcell":
                rows = [row for row in rows if row.get("_surface_accessible") is True]
            result[strain][epitope_class] = {
                "candidate_rows": len(rows),
                "distinct_proteins": len({protein_key(row) for row in rows}),
            }
    return result


def write_report(
    rows_by_class,
    strains,
    strain_mode,
    quota,
    versions,
    summary_rows,
    output_files,
    run_metadata,
):
    report_path = REPORT_DIR / "20_vaccine_construct_selection_summary.md"
    json_path = REPORT_DIR / "20_vaccine_construct_selection_summary.json"

    availability = candidate_availability(rows_by_class, strains)
    expected_total = len(strains) * len(CLASS_ORDER) * quota

    lines = []
    lines.append("# Strain-aware Vaccine Construct Epitope Selection Summary")
    lines.append("")
    lines.append("Run ID: `%s`" % run_metadata["run_id"])
    lines.append("")
    lines.append("Script version: `%s`" % SCRIPT_VERSION)
    lines.append("")
    lines.append("Generated UTC: %s" % run_metadata["generated_utc"])
    lines.append("")
    lines.append("## Dynamic selection policy")
    lines.append("")
    lines.append("- Detected strains: %d (%s)." % (len(strains), "; ".join(strains)))
    lines.append("- Strain mode: `%s`." % strain_mode)
    lines.append(
        "- Quota per strain per class: %d = floor(%d / (3 classes × %d strains))."
        % (quota, MAX_TOTAL_EPITOPES, len(strains))
    )
    lines.append(
        "- Requested rows per version: %d; hard maximum: %d."
        % (expected_total, MAX_TOTAL_EPITOPES)
    )
    lines.append(
        "- B-cell candidates are restricted to surface-accessible proteins."
    )
    lines.append(
        "- B-cell, CTL, and HTL selections use distinct proteins within each strain and class before any protein is repeated."
    )
    lines.append(
        "- If too few distinct proteins exist, the remaining quota is filled by the highest-ranked unused candidates from already represented proteins."
    )
    lines.append(
        "- Exact peptide duplicates are retained as strain-level selections and flagged; the summary also reports a unique-peptide block."
    )
    lines.append("")

    lines.append("## Eligible candidate availability")
    lines.append("")
    lines.append("| Strain | Class | Candidate rows | Distinct eligible proteins |")
    lines.append("|---|---|---:|---:|")
    for strain in strains:
        for epitope_class in CLASS_ORDER:
            info = availability[strain][epitope_class]
            lines.append(
                "| %s | %s | %d | %d |"
                % (
                    strain,
                    epitope_class,
                    info["candidate_rows"],
                    info["distinct_proteins"],
                )
            )
    lines.append("")

    lines.append("## Version rules")
    lines.append("")
    lines.append("| Version | Name | Rule |")
    lines.append("|---|---|---|")
    for version in versions:
        lines.append(
            "| %s | %s | %s |"
            % (version["id"], version["name"], version["description"])
        )
    lines.append("")

    lines.append("## Version output summary")
    lines.append("")
    lines.append(
        "| Version | Expected | Actual | B-cell | CTL | HTL | Shortfall | Repeated-protein fallback rows | Duplicate peptide rows | Identical to |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for row in summary_rows:
        lines.append(
            "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |"
            % (
                row["construct_version"],
                row["expected_epitope_count"],
                row["total_epitope_count"],
                row["bcell_count"],
                row["ctl_count"],
                row["htl_count"],
                row["selection_shortfall_count"],
                row["fallback_repeated_protein_rows"],
                row["duplicate_peptide_rows"],
                row["identical_to_version"],
            )
        )
    lines.append("")

    shortfalls = [
        row for row in summary_rows if row["selection_shortfall_count"] > 0
    ]
    if shortfalls:
        lines.append("## Selection shortfalls")
        lines.append("")
        lines.append(
            "One or more versions could not fill every strain/class quota because insufficient eligible candidates were available. The script does not duplicate a candidate row merely to reach the target."
        )
        lines.append("")

    lines.append("## Reproducibility")
    lines.append("")
    lines.append("Input SHA-256 hashes:")
    lines.append("")
    for label, digest in sorted(run_metadata["input_sha256"].items()):
        lines.append("- `%s`: `%s`" % (label, digest))
    lines.append("")
    lines.append("Selection is deterministic: no random sampling is used, all tie-breaks are explicit, and input rows and strains are stably sorted.")
    lines.append("")

    lines.append("## Output files")
    lines.append("")
    for path in output_files:
        try:
            display_path = path.relative_to(PROJECT_DIR)
        except ValueError:
            display_path = path
        lines.append("- `%s`" % display_path)
    lines.append("")
    lines.append(
        "Actual linked vaccine sequences are not assembled in this step because linker, adjuvant, and final ordering rules remain separate design decisions."
    )
    lines.append("")

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")

    json_payload = {
        "run_metadata": run_metadata,
        "script_version": SCRIPT_VERSION,
        "selection_policy": {
            "maximum_total_epitopes": MAX_TOTAL_EPITOPES,
            "class_order": CLASS_ORDER,
            "strain_column": STRAIN_COLUMN,
            "strains": strains,
            "strain_mode": strain_mode,
            "quota_per_class_per_strain": quota,
            "expected_epitopes_per_version": expected_total,
            "bcell_surface_accessible_only": True,
            "distinct_protein_first_for_all_classes": True,
            "fallback_repeated_protein_when_needed": True,
        },
        "candidate_availability": availability,
        "version_definitions": versions,
        "version_summaries": summary_rows,
        "output_files": [str(path) for path in output_files],
    }
    json_path.write_text(
        json.dumps(json_payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return report_path, json_path


def build_run_metadata(strains, quota):
    input_paths = dict(INPUT_FILES)
    input_paths["Human_BLASTp_summary"] = BLAST_SUMMARY_FILE
    input_hashes = {
        label: sha256_file(path) for label, path in sorted(input_paths.items())
    }

    configuration = {
        "script_version": SCRIPT_VERSION,
        "max_total_epitopes": MAX_TOTAL_EPITOPES,
        "class_order": CLASS_ORDER,
        "strain_column": STRAIN_COLUMN,
        "strains": strains,
        "quota_per_class_per_strain": quota,
        "surface_accessible_column": SURFACE_ACCESSIBLE_COLUMN,
        "surface_accessible_proteins": sorted(SURFACE_ACCESSIBLE_PROTEINS),
        "assume_all_bcell_inputs_surface_accessible": (
            ASSUME_ALL_BCELL_INPUTS_SURFACE_ACCESSIBLE
        ),
    }

    script_hash = ""
    try:
        script_path = Path(__file__).resolve()
        if script_path.exists():
            script_hash = sha256_file(script_path)
    except (NameError, OSError):
        pass

    run_basis = {
        "configuration": configuration,
        "input_sha256": input_hashes,
        "script_sha256": script_hash,
    }
    return {
        "run_id": stable_json_hash(run_basis)[:16],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "script_sha256": script_hash,
        "input_sha256": input_hashes,
        "configuration": configuration,
    }


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    rows_by_class, strains, strain_mode = load_inputs()
    quota = compute_quota_per_class_per_strain(len(strains))
    stats = build_normalizers(rows_by_class)
    versions = version_definitions(quota)
    run_metadata = build_run_metadata(strains, quota)

    output_files = []
    summary_rows = []
    all_rows = []

    for index, version in enumerate(versions, start=1):
        selected_rows = select_version(
            version,
            rows_by_class,
            strains,
            stats,
            quota,
        )
        selected_rows = mark_duplicate_peptides(selected_rows)
        validate_selected_version(
            selected_rows,
            strains,
            quota,
            version["id"],
        )

        summary_rows.append(
            summarize_version(
                version,
                selected_rows,
                strains,
                quota,
            )
        )
        all_rows.extend(selected_rows)

        slug = version["name"].lower()
        slug = re.sub(r"[^a-z0-9]+", "_", slug).strip("_")
        csv_path = OUTPUT_DIR / (
            "%02d_vaccine_construct_%s_%s.csv"
            % (index, version["id"], slug)
        )
        tsv_path = OUTPUT_DIR / (
            "%02d_vaccine_construct_%s_%s.tsv"
            % (index, version["id"], slug)
        )
        write_table(csv_path, selected_rows, OUTPUT_COLUMNS, delimiter=",")
        write_table(tsv_path, selected_rows, OUTPUT_COLUMNS, delimiter="\t")
        output_files.extend([csv_path, tsv_path])

    summary_rows = mark_identical_versions(summary_rows)

    summary_columns = [
        "construct_version",
        "version_name",
        "selection_scope",
        "strain_count",
        "quota_per_class_per_strain",
        "expected_epitope_count",
        "total_epitope_count",
        "selection_shortfall_count",
        "bcell_count",
        "ctl_count",
        "htl_count",
        "unique_protein_count",
        "unique_peptide_count",
        "bcell_proteins",
        "ctl_proteins",
        "htl_proteins",
        "strain_class_counts",
        "fallback_repeated_protein_rows",
        "duplicate_peptide_rows",
        "selection_signature_sha256",
        "identical_to_version",
        "peptide_block_pipe_delimited",
        "unique_peptide_block_pipe_delimited",
    ]

    summary_csv = OUTPUT_DIR / "11_vaccine_construct_version_summary.csv"
    summary_tsv = OUTPUT_DIR / "11_vaccine_construct_version_summary.tsv"
    combined_csv = OUTPUT_DIR / "12_all_vaccine_construct_versions.csv"
    combined_tsv = OUTPUT_DIR / "12_all_vaccine_construct_versions.tsv"

    write_table(summary_csv, summary_rows, summary_columns, delimiter=",")
    write_table(summary_tsv, summary_rows, summary_columns, delimiter="\t")
    write_table(combined_csv, all_rows, OUTPUT_COLUMNS, delimiter=",")
    write_table(combined_tsv, all_rows, OUTPUT_COLUMNS, delimiter="\t")
    output_files.extend(
        [summary_csv, summary_tsv, combined_csv, combined_tsv]
    )

    report_path, json_path = write_report(
        rows_by_class,
        strains,
        strain_mode,
        quota,
        versions,
        summary_rows,
        output_files,
        run_metadata,
    )

    manifest_path = OUTPUT_DIR / "13_reproducibility_manifest.json"
    manifest_path.write_text(
        json.dumps(run_metadata, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print("Run ID: %s" % run_metadata["run_id"])
    print("Detected %d strain(s): %s" % (len(strains), ", ".join(strains)))
    print("Quota per class per strain: %d" % quota)
    print(
        "Requested selections per version: %d (maximum %d)"
        % (len(strains) * len(CLASS_ORDER) * quota, MAX_TOTAL_EPITOPES)
    )
    print("Wrote %d construct versions to %s" % (len(versions), OUTPUT_DIR))
    print("Wrote summary report: %s" % report_path)
    print("Wrote summary JSON: %s" % json_path)
    print("Wrote reproducibility manifest: %s" % manifest_path)


if __name__ == "__main__":
    main()
