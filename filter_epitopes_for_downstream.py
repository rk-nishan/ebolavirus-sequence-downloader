import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path


DEFAULT_BASE = Path(r"F:\3.Project_Drive\11. Vaccine\000. Version 1.0\10. Epitope prediction")
BASE = DEFAULT_BASE if DEFAULT_BASE.exists() else Path("./10_epitope_prediction")
INPUT = BASE / "merged_deduplicated" / "all_predictions_exact_deduplicated.tsv"
OUT = BASE / "filtered_cleaned"

CTL_MAX_PERCENTILE = 1.0
HTL_MAX_PERCENTILE = 10.0
BCELL_MIN_LENGTH = 6

FIELDS = [
    "epitope_class",
    "protein",
    "tool",
    "peptide",
    "start",
    "end",
    "length",
    "allele",
    "score",
    "rank",
    "bind_level",
    "source_table",
    "source_file",
    "extra_json",
]


def parse_float(value):
    try:
        text = str(value).strip().replace("%", "")
        if not text:
            return None
        return float(text)
    except ValueError:
        return None


def parse_int(value):
    try:
        return int(float(str(value).strip()))
    except ValueError:
        return None


def keep_row(row, warnings):
    epitope_class = row["epitope_class"]
    if epitope_class == "CTL":
        percentile = parse_float(row.get("rank", ""))
        if percentile is None:
            warnings.append(f"Excluded CTL row with nonnumeric percentile: {row.get('protein')} {row.get('tool')} {row.get('peptide')}")
            return False
        return percentile <= CTL_MAX_PERCENTILE
    if epitope_class == "HTL":
        percentile = parse_float(row.get("rank", ""))
        if percentile is None:
            warnings.append(f"Excluded HTL row with nonnumeric percentile: {row.get('protein')} {row.get('tool')} {row.get('peptide')}")
            return False
        return percentile <= HTL_MAX_PERCENTILE
    if epitope_class == "Bcell":
        length = parse_int(row.get("length", ""))
        if length is None:
            length = len(row.get("peptide", ""))
        return length >= BCELL_MIN_LENGTH
    warnings.append(f"Excluded row with unknown epitope class: {epitope_class}")
    return False


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sort_number_set(values):
    return ";".join(sorted(values, key=lambda x: int(x) if str(x).isdigit() else 0))


def aggregate_unique(rows, key_fields):
    groups = {}
    for row in rows:
        key = tuple(row.get(field, "") for field in key_fields)
        item = groups.setdefault(
            key,
            {
                **{field: row.get(field, "") for field in key_fields},
                "length": row.get("length", ""),
                "tools": set(),
                "alleles": set(),
                "starts": set(),
                "ends": set(),
                "evidence_count": 0,
                "source_files": set(),
            },
        )
        item["tools"].add(row.get("tool", ""))
        if row.get("allele"):
            item["alleles"].add(row["allele"])
        if row.get("start"):
            item["starts"].add(row["start"])
        if row.get("end"):
            item["ends"].add(row["end"])
        if row.get("source_file"):
            item["source_files"].add(row["source_file"])
        item["evidence_count"] += 1

    output = []
    for item in groups.values():
        row = dict(item)
        row["tools"] = ";".join(sorted(item["tools"]))
        row["alleles"] = ";".join(sorted(item["alleles"]))
        row["starts"] = sort_number_set(item["starts"])
        row["ends"] = sort_number_set(item["ends"])
        row["source_files"] = ";".join(sorted(item["source_files"]))
        output.append(row)
    return sorted(output, key=lambda r: (r.get("epitope_class", ""), r.get("protein", ""), r.get("peptide", ""), r.get("start", "")))


def count_by_class_tool(rows):
    counts = defaultdict(lambda: {"rows": 0, "unique_peptides": set()})
    for row in rows:
        key = (row["epitope_class"], row["protein"], row["tool"])
        counts[key]["rows"] += 1
        counts[key]["unique_peptides"].add(row["peptide"])
    out = []
    for (epitope_class, protein, tool), item in sorted(counts.items()):
        out.append(
            {
                "epitope_class": epitope_class,
                "protein": protein,
                "tool": tool,
                "rows": item["rows"],
                "unique_peptide_sequences": len(item["unique_peptides"]),
            }
        )
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    warnings = []
    all_rows = []
    kept_rows = []

    with INPUT.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            all_rows.append(row)
            if keep_row(row, warnings):
                kept_rows.append(row)

    write_tsv(OUT / "all_predictions_filtered.tsv", kept_rows, FIELDS)
    for epitope_class, filename in [
        ("CTL", "ctl_epitopes_filtered_netmhcpan_el_percentile_le1.tsv"),
        ("HTL", "htl_epitopes_filtered_netmhciipan_el_percentile_le10.tsv"),
        ("Bcell", "bcell_epitopes_filtered_min_length6.tsv"),
    ]:
        write_tsv(OUT / filename, [row for row in kept_rows if row["epitope_class"] == epitope_class], FIELDS)

    by_position = aggregate_unique(kept_rows, ["epitope_class", "protein", "peptide", "start", "end"])
    write_tsv(OUT / "unique_filtered_epitopes_by_position.tsv", by_position, ["epitope_class", "protein", "peptide", "start", "end", "length", "tools", "alleles", "evidence_count", "source_files"])

    by_sequence = aggregate_unique(kept_rows, ["epitope_class", "protein", "peptide"])
    write_tsv(OUT / "unique_filtered_epitopes_by_sequence.tsv", by_sequence, ["epitope_class", "protein", "peptide", "length", "tools", "alleles", "starts", "ends", "evidence_count", "source_files"])

    write_tsv(OUT / "filtering_counts_before.tsv", count_by_class_tool(all_rows), ["epitope_class", "protein", "tool", "rows", "unique_peptide_sequences"])
    write_tsv(OUT / "filtering_counts_after.tsv", count_by_class_tool(kept_rows), ["epitope_class", "protein", "tool", "rows", "unique_peptide_sequences"])
    (OUT / "filtering_warnings.txt").write_text("No filtering warnings.\n" if not warnings else "\n".join(warnings) + "\n", encoding="utf-8")

    totals_before = defaultdict(int)
    totals_after = defaultdict(int)
    for row in all_rows:
        totals_before[row["epitope_class"]] += 1
    for row in kept_rows:
        totals_after[row["epitope_class"]] += 1

    report = [
        "# Downstream Epitope Filtering Report",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Filters Applied",
        "",
        f"- CTL: retained rows with NetMHCpan EL percentile/rank <= {CTL_MAX_PERCENTILE:g}.",
        f"- HTL: retained rows with NetMHCIIpan EL percentile/rank <= {HTL_MAX_PERCENTILE:g}.",
        f"- B-cell: retained rows with peptide length >= {BCELL_MIN_LENGTH}.",
        "",
        "## Totals",
        "",
        "| Class | Before | After |",
        "|---|---:|---:|",
    ]
    for epitope_class in ["CTL", "HTL", "Bcell"]:
        report.append(f"| {epitope_class} | {totals_before[epitope_class]} | {totals_after[epitope_class]} |")
    report.extend(["", "## Validation", "", "No filtering warnings." if not warnings else f"{len(warnings)} filtering warnings; see filtering_warnings.txt."])
    (OUT / "downstream_filtering_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Input rows: {len(all_rows)}")
    print(f"Filtered rows: {len(kept_rows)}")
    print(f"Warnings: {len(warnings)}")
    print(f"Output folder: {OUT}")


if __name__ == "__main__":
    main()
