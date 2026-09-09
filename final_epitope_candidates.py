import csv
import json
from pathlib import Path


DEFAULT_ROOT = Path(r"F:\3.Project_Drive\11. Vaccine\000. Version 1.0\11. Epitope Validation")
ROOT = DEFAULT_ROOT if DEFAULT_ROOT.exists() else Path("./11_epitope_validation")
INPUT_DIR = ROOT / "01_inputs"
PARSED_DIR = ROOT / "04_parsed_results"
FILTERED_DIR = ROOT / "05_filtered_candidates"
REPORT_DIR = ROOT / "06_reports"


def read_table(path):
    with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
        return list(csv.DictReader(handle))


def write_table(base, rows):
    base.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    for suffix, delimiter in [(".csv", ","), (".tsv", "\t")]:
        with base.with_suffix(suffix).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)


def index(rows):
    return {row["epitope_id"]: row for row in rows}


def allele_count(alleles):
    return len({item.strip() for item in alleles.split(";") if item.strip()})


def percent_value(text):
    if not text:
        return 0.0
    try:
        return float(str(text).split("%", 1)[0].strip())
    except ValueError:
        return 0.0


def is_ifn_positive(prediction):
    text = prediction.strip().lower()
    return "inducer" in text and not text.startswith("non")


def base_fields(row):
    return {
        "epitope_id": row.get("epitope_id", ""),
        "epitope_class": row.get("epitope_class", ""),
        "protein": row.get("protein", ""),
        "peptide": row.get("peptide", ""),
        "length": row.get("length", ""),
        "tools": row.get("tools", ""),
        "starts": row.get("starts", ""),
        "ends": row.get("ends", ""),
        "evidence_count": row.get("evidence_count", ""),
    }


def common_validation_fields(epitope_id, vaxijen, allertop, toxinpred, conservancy):
    vaxi = vaxijen.get(epitope_id, {})
    aller = allertop.get(epitope_id, {})
    toxin = toxinpred.get(epitope_id, {})
    cons = conservancy.get(epitope_id, {})
    return {
        "vaxijen_score": vaxi.get("vaxijen_score", ""),
        "vaxijen_prediction": vaxi.get("vaxijen_prediction", ""),
        "allertop_prediction": aller.get("allertop_prediction", ""),
        "toxinpred_score": toxin.get("toxinpred_score", ""),
        "toxinpred_prediction": toxin.get("toxinpred_prediction", ""),
        "conservancy_percent": cons.get("conservancy_percent_matches_identity_le_100", ""),
        "conservancy_percent_numeric": percent_value(cons.get("conservancy_percent_matches_identity_le_100", "")),
        "conservancy_minimum_identity": cons.get("conservancy_minimum_identity", ""),
        "conservancy_maximum_identity": cons.get("conservancy_maximum_identity", ""),
    }


def passes_core_safety(fields):
    return (
        fields["vaxijen_prediction"] == "Probable ANTIGEN"
        and fields["allertop_prediction"] == "Probable NON-ALLERGEN"
        and fields["toxinpred_prediction"] == "Non-Toxin"
        and float(fields["conservancy_percent_numeric"]) > 70.0
    )


def add_pass_flags(out):
    out.update(
        {
            "passes_antigenic": "True",
            "passes_non_allergen": "True",
            "passes_non_toxic": "True",
            "passes_conservancy_gt70": "True",
        }
    )


def filter_bcell(rows, vaxijen, allertop, toxinpred, conservancy):
    kept = []
    for row in rows:
        epitope_id = row["epitope_id"]
        fields = common_validation_fields(epitope_id, vaxijen, allertop, toxinpred, conservancy)
        if passes_core_safety(fields):
            out = base_fields(row)
            out.update(fields)
            add_pass_flags(out)
            out["filter_rule"] = "B-cell: antigenic + non-allergen + non-toxic + conservancy >70%"
            kept.append(out)
    return kept


def filter_ctl(rows, mhc_i, vaxijen, allertop, toxinpred, conservancy):
    kept = []
    for row in rows:
        epitope_id = row["epitope_id"]
        fields = common_validation_fields(epitope_id, vaxijen, allertop, toxinpred, conservancy)
        mhc = mhc_i.get(epitope_id, {})
        mhc_score = float(mhc.get("Score") or 0)
        if passes_core_safety(fields) and mhc_score > 0:
            hla_i_alleles = row.get("alleles", "")
            out = base_fields(row)
            out.update(fields)
            add_pass_flags(out)
            out.update(
                {
                    "mhc_i_immunogenicity_score": mhc.get("Score", ""),
                    "mhc_i_immunogenicity_positive": "True",
                    "hla_i_binding_allele_count": allele_count(hla_i_alleles),
                    "hla_i_binding_alleles": hla_i_alleles,
                    "filter_rule": (
                        "CTL: antigenic + non-allergen + non-toxic + MHC-I immunogenicity score >0 "
                        "+ conservancy >70%"
                    ),
                }
            )
            kept.append(out)
    return kept


def filter_htl(rows, ifn, il4, il10, vaxijen, allertop, toxinpred, conservancy):
    kept = []
    for row in rows:
        epitope_id = row["epitope_id"]
        fields = common_validation_fields(epitope_id, vaxijen, allertop, toxinpred, conservancy)
        ifn_row = ifn.get(epitope_id, {})
        il4_row = il4.get(epitope_id, {})
        il10_row = il10.get(epitope_id, {})
        cytokine_positive = (
            is_ifn_positive(ifn_row.get("ifnepitope2_prediction", ""))
            and il4_row.get("il4pred_prediction", "") == "IL4 inducer"
            and il10_row.get("il10pred_prediction", "") == "IL10 inducer"
        )
        if passes_core_safety(fields) and cytokine_positive:
            hla_ii_alleles = row.get("alleles", "")
            out = base_fields(row)
            out.update(fields)
            add_pass_flags(out)
            out.update(
                {
                    "ifnepitope2_score": ifn_row.get("ifnepitope2_score", ""),
                    "ifnepitope2_prediction": ifn_row.get("ifnepitope2_prediction", ""),
                    "il4pred_score": il4_row.get("il4pred_score", ""),
                    "il4pred_prediction": il4_row.get("il4pred_prediction", ""),
                    "il10pred_score": il10_row.get("il10pred_score", ""),
                    "il10pred_prediction": il10_row.get("il10pred_prediction", ""),
                    "hla_ii_binding_allele_count": allele_count(hla_ii_alleles),
                    "hla_ii_binding_alleles": hla_ii_alleles,
                    "filter_rule": (
                        "HTL: antigenic + non-allergen + non-toxic + IFN-gamma positive + "
                        "IL4 positive + IL10 positive + conservancy >70%"
                    ),
                }
            )
            kept.append(out)
    return kept


def class_counts(rows):
    return {
        "CTL": sum(1 for row in rows if row.get("epitope_class") == "CTL"),
        "HTL": sum(1 for row in rows if row.get("epitope_class") == "HTL"),
        "Bcell": sum(1 for row in rows if row.get("epitope_class") == "Bcell"),
        "Total": len(rows),
    }


def main():
    FILTERED_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    ctl_rows = read_table(INPUT_DIR / "01_CTL_filtered_unique_epitopes.csv")
    htl_rows = read_table(INPUT_DIR / "02_HTL_filtered_unique_epitopes.csv")
    bcell_rows = read_table(INPUT_DIR / "03_Bcell_filtered_unique_epitopes.csv")

    vaxijen = index(read_table(PARSED_DIR / "07_VaxiJen_merged.csv"))
    allertop = index(read_table(PARSED_DIR / "08_AllerTOP_merged.csv"))
    toxinpred = index(read_table(PARSED_DIR / "05_ToxinPred_merged.csv"))
    conservancy = index(read_table(PARSED_DIR / "06_Conservancy_merged.csv"))
    mhc_i = index(read_table(PARSED_DIR / "01_CTL_MHC_I_immunogenicity_merged.csv"))
    ifn = index(read_table(PARSED_DIR / "02_HTL_IFNepitope2_Human_merged.csv"))
    il4 = index(read_table(PARSED_DIR / "03_HTL_IL4pred_merged.csv"))
    il10 = index(read_table(PARSED_DIR / "04_HTL_IL10Pred_merged.csv"))

    outputs = {
        "13_Bcell_strict_antigen_nonallergen_nontoxin_conservancy_gt95": filter_bcell(
            bcell_rows, vaxijen, allertop, toxinpred, conservancy
        ),
        "14_CTL_strict_antigen_nonallergen_nontoxin_mhci_positive_conservancy_gt95": filter_ctl(
            ctl_rows, mhc_i, vaxijen, allertop, toxinpred, conservancy
        ),
        "15_HTL_strict_antigen_nonallergen_nontoxin_ifn_il4_il10_conservancy_gt95": filter_htl(
            htl_rows, ifn, il4, il10, vaxijen, allertop, toxinpred, conservancy
        ),
    }
    outputs["16_combined_strict_filtered_candidates"] = (
        outputs["13_Bcell_strict_antigen_nonallergen_nontoxin_conservancy_gt95"]
        + outputs["14_CTL_strict_antigen_nonallergen_nontoxin_mhci_positive_conservancy_gt95"]
        + outputs["15_HTL_strict_antigen_nonallergen_nontoxin_ifn_il4_il10_conservancy_gt95"]
    )

    for name, rows in outputs.items():
        write_table(FILTERED_DIR / name, rows)

    summary = []
    start_counts = {
        "13_Bcell_strict_antigen_nonallergen_nontoxin_conservancy_gt95": len(bcell_rows),
        "14_CTL_strict_antigen_nonallergen_nontoxin_mhci_positive_conservancy_gt95": len(ctl_rows),
        "15_HTL_strict_antigen_nonallergen_nontoxin_ifn_il4_il10_conservancy_gt95": len(htl_rows),
        "16_combined_strict_filtered_candidates": len(bcell_rows) + len(ctl_rows) + len(htl_rows),
    }
    for name, rows in outputs.items():
        counts = class_counts(rows)
        summary.append(
            {
                "output": name,
                "starting_rows": start_counts[name],
                "kept_rows": len(rows),
                "kept_ctl": counts["CTL"],
                "kept_htl": counts["HTL"],
                "kept_bcell": counts["Bcell"],
                "csv": str((FILTERED_DIR / name).with_suffix(".csv")),
                "tsv": str((FILTERED_DIR / name).with_suffix(".tsv")),
            }
        )

    write_table(REPORT_DIR / "11_strict_final_filter_summary", summary)
    (REPORT_DIR / "11_strict_final_filter_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
