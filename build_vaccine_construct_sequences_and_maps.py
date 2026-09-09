import csv
import hashlib
import html
import json
import math
import re
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


DEFAULT_ROOT = Path(r"F:\3.Project_Drive\11. Vaccine\000. Version 1.0")
ROOT = DEFAULT_ROOT if DEFAULT_ROOT.exists() else Path(".")
VERSION_DIR = ROOT / "11. Epitope Validation" / "12. vaccine_constructs"
OUTPUT_DIR = ROOT / "12. Vaccine_construct_designs"
FIGURE_DIR = OUTPUT_DIR / "figures"
MASTER_FIGURE_DIR = FIGURE_DIR / "master_architecture_maps"
REPORT_DIR = ROOT / "11. Epitope Validation" / "06_reports"

ADJUVANTS = [
    {
        "code": "BD",
        "name": "Beta_defensin",
        "sequence": "GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK",
        "description": "Human beta-defensin (hBD) adjuvant sequence.",
    },
    {
        "code": "RS09",
        "name": "RS09",
        "sequence": "APPHALS",
        "description": "RS09 TLR4 agonist peptide sequence.",
    },
]

LINKERS = {
    "adjuvant_linker": "EAAAK",
    "block_linker": "HEYGAEALERAG",
    "CTL": "AAY",
    "HTL": "GPGPG",
    "Bcell": "KK",
}

ARCHITECTURES = [
    {
        "id": "L01",
        "name": "Direct CTL-HTL-Bcell blocks",
        "short_name": "Direct_CTL_HTL_Bcell",
        "block_order": ["CTL", "HTL", "Bcell"],
        "internal_linkers": {"CTL": "", "HTL": "", "Bcell": ""},
        "formula": "Adjuvant-EAAAK-[CTL direct block]-HEYGAEALERAG-[HTL direct block]-HEYGAEALERAG-[B-cell direct block]",
    },
    {
        "id": "L02",
        "name": "Bcell-KK then CTL-AAY then HTL-GPGPG",
        "short_name": "Bcell_CTL_HTL_class_linkers",
        "block_order": ["Bcell", "CTL", "HTL"],
        "internal_linkers": {"CTL": "AAY", "HTL": "GPGPG", "Bcell": "KK"},
        "formula": "Adjuvant-EAAAK-[B-cell-KK block]-HEYGAEALERAG-[CTL-AAY block]-HEYGAEALERAG-[HTL-GPGPG block]",
    },
    {
        "id": "L03",
        "name": "CTL-AAY then HTL-GPGPG then Bcell-KK",
        "short_name": "CTL_HTL_Bcell_class_linkers",
        "block_order": ["CTL", "HTL", "Bcell"],
        "internal_linkers": {"CTL": "AAY", "HTL": "GPGPG", "Bcell": "KK"},
        "formula": "Adjuvant-EAAAK-[CTL-AAY block]-HEYGAEALERAG-[HTL-GPGPG block]-HEYGAEALERAG-[B-cell-KK block]",
    },
]

AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
CLASS_ORDER = ["CTL", "HTL", "Bcell"]

COLORS = {
    "background": "#F7F9FC",
    "ink": "#172033",
    "muted": "#64748B",
    "line": "#CBD5E1",
    "adjuvant": "#C64636",
    "EAAAK": "#7B2CBF",
    "HEYGAEALERAG": "#3446B3",
    "AAY": "#F27E2E",
    "GPGPG": "#45A763",
    "KK": "#D5A300",
    "CTL": "#2577B8",
    "HTL": "#52AD5E",
    "Bcell": "#E7B93E",
    "epitope_outline": "#FFFFFF",
}


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
        return list(csv.DictReader(handle))


def write_table(path, rows, columns, delimiter=","):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def wrap_sequence(sequence, width=80):
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def numeric(row, column):
    value = row.get(column, "")
    if value in ("", None):
        return None
    try:
        return float(value)
    except ValueError:
        match = re.search(r"-?\d+(?:\.\d+)?", str(value))
        if not match:
            return None
        return float(match.group(0))


def mean(values):
    values = [value for value in values if value is not None]
    if not values:
        return ""
    return round(sum(values) / len(values), 6)


def fmt(value):
    if value == "" or value is None:
        return ""
    return ("%.6f" % float(value)).rstrip("0").rstrip(".")


def allele_union(rows, column):
    alleles = set()
    for row in rows:
        for item in row.get(column, "").split(";"):
            item = item.strip()
            if item:
                alleles.add(item)
    return sorted(alleles)


def load_epitope_versions():
    paths = sorted(VERSION_DIR.glob("[0-9][0-9]_vaccine_construct_V*.csv"))
    paths = [path for path in paths if "version_summary" not in path.name and "all_vaccine" not in path.name]
    versions = []
    for path in paths:
        rows = read_csv(path)
        if not rows:
            raise ValueError("%s contains no epitope rows" % path)

        if len(rows) > 27:
            raise ValueError(
                "%s has %d rows; maximum allowed is 27"
                % (path, len(rows))
            )
        version = rows[0]["construct_version"]
        versions.append(
            {
                "version": version,
                "version_name": rows[0]["version_name"],
                "version_description": rows[0]["version_description"],
                "source_file": str(path),
                "rows": sorted(rows, key=lambda row: int(row["construct_order"])),
            }
        )
    if len(versions) != 10:
        raise ValueError("Expected 10 epitope versions, found %d" % len(versions))
    return versions


def rows_by_class(rows):
    grouped = {epitope_class: [] for epitope_class in CLASS_ORDER}
    for row in rows:
        grouped[row["epitope_class"]].append(row)
    for epitope_class in CLASS_ORDER:
        grouped[epitope_class] = sorted(grouped[epitope_class], key=lambda row: int(row["construct_order"]))
        if not grouped[epitope_class]:
            raise ValueError(
                "No %s epitopes found in the construct version"
                % epitope_class
            )
    return grouped


def make_segment(segment_index, segment_type, segment_subtype, block_class, block_order, epitope_in_block, sequence, label, epitope_row=None):
    row = {
        "segment_index": segment_index,
        "segment_type": segment_type,
        "segment_subtype": segment_subtype,
        "block_class": block_class,
        "block_order": block_order,
        "epitope_in_block": epitope_in_block,
        "sequence": sequence,
        "segment_length": len(sequence),
        "label": label,
        "epitope_id": "",
        "epitope_class": "",
        "protein": "",
        "peptide": "",
        "source_start": "",
        "source_end": "",
        "vaxijen_score": "",
        "vaxijen_prediction": "",
        "allertop_prediction": "",
        "toxinpred_score": "",
        "toxinpred_prediction": "",
        "conservancy_percent": "",
        "mhc_i_immunogenicity_score": "",
        "hla_i_binding_allele_count": "",
        "hla_i_binding_alleles": "",
        "ifnepitope2_score": "",
        "ifnepitope2_prediction": "",
        "il4pred_score": "",
        "il4pred_prediction": "",
        "il10pred_score": "",
        "il10pred_prediction": "",
        "hla_ii_binding_allele_count": "",
        "hla_ii_binding_alleles": "",
        "human_blastp_full_length_exact_hit": "",
        "human_blastp_any_80pct_identity_80pct_query_span_hit": "",
    }
    if epitope_row:
        row.update(
            {
                "epitope_id": epitope_row.get("epitope_id", ""),
                "epitope_class": epitope_row.get("epitope_class", ""),
                "protein": epitope_row.get("protein", ""),
                "peptide": epitope_row.get("peptide", ""),
                "source_start": epitope_row.get("starts", ""),
                "source_end": epitope_row.get("ends", ""),
                "vaxijen_score": epitope_row.get("vaxijen_score", ""),
                "vaxijen_prediction": epitope_row.get("vaxijen_prediction", ""),
                "allertop_prediction": epitope_row.get("allertop_prediction", ""),
                "toxinpred_score": epitope_row.get("toxinpred_score", ""),
                "toxinpred_prediction": epitope_row.get("toxinpred_prediction", ""),
                "conservancy_percent": epitope_row.get("conservancy_percent", ""),
                "mhc_i_immunogenicity_score": epitope_row.get("mhc_i_immunogenicity_score", ""),
                "hla_i_binding_allele_count": epitope_row.get("hla_i_binding_allele_count", ""),
                "hla_i_binding_alleles": epitope_row.get("hla_i_binding_alleles", ""),
                "ifnepitope2_score": epitope_row.get("ifnepitope2_score", ""),
                "ifnepitope2_prediction": epitope_row.get("ifnepitope2_prediction", ""),
                "il4pred_score": epitope_row.get("il4pred_score", ""),
                "il4pred_prediction": epitope_row.get("il4pred_prediction", ""),
                "il10pred_score": epitope_row.get("il10pred_score", ""),
                "il10pred_prediction": epitope_row.get("il10pred_prediction", ""),
                "hla_ii_binding_allele_count": epitope_row.get("hla_ii_binding_allele_count", ""),
                "hla_ii_binding_alleles": epitope_row.get("hla_ii_binding_alleles", ""),
                "human_blastp_full_length_exact_hit": epitope_row.get("human_blastp_full_length_exact_hit", ""),
                "human_blastp_any_80pct_identity_80pct_query_span_hit": epitope_row.get(
                    "human_blastp_any_80pct_identity_80pct_query_span_hit", ""
                ),
            }
        )
    return row


def assemble_construct(version_info, adjuvant, architecture, construct_number):
    grouped = rows_by_class(version_info["rows"])
    segments = []
    segment_index = 1

    segments.append(
        make_segment(
            segment_index,
            "adjuvant",
            adjuvant["name"],
            "",
            "",
            "",
            adjuvant["sequence"],
            adjuvant["name"].replace("_", "-"),
        )
    )
    segment_index += 1
    segments.append(
        make_segment(
            segment_index,
            "linker",
            "EAAAK_adjuvant_linker",
            "",
            "",
            "",
            LINKERS["adjuvant_linker"],
            "EAAAK",
        )
    )
    segment_index += 1

    for block_number, epitope_class in enumerate(architecture["block_order"], start=1):
        if block_number > 1:
            segments.append(
                make_segment(
                    segment_index,
                    "linker",
                    "HEYGAEALERAG_block_linker",
                    "",
                    block_number - 1,
                    "",
                    LINKERS["block_linker"],
                    "HEYGAEALERAG",
                )
            )
            segment_index += 1
        internal_linker = architecture["internal_linkers"].get(epitope_class, "")
        for epitope_number, epitope_row in enumerate(grouped[epitope_class], start=1):
            if epitope_number > 1 and internal_linker:
                segments.append(
                    make_segment(
                        segment_index,
                        "linker",
                        "%s_internal_linker" % internal_linker,
                        epitope_class,
                        block_number,
                        "",
                        internal_linker,
                        internal_linker,
                    )
                )
                segment_index += 1
            label_prefix = {"CTL": "C", "HTL": "H", "Bcell": "B"}[epitope_class]
            segments.append(
                make_segment(
                    segment_index,
                    "epitope",
                    "%s_epitope" % epitope_class,
                    epitope_class,
                    block_number,
                    epitope_number,
                    epitope_row["peptide"],
                    "%s%d" % (label_prefix, epitope_number),
                    epitope_row,
                )
            )
            segment_index += 1

    sequence = "".join(segment["sequence"] for segment in segments)
    position = 1
    for segment in segments:
        segment["start_aa"] = position
        segment["end_aa"] = position + segment["segment_length"] - 1
        position = segment["end_aa"] + 1

    if not AA_PATTERN.match(sequence):
        raise ValueError("Construct contains non-standard amino-acid characters: %s" % sequence)
    if position - 1 != len(sequence):
        raise ValueError("Position accounting failed")

    construct_id = "VC%03d_%s_%s_%s" % (
        construct_number,
        version_info["version"],
        adjuvant["code"],
        architecture["id"],
    )
    for segment in segments:
        segment.update(
            {
                "construct_id": construct_id,
                "epitope_version": version_info["version"],
                "adjuvant_code": adjuvant["code"],
                "adjuvant_name": adjuvant["name"],
                "architecture_id": architecture["id"],
                "architecture_name": architecture["name"],
            }
        )

    metadata = build_metadata(construct_id, version_info, adjuvant, architecture, sequence, segments)
    return metadata, segments


def build_metadata(construct_id, version_info, adjuvant, architecture, sequence, segments):
    epitope_segments = [segment for segment in segments if segment["segment_type"] == "epitope"]
    ctl = [segment for segment in epitope_segments if segment["epitope_class"] == "CTL"]
    htl = [segment for segment in epitope_segments if segment["epitope_class"] == "HTL"]
    bcell = [segment for segment in epitope_segments if segment["epitope_class"] == "Bcell"]
    hla_i = allele_union(ctl, "hla_i_binding_alleles")
    hla_ii = allele_union(htl, "hla_ii_binding_alleles")

    def class_field(rows, field):
        return "|".join(row.get(field, "") for row in rows)

    linker_counts = {
        "EAAAK": sum(1 for segment in segments if segment["sequence"] == LINKERS["adjuvant_linker"]),
        "HEYGAEALERAG": sum(1 for segment in segments if segment["sequence"] == LINKERS["block_linker"]),
        "AAY": sum(1 for segment in segments if segment["sequence"] == LINKERS["CTL"]),
        "GPGPG": sum(1 for segment in segments if segment["sequence"] == LINKERS["HTL"]),
        "KK": sum(1 for segment in segments if segment["sequence"] == LINKERS["Bcell"]),
    }

    return {
        "construct_id": construct_id,
        "epitope_version": version_info["version"],
        "version_name": version_info["version_name"],
        "version_description": version_info["version_description"],
        "adjuvant_code": adjuvant["code"],
        "adjuvant_name": adjuvant["name"],
        "adjuvant_sequence": adjuvant["sequence"],
        "adjuvant_length": len(adjuvant["sequence"]),
        "architecture_id": architecture["id"],
        "architecture_name": architecture["name"],
        "architecture_formula": architecture["formula"],
        "block_order": ">".join(architecture["block_order"]),
        "internal_linker_CTL": architecture["internal_linkers"].get("CTL", ""),
        "internal_linker_HTL": architecture["internal_linkers"].get("HTL", ""),
        "internal_linker_Bcell": architecture["internal_linkers"].get("Bcell", ""),
        "sequence_length": len(sequence),
        "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        "vaccine_sequence": sequence,
        "ctl_epitope_count": len(ctl),
        "htl_epitope_count": len(htl),
        "bcell_epitope_count": len(bcell),
        "total_epitope_count": len(epitope_segments),
        "ctl_epitope_ids": class_field(ctl, "epitope_id"),
        "htl_epitope_ids": class_field(htl, "epitope_id"),
        "bcell_epitope_ids": class_field(bcell, "epitope_id"),
        "ctl_peptides": class_field(ctl, "peptide"),
        "htl_peptides": class_field(htl, "peptide"),
        "bcell_peptides": class_field(bcell, "peptide"),
        "ctl_proteins": class_field(ctl, "protein"),
        "htl_proteins": class_field(htl, "protein"),
        "bcell_proteins": class_field(bcell, "protein"),
        "unique_protein_count": len({segment["protein"] for segment in epitope_segments if segment["protein"]}),
        "unique_proteins": "|".join(sorted({segment["protein"] for segment in epitope_segments if segment["protein"]})),
        "mean_vaxijen_all_epitopes": fmt(mean([numeric(segment, "vaxijen_score") for segment in epitope_segments])),
        "mean_vaxijen_ctl": fmt(mean([numeric(segment, "vaxijen_score") for segment in ctl])),
        "mean_vaxijen_htl": fmt(mean([numeric(segment, "vaxijen_score") for segment in htl])),
        "mean_vaxijen_bcell": fmt(mean([numeric(segment, "vaxijen_score") for segment in bcell])),
        "mean_ctl_mhci_immunogenicity": fmt(mean([numeric(segment, "mhc_i_immunogenicity_score") for segment in ctl])),
        "sum_hla_i_binding_allele_counts": int(sum(numeric(segment, "hla_i_binding_allele_count") or 0 for segment in ctl)),
        "unique_hla_i_allele_count": len(hla_i),
        "unique_hla_i_alleles": "|".join(hla_i),
        "mean_htl_ifn_score": fmt(mean([numeric(segment, "ifnepitope2_score") for segment in htl])),
        "mean_htl_il4_score": fmt(mean([numeric(segment, "il4pred_score") for segment in htl])),
        "mean_htl_il10_score": fmt(mean([numeric(segment, "il10pred_score") for segment in htl])),
        "sum_hla_ii_binding_allele_counts": int(sum(numeric(segment, "hla_ii_binding_allele_count") or 0 for segment in htl)),
        "unique_hla_ii_allele_count": len(hla_ii),
        "unique_hla_ii_alleles": "|".join(hla_ii),
        "linker_EAAAK_count": linker_counts["EAAAK"],
        "linker_HEYGAEALERAG_count": linker_counts["HEYGAEALERAG"],
        "linker_AAY_count": linker_counts["AAY"],
        "linker_GPGPG_count": linker_counts["GPGPG"],
        "linker_KK_count": linker_counts["KK"],
        "full_length_exact_human_hit_epitope_count": sum(
            1 for segment in epitope_segments if segment.get("human_blastp_full_length_exact_hit") == "True"
        ),
    }


def load_font(size, bold=False, mono=False):
    candidates = []
    if mono:
        candidates.extend(
            [
                r"C:\Windows\Fonts\consolab.ttf" if bold else r"C:\Windows\Fonts\consola.ttf",
                r"C:\Windows\Fonts\courbd.ttf" if bold else r"C:\Windows\Fonts\cour.ttf",
            ]
        )
    candidates.extend(
        [
            r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
            r"C:\Windows\Fonts\calibrib.ttf" if bold else r"C:\Windows\Fonts\calibri.ttf",
        ]
    )
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def text_size(draw, text, font):
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0], bbox[3] - bbox[1]


def draw_centered(draw, box, text, font, fill):
    x0, y0, x1, y1 = box
    width, height = text_size(draw, text, font)
    if width > (x1 - x0 - 6):
        return
    draw.text((x0 + (x1 - x0 - width) / 2, y0 + (y1 - y0 - height) / 2), text, font=font, fill=fill)


def segment_color(segment):
    if segment["segment_type"] == "adjuvant":
        return COLORS["adjuvant"]
    if segment["segment_type"] == "epitope":
        return COLORS.get(segment["epitope_class"], "#94A3B8")
    sequence = segment["sequence"]
    return COLORS.get(sequence, "#94A3B8")


def layout_segments(segments, width, left_margin=64, right_margin=64):
    available = width - left_margin - right_margin
    total_length = sum(int(segment["segment_length"]) for segment in segments)
    scale = available / float(total_length)
    widths = []
    for segment in segments:
        base = int(segment["segment_length"]) * scale
        if segment["segment_type"] == "linker":
            minimum = 16 if segment["sequence"] in {"AAY", "KK", "GPGPG"} else 34
        elif segment["segment_type"] == "adjuvant":
            minimum = 70
        else:
            minimum = 20
        widths.append(max(base, minimum))
    total_width = sum(widths)
    if total_width > available:
        factor = available / total_width
        widths = [max(4, value * factor) for value in widths]
    coords = []
    x = left_margin
    for segment, segment_width in zip(segments, widths):
        x0 = x
        x1 = x + segment_width
        coords.append((segment, x0, x1))
        x = x1
    return coords


def block_ranges(coords):
    ranges = []
    current_class = None
    start_x = None
    end_x = None
    count = 0
    for segment, x0, x1 in coords:
        block_class = segment.get("block_class", "")
        if segment["segment_type"] == "epitope":
            if block_class != current_class:
                if current_class:
                    ranges.append((current_class, start_x, end_x, count))
                current_class = block_class
                start_x = x0
                count = 0
            end_x = x1
            count += 1
        elif current_class and segment["segment_type"] == "linker" and segment.get("block_class") == current_class:
            end_x = x1
    if current_class:
        ranges.append((current_class, start_x, end_x, count))
    return ranges


def render_master_png(metadata, segments, out_path):
    width, height = 1800, 720
    image = Image.new("RGB", (width, height), COLORS["background"])
    draw = ImageDraw.Draw(image)
    title_font = load_font(34, bold=True)
    subtitle_font = load_font(20)
    label_font = load_font(17, bold=True)
    small_font = load_font(13)
    tiny_font = load_font(11, bold=True)
    mono_font = load_font(13, mono=True)

    display_adjuvant = metadata["adjuvant_name"].replace("_", "-")
    display_architecture = metadata["architecture_name"].replace("Bcell", "B-cell")

    draw.rectangle((0, 0, width, 88), fill="#111827")
    draw.text((64, 24), "%s | %s | %s" % (metadata["architecture_id"], display_adjuvant, display_architecture), font=title_font, fill="#FFFFFF")
    draw.text((64, 96), "Representative construct: %s  |  Epitope set: %s  |  Length: %s aa" % (metadata["construct_id"], metadata["epitope_version"], metadata["sequence_length"]), font=subtitle_font, fill=COLORS["ink"])

    map_y = 250
    map_h = 112
    coords = layout_segments(segments, width)

    for block_class, x0, x1, count in block_ranges(coords):
        draw.line((x0, map_y - 38, x1, map_y - 38), fill=COLORS[block_class], width=4)
        draw.line((x0, map_y - 38, x0, map_y - 26), fill=COLORS[block_class], width=4)
        draw.line((x1, map_y - 38, x1, map_y - 26), fill=COLORS[block_class], width=4)
        label = "%s block (%d epitopes)" % (block_class.replace("Bcell", "B-cell"), count)
        tw, th = text_size(draw, label, label_font)
        draw.rectangle((x0 + (x1 - x0 - tw) / 2 - 12, map_y - 76, x0 + (x1 - x0 + tw) / 2 + 12, map_y - 48), fill="#FFFFFF", outline=COLORS[block_class], width=2)
        draw.text((x0 + (x1 - x0 - tw) / 2, map_y - 72), label, font=label_font, fill=COLORS["ink"])

    for segment, x0, x1 in coords:
        color = segment_color(segment)
        draw.rounded_rectangle((x0, map_y, x1, map_y + map_h), radius=5, fill=color, outline="#FFFFFF", width=2)
        label = segment["label"]
        if segment["segment_type"] == "linker" and segment["sequence"] in {"AAY", "KK", "GPGPG"}:
            label = segment["sequence"]
        fill = "#FFFFFF" if segment["segment_type"] != "epitope" or segment["epitope_class"] == "CTL" else "#172033"
        draw_centered(draw, (x0, map_y, x1, map_y + map_h), label, tiny_font, fill)

    y = 410
    legend = [
        ("Adjuvant", COLORS["adjuvant"]),
        ("EAAAK", COLORS["EAAAK"]),
        ("HEYGAEALERAG", COLORS["HEYGAEALERAG"]),
        ("CTL epitope", COLORS["CTL"]),
        ("HTL epitope", COLORS["HTL"]),
        ("B-cell epitope", COLORS["Bcell"]),
        ("AAY", COLORS["AAY"]),
        ("GPGPG", COLORS["GPGPG"]),
        ("KK", COLORS["KK"]),
    ]
    x = 64
    for label, color in legend:
        draw.rounded_rectangle((x, y, x + 26, y + 18), radius=3, fill=color)
        draw.text((x + 34, y - 1), label, font=small_font, fill=COLORS["ink"])
        x += 160 if len(label) < 8 else 210
        if x > width - 220:
            x = 64
            y += 30

    summary_lines = [
        "Formula: %s" % metadata["architecture_formula"],
        "Mean VaxiJen: all %s | CTL %s | HTL %s | B-cell %s"
        % (
            metadata["mean_vaxijen_all_epitopes"],
            metadata["mean_vaxijen_ctl"],
            metadata["mean_vaxijen_htl"],
            metadata["mean_vaxijen_bcell"],
        ),
        "HLA breadth: %s unique HLA-I alleles, %s unique HLA-II alleles"
        % (metadata["unique_hla_i_allele_count"], metadata["unique_hla_ii_allele_count"]),
        "Linkers: EAAAK x%s, HEYGAEALERAG x%s, AAY x%s, GPGPG x%s, KK x%s"
        % (
            metadata["linker_EAAAK_count"],
            metadata["linker_HEYGAEALERAG_count"],
            metadata["linker_AAY_count"],
            metadata["linker_GPGPG_count"],
            metadata["linker_KK_count"],
        ),
    ]
    y = 505
    for line in summary_lines:
        draw.text((64, y), line, font=subtitle_font if line.startswith("Formula") else small_font, fill=COLORS["ink"])
        y += 34 if line.startswith("Formula") else 24

    sequence_preview = metadata["vaccine_sequence"][:145] + ("..." if len(metadata["vaccine_sequence"]) > 145 else "")
    draw.rectangle((64, 642, width - 64, 687), fill="#FFFFFF", outline=COLORS["line"], width=1)
    draw.text((80, 655), sequence_preview, font=mono_font, fill="#111827")

    image.save(out_path)


def svg_text(x, y, text, size=14, weight="400", fill="#172033", anchor="start", family="Arial"):
    return '<text x="%s" y="%s" font-family="%s" font-size="%s" font-weight="%s" fill="%s" text-anchor="%s">%s</text>' % (
        round(x, 2),
        round(y, 2),
        family,
        size,
        weight,
        fill,
        anchor,
        html.escape(str(text)),
    )


def render_master_svg(metadata, segments, out_path):
    width, height = 1800, 720
    map_y = 250
    map_h = 112
    coords = layout_segments(segments, width)
    display_adjuvant = metadata["adjuvant_name"].replace("_", "-")
    display_architecture = metadata["architecture_name"].replace("Bcell", "B-cell")
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (width, height, width, height),
        '<rect width="100%%" height="100%%" fill="%s"/>' % COLORS["background"],
        '<rect x="0" y="0" width="%d" height="88" fill="#111827"/>' % width,
        svg_text(64, 54, "%s | %s | %s" % (metadata["architecture_id"], display_adjuvant, display_architecture), 34, "700", "#FFFFFF"),
        svg_text(64, 118, "Representative construct: %s | Epitope set: %s | Length: %s aa" % (metadata["construct_id"], metadata["epitope_version"], metadata["sequence_length"]), 20),
    ]
    for block_class, x0, x1, count in block_ranges(coords):
        color = COLORS[block_class]
        lines.append('<line x1="%.2f" y1="%d" x2="%.2f" y2="%d" stroke="%s" stroke-width="4"/>' % (x0, map_y - 38, x1, map_y - 38, color))
        lines.append('<line x1="%.2f" y1="%d" x2="%.2f" y2="%d" stroke="%s" stroke-width="4"/>' % (x0, map_y - 38, x0, map_y - 26, color))
        lines.append('<line x1="%.2f" y1="%d" x2="%.2f" y2="%d" stroke="%s" stroke-width="4"/>' % (x1, map_y - 38, x1, map_y - 26, color))
        label = "%s block (%d epitopes)" % (block_class.replace("Bcell", "B-cell"), count)
        cx = x0 + (x1 - x0) / 2
        label_w = max(180, len(label) * 10)
        lines.append('<rect x="%.2f" y="%d" width="%d" height="30" rx="5" fill="#FFFFFF" stroke="%s" stroke-width="2"/>' % (cx - label_w / 2, map_y - 78, label_w, color))
        lines.append(svg_text(cx, map_y - 57, label, 17, "700", COLORS["ink"], "middle"))

    for segment, x0, x1 in coords:
        color = segment_color(segment)
        label = segment["label"]
        if segment["segment_type"] == "linker" and segment["sequence"] in {"AAY", "KK", "GPGPG"}:
            label = segment["sequence"]
        lines.append('<rect x="%.2f" y="%d" width="%.2f" height="%d" rx="5" fill="%s" stroke="#FFFFFF" stroke-width="2"/>' % (x0, map_y, x1 - x0, map_h, color))
        if x1 - x0 > 26:
            fill = "#FFFFFF" if segment["segment_type"] != "epitope" or segment["epitope_class"] == "CTL" else "#172033"
            lines.append(svg_text(x0 + (x1 - x0) / 2, map_y + 62, label, 11, "700", fill, "middle"))

    y = 410
    x = 64
    legend = [
        ("Adjuvant", COLORS["adjuvant"]),
        ("EAAAK", COLORS["EAAAK"]),
        ("HEYGAEALERAG", COLORS["HEYGAEALERAG"]),
        ("CTL epitope", COLORS["CTL"]),
        ("HTL epitope", COLORS["HTL"]),
        ("B-cell epitope", COLORS["Bcell"]),
        ("AAY", COLORS["AAY"]),
        ("GPGPG", COLORS["GPGPG"]),
        ("KK", COLORS["KK"]),
    ]
    for label, color in legend:
        lines.append('<rect x="%d" y="%d" width="26" height="18" rx="3" fill="%s"/>' % (x, y, color))
        lines.append(svg_text(x + 34, y + 14, label, 13))
        x += 160 if len(label) < 8 else 210
        if x > width - 220:
            x = 64
            y += 30
    summary_lines = [
        "Formula: %s" % metadata["architecture_formula"],
        "Mean VaxiJen: all %s | CTL %s | HTL %s | B-cell %s"
        % (
            metadata["mean_vaxijen_all_epitopes"],
            metadata["mean_vaxijen_ctl"],
            metadata["mean_vaxijen_htl"],
            metadata["mean_vaxijen_bcell"],
        ),
        "HLA breadth: %s unique HLA-I alleles, %s unique HLA-II alleles"
        % (metadata["unique_hla_i_allele_count"], metadata["unique_hla_ii_allele_count"]),
        "Linkers: EAAAK x%s, HEYGAEALERAG x%s, AAY x%s, GPGPG x%s, KK x%s"
        % (
            metadata["linker_EAAAK_count"],
            metadata["linker_HEYGAEALERAG_count"],
            metadata["linker_AAY_count"],
            metadata["linker_GPGPG_count"],
            metadata["linker_KK_count"],
        ),
    ]
    y = 505
    for index, line in enumerate(summary_lines):
        lines.append(svg_text(64, y, line, 20 if index == 0 else 13, "400"))
        y += 34 if index == 0 else 24
    sequence_preview = metadata["vaccine_sequence"][:145] + ("..." if len(metadata["vaccine_sequence"]) > 145 else "")
    lines.append('<rect x="64" y="642" width="%d" height="45" fill="#FFFFFF" stroke="%s"/>' % (width - 128, COLORS["line"]))
    lines.append(svg_text(80, 672, sequence_preview, 13, "400", "#111827", "start", "Consolas"))
    lines.append("</svg>")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def architecture_summary_rows(construct_rows):
    grouped = {}
    for row in construct_rows:
        key = (row["adjuvant_code"], row["architecture_id"])
        grouped.setdefault(key, []).append(row)
    rows = []
    for (adjuvant_code, architecture_id), items in sorted(grouped.items()):
        first = items[0]
        rows.append(
            {
                "adjuvant_code": adjuvant_code,
                "adjuvant_name": first["adjuvant_name"],
                "architecture_id": architecture_id,
                "architecture_name": first["architecture_name"],
                "architecture_formula": first["architecture_formula"],
                "construct_count": len(items),
                "min_sequence_length": min(int(item["sequence_length"]) for item in items),
                "max_sequence_length": max(int(item["sequence_length"]) for item in items),
                "mean_sequence_length": fmt(mean([numeric(item, "sequence_length") for item in items])),
                "linker_EAAAK_count": first["linker_EAAAK_count"],
                "linker_HEYGAEALERAG_count": first["linker_HEYGAEALERAG_count"],
                "linker_AAY_count": first["linker_AAY_count"],
                "linker_GPGPG_count": first["linker_GPGPG_count"],
                "linker_KK_count": first["linker_KK_count"],
            }
        )
    return rows


def write_architecture_markdown(rows):
    path = OUTPUT_DIR / "04_vaccine_construct_architecture_summary.md"
    lines = [
        "# Vaccine Construct Architecture Summary",
        "",
        "Generated: %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "",
        "| Adjuvant | Architecture | Constructs | Length range | Formula |",
        "|---|---|---:|---|---|",
    ]
    for row in rows:
        lines.append(
            "| %s | %s | %s | %s-%s aa | %s |"
            % (
                row["adjuvant_name"],
                row["architecture_id"],
                row["construct_count"],
                row["min_sequence_length"],
                row["max_sequence_length"],
                row["architecture_formula"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(construct_rows, figure_rows):
    path = REPORT_DIR / "22_vaccine_construct_design_summary.md"
    lengths = [int(row["sequence_length"]) for row in construct_rows]
    lines = [
        "# Vaccine Construct Design Summary",
        "",
        "Generated: %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "",
        "A total of 60 vaccine constructs were generated from 10 epitope-set versions, 2 adjuvants, and 3 architecture layouts.",
        "",
        "Adjuvants:",
        "",
        "- Beta_defensin: `%s`" % ADJUVANTS[0]["sequence"],
        "- RS09: `%s`" % ADJUVANTS[1]["sequence"],
        "",
        "Core linkers:",
        "",
        "- Adjuvant linker: `EAAAK`",
        "- Block linker: `HEYGAEALERAG`",
        "- CTL internal linker: `AAY`",
        "- HTL internal linker: `GPGPG`",
        "- B-cell internal linker: `KK`",
        "",
        "Output summary:",
        "",
        "- Construct count: `%d`" % len(construct_rows),
        "- Sequence length range: `%d-%d aa`" % (min(lengths), max(lengths)),
        "- Epitope counts are determined dynamically from each strain-aware epitope version.",
        "- No C_Protein HTL rescue exception was added; strictly adheres to baseline epitope filtering thresholds.",
        "",
        "Main output files:",
        "",
        "- `08_vaccine_construct_designs/01_60_vaccine_constructs_metadata.csv`",
        "- `08_vaccine_construct_designs/02_60_vaccine_constructs_sequences.fasta`",
        "- `08_vaccine_construct_designs/03_60_vaccine_construct_segment_map.csv`",
        "- `08_vaccine_construct_designs/04_vaccine_construct_architecture_summary.md`",
        "- `08_vaccine_construct_designs/05_master_visual_design_index.csv`",
        "",
        "Master figure files:",
        "",
    ]
    for row in figure_rows:
        lines.append("- `%s`" % Path(row["png_file"]).relative_to(ROOT))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    MASTER_FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    versions = load_epitope_versions()
    construct_rows = []
    all_segments = []
    construct_segments = {}
    construct_number = 1

    for version in versions:
        for adjuvant in ADJUVANTS:
            for architecture in ARCHITECTURES:
                metadata, segments = assemble_construct(version, adjuvant, architecture, construct_number)
                construct_rows.append(metadata)
                all_segments.extend(segments)
                construct_segments[metadata["construct_id"]] = segments
                construct_number += 1

    construct_columns = list(construct_rows[0].keys())
    segment_columns = [
        "construct_id",
        "epitope_version",
        "adjuvant_code",
        "adjuvant_name",
        "architecture_id",
        "architecture_name",
        "segment_index",
        "segment_type",
        "segment_subtype",
        "block_class",
        "block_order",
        "epitope_in_block",
        "start_aa",
        "end_aa",
        "segment_length",
        "label",
        "sequence",
        "epitope_id",
        "epitope_class",
        "protein",
        "peptide",
        "source_start",
        "source_end",
        "vaxijen_score",
        "vaxijen_prediction",
        "allertop_prediction",
        "toxinpred_score",
        "toxinpred_prediction",
        "conservancy_percent",
        "mhc_i_immunogenicity_score",
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
        "human_blastp_full_length_exact_hit",
        "human_blastp_any_80pct_identity_80pct_query_span_hit",
    ]

    write_table(OUTPUT_DIR / "01_60_vaccine_constructs_metadata.csv", construct_rows, construct_columns, ",")
    write_table(OUTPUT_DIR / "01_60_vaccine_constructs_metadata.tsv", construct_rows, construct_columns, "\t")
    write_json(OUTPUT_DIR / "01_60_vaccine_constructs_metadata.json", construct_rows)

    fasta_lines = []
    txt_lines = []
    for row in construct_rows:
        header = "%s|%s|%s|%s|length=%s" % (
            row["construct_id"],
            row["epitope_version"],
            row["adjuvant_name"],
            row["architecture_id"],
            row["sequence_length"],
        )
        fasta_lines.append(">" + header)
        fasta_lines.append(wrap_sequence(row["vaccine_sequence"]))
        txt_lines.append(header)
        txt_lines.append(row["vaccine_sequence"])
        txt_lines.append("")
    (OUTPUT_DIR / "02_60_vaccine_constructs_sequences.fasta").write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    (OUTPUT_DIR / "02_60_vaccine_constructs_sequences.txt").write_text("\n".join(txt_lines), encoding="utf-8")

    write_table(OUTPUT_DIR / "03_60_vaccine_construct_segment_map.csv", all_segments, segment_columns, ",")
    write_table(OUTPUT_DIR / "03_60_vaccine_construct_segment_map.tsv", all_segments, segment_columns, "\t")
    write_json(OUTPUT_DIR / "03_60_vaccine_construct_segment_map.json", all_segments)

    arch_rows = architecture_summary_rows(construct_rows)
    arch_columns = list(arch_rows[0].keys())
    write_table(OUTPUT_DIR / "04_vaccine_construct_architecture_summary.csv", arch_rows, arch_columns, ",")
    write_table(OUTPUT_DIR / "04_vaccine_construct_architecture_summary.tsv", arch_rows, arch_columns, "\t")
    write_json(OUTPUT_DIR / "04_vaccine_construct_architecture_summary.json", arch_rows)
    write_architecture_markdown(arch_rows)

    figure_rows = []
    representative_rows = [row for row in construct_rows if row["epitope_version"] == "V1"]
    for row in representative_rows:
        basename = "%s_%s_%s_master_map" % (row["adjuvant_code"], row["architecture_id"], row["architecture_name"].replace(" ", "_").replace("-", "_"))
        basename = re.sub(r"[^A-Za-z0-9_]+", "_", basename).strip("_")
        png_path = MASTER_FIGURE_DIR / ("%s.png" % basename)
        svg_path = MASTER_FIGURE_DIR / ("%s.svg" % basename)
        render_master_png(row, construct_segments[row["construct_id"]], png_path)
        render_master_svg(row, construct_segments[row["construct_id"]], svg_path)
        figure_rows.append(
            {
                "adjuvant_code": row["adjuvant_code"],
                "adjuvant_name": row["adjuvant_name"],
                "architecture_id": row["architecture_id"],
                "architecture_name": row["architecture_name"],
                "representative_construct_id": row["construct_id"],
                "png_file": str(png_path),
                "svg_file": str(svg_path),
            }
        )

    figure_columns = list(figure_rows[0].keys())
    write_table(OUTPUT_DIR / "05_master_visual_design_index.csv", figure_rows, figure_columns, ",")
    write_table(OUTPUT_DIR / "05_master_visual_design_index.tsv", figure_rows, figure_columns, "\t")
    write_json(OUTPUT_DIR / "05_master_visual_design_index.json", figure_rows)
    write_report(construct_rows, figure_rows)

    print("Generated %d vaccine constructs" % len(construct_rows))
    print("Generated %d segment rows" % len(all_segments))
    print("Generated %d master architecture maps" % len(figure_rows))
    print("Output folder: %s" % OUTPUT_DIR)


if __name__ == "__main__":
    main()
