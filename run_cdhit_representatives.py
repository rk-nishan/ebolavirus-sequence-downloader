#!/usr/bin/env python3
"""Run CD-HIT at 99% identity with 99% short/long coverage."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_DIR = SCRIPT_DIR.parent
# Default directories: use project path if present, otherwise default to local directories
PROJECT_INPUT = Path(r"F:\3.Project_Drive\11. Vaccine\000. Version 1.0\2. Cleaned Protein Fasta\cleaned_fasta")
PROJECT_OUTPUT = Path(r"F:\3.Project_Drive\11. Vaccine\000. Version 1.0\4. CD-HIT")
INPUT_DIR = PROJECT_INPUT if PROJECT_INPUT.exists() else SCRIPT_DIR / "cleaned_fasta"
OUTPUT_DIR = PROJECT_OUTPUT if PROJECT_OUTPUT.exists() else SCRIPT_DIR / "cdhit_output"
REPORTS_DIR = OUTPUT_DIR
LOGS_DIR = OUTPUT_DIR


def numeric_prefix(path: Path) -> int:
    match = re.match(r"^(\d+)\.", path.parent.name)
    return int(match.group(1)) if match else 10_000


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")


def count_fasta_records(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def count_clusters(path: Path) -> int:
    if not path.exists():
        return 0
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">Cluster"):
                count += 1
    return count


def find_cdhit_executable(user_path: str | None) -> str:
    if user_path:
        path = Path(user_path)
        if not path.exists():
            raise FileNotFoundError(f"CD-HIT executable not found: {path}")
        return str(path)

    found = shutil.which("cd-hit") or shutil.which("cd-hit.exe")
    if not found:
        raise FileNotFoundError(
            "CD-HIT executable was not found on PATH. Install CD-HIT or pass --cd-hit C:\\path\\to\\cd-hit.exe"
        )
    return found


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run CD-HIT representatives for each cleaned protein FASTA.")
    parser.add_argument("--cd-hit", default=None, help="Path to cd-hit executable.")
    parser.add_argument("--input-dir", type=Path, default=INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--identity", default="0.99")
    parser.add_argument("--short-coverage", default="0.99")
    parser.add_argument("--long-coverage", default="0.99")
    parser.add_argument("--word-size", default="5")
    parser.add_argument("--threads", default="0")
    parser.add_argument("--memory", default="0")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cdhit = find_cdhit_executable(args.cd_hit)

    fasta_files = sorted(args.input_dir.glob("*/*.fasta"), key=lambda path: (numeric_prefix(path), str(path)))
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found under {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []
    for input_fasta in fasta_files:
        protein_folder = input_fasta.parent.name
        output_subdir = args.output_dir / protein_folder
        output_subdir.mkdir(parents=True, exist_ok=True)
        output_fasta = output_subdir / f"{safe_name(protein_folder)}_{input_fasta.stem}_cdhit099.fasta"
        log_path = LOGS_DIR / f"{safe_name(protein_folder)}_cdhit099.log"

        command = [
            cdhit,
            "-i",
            str(input_fasta),
            "-o",
            str(output_fasta),
            "-c",
            args.identity,
            "-n",
            args.word_size,
            "-G",
            "1",
            "-aS",
            args.short_coverage,
            "-aL",
            args.long_coverage,
            "-g",
            "1",
            "-d",
            "0",
            "-T",
            args.threads,
            "-M",
            args.memory,
        ]

        status = "pending"
        if output_fasta.exists() and not args.overwrite:
            status = "skipped_existing"
        elif args.dry_run:
            status = "dry_run"
            print(" ".join(command))
        else:
            with log_path.open("w", encoding="utf-8", errors="replace") as log_handle:
                completed = subprocess.run(
                    command,
                    check=False,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
            status = "completed" if completed.returncode == 0 else f"failed_return_{completed.returncode}"
            if completed.returncode != 0:
                raise RuntimeError(f"CD-HIT failed for {input_fasta}. See {log_path}")

        summary_rows.append(
            {
                "protein_folder": protein_folder,
                "input_fasta": str(input_fasta),
                "output_fasta": str(output_fasta),
                "cluster_file": str(output_fasta) + ".clstr",
                "input_sequences": count_fasta_records(input_fasta),
                "representative_sequences": count_fasta_records(output_fasta) if output_fasta.exists() else "",
                "clusters": count_clusters(Path(str(output_fasta) + ".clstr")),
                "identity": args.identity,
                "short_coverage_aS": args.short_coverage,
                "long_coverage_aL": args.long_coverage,
                "global_identity_G": "1",
                "accurate_mode_g": "1",
                "word_size_n": args.word_size,
                "status": status,
            }
        )

    report_path = REPORTS_DIR / "cdhit099_summary.csv"
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()

