#!/usr/bin/env python3
"""
===============================================================================
Pan-ebolavirus Vaccine Design Study | Automated NCBI Entrez Sequence Downloader
Species: Zaire ebolavirus (EBOV)
NCBI Taxonomy ID: txid186538
Reference Assembly: NC_002549.1 (AF086833.2)
Expected Complete Protein Dataset: 19,286 sequences
===============================================================================
Description:
    This script automates the retrieval of full-length protein sequences for all
    seven structural and non-structural proteins of Zaire ebolavirus
    from the NCBI Entrez Protein database in accordance with the study methodology
    and Supplementary File 1.

Requirements:
    - Python >= 3.8
    - biopython >= 1.80

Usage:
    python Zaire_sequence_downloader.py
"""

import os
import time
from Bio import Entrez

# =============================================================================
# 1. NCBI API CONFIGURATION
# =============================================================================
# NCBI requires a contact email address to track automated query traffic.
Entrez.email = "your.email@example.com"

# Optional: Add your NCBI personal API key to elevate request limits to 10 req/sec
# Obtain a key at: https://www.ncbi.nlm.nih.gov/account/settings/
Entrez.api_key = None

# =============================================================================
# 2. DEFINITIVE SPECIES-SPECIFIC PROTEIN QUERIES
# =============================================================================
# Exact queries as reported in Supplementary File 1 and Manuscript Table 1
protein_tasks = {
    "1. Nucleoprotein": (
        'txid186538[Organism] AND ("Nucleoprotein"[Protein Name] OR "NP"[Protein Name] OR '
        '"nucleocapsid protein"[Protein Name]) NOT (partial[Title] OR fragment[Title] OR '
        'chain[Title] OR mutant[Title] OR synthetic[Title] OR construct[Title]) AND (700:800[SLEN])'
    ),

    "2. Polymerase cofactor VP35": (
        'txid186538[Organism] AND ("VP35"[Protein Name] OR "polymerase cofactor"[Protein Name] OR '
        '"polymerase complex protein"[Protein Name]) NOT (partial[Title] OR fragment[Title] OR '
        'chain[Title] OR mutant[Title] OR synthetic[Title] OR construct[Title]) AND (300:400[SLEN])'
    ),

    "3. Matrix protein VP40": (
        'txid186538[Organism] AND ("VP40"[Protein Name] OR "matrix protein"[Protein Name]) '
        'NOT (partial[Title] OR fragment[Title] OR chain[Title] OR mutant[Title] OR '
        'synthetic[Title] OR construct[Title]) AND (290:360[SLEN])'
    ),

    "4. Glycoprotein GP": (
        'txid186538[Organism] AND ("glycoprotein"[Title] OR "GP"[Title] OR "spike protein"[Title] OR '
        '"envelope glycoprotein"[Title]) NOT (partial[Title] OR fragment[Title] OR '
        'mutant[Title] OR chain[Title] OR synthetic[Title] OR construct[Title]) AND (250:750[SLEN])'
    ),

    "5. Minor nucleoprotein VP30": (
        'txid186538[Organism] AND ("VP30"[Protein Name] OR "transcription factor"[Protein Name] OR '
        '"minor nucleoprotein"[Protein Name] OR "transcription activator"[Protein Name]) '
        'NOT (partial[Title] OR fragment[Title] OR chain[Title] OR mutant[Title] OR '
        'synthetic[Title] OR construct[Title]) AND (250:330[SLEN])'
    ),

    "6. Membrane-associated protein VP24": (
        'txid186538[Organism] AND ("VP24"[Protein Name] OR "minor matrix protein"[Protein Name] OR '
        '"membrane-associated protein"[Protein Name]) NOT (partial[Title] OR fragment[Title] OR '
        'chain[Title] OR mutant[Title] OR synthetic[Title] OR construct[Title]) AND (220:280[SLEN])'
    ),

    "7. RNA polymerase L": (
        'txid186538[Organism] AND ("polymerase"[Title] OR "L protein"[Title] OR "large protein"[Title] OR '
        '"RNA-directed RNA polymerase"[Title] OR "RNA-dependent RNA polymerase"[Title]) '
        'NOT (partial[Title] OR fragment[Title] OR chain[Title] OR mutant[Title] OR synthetic[Title] OR '
        'construct[Title]) AND 1900:2400[SLEN]'
    )
}

# =============================================================================
# 3. DOWNLOAD ENGINE EXECUTION
# =============================================================================
def run_downloader():
    print("=" * 70)
    print("Pan-ebolavirus Sequence Retrieval Engine")
    print(f"Target Species: Zaire ebolavirus (EBOV, txid186538)")
    print("=" * 70)
    
    download_report = {}
    total_downloaded = 0

    for folder_name, query_string in protein_tasks.items():
        print(f"\n[Processing Target] {folder_name}")
        os.makedirs(folder_name, exist_ok=True)

        # File naming convention
        protein_filename = folder_name.split(". ", 1)[-1].replace(" ", "_")
        query_out_path = os.path.join(folder_name, f"{protein_filename}_query.txt")
        fasta_out_path = os.path.join(folder_name, f"{protein_filename}.fasta")

        # Save query text for record keeping
        with open(query_out_path, "w", encoding="utf-8") as q_file:
            q_file.write(query_string)

        try:
            # 1. ESearch: query NCBI Entrez with History Server enabled
            search_handle = Entrez.esearch(db="protein", term=query_string, usehistory="y")
            search_results = Entrez.read(search_handle)
            search_handle.close()

            count = int(search_results["Count"])
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]

            print(f"  -> Matching NCBI records: {count}")
            download_report[folder_name] = count

            if count == 0:
                print("  -> Skipping download (0 records found).")
                continue

            total_downloaded += count

            # 2. EFetch: batch retrieve sequences in fasta format
            batch_size = 200
            with open(fasta_out_path, "w", encoding="utf-8") as fasta_file:
                for start in range(0, count, batch_size):
                    end_idx = min(start + batch_size, count)
                    print(f"  -> Streaming records {start + 1} to {end_idx} of {count}...")

                    fetch_handle = Entrez.efetch(
                        db="protein",
                        rettype="fasta",
                        retmode="text",
                        retstart=start,
                        retmax=batch_size,
                        webenv=webenv,
                        query_key=query_key
                    )
                    data = fetch_handle.read()
                    fetch_handle.close()

                    fasta_file.write(data)
                    # Respect NCBI Entrez rate limits (max 3 req/sec without API key)
                    delay = 0.15 if Entrez.api_key else 0.40
                    time.sleep(delay)

            print(f"  -> Successfully saved: {fasta_out_path}")

        except Exception as err:
            print(f"  [ERROR] Failed to download {folder_name}: {err}")
            download_report[folder_name] = f"FAILED: {err}"

    # =========================================================================
    # 4. SUMMARY REPORT GENERATION
    # =========================================================================
    report_file = "download_report.txt"
    print("\n" + "=" * 70)
    print("NCBI ENTREZ PROTEIN DOWNLOAD SUMMARY REPORT")
    print(f"Species: Zaire ebolavirus (EBOV)")
    print("=" * 70)
    
    with open(report_file, "w", encoding="utf-8") as rep:
        rep.write("=" * 60 + "\n")
        rep.write(f"NCBI ENTREZ PROTEIN DOWNLOAD REPORT: Zaire ebolavirus (EBOV)\n")
        rep.write(f"Taxonomy ID: txid186538\n")
        rep.write("=" * 60 + "\n\n")
        for target, res in download_report.items():
            line = f"{target}: {res} sequences successfully written"
            print(f"  {line}")
            rep.write(line + "\n")
        rep.write(f"\nTotal sequences downloaded across 7 proteins: {total_downloaded}\n")
    
    print(f"\nReport written to {report_file}.")
    print(f"Total sequences retrieved: {total_downloaded} (Expected: 19,286)")

if __name__ == "__main__":
    run_downloader()
