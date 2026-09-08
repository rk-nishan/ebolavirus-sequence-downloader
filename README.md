# Pan-ebolavirus Vaccine Design: Automated NCBI Sequence Retrieval Scripts

This repository contains the official, reproducible Python sequence retrieval pipeline supporting the study:

> **"Computationally Designed Pan-Ebolavirus Multi-Epitope Vaccine Candidates Targeting Conserved Structural and Non-Structural Viral Proteomes"**

## Overview

The scripts automate the complete retrieval of structural and non-structural protein sequences from the **NCBI Entrez Protein database** for all four human-pathogenic species of the genus *Orthoebolavirus*:

| Script Name | Target Species | NCBI Taxonomy ID | Reference Genome | Target Proteins | Retrieved Count |
| :--- | :--- | :---: | :---: | :---: | :---: |
| `Bundibugyo_sequence_downloader.py` | *Bundibugyo ebolavirus* (BDBV) | `txid565995` | `NC_014373.1` | 7 Proteins | 269 |
| `Sudan_sequence_downloader.py` | *Sudan ebolavirus* (SUDV) | `txid186540` | `NC_006432.1` | 7 Proteins | 1,132 |
| `Tai_forest_sequence_downloader.py` | *Taï Forest ebolavirus* (TAFV) | `txid186541` | `NC_014372.1` | 7 Proteins | 36 |
| `Zaire_sequence_downloader.py` | *Zaire ebolavirus* (EBOV) | `txid186538` | `NC_002549.1` | 7 Proteins | 19,286 |
| **Combined Proteomic Set** | **All 4 Species** | - | - | **28 Targets** | **20,723** |

### Target Viral Proteins
For each species, the seven canonical viral targets are systematically retrieved:
1. **Nucleoprotein (NP)**: Length window `700:800` aa
2. **Polymerase cofactor (VP35)**: Length window `300:400` aa
3. **Matrix protein (VP40)**: Length window `290:360` aa
4. **Surface glycoprotein (GP)**: Length window `250:750` aa
5. **Minor nucleoprotein (VP30)**: Length window `250:330` aa
6. **Membrane-associated protein (VP24)**: Length window `220:280` aa
7. **RNA-dependent RNA polymerase (L)**: Length window `1900:2400` aa

### Quality-Filtering at Query Level
All queries enforce the universal negative-selection quality exclusion:
```sql
NOT (partial[Title] OR fragment[Title] OR chain[Title] OR mutant[Title] OR synthetic[Title] OR construct[Title])
```
to exclude incomplete fragments, laboratory mutants, crystallographic chains, and artificial constructs.

---

## Installation & Requirements

### Dependencies
- Python >= 3.8
- Biopython >= 1.80

Install dependencies via:
```bash
pip install -r requirements.txt
```

---

## Usage Instructions

1. **Configure your email** (required by NCBI Entrez API):
   Open the respective script and set your email address:
   ```python
   Entrez.email = "your.email@institution.edu"
   ```
   *(Optional)* If you have an NCBI personal API key, set `Entrez.api_key = "your_key"` to enable up to 10 requests/second.

2. **Execute retrieval for a target species:**
   ```bash
   # For Bundibugyo ebolavirus
   python Bundibugyo_sequence_downloader.py

   # For Sudan ebolavirus
   python Sudan_sequence_downloader.py

   # For Taï Forest ebolavirus
   python Tai_forest_sequence_downloader.py

   # For Zaire ebolavirus
   python Zaire_sequence_downloader.py
   ```

3. **Output Structure:**
   Each script generates dedicated subdirectories for each protein containing:
   - `<Protein_Name>.fasta`: Full FASTA file containing all retrieved sequences.
   - `<Protein_Name>_query.txt`: Verbatim text of the Entrez search query executed.
   - `download_report.txt`: Summary report of record counts written.

---

## Reproducibility & Data Integrity
All queries, parameter windows, and canonical yields in these scripts match **Supplementary File 1**, **Supplementary Table S1**, and **Table 1** of the manuscript with 100% precision.

## License
MIT License. Openly available for research and verification under FAIR data principles.
