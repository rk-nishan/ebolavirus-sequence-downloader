# Orthoebolavirus Proteome Curation & Multi-Epitope Vaccine Engineering Pipeline

[![DOI](https://img.shields.io/badge/DOI-10.6084%2Fm9.figshare.33471685-blue.svg)](https://doi.org/10.6084/m9.figshare.33471685)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.8+](https://img.shields.io/badge/Python-3.8%2B-brightgreen.svg)](https://www.python.org/)
[![R: 4.0+](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://www.r-project.org/)

An automated computational pipeline for batch sequence retrieval, quality cleaning, redundancy reduction, immunological filtering, and multi-epitope construct engineering across four human-pathogenic *Orthoebolavirus* species:
- **Bundibugyo ebolavirus (BDBV)** (NCBI Taxonomy ID: 565995)
- **Sudan ebolavirus (SUDV)** (NCBI Taxonomy ID: 186540)
- **Taï Forest ebolavirus (TAFV)** (NCBI Taxonomy ID: 186541)
- **Zaire ebolavirus (EBOV)** (NCBI Taxonomy ID: 186538)

---

## Permanent Archive & Citation
This pipeline is permanently archived on Figshare:
- **DOI:** [10.6084/m9.figshare.33471685](https://doi.org/10.6084/m9.figshare.33471685)
- **Public Record:** [https://figshare.com/articles/software/Automated_NCBI_Sequence_Retrieval_Pipeline_for_Orthoebolavirus_Proteomes/33471685](https://figshare.com/articles/software/Automated_NCBI_Sequence_Retrieval_Pipeline_for_Orthoebolavirus_Proteomes/33471685)

---

## Pipeline Architecture & Included Scripts

```
[NCBI Entrez Database]
         │
         ▼ (Stage 1: Automated Retrieval via Python)
Raw FASTA Sequences (7 Viral Proteins × 4 Species)
         │
         ▼ (Stage 2: Canonical Filtering via R)
Cleaned FASTA Files (Standard 20 Amino Acids Only)
         │
         ▼ (Stage 3: Redundancy Reduction via CD-HIT & Python)
Non-Redundant Clusters (99% Identity & Coverage)
         │
         ▼ (Stage 4 & 5: Epitope Filtering & Selection via Python)
High-Affinity, Non-Allergenic, Non-Toxic CTL/HTL/LBL Epitopes
         │
         ▼ (Stage 6 & 7: Combinatorial Construct Assembly via Python)
Multi-Epitope Construct Sequences & Architecture Maps
         │
         ▼ (Stage 8: Physicochemical & Immunological Ranking via Python)
Prioritized Vaccine Constructs
```

### 1. Sequence Retrieval Downloaders (`*_sequence_downloader.py`)
- `Bundibugyo_sequence_downloader.py`
- `Sudan_sequence_downloader.py`
- `Tai_forest_sequence_downloader.py`
- `Zaire_sequence_downloader.py`
*Features:* Automated batch querying of the NCBI Entrez Protein database using Biopython E-utilities for seven target proteins (NP, VP35, VP40, GP, VP30, VP24, L) with negative selection filters and strict length constraints.

### 2. Canonical Amino Acid Quality Cleaning (`clean_standard_amino_acids.R`)
*Features:* Validates each sequence strictly against the 20 standard canonical amino acids (`ACDEFGHIKLMNPQRSTVWY`), filtering out non-standard or ambiguous residues (`X`, `B`, `Z`, `J`) with automated CSV audit trails.

### 3. High-Stringency Redundancy Reduction (`run_cdhit_representatives.py`)
*Features:* Automates CD-HIT representative clustering at `-c 0.99 -aS 0.99 -aL 0.99 -G 1 -g 1 -n 5` across all cleaned datasets.

### 4. Multi-Criteria Epitope Filtering (`filter_epitopes_for_downstream.py`)
*Features:* Filters CTL (percentile rank $\le 1.0$), HTL (percentile rank $\le 10.0$), and linear B-cell epitopes ($\ge 6$ aa) across prediction tools.

### 5. Final Candidate Epitope Selection (`final_epitope_candidates.py`)
*Features:* Cross-evaluates predicted epitopes against allergenicity (AllerTOP), toxicity (ToxinPred2), antigenicity (VaxiJen), cytokine-inducing capacity (IFN-γ, IL-4, IL-10), and proteomic conservancy ($\ge 95\%$).

### 6. Strain-Aware Vaccine Construct Selector (`strain_aware_vaccine_construct_selector.py`)
*Features:* Combinatorially optimizes epitope composition to maximize balanced cross-strain coverage across BDBV, SUDV, TAFV, and EBOV.

### 7. Construct Sequence & Architecture Generator (`build_vaccine_construct_sequences_and_maps.py`)
*Features:* Assembles full multi-epitope constructs incorporating adjuvant sequences (e.g., human $\beta$-defensin), linkers (EAAAK, AAY, GPGPG, KK), CTL/HTL/LBL epitopes, and terminal tags; exports FASTA sequences and visual architecture diagrams.

### 8. Multi-Parametric Construct Screening & Ranking (`filter_and_rank_final_constructs.py`)
*Features:* Filters and ranks candidate vaccine constructs based on global antigenicity, non-allergenicity, non-toxicity, ProtParam stability metrics (instability index, GRAVY, aliphatic index), and solubility.

---

## Installation & Requirements

### 1. Python Environment (>= 3.8)
Install required Python dependencies:
```bash
pip install -r requirements.txt
```
*Required packages:* `biopython>=1.80`, `Pillow>=9.0.0`

### 2. R Environment (>= 4.0)
Base R is required for `clean_standard_amino_acids.R`. No external CRAN packages required.

### 3. CD-HIT (>= 4.8.1)
Ensure `cd-hit` (or `cd-hit.exe`) is available in system `PATH` or passed via `--cd-hit`.

---

## Usage Guide

1. **Retrieve Sequences:**
   ```bash
   python Bundibugyo_sequence_downloader.py
   python Sudan_sequence_downloader.py
   python Tai_forest_sequence_downloader.py
   python Zaire_sequence_downloader.py
   ```

2. **Clean Sequences:**
   ```bash
   Rscript clean_standard_amino_acids.R --input-root=path/to/raw_fasta --output-dir=path/to/cleaned_output
   ```

3. **Reduce Redundancy:**
   ```bash
   python run_cdhit_representatives.py --input-dir=path/to/cleaned_output --output-dir=path/to/cdhit_output
   ```

4. **Filter & Select Epitopes:**
   ```bash
   python filter_epitopes_for_downstream.py
   python final_epitope_candidates.py
   ```

5. **Assemble & Rank Vaccine Constructs:**
   ```bash
   python strain_aware_vaccine_construct_selector.py
   python build_vaccine_construct_sequences_and_maps.py
   python filter_and_rank_final_constructs.py
   ```

---

## License
Distributed under the open-source **MIT License**.
