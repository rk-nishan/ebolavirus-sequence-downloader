# Orthoebolavirus Proteome Curation & Redundancy Reduction Pipeline

[![DOI](https://img.shields.io/badge/DOI-10.6084%2Fm9.figshare.33471685-blue.svg)](https://doi.org/10.6084/m9.figshare.33471685)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.8+](https://img.shields.io/badge/Python-3.8%2B-brightgreen.svg)](https://www.python.org/)
[![R: 4.0+](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://www.r-project.org/)

An automated computational pipeline for batch sequence retrieval, canonical amino acid quality filtering, and high-stringency redundancy reduction of viral proteomes across four human-pathogenic *Orthoebolavirus* species:
- **Bundibugyo ebolavirus (BDBV)** (NCBI Taxonomy ID: 565995)
- **Sudan ebolavirus (SUDV)** (NCBI Taxonomy ID: 186540)
- **Taï Forest ebolavirus (TAFV)** (NCBI Taxonomy ID: 186541)
- **Zaire ebolavirus (EBOV)** (NCBI Taxonomy ID: 186538)

---

## Permanent Archive & Citation
This pipeline is permanently archived on Figshare:
- **DOI:** [10.6084/m9.figshare.33471685](https://doi.org/10.6084/m9.figshare.33471685)
- **Direct Record:** [https://figshare.com/articles/software/Automated_NCBI_Sequence_Retrieval_Pipeline_for_Orthoebolavirus_Proteomes/33471685](https://figshare.com/articles/software/Automated_NCBI_Sequence_Retrieval_Pipeline_for_Orthoebolavirus_Proteomes/33471685)

---

## Pipeline Architecture & Workflow

```
[NCBI Entrez Database]
         │
         ▼ (Stage 1: Automated Retrieval via Python)
Raw FASTA Files (7 Proteins × 4 Species)
         │
         ▼ (Stage 2: Canonical Filtering via R)
Cleaned FASTA Files (Standard 20 Amino Acids Only, ACDEFGHIKLMNPQRSTVWY)
         │
         ▼ (Stage 3: Redundancy Reduction via CD-HIT & Python)
Representative Non-Redundant Clusters (99% Identity & 99% Coverage)
```

### Stage 1: Automated NCBI Entrez Sequence Retrieval
Dedicated Python scripts query the NCBI Entrez Protein database using Biopython E-utilities:
- `Bundibugyo_sequence_downloader.py`
- `Sudan_sequence_downloader.py`
- `Tai_forest_sequence_downloader.py`
- `Zaire_sequence_downloader.py`

**Target Proteins:**
1. Nucleoprotein (NP)
2. Polymerase cofactor (VP35)
3. Matrix protein (VP40)
4. Surface glycoprotein (GP)
5. Minor nucleoprotein (VP30)
6. Membrane-associated protein (VP24)
7. RNA-dependent RNA polymerase (L)

Each query enforces strict negative-selection exclusion filters (`NOT (partial OR fragment OR chain OR mutant OR synthetic OR construct)`) and species-specific coding length constraints (`[SLEN]`).

### Stage 2: Canonical Amino Acid Quality Filtering
Implemented in `clean_standard_amino_acids.R`:
- Validates each retrieved sequence against the 20 standard canonical amino acids (`ACDEFGHIKLMNPQRSTVWY`).
- Excludes records containing non-standard, ambiguous, or undetermined characters (e.g., `X`, `B`, `Z`, `J`).
- Generates structured CSV audit trails:
  - `fasta_cleaning_summary.csv`: Per-protein counts of input, retained, and removed records.
  - `fasta_cleaning_record_log.csv`: Full accession-level audit log.
  - `fasta_cleaning_removed_sequences.csv`: List of excluded accessions with detected invalid characters.

### Stage 3: High-Stringency Redundancy Reduction (CD-HIT)
Implemented in `run_cdhit_representatives.py`:
- Coordinates CD-HIT execution across all cleaned protein datasets.
- Parameters enforced:
  - Sequence identity: `-c 0.99` (99% identity)
  - Short-sequence coverage: `-aS 0.99` (99% alignment coverage)
  - Long-sequence coverage: `-aL 0.99` (99% alignment coverage)
  - Global sequence identity: `-G 1`
  - Accurate clustering mode: `-g 1`
  - Word length: `-n 5`
- Preserves natural strain diversity while removing duplicate isolates.
- Generates representative FASTA files (`*_cdhit099.fasta`), cluster files (`.clstr`), and an execution summary report (`cdhit099_summary.csv`).

---

## Installation & Dependencies

### 1. Python Environment
```bash
pip install -r requirements.txt
```
*Required packages:*
- `biopython>=1.80`

### 2. R Environment
Requires base R (>= 4.0). No external CRAN packages required.

### 3. CD-HIT
Download and install CD-HIT from [CD-HIT GitHub releases](https://github.com/weizhongli/cdhit/releases) and ensure `cd-hit` (or `cd-hit.exe`) is available in your system `PATH`, or specify its location via `--cd-hit`.

---

## Usage Instructions

### Step 1: Sequence Retrieval
Configure your email in the target script (`Entrez.email = "your.email@example.com"`) and execute:
```bash
python Bundibugyo_sequence_downloader.py
python Sudan_sequence_downloader.py
python Tai_forest_sequence_downloader.py
python Zaire_sequence_downloader.py
```

### Step 2: Quality Cleaning
```bash
Rscript clean_standard_amino_acids.R --input-root=path/to/downloaded_fasta --output-dir=path/to/cleaned_output
```

### Step 3: CD-HIT Clustering
```bash
python run_cdhit_representatives.py --input-dir=path/to/cleaned_output --output-dir=path/to/cdhit_output
```

---

## License
Distributed under the open-source **MIT License**.
