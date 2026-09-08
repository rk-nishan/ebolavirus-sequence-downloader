# Orthoebolavirus NCBI Sequence Downloader

[![DOI](https://img.shields.io/badge/DOI-10.6084%2Fm9.figshare.33471685-blue.svg)](https://doi.org/10.6084/m9.figshare.33471685)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.8+](https://img.shields.io/badge/Python-3.8%2B-brightgreen.svg)](https://www.python.org/)

Automated Python pipeline for batch retrieval of full-length structural and non-structural protein sequences from the NCBI Entrez Protein database across human-pathogenic *Orthoebolavirus* species.

## Permanent Archive & Citation
This pipeline is permanently archived on Figshare:
- **DOI:** [10.6084/m9.figshare.33471685](https://doi.org/10.6084/m9.figshare.33471685)
- **Direct Figshare Record:** https://figshare.com/articles/software/Automated_NCBI_Sequence_Retrieval_Pipeline_for_Orthoebolavirus_Proteomes/33471685

## Included Scripts

- `Bundibugyo_sequence_downloader.py`: Retrieval pipeline for *Bundibugyo ebolavirus* (NCBI Taxonomy ID: 565995).
- `Sudan_sequence_downloader.py`: Retrieval pipeline for *Sudan ebolavirus* (NCBI Taxonomy ID: 186540).
- `Tai_forest_sequence_downloader.py`: Retrieval pipeline for *Taï Forest ebolavirus* (NCBI Taxonomy ID: 186541).
- `Zaire_sequence_downloader.py`: Retrieval pipeline for *Zaire ebolavirus* (NCBI Taxonomy ID: 186538).

## Target Proteins & Curation
For each species, the pipeline systematically retrieves seven viral targets with predefined sequence length boundaries:
1. Nucleoprotein (NP)
2. Polymerase cofactor (VP35)
3. Matrix protein (VP40)
4. Surface glycoprotein (GP)
5. Minor nucleoprotein (VP30)
6. Membrane-associated protein (VP24)
7. RNA-dependent RNA polymerase (L)

All query strings incorporate negative-selection filtering to exclude partial, fragment, chain, mutant, synthetic, and construct entries directly at the database query stage.

## Requirements
- Python >= 3.8
- Biopython >= 1.80

Install dependencies via:
```bash
pip install -r requirements.txt
```

## Usage

1. Open the target script and configure your contact email (required by NCBI Entrez API policies):
   ```python
   Entrez.email = "your.email@example.com"
   ```
   *(Optional)* If you have an NCBI personal API key, specify `Entrez.api_key = "your_key"`.

2. Execute the downloader for the target species:
   ```bash
   python Bundibugyo_sequence_downloader.py
   python Sudan_sequence_downloader.py
   python Tai_forest_sequence_downloader.py
   python Zaire_sequence_downloader.py
   ```

3. Output Structure:
   Each script generates dedicated folders per protein containing:
   - `<Protein_Name>.fasta`: Retrieved protein sequences in FASTA format.
   - `<Protein_Name>_query.txt`: Verbatim Entrez search query text.
   - `download_report.txt`: Summary of records retrieved and written.

## License
Distributed under the MIT License.
