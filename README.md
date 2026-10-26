# Glioblastoma Temozolomide Resistance: Gene-Expression Analysis

## Project Overview

This project explores gene-expression differences between temozolomide (TMZ)-sensitive and TMZ-resistant glioblastoma cell models using publicly available RNA-sequencing data.

The analysis was developed in Python to practice computational biology, gene-expression analysis, scientific visualization, and reproducible research.

## Research Question

Which genes show consistent expression differences between TMZ-sensitive and TMZ-resistant glioblastoma cell models?

## Dataset

**Source:** NCBI Gene Expression Omnibus (GEO)

**Accession:** GSE327029

**Dataset:** GSE327029_mRNA-12.TPM.genes.tsv

The dataset contains:

- 60,498 gene entries
- 12 RNA-sequencing samples
- Two glioblastoma cell models: SF126 and U118MG
- Sensitive and resistant conditions, with three samples per group

Data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327029

## Technologies Used

- Python
- pandas
- NumPy
- Matplotlib
- PyCharm
- Git and GitHub

## Methodology

1. Imported the TPM-normalized gene-expression matrix.
2. Organized the 12 samples into four groups.
3. Calculated mean expression for each gene in each group.
4. Calculated exploratory log2 fold changes between resistant and sensitive cells.
5. Selected candidate genes with changes in the same direction in both models, absolute log2 fold changes of at least 1, and minimum expression criteria.
6. Created a scatter plot comparing expression changes.
7. Examined KPNA3, MGMT, and STAT3.
8. Checked missing values, duplicate gene IDs, and within-group sample correlations.

## Key Results

The exploratory analysis identified **44 candidate gene entries**:

- 37 showed lower expression in both resistant models.
- 7 showed higher expression in both resistant models.

Examples of candidates with lower expression included CXCL1, FABP5, IER3, CYP24A1, and OCIAD2.

KPNA3 and MGMT showed higher average RNA expression in both resistant cell models.

Quality-control checks found zero missing values and zero duplicate gene identifiers. Within-group sample correlations ranged from 0.972 to 0.997.

## Research Figure

The file `figures/candidate_genes.png` displays the log2 fold-change comparison between the two cell models for the 44 selected candidates.

## Project Files

- `inspect_data.py` — Python analysis script
- `research_report.md` — Scientific report
- `results/all_gene_changes.csv` — Exploratory expression-change measurements
- `results/candidate_genes.csv` — Selected candidate genes
- `results/important_genes.csv` — KPNA3, MGMT, and STAT3 results
- `figures/candidate_genes.png` — Scientific scatter plot

## Reproducing the Analysis

Download the gene-expression dataset from GEO accession GSE327029.

Place the decompressed file `GSE327029_mRNA-12.TPM.genes.tsv` inside a folder named `data`.

Install the required Python libraries:

```bash
python -m pip install pandas numpy matplotlib
```

Run the analysis:

```bash
python inspect_data.py
```

The script generates CSV result files and a scientific figure. Close the displayed graph window if necessary to allow the script to finish.

## Limitations

This is an exploratory reanalysis of previously published data, not an independent experimental study.

The project uses TPM-normalized values, and formal differential-expression statistical testing was not performed. The fold-change calculations use a pseudocount of 1, and candidate-selection thresholds are descriptive rather than statistical.

The results do not establish causation, clinical utility, or experimentally validated biomarkers.

## Future Research

Future work could include:

- Formal differential-expression analysis using raw RNA-sequencing counts
- Gene-set and pathway enrichment
- Validation using independent datasets
- Investigation of potential molecular mechanisms

## Reference

Gao et al. (2026). *KPNA3 drives temozolomide resistance in glioblastoma by upregulating MGMT and activating STAT3 to sustain glioma stem cells*. Cell Biology International.

NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327029
