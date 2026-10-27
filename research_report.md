# Exploratory Gene-Expression Patterns in Temozolomide-Resistant Glioblastoma Cell Models

**Computational Research Portfolio Report**  
**Salma Mifdal**  
**October 2026**

*This project is a computational reanalysis of publicly available data, not a new laboratory experiment or a peer-reviewed publication.*

## Abstract

Temozolomide (TMZ) resistance is a major challenge in the treatment of glioblastoma. This project explored RNA-expression patterns associated with TMZ resistance in two glioblastoma cell models using the publicly available Gene Expression Omnibus (GEO) dataset GSE327029.

The dataset contained 60,498 gene entries and 12 samples divided into sensitive and resistant SF126 and U118MG groups. Using Python, pandas, NumPy, and Matplotlib, I calculated average gene expression, compared expression changes, and identified exploratory candidate genes.

The analysis identified 44 candidate gene entries showing consistent directional changes across both models. Of these, 37 showed lower expression and 7 showed higher expression in resistant cells.

KPNA3 and MGMT showed higher average expression in resistant cells in both models, while STAT3 exhibited different expression patterns between the two models. Data-quality checks identified no missing measurements or duplicate gene identifiers, and within-group sample correlations ranged from 0.972 to 0.997.

These findings are exploratory and require further statistical and biological validation.

## 1. Introduction

Glioblastoma is an aggressive form of brain cancer. Temozolomide is a chemotherapy drug commonly used to treat glioblastoma, but cancer cells can become resistant to treatment.

Understanding the biological differences between sensitive and resistant cancer cells may help researchers investigate why resistance develops.

Gene-expression analysis provides a way to examine how genes behave under different biological conditions.

In this project, I used publicly available RNA-sequencing data to explore expression differences between TMZ-sensitive and TMZ-resistant glioblastoma cell models.

### Research Question

**Which gene-expression patterns differ between TMZ-sensitive and TMZ-resistant glioblastoma cell models, and which changes are shared by SF126 and U118MG?**

The goal was to build a reproducible computational analysis and identify patterns for further investigation.

## 2. Dataset and Experimental Design

The dataset was obtained from the National Center for Biotechnology Information (NCBI) Gene Expression Omnibus.

**GEO accession:** GSE327029

**Data file:** `GSE327029_mRNA-12.TPM.genes.tsv`

The dataset contained 60,498 gene entries and 14 columns, including gene identifiers, gene names, and expression measurements for 12 samples.

The samples were divided into four groups:

| Group | Samples |
|---|---|
| SF126 sensitive | SF126-1, SF126-2, SF126-3 |
| SF126 resistant | SF126R-1, SF126R-2, SF126R-3 |
| U118MG sensitive | U118MG-1, U118MG-2, U118MG-3 |
| U118MG resistant | U118MGR-1, U118MGR-2, U118MGR-3 |

Each group contained three labeled samples.

The two cell models were analyzed separately before comparing their shared expression patterns. The biological independence of the samples was not independently verified.

## 3. Materials and Methods

### 3.1 Programming Tools

I used the following tools to complete the analysis:

- Python for data processing and calculations.
- pandas for reading and organizing the gene-expression data.
- NumPy for numerical calculations.
- Matplotlib for scientific visualization.
- PyCharm for writing and executing the code.
- GitHub for documenting and sharing the project.

### 3.2 Data Processing

The RNA-expression dataset was imported into Python using pandas.

For every gene entry, I calculated the average expression across the three samples in each experimental group.

This produced four average-expression measurements per gene entry: SF126 sensitive, SF126 resistant, U118MG sensitive, and U118MG resistant.

### 3.3 Expression-Change Calculations

I used a pseudocount-adjusted log2 fold-change calculation to compare resistant and sensitive conditions.

**log2FC = log2[(mean resistant TPM + 1) / (mean sensitive TPM + 1)]**

The addition of 1 prevented division by zero.

A positive log2FC indicated higher average expression in resistant cells, while a negative log2FC indicated lower average expression in resistant cells.

Because this calculation used TPM-normalized measurements rather than raw sequencing counts, it was used for exploratory comparisons rather than formal statistical differential-expression testing.

### 3.4 Candidate-Gene Selection

Candidate genes were selected using three criteria:

1. Expression had to change in the same direction in both cell models.
2. The absolute log2 fold change had to be at least 1 in both models.
3. At least one group mean had to be 5 TPM or greater within each model.

Candidates were ranked according to the smaller absolute log2 fold change across the two models.

These were descriptive filtering rules and did not establish statistical significance.

### 3.5 Data-Quality Checks

I checked the dataset for missing expression values and duplicate gene identifiers.

I also calculated pairwise correlations among samples within each group using log2(TPM + 1) values.

These checks helped evaluate sample consistency but did not replace comprehensive RNA-sequencing quality control.

## 4. Results

### 4.1 Dataset Quality

The initial analysis produced the following results:

| Quality Measure | Result |
|---|---|
| Gene entries | 60,498 |
| Samples | 12 |
| Experimental groups | 4 |
| Missing expression values | 0 |
| Duplicate gene IDs | 0 |
| Within-group correlations | 0.972–0.997 |

Most replicate correlations were between 0.996 and 0.997.

The correlations involving U118MG-1 and the other two U118MG sensitive samples were slightly lower, at 0.972.

These observations suggest generally similar expression profiles within groups, although additional quality-control analysis would be needed.

### 4.2 Candidate-Gene Identification

Using the exploratory filtering criteria, the analysis identified **44 candidate gene entries**.

Among these:

- 37 showed lower expression in resistant cells in both models.
- 7 showed higher expression in resistant cells in both models.

The following table presents ten selected candidates:

| Gene | SF126 log2FC | U118MG log2FC |
|---|---:|---:|
| CXCL1 | -7.127 | -5.281 |
| FABP5 | -4.318 | -4.984 |
| IER3 | -7.685 | -3.244 |
| CYP24A1 | -3.165 | -4.063 |
| OCIAD2 | -3.152 | -3.709 |
| CD24 | -3.093 | -3.905 |
| NETO2 | -5.459 | -2.897 |
| SNRPGP10 | -3.521 | -2.828 |
| TSPAN15 | -4.079 | -2.747 |
| KRT81 | -2.657 | -5.091 |

These genes showed substantial descriptive decreases in expression under the selected comparison method.

The higher number of negatively changing candidates applies only to the selected gene set and should not be interpreted as a global decrease in gene activity.

### 4.3 Scientific Visualization

The scatter plot below compares gene-expression changes between the two cell models.

![Scatter plot of shared gene-expression changes](figures/candidate_genes.png)

**Figure 1.** Scatter plot of pseudocount-adjusted log2 fold changes for 44 exploratory candidate gene entries.

The horizontal axis represents expression changes in SF126 cells, while the vertical axis represents expression changes in U118MG cells.

Genes in the lower-left quadrant have lower average expression in resistant cells in both models.

Genes in the upper-right quadrant have higher average expression in resistant cells in both models.

This figure visualizes expression patterns but does not represent statistical significance.

### 4.4 Investigation of KPNA3, MGMT, and STAT3

I separately examined three genes highlighted in the original research study.

| Gene | SF126 log2FC | U118MG log2FC |
|---|---:|---:|
| KPNA3 | +5.432 | +0.449 |
| MGMT | +4.401 | +0.735 |
| STAT3 | -0.619 | +0.341 |

**KPNA3:** Average RNA expression was higher in resistant cells in both models. The increase was considerably larger in SF126.

**MGMT:** Average RNA expression was also higher in both resistant models. MGMT is relevant to the repair of certain types of DNA damage and has an established connection to temozolomide response.

**STAT3:** RNA-expression changes differed between the two cell models. Its expression decreased in resistant SF126 cells but increased in resistant U118MG cells.

These differences emphasize that gene-expression behavior can vary across cellular backgrounds.

Importantly, STAT3 RNA-expression measurements cannot independently establish whether STAT3 signaling was activated.

## 5. Discussion

The main finding of this project was the identification of 44 exploratory gene entries showing expression changes in the same direction across two glioblastoma cell models.

Several genes, including CXCL1, FABP5, IER3, CYP24A1, and OCIAD2, showed relatively large decreases in expression in resistant cells under the chosen method.

The results also highlighted KPNA3 and MGMT, which had higher average RNA expression in resistant cells in both models.

These patterns are relevant to the original study's investigation of TMZ resistance, but they do not prove that these genes cause resistance.

The analysis also showed that STAT3 RNA expression was not consistent across both models, demonstrating the importance of considering biological differences between cell lines.

This project helped demonstrate how computational approaches can be used to organize biological data, investigate patterns, create scientific visualizations, and generate testable research questions.

## 6. Limitations

Several limitations should be considered when interpreting the results.

**Use of TPM data:** The analysis used normalized TPM expression values rather than raw sequencing counts. Formal differential-expression testing was not performed.

**Pseudocount adjustment:** Adding 1 to the expression means influenced the calculated fold changes, particularly for genes with very low expression.

**Selection criteria:** The 44 candidates depended on the chosen expression thresholds. Different thresholds could produce different results.

**Limited datasets:** Both cell models were included in a single published dataset. Independent datasets were not used for external validation.

**Replicate information:** The biological independence of the labeled replicates was not independently confirmed.

**Limited quality control:** Correlation checks cannot rule out all technical differences, batch effects, or sample-quality issues.

**No experimental validation:** The project did not experimentally test gene functions or establish causal relationships between gene expression and TMZ resistance.

Therefore, these results should be interpreted as exploratory expression patterns rather than confirmed biomarkers or therapeutic targets.

## 7. Future Research Directions

Future work could improve this analysis by obtaining raw RNA-sequencing counts and performing appropriate statistical differential-expression testing.

Additional research could also investigate biological pathways associated with candidate genes, examine possible batch effects, and compare the results with independent datasets.

Experimental studies would be required to investigate whether the identified genes directly contribute to chemotherapy resistance.

## 8. Conclusion

This project used Python to explore RNA-expression differences between TMZ-sensitive and TMZ-resistant glioblastoma cell models.

The analysis examined 60,498 gene entries and identified 44 exploratory candidates showing consistent directional changes across both models.

Of these candidates, 37 had lower expression and 7 had higher expression in resistant cells.

The project also examined KPNA3, MGMT, and STAT3 and evaluated initial sample-quality measurements.

Although the findings are descriptive and require further validation, the project demonstrates practical skills in Python programming, computational biology, data processing, visualization, and scientific reporting.

## Reproducibility

The Python analysis, scientific figure, and generated results are available in the public GitHub repository:

https://github.com/salmamifdal/GBM-TMZ-Resistance-Analysis

The project includes the analysis script, results tables, scientific visualization, and Python dependencies.

The original gene-expression dataset can be obtained from GEO accession GSE327029.

The downloaded dataset should be placed at `data/GSE327029_mRNA-12.TPM.genes.tsv` before running `inspect_data.py`.

## References

1. NCBI Gene Expression Omnibus. GEO Series GSE327029: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327029

2. Gao, Y., et al. (2026). *KPNA3 drives temozolomide resistance in glioblastoma by upregulating MGMT and activating STAT3 to sustain glioma stem cells*. *Cell Biology International*. Full publication metadata should be verified before formal academic submission.

*This report presents an exploratory computational reanalysis of publicly available gene-expression data.*
