# Exploratory Analysis of Gene-Expression Patterns Associated with Temozolomide Resistance in Glioblastoma Cell Models

## Abstract

Glioblastoma is an aggressive brain cancer, and resistance to temozolomide (TMZ) presents a major treatment challenge. This project explored gene-expression differences between TMZ-sensitive and TMZ-resistant glioblastoma cell models using publicly available RNA-sequencing data from NCBI GEO accession GSE327029.

Using Python, pandas, and NumPy, I analyzed a dataset containing 60,498 gene entries across 12 samples. Genes were compared separately in the SF126 and U118MG cell models. An exploratory filtering procedure identified 44 candidate genes showing changes in the same direction in both models, with 37 showing decreased expression and 7 showing increased expression in resistant cells.

The results highlight potentially interesting expression patterns for further investigation, but they do not establish statistical significance or causal relationships.

## 1. Introduction

Temozolomide is a chemotherapy drug used to treat glioblastoma. However, cancer cells may develop resistance, reducing the effectiveness of treatment.

Gene-expression analysis can help researchers explore biological differences between treatment-sensitive and treatment-resistant cells.

The objective of this project was to investigate which genes show consistent expression differences between TMZ-sensitive and TMZ-resistant glioblastoma cell models.

## 2. Materials and Methods

### Dataset

The study used publicly available RNA-sequencing expression data from GEO accession GSE327029.

The dataset contained 12 samples divided into four groups:

- SF126 sensitive: 3 samples
- SF126 resistant: 3 samples
- U118MG sensitive: 3 samples
- U118MG resistant: 3 samples

The dataset included 60,498 gene entries and contained TPM-normalized expression measurements.

### Computational analysis

Python was used to load and analyze the expression matrix.

For each gene, mean expression was calculated separately for the sensitive and resistant groups in each cell model.

An exploratory log2 fold-change measure was calculated:

**log2FC = log2((mean resistant TPM + 1) / (mean sensitive TPM + 1))**

The addition of 1 allowed calculations involving zero-expression values.

Candidate genes were selected when they met all of the following criteria:

1. Expression changed in the same direction in both models.
2. The absolute log2 fold-change was at least 1 in both models.
3. At least one condition had a mean expression of 5 TPM or greater within each model.

These thresholds were chosen for descriptive exploration, not formal statistical differential-expression testing.

### Data-quality checks

Missing measurements and duplicate gene identifiers were checked. Within-group sample correlations were calculated using log2(TPM + 1) measurements.

## 3. Results

The exploratory filtering process identified **44 candidate gene entries**.

Among these:

- 37 showed decreased mean expression in resistant cells in both models.
- 7 showed increased mean expression in resistant cells in both models.

Examples of genes showing substantial decreases included CXCL1, FABP5, IER3, CYP24A1, and OCIAD2.

Three genes discussed in the original publication were also examined:

| Gene | SF126 log2FC | U118MG log2FC |
|---|---:|---:|
| KPNA3 | +5.43 | +0.45 |
| MGMT | +4.40 | +0.73 |
| STAT3 | −0.62 | +0.34 |

KPNA3 and MGMT showed higher average RNA expression in resistant cells in both models. STAT3 showed changes in opposite directions between the models.

The quality-control analysis found no missing values or duplicate gene identifiers. Within-group sample correlations ranged from 0.972 to 0.997.

**Figure 1:** Scatter plot of exploratory log2 fold changes for the 44 selected candidate genes in the SF126 and U118MG cell models.

## 4. Discussion

The findings suggest that some gene-expression differences are shared between the two TMZ-resistant glioblastoma cell models.

The exploratory results showed a greater number of candidate genes with decreased expression than increased expression in resistant cells.

KPNA3 and MGMT also displayed higher mean RNA expression in both resistant groups. These observations are relevant to the biological questions investigated by the original study, although they do not independently establish the mechanisms responsible for resistance.

STAT3 expression changed in different directions across the two models, demonstrating that gene-expression patterns are not necessarily consistent between different cellular backgrounds. STAT3 signaling activity cannot be determined from RNA abundance alone.

## 5. Limitations

This project was an exploratory reanalysis of previously published data.

The following limitations should be considered:

- The dataset contained TPM-normalized measurements rather than raw sequencing counts.
- Formal differential-expression significance tests were not performed.
- The calculations used a pseudocount of 1, which can influence fold changes for genes with very low expression.
- Only two cell models from one published study were investigated.
- The independence of the reported replicates was not independently verified.
- The candidate genes have not been experimentally validated in this project.

Therefore, the findings should be interpreted as descriptive gene-expression patterns rather than validated biomarkers or causal mechanisms.

## 6. Conclusion

This project successfully used Python to explore RNA-sequencing gene-expression data associated with temozolomide resistance in glioblastoma cell models.

The analysis identified 44 candidate gene entries with consistent directional expression changes across two cell-model comparisons. It also demonstrated the application of Python programming, biological data processing, visualization, and basic quality-control methods.

Future research could include formal differential-expression analysis using raw count data, pathway enrichment, independent dataset validation, and experimental investigations.

## Reference

Gao et al. (2026). *KPNA3 drives temozolomide resistance in glioblastoma by upregulating MGMT and activating STAT3 to sustain glioma stem cells*. Cell Biology International.

NCBI Gene Expression Omnibus: GSE327029  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327029
