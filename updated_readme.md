# Statistical Analysis of Gene Expression in Breast Cancer (TCGA-BRCA)

Course project for *Statistics* (MSc in Bioinformatics for Computational Genomics, A.Y. 2026/27).
Team: [add team members' names].

## Overview

This project uses real breast cancer gene-expression data from the TCGA-BRCA cohort to apply the statistical methods learned in the course: descriptive statistics, graphical analysis, confidence intervals, t-tests, correlation and multiple linear regression. The biological context is breast cancer, but the focus is on statistical analysis, interpretation and visualization. All analyses are done in R.

We work on a selected set of ten genes instead of the whole transcriptome, and we study how their expression differs between tumor and normal tissue, how the genes relate to each other and to patient age, and whether the expression of the proliferation marker MKI67 can be explained statistically by the other genes and by age.

## Research questions

| # | Question | Status |
|---|---|---|
| 1 | **Tumor vs normal:** is the expression of the selected genes different between tumor and normal breast tissue? | Exploratory analysis done; formal tests planned |
| 2 | **Subtypes:** do the genes show different expression between molecular subtypes (e.g. Basal-like vs Luminal A)? | Planned extension (see *Future work*) |
| 3 | **Relationships:** is the expression of the genes related to each other and to patient age? | Exploratory analysis done |
| 4 | **MKI67 model:** can MKI67 expression be explained by the other genes and age? | Exploratory plots done; regression planned |

## Data

Data come from the TCGA-BRCA project (Breast Invasive Carcinoma), downloaded from the UCSC Xena browser (GDC hub):

- Gene expression (RNA-seq, STAR, FPKM, log2 scale): [TCGA-BRCA.star_fpkm.tsv](https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.star_fpkm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
- Clinical / sample information: [TCGA-BRCA.clinical.tsv](https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.clinical.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
- Gene annotation (`gencode.v36.annotation.gtf.gene.probemap`): maps Ensembl IDs to gene symbols.

### Selected genes

| Role | Genes |
|---|---|
| Cell proliferation | MKI67, CCNB1, CDK1 |
| Hormone receptors / luminal markers | ESR1, PGR, FOXA1 |
| Growth factor receptors | ERBB2, EGFR |
| DNA repair / tumor suppression | BRCA1, TP53 |

### Data preparation

1. Read the expression matrix (60,660 genes x 1,226 samples), the clinical table (1,255 records) and the gene map.
2. Keep the rows of the ten selected genes (Ensembl IDs matched to gene symbols) and transpose: one row per sample.
3. Join the clinical table by sample barcode (all 1,226 samples were found) to add sample type (`sample_type.samples`) and age (`age_at_index.demographic`).
4. Remove 7 metastatic samples, keeping only *Primary Tumor* and *Solid Tissue Normal* (1,219 samples).
5. Remove 11 duplicate tumor samples (same patient, same tissue type), keeping the first one per patient.

Final dataset: **1,208 samples = 1,095 primary tumors + 113 normal tissue samples**. All 113 normal samples come from patients who also have a tumor sample, so there are **113 matched tumor/normal pairs** for paired tests. Age is missing for one record in the clinical table.

## Repository contents (Flat directory structure)

All project files (scripts, raw data, processed data, tables, figures, and presentation materials) are stored together in the main repository directory.

### Scripts & Source Code

| File | Description |
|---|---|
| `work-in-progres_R-scripts.R` (or `work-in-progress_R-scripts_en.R`) | R script containing all data cleaning, descriptive statistics, and plotting steps (with English comments). |

### Downloaded data (raw)

| File | Description |
|---|---|
| `TCGA-BRCA.star_fpkm.tsv.gz` | Gene-expression matrix from Xena (rows = genes, columns = samples). *Note: Large file.* |
| `TCGA-BRCA.clinical.tsv.gz` | Clinical and sample information (age, sample type, stage, etc.). |
| `gencode.v36.annotation.gtf.gene.probemap` | Mapping from Ensembl gene IDs to gene symbols. |

### Processed data

| File | Description |
|---|---|
| `brca_clean.csv` | Clean dataset used for analyses. One row per sample, standard CSV format. |
| `brca_clean_excel.csv` | Same data, saved with semicolon separator and decimal comma for Excel (Italian locale). |

### Tables

| File | Description |
|---|---|
| `table_descriptive.csv` | Mean (standard deviation) of each gene in normal and tumor samples. |
| `table_descriptive_excel.csv` | Same table in Excel-friendly format. |

### Figures

| File | Description |
|---|---|
| `hist_genes.png` | Histograms of the ten genes (all samples, tumor and normal together). |
| `box_tumor_normal.png` | Boxplots of each gene in normal vs tumor tissue. |
| `corr_heatmap.png` | Correlation heatmap of the ten genes (tumor samples only). |
| `scatter_mki67.png` | Scatter plots of MKI67 against predictors and age (tumor samples only). |

### Presentation & Documentation Materials

| File | Description |
|---|---|
| `Statistics_ppt_bioinformatics.pptx` | 8-slide work-in-progress presentation. |
| `Speech_Work_in_Progress_Presentation.docx` | Timed spoken script (~8 minutes) mapped slide by slide. |
| `Gene_Expression_Analysis_Summary_EN.docx` | Detailed exploratory data analysis summary and interpretation guide in English. |

## Preliminary results (exploratory, no formal tests yet)

- The proliferation genes MKI67, CCNB1 and CDK1 are clearly higher in tumors (mean log2 expression 0.89 to 2.96, 2.26 to 4.38 and 1.23 to 3.50, respectively), with non-overlapping boxes.
- EGFR is lower in tumors (4.17 to 1.80); TP53 is almost unchanged (3.77 vs 3.88).
- ESR1, FOXA1 and PGR are much more variable in tumors. ESR1 and FOXA1 have two peaks in their distribution and PGR has a large group of samples near zero, which suggests hidden tumor subgroups.
- Within tumors, genes cluster by function: the proliferation genes are strongly correlated (r between 0.76 and 0.82), as are ESR1 and FOXA1 (r = 0.78); the two groups are negatively correlated (MKI67 with ESR1, FOXA1 and PGR: r between -0.31 and -0.35).
- Correlations with age are weak (the largest is ESR1, r = 0.26).
- MKI67 increases with BRCA1 (r = 0.47) and shows no visible relation with TP53, ERBB2 or age.

## How to reproduce

1. Install R and the required packages (`data.table` and `R.utils`).
2. Download the three raw files (`TCGA-BRCA.star_fpkm.tsv.gz`, `TCGA-BRCA.clinical.tsv.gz`, and `gencode.v36.annotation.gtf.gene.probemap`) into your working directory.
3. Run `work-in-progress_R-scripts_en.R` to perform data processing, generate descriptive tables, and export all figures.

## Future work

**Inference (until 30 December)**

- Confidence intervals and t-tests for tumor vs normal: a paired t-test on the 113 matched pairs, and a two-sample Welch test on all samples without assuming equal variances.
- Multiple-testing correction across genes.
- Multiple linear regression with MKI67 as the outcome and BRCA1, TP53, ERBB2, ESR1, FOXA1, EGFR, PGR and age as predictors (tumor samples). Check for multicollinearity (VIF) and evaluate model residuals and R².

**Molecular subtypes (optional extension)**

- Motivation: bimodal distributions of ESR1, FOXA1 and PGR point to tumor mixtures that PAM50 molecular subtypes can clarify.
- Timeline: work in progress presented on 16 November; inferential analyses by 30 December; final presentation in January.

## Notes and limitations

- Expression values are on a log2 scale ($\log_2(\text{FPKM}+1)$).
- Normal samples are fewer than tumor samples (113 vs 1,095).
- Outliers are retained as true biological signals (e.g., ERBB2 amplification).

## Data source

Data generated by the TCGA Research Network (https://www.cancer.gov/tcga), distributed through the NCI Genomic Data Commons and the UCSC Xena platform.