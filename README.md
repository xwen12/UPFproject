# Longitudinal plasma lipidomic signatures of ultra-processed food intake and incident cardiovascular disease: the Strong Heart Family Study

Analysis code for the manuscript above.

## Data

This repository contains code only. Strong Heart Study (SHS) data are not public. They are available from the
SHS upon approval, subject to SHS and Tribal data-sharing policies. No individual-level data, results or
participant identifiers are included in this repository.

## Setup

Each script starts with two folder paths. Set them before running:

```r
proj_dir <- "path/to/proj16_UPF_rework"   # project folder with data/, result/, result/tabsFigs/
raw_dir  <- "path/to/SHS_raw_data"        # SHS raw data (incl. "lipidomics data/")
```

Analyses were run in R 4.5.1. R packages: dplyr, tidyr, tidyverse, stringr, tibble, lubridate, haven, readxl,
openxlsx, officer, flextable, lme4, lmerTest, r2glmm, survival, mediation, sandwich, geepack, multgee, mice,
gmodels, ggplot2, ggrepel, ggtext, ggbeeswarm, patchwork, scales, gtable, ragg.
Figure 1 also uses Node.js (pptxgenjs, react, react-dom, react-icons, sharp) and Python (numpy, matplotlib).

## Run order and outputs

Scripts are numbered; run them in this order.

| Step | Script | Manuscript item |
|---|---|---|
| Data cleaning | 1.1, 1.2, 1.3 | Analytic datasets |
| Sample selection | 1.4 (run in the same R session right after 1.3) | Figure 1a, Supplemental Figure S1 |
| Food items | 1.5 | Supplemental Table S1 |
| Study overview | 13Fig1_lipidNumbers.R, 1.7, 1.6 | Figure 1 |
| Table 1 | 4.1 table1_baseline_followup.R | Table 1 |
| Aim 1: UPF and incident CVD outcomes | 2.0.1, 2.0.2 | Figure 2b-c |
| | 2.1.3, 2.1.4, 2.1.5 | Figure 2 |
| | 2.0.2.1 | Supplemental Table S3 |
| Aim 1 sensitivity: diet quality (AHEI) | 2.2, 2.2.1, 2.3 | Supplemental Table S2 |
| Aim 2: UPF and lipids | 2.1 | Figure 3, Supplemental Table S4 |
| | 6.2.1, 6.1.1, 6.3.1, 6.3.2, 6.4 | Figure 3a-c |
| | 6.2.2 | Supplemental Table S4 |
| Aim 2 sensitivity: center split (AZ + SD discovery, OK replication) | 2.4, 2.4.1, 2.5 | Supplemental Table S5 |
| Aim 2 pathways | 12.pathway_metaboanalyst.R (lists uploaded to MetaboAnalyst 6.0 and LIPEA), 12.2 | Figure 5 |
| Aim 3: baseline lipids and incident outcomes | 3.1, 3.2, 3.3 | Figure 4a, Supplemental Table S6 |
| Aim 3: 5-year lipid change and incident outcomes | 4.0, 4.1 lipidChange_cvd_cox.R, 4.2, 4.3 | Figure 4a, Supplemental Table S7 |
| | 8.1, 8.1.1 | Figure 4a, Supplemental Tables S6-S7 |
| Aim 3: mediation | 5.6 | Figure 4b, Supplemental Table S8 |
| | 8.5, 8.5.1, 8.4 | Figure 4b, Supplemental Table S8, Figure 4 |

## Contact

Xiaoxiao Wen, University of South Florida.
