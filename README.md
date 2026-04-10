# Spatial Multiomics

## Overview

This repository contains the full analysis pipeline for the Spatial Multiomics study, including:

* data preprocessing
* batch correction
* sample selection
* clustering analysis
* spatial sparse CCA (spCCA)
* figure generation for the manuscript

⚠️ **Important:** Two large raw proteomics input files are **not included in this repository** due to file size limitations. They are deposited separately in the Duke Library:

* `10788_SNE_combine_final_012026_Report_protein_export (Normal).tsv`
* `Bonnie Proteomics with QC variables.txt`

These files are **only required for the initial reformatting step** (P1 below).

---

## Full Pipeline and Figure Reproduction

To fully reproduce the results and manuscript figures, users should follow the complete pipeline:

➡️ Full workflow diagram:
https://github.com/QiuyiWu/Spatial-Multiomics/blob/main/DIAGRAM.md

➡️ Figures and outputs:
https://github.com/QiuyiWu/Spatial-Multiomics/blob/main/README_Figures.md

### Key dependency structure

* Raw data (Duke Library) → required for **P1 reformatting**
* Pipeline (`DIAGRAM.md`) → defines full workflow
* This repository → contains all downstream analysis and figure generation

Users should run the pipeline **sequentially** to ensure all intermediate data objects are correctly generated.

---

# A26 Analysis Pipeline

## Metabolomics

### M1. Data preprocessing and Batch Correction

* Input: `10412-Q500 Data.xlsx`
* Output: `PostCombat_A26_Metabolite_beforeHeldout.xlsx`
  (403 samples × 328 variables)
* Code: `Analysis_BatchEffect_A26_met_adj.Rmd`

Steps:

* Remove samples 244, 96, 188
* Missingness filtering (TG >70%, others >85%)
* LOD/2 imputation
* Dilution factor adjustment
* Remove 3 SPQC outliers
* Log transform
* Subtract SPQC
* ComBat batch correction
* PCA analysis (exclude HeldOut_AnalysisIDs)

---

## Proteomics

### P1. Data Reformatting (requires Duke Library data)

* Input: `Bonnie Proteomics with QC variables.txt`
* Output: `Adjusted_Protein_A26.RData`
  (428 samples × 9990 variables)
* Code: `Reformatting_A26_pro.R`

This step converts **log2-transformed raw long-format proteomics data** into a **wide-format dataset** used for downstream analysis.

---

### P2. Data preprocessing and Batch Correction

* Input: `Adjusted_Protein_A26.RData`
* Output: `PostCombat_A26_Protein.xlsx`
  (427 samples × 9990 variables)
* Code: `Analysis_BatchEffect_A26_pro_adj.Rmd`

Steps:

* Remove outlier `ID117530`
* LOD/2 imputation
* Dilution factor
* ComBat batch correction
* PCA analysis

---

## Combined Analysis

### C3. Heldout Preprocessing

* Output: `Heldout_Updated.xlsx` (267 samples selected)
* Code: `HeldoutPreprocessing.R`

Key idea:

* Robust PCA-based selection using MAD-normalized distance
* Select samples closest to robust center

---

### C4. Sample Selection

* Output: `a26.RData`

  * metabolomics: 267 × 328
  * proteomics: 267 × 9990
* Code: `CombineClustering.R`

---

### C5. Clustering Analysis

* Output: `A26_PChclust.RData`
* Code: `CombineClustering.R`

Contains:

* group
* score
* loading
* cluster summary
* variance explained

---

### C6. Spatial Sparse CCA

* Code: `A26_analysis_sCCA_laplacians.Rmd`
* Result: 233 matched samples

---

# A27 Analysis Pipeline

## Metabolomics

### M1. Data preprocessing and Batch Correction

* Input: `10788-Q500 Data.xlsx`
* Output: `PostCombat_A27_Metabolite_beforeHeldout.xlsx`
  (548 × 343)
* Code: `Analysis_BatchEffect_A27_met_adj.Rmd`

Steps:

* Missingness filtering
* Remove outlier 212
* LOD/2 imputation
* Dilution factor
* Log transform
* Subtract SPQC
* ComBat
* PCA

---

## Proteomics

### P1. Data Reformatting (requires Duke Library data)

* Input:
  `10788_SNE_combine_final_012026_Report_protein_export (Normal).tsv`

* Output:
  `Adjusted_Protein_A27.RData`
  (566 samples × 10,253 variables)

* Code: `Reformatting_A27_pro.R`

This step converts **log2-transformed long-format proteomics data** into a **wide-format dataset**, which serves as the starting point for all downstream proteomics analyses.

---

### P2. Data preprocessing and Batch Correction

* Input: `Adjusted_Protein_A27.RData`
* Output: `PostCombat_A27_Protein.xlsx`
  (566 × 10,253)
* Code: `Analysis_BatchEffect_A27_pro_adj.Rmd`

Steps:

* LOD/2 imputation
* Dilution factor
* ComBat
* PCA

---

## Combined Analysis

### C3. Heldout Preprocessing

* Output: `Heldout_Updated_a27.xlsx` (527 samples)
* Code: `HeldoutPreprocessing.R`

---

### C4. Sample Selection

* Output: `a27.RData`

  * metabolomics: 527 × 342
  * proteomics: 527 × 10,253
* Code: `CombineClustering.R`

---

### C5. Clustering Analysis

* Output: `A27_PChclust.RData`
* Code: `CombineClustering.R`

---

### C6. Spatial Sparse CCA

* Code: `A27_analysis_sCCA_laplacians.Rmd`
* Result: 508 matched samples

---

## Reproducibility Summary

To fully reproduce the analysis:

1. Obtain raw proteomics files from Duke Library
2. Run P1 (reformatting scripts)
3. Run preprocessing (M1, P2)
4. Run combined analysis (C3–C6)
5. Follow `DIAGRAM.md` for full pipeline
6. Generate figures using scripts referenced in repository `README_Figures.md`

Users interested only in reproducing figures may use processed data included in this repository without rerunning the raw data reformatting step.
