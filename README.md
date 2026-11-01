# Spatial Multiomics

This README describes the **main** branch and the **v1.1** release. For environment setup, execution commands, reproduction/update modes, and validation, see [REPRODUCIBILITY.md](REPRODUCIBILITY.md). For manuscript figure exports, see [README_Figures.md](README_Figures.md).

## Overview

This repository contains the full analysis pipeline for the Spatial Multiomics study, including:

* data preprocessing
* batch correction
* sample selection
* clustering analysis
* spatial sparse CCA (spCCA)
* figure generation for the manuscript

## Data

External inputs are available from their own citable data records:

* **Raw metabolomics data and spatial flatmap grids:** Zenodo, [doi:10.5281/zenodo.23857927](https://doi.org/10.5281/zenodo.23857927). From the repository root, run:

  ```sh
  Rscript --vanilla Code/get_raw_data.R
  ```

  The script downloads the A26/A27 Q500 workbooks, the 5583 FIA/UPLC Q500 workbooks, and both flatmap grids into `Data/`, verifies their checksums, and restores the filenames expected by the analysis scripts. The grids are deposited without animal names and restored locally as `A26_Bonnie_Flatmap_Grid_20260315.xlsx` and `A27_Frizzle_Flatmap_Grid_20260317.xlsx`. These downloaded raw inputs are ignored by Git.

* **Proteomics quantification files:** Duke Research Data Repository, [doi:10.7924/r4r506](https://doi.org/10.7924/r4r506). The external files `A26 Proteomics with QC variables.txt` and `A27_SNE_combine_protein.tsv` are required for the initial proteomics reformatting step (P1). See [README_DukeLibrary.md](README_DukeLibrary.md) for input details, including the additional A27 annotation input and legacy input paths.

Processed tables and analysis objects remain in `Data/`.

The license associated with the code record covers the code only. Data are distributed under the terms of their respective data records and must be cited using those records, separately from the code record.

---

## Full Pipeline and Figure Reproduction

To fully reproduce the results and manuscript figures, users should follow the complete pipeline:

➡️ Full workflow diagram:
[DIAGRAM.md](DIAGRAM.md)

➡️ Figures and outputs:
[README_Figures.md](README_Figures.md)

### Key dependency structure

* Large protein level quantification proteomics data (Duke Library) → required for **P1 reformatting**
* Pipeline (`DIAGRAM.md`) → defines full workflow
* This repository → contains all downstream analysis and figure generation

Users should run the pipeline **sequentially** to ensure all intermediate data objects are correctly generated.

---

# A26 Analysis Pipeline

## Metabolomics

### M1. Data preprocessing and Batch Correction

* Input: `A26-Q500 Data.xlsx`
* Output: `PostCombat_A26_Metabolite_beforeHeldout.xlsx`
  (403 samples × 328 variables)
* Code: `Analysis_BatchEffect_A26_met_adj.Rmd`

Steps:

* Remove samples 244, 96, 188
* Missingness filtering (TG >70%, others >85%)
* Imputation: half the minimum observed intensity
* Dilution factor adjustment
* Remove 3 SPQC outliers
* Log transform
* Subtract SPQC
* ComBat batch correction
* PCA analysis (before heldout selection)

---

## Proteomics

### P1. Data Reformatting (requires Duke Library data)

* Input: `A26 Proteomics with QC variables.txt`
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
* Imputation: for each protein, replace missing log2 values with the minimum observed log2 value minus 1 (`min(x, na.rm = TRUE) - 1`), equivalent to half the minimum observed intensity on the original scale
* Dilution factor
* ComBat batch correction
* PCA analysis (before heldout selection)

---

## Combined Analysis

### C3. Heldout Preprocessing

Compare eligible mapped and heldout candidates and record the selected `finalID` for each spatial location.

* Output: `Data/Heldout_Updated_a26.xlsx` (267 samples selected)
* Code: `Code/A26/HeldoutPreprocessing.R`
* Run from the repository root: `source("Code/A26/HeldoutPreprocessing.R")`

Key idea:

* Robust PCA-based selection using MAD-normalized distance
* Select samples closest to robust center

---

### C4. Sample Selection

Filter both metabolomics and proteomics to the `finalID` selections recorded in C3.

* Output: `Data/a26.RData`

  * metabolomics: 267 × 328
  * proteomics: 267 × 9990
* Code: `CombineClustering.R`

---

### C5. Clustering Analysis

* Output: `Data/A26_PChclust.RData`
* Code: C5 section of `Code/A26/CombineClustering.R`
* Input: `Data/a26.RData` from C4, with aligned samples across modalities.
* Settings: Spearman correlations, average linkage, and a 70% PC1 variance threshold. Constant protein variables are excluded before clustering.

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

* Input: `A27-Q500 Data.xlsx`
* Output: `PostCombat_a27_Metabolite_beforeHeldout.xlsx`
  (548 × 300)
* Code: `Analysis_BatchEffect_A27_met_adj.Rmd`

Steps:

* Remove the 63 commercial plasma QC Level 1–3 rows before feature filtering and imputation; retain animal-specific SPQC for normalization
* Missingness filtering
* Remove outlier 212
* Imputation: half the minimum observed intensity
* Dilution factor
* Log transform
* Subtract SPQC
* ComBat
* PCA analysis (before heldout selection)

---

## Proteomics

### P1. Data Reformatting (requires Duke Library data)

* Input:
  `A27_SNE_combine_protein.tsv`

* Output:
  `Adjusted_Protein_A27.RData`
  (566 samples × 10,253 variables)

* Code: `Reformatting_A27_pro.R`
* Additional required input: `20260120_122736_Identifications.tsv` for plate and SPQC annotations; see [README_DukeLibrary.md](README_DukeLibrary.md).

This step converts **long-format proteomics data** into a **log2-transformed wide-format dataset**, which serves as the starting point for all downstream proteomics analyses.

---

### P2. Data preprocessing and Batch Correction

* Input: `Adjusted_Protein_A27.RData`
* Output: `PostCombat_a27_Protein.xlsx`
  (566 × 10,253)
* Code: `Analysis_BatchEffect_A27_pro_adj.Rmd`

Steps:

* Imputation: for each protein, replace missing log2 values with the minimum observed log2 value minus 1 (`min(x, na.rm = TRUE) - 1`), equivalent to half the minimum observed intensity on the original scale
* Dilution factor
* ComBat
* PCA analysis (before heldout selection)

---

## Combined Analysis

### C3. Heldout Preprocessing

Compare eligible mapped and heldout candidates and record the selected `finalID` for each spatial location.

* Output: `Data/Heldout_Updated_a27.xlsx` (527 samples)
* Code: `Code/A27/HeldoutPreprocessing.R`
* Run from the repository root: `source("Code/A27/HeldoutPreprocessing.R")`
* Metabolomics PCA uses all 300 retained analytes, excluding only the two metadata columns.

---

### C4. Sample Selection

Filter both metabolomics and proteomics to the `finalID` selections recorded in C3.

* Output: `Data/a27.RData`

  * metabolomics: 527 × 300
  * proteomics: 527 × 10,253
* Code: `CombineClustering.R`

---

### C5. Clustering Analysis

* Output: `Data/A27_PChclust.RData`
* Code: C5 section of `Code/A27/CombineClustering.R`
* Input: `Data/a27.RData` from C4, with aligned samples across modalities.
* Settings: Spearman correlations, average linkage, and a 70% PC1 variance threshold. Constant protein variables are excluded before clustering.

---

### C6. Spatial Sparse CCA

* Code: `A27_analysis_sCCA_laplacians.Rmd`
* Result: 508 matched samples

---

## Figure exports and QCFix corrections

See [QCFIX_VALIDATION.md](QCFIX_VALIDATION.md) for verified before/after counts and the correction record.

See [README_Figures.md](README_Figures.md) for the final TIFF workflow and [REPRODUCIBILITY.md](REPRODUCIBILITY.md#regenerate-the-coordinated-qcfix-corrections) for the coordinated rerun command. Figures 3, S3, and S4 are assembled from component PNGs; Figures 4, 5, and S5 are exported directly as TIFFs. Figures 4 and S5 both use loading-weighted protein clusters; Figure 4 maps sample-level neighbor means, and S5 displays the distribution of those same per-sample means.

Figure 5 now removes only flagged outlier rows, using row identity rather than sample numbers shared across animals. Its later PCA/sCCA loading comparisons harmonize metabolite names before matching. QC3 participation in Figure 5 filtering, imputation, and normalization is retained.

## Reproducibility Summary

To fully reproduce the analysis:

1. Obtain raw proteomics files from Duke Library
2. Run P1 (reformatting scripts)
3. Run preprocessing (M1, P2)
4. Run combined analysis (C3–C6)
5. Follow `DIAGRAM.md` for full pipeline
6. Generate figures using scripts referenced in repository `README_Figures.md`

Users interested only in reproducing figures may use processed data included in this repository without rerunning the raw data reformatting step.
