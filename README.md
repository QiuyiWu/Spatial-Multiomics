# Spatial Multiomics — ProteinFix

Analysis code for *Integrated metabolomics and proteomics from voxelated cortical hemispheres of adult rhesus monkeys*, by Qiuyi Wu, Alev M. Brigande, Michael W. Lutz, Pixu Shi, and Anita A. Disney.

This branch contains corrected proteomics imputation and downstream A26/A27 analyses. The study integrates metabolomics and proteomics, selects samples for spatial mapping, clusters analytes, and fits spatially regularized sparse canonical correlation models (sr-sCCA).

## Start here

- [Workflow and dependencies](DIAGRAM.md)
- [Figure and supplementary-table reproduction](README_Figures.md)
- [External proteomics inputs](README_DukeLibrary.md)

Run analyses from the repository root. Files under `Data/` are tracked with Git LFS; a pointer file is not an analysis-ready input.

```sh
git lfs install
git clone --branch ProteinFix https://github.com/QiuyiWu/Spatial-Multiomics.git
cd Spatial-Multiomics
git lfs pull
```

### R environment

`renv.lock` records R 4.5.2 and package versions. Start R in the repository root and restore the project library:

```r
renv::restore()
```

The repository's `.Rprofile` activates `renv`. Install `renv` first if it is unavailable. Use the recorded R version where possible; operating-system dependencies for spatial and compiled packages may also be required. The lockfile includes `adegenet`, the runtime-required `adespatial`, and their dependencies, needed to rebuild the C6 sPCA graphs.

The commands below assume the restored project environment. `Rscript --vanilla` bypasses `.Rprofile`; use it only if the intended packages have been loaded by another explicit setup.

## What ProteinFix changes

Both reformatted proteomics inputs are already on the log2 scale when P2 starts. A26's input supplies log2 values; A27's P1 script applies `log2()` to its intensity column.

P2 replaces missing values within each protein with:

```r
min(x, na.rm = TRUE) - 1
```

Subtracting one on the log2 scale corresponds to **half the minimum observed intensity** on the original intensity scale. Dividing the log2 value by two does not implement that rule. The dilution adjustment follows imputation, then ordinary ComBat is applied without SPQC subtraction for the protein outputs used downstream.

Metabolomics imputation occurs before its log2 transformation and uses half the minimum observed intensity on the original scale. This is an observed-minimum rule, not a measured assay limit of detection.

## Inputs and execution order

Processed inputs are included under `Data/` through Git LFS. To reproduce downstream analyses, start at the appropriate saved stage. To reproduce proteomics from the long-format source files, complete P1 first; see the external-input notes below.

### P1. Reformat proteomics, if starting from long-format inputs

| Animal | Script | Input | Output object |
| --- | --- | --- | --- |
| A26 | [Reformatting_A26_pro.R](Code/A26/Reformatting_A26_pro.R) | `A26 Proteomics with QC variables.txt` | `pro.a26.adj` in `Data/Adjusted_Protein_A26.RData` |
| A27 | [Reformatting_A27_pro.R](Code/A27/Reformatting_A27_pro.R) | `A27_SNE_combine_protein.tsv` and `20260120_122736_Identifications.tsv` | `pro.a27.adj` in `Data/Adjusted_Protein_A27.RData` |

The P1 scripts retain legacy relative input/output paths. Set these paths to your input locations and the repository's `Data/` directory before running them. They also assume some functions are already attached; load `readr`, `dplyr`, and `tidyr` in a clean session. The A27 identification file supplies plate and SPQC annotations and is an additional required input; its availability must be checked separately from the two large quantification files.

A26 P1 produces 428 samples and 9,990 protein variables. A27 P1 produces 566 samples and 10,253 protein variables. Metadata columns are additional.

### M1 and P2. Preprocess and correct batches

| Animal / modality | Rmd under `Code/` | Post-ComBat workbook | Samples × analytes |
| --- | --- | --- | --- |
| A26 metabolomics | `A26/Analysis_BatchEffect_A26_met_adj.Rmd` | `Data/PostCombat_A26_Metabolite_beforeHeldout.xlsx` | 403 × 328 |
| A26 proteomics | `A26/Analysis_BatchEffect_A26_pro_adj.Rmd` | `Data/PostCombat_A26_Protein.xlsx` | 427 × 9,990 |
| A27 metabolomics | `A27/Analysis_BatchEffect_A27_met_adj.Rmd` | `Data/PostCombat_a27_Metabolite_beforeHeldout.xlsx` | 548 × 343 |
| A27 proteomics | `A27/Analysis_BatchEffect_A27_pro_adj.Rmd` | `Data/PostCombat_a27_Protein.xlsx` | 566 × 10,253 |

File-name capitalization in this table matches the repository. Dimensions exclude metadata columns. Downstream protein analyses read the `combat_corrected_unscaled` worksheet.

For metabolomics, the scripts apply missingness thresholds of 70% for TG analytes and 85% for other analytes, imputation, dilution adjustment, log2 transformation, SPQC adjustment, and ComBat. A26 excludes samples 244, 96, and 188; A27 excludes sample 212. QC exclusions account for differences between pre- and post-ComBat sample counts. A26 proteomics removes outlier `ID117530`.

PCA exports at this stage precede heldout selection. Each Rmd saves a reproducible pre/post matrix pair with matching sample metadata:

| Saved PCA object | Pre-ComBat | Post-ComBat |
| --- | --- | --- |
| `Data/met_a26_PCAdata.RData` | 433 × 328 | 403 × 328 |
| `Data/pro_a26_PCAdata.RData` | 427 × 9,990 | 427 × 9,990 |
| `Data/met_a27_PCAdata.RData` | 644 × 343 | 548 × 343 |
| `Data/pro_a27_PCAdata.RData` | 566 × 10,253 | 566 × 10,253 |

For example, render an analysis from R in the repository root:

```r
rmarkdown::render(
  "Code/A26/Analysis_BatchEffect_A26_met_adj.Rmd",
  knit_root_dir = getwd()
)
```

### C3. Select mapped or heldout candidates

```r
source("Code/A26/HeldoutPreprocessing.R")
source("Code/A27/HeldoutPreprocessing.R")
```

The scripts compare candidate samples using robust distances in the first two PCA dimensions for each modality. The maximum normalized distance across available modalities is used to select the candidate with the smallest value.

Outputs are `Data/Heldout_Updated_a26.xlsx` (267 selected samples) and `Data/Heldout_Updated_a27.xlsx` (527). A27 metabolomics uses all 343 analytes, excluding its two metadata columns.

### C4–C5. Select samples and cluster analytes

```r
source("Code/A26/CombineClustering.R")
source("Code/A27/CombineClustering.R")
```

Each script performs C4 and C5 sequentially.

| Animal | C4 output | Selected metabolomics | Selected proteomics | C5 output |
| --- | --- | --- | --- | --- |
| A26 | `Data/a26.RData` | 267 × 328 | 267 × 9,990 | `Data/A26_PChclust.RData` |
| A27 | `Data/a27.RData` | 527 × 343 | 527 × 10,253 | `Data/A27_PChclust.RData` |

C5 uses Spearman correlations, average linkage, and a 70% PC1 variance threshold. Constant protein columns are excluded before clustering. The saved cluster objects contain group assignments, scores, loadings, cluster summaries, and variance explained.

### C6. Fit spatial models and export component maps

```sh
Rscript Code/regenerate_C6_outputs.R A26 /absolute/path/to/CorrectedFigure
Rscript Code/regenerate_C6_outputs.R A27 /absolute/path/to/CorrectedFigure
```

The driver runs the spatial-analysis chunks of the corresponding `*_analysis_sCCA_laplacians.Rmd`, rebuilding the Delaunay graph with `adegenet::spca()`. It uses seed 123, 25 tuning permutations, and 10 components. Spatially matched sample counts are 233 for A26 and 508 for A27.

The driver saves the current graphs and spatial-model weights under `Data/`, exports paired component maps into `A26_components/` and `A27_components/`, and writes diagnostics under `Figure/C6_diagnostics/`. The A26 run also exports Table S3. It does not run the separate nonspatial CCA comparison. See [figure reproduction](README_Figures.md) for output names, presentation conventions, and supplementary-table annotations.

## Manuscript figures

Figure 2 and S2 use the four pre/post PCA objects. Figure 3 and the additional component maps use C6 and the cross-animal projection. Figure 4 and S5 use corrected exact-distance neighbor groups. Figure 5 is a separate **metabolomics-only** U5/A26/A27 PCA; it does not use proteomics imputation.

`Figure/` and exported `Supp_Table_S3*.xlsx` workbooks are ignored by Git. Generated files in an external output folder are not automatically published to this branch. Follow [README_Figures.md](README_Figures.md) for the current reproducible exports and the remaining presentation steps.
