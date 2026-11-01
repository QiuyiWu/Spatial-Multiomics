# Reproducibility guide — QCFix

Analysis code for *Integrated metabolomics and proteomics from voxelated cortical hemispheres of adult rhesus monkeys*, by Qiuyi Wu, Alev M. Brigande, Michael W. Lutz, Pixu Shi, and Anita A. Disney.

This branch contains corrected proteomics imputation and downstream A26/A27 analyses. The study integrates metabolomics and proteomics, selects samples for spatial mapping, clusters analytes, and fits spatially regularized sparse canonical correlation models (sr-sCCA).

For the study overview and animal-by-animal pipeline, see [README.md](README.md).

## Start here

- [Workflow and dependencies](DIAGRAM.md)
- [Figure and supplementary-table reproduction](README_Figures.md)
- [External proteomics inputs](README_DukeLibrary.md)

Run analyses from the repository root. Files under `Data/` are tracked with Git LFS; a pointer file is not an analysis-ready input.

```sh
git lfs install
git clone --branch QCFix https://github.com/QiuyiWu/Spatial-Multiomics.git
cd Spatial-Multiomics
git lfs pull
```

### R environment

`renv.lock` records R 4.5.2 and package versions. Start R in the repository root and restore the project library:

```r
renv::restore()
```

The repository's `.Rprofile` activates `renv`. Install `renv` first if it is unavailable. Use the recorded R version where possible; operating-system dependencies for spatial and compiled packages may also be required. The lockfile includes `adegenet`, the runtime-required `adespatial`, and their dependencies for optional C6 graph rebuilding.

The commands below assume the restored project environment. `Rscript --vanilla` bypasses `.Rprofile`; use it only if the intended packages have been loaded by another explicit setup.

## Regenerate the coordinated QCFix corrections

After restoring the R environment, run from the repository root with fresh, empty run and figure directories:

```sh
Rscript Code/regenerate_QCFix.R /absolute/path/to/fresh-run /absolute/path/to/CorrectedFigure2
```

This reruns A27 metabolomics preprocessing, heldout selection, aligned samples, and clustering; refits the A26/A27 spatial models; regenerates the projection, PCA figures, Figures 4/5/S5, and final component TIFF composites. It writes regenerated data under the run directory's `Data/`, stage logs under `logs/`, and the complete figure set to the supplied figure directory. Review `comparison.csv` before replacing tracked references. Saved graphs are validated; explicitly request rebuilding only if sample/coordinate changes require it.

The correction excludes 63 A27 commercial plasma QC Level 1–3 rows before filtering/imputation in the SPQC/ComBat workflow. Figure 5 removes exactly its flagged rows, matches metabolite names consistently in the loading comparisons, and **retains QC3 in its own filtering, imputation, and normalization**. Original labels and matching/outlier audits are exported.

The recorded local validation uses R 4.4.3 with an independently configured system library, invoked explicitly with `--vanilla`; `renv.lock` remains the R 4.5.2 reference environment. For that explicit alternative, set `PROTEINFIX_RSCRIPT_FLAGS=--vanilla` when calling the coordinated runner. Session information accompanies the outputs.

## Reproduction mode and reference updates

The pipeline defaults to **reproduction mode**. Generated data artifacts go to
`outputs/reproduction/Data/`, leaving committed `Data/` files intact. This applies
to P1 output objects, M1/P2 workbooks and PCA objects, C3 selections, C4/C5 objects,
C6 fitted weights, optional rebuilt graphs, and cross-animal projection results.
Figures retain their documented figure/output directories.

`Code/pipeline_paths.R` routes downstream reads to a regenerated file in the same
run folder when present, otherwise to the committed `Data/` reference. Every
resolved input and output path is printed. This permits starting from any saved
stage; for an independent full rerun, use a new empty run folder and execute the
stages in order so old run files cannot be reused accidentally.

Configure one R session with:

```r
options(proteinfix.mode = "reproduce",
        proteinfix.run_dir = here::here("outputs", "reproduction"))
```

For commands in separate R processes, share the configuration through environment
variables (use an absolute run path):

```sh
export PROTEINFIX_MODE=reproduce
export PROTEINFIX_RUN_DIR=/absolute/path/to/reproduction-run
```

C6 **loads saved spatial graphs by default** and checks sample IDs, order,
coordinates, and adjacency validity. An incompatible graph stops the analysis;
it is not silently reused or rebuilt. To deliberately rebuild:

```r
options(proteinfix.rebuild_graph = TRUE)
```

or set `PROTEINFIX_REBUILD_GRAPH=true` for command-line runs. Rebuilt graphs follow
the same output routing. `adegenet` and `adespatial` are needed for this optional
rebuild, not the default graph-loading step.

After regeneration, compare data contents with committed references:

```sh
Rscript Code/compare_reproduction.R /absolute/path/to/reproduction-run
```

The report is written to the run folder's `comparison.csv`. It compares workbook
sheets and values, R data objects, and graph IDs/coordinates/neighbors with numeric
tolerance 1e-8. Binary serialization and Excel timestamps are not used as evidence
of a numerical change. The report covers files present in that run folder; it does
not certify that every pipeline stage was run. Investigate differences before
updating reference data; fitted-component sign changes can require interpretation.

**Update mode is explicit and does replace reference files.** Set
`options(proteinfix.mode = "update")` or `PROTEINFIX_MODE=update`, then rerun the
intended stage. Update-mode reads and writes use `Data/`. Reset to reproduction
mode afterward and review the resulting Git diff before committing. No existing
references are promoted automatically.

To check routing, C3/C4 results, graph validation/rebuilding, and reference-file
preservation in a temporary folder:

```sh
Rscript Code/tests/test_reproduction_mode.R
```

## Inputs and execution order

Processed inputs are included under `Data/` through Git LFS. To reproduce downstream analyses, start at the appropriate saved stage. To reproduce proteomics from the long-format source files, complete P1 first; see the external-input notes below.

### P1. Reformat proteomics, if starting from long-format inputs

| Animal | Script | Input | Output object |
| --- | --- | --- | --- |
| A26 | [Reformatting_A26_pro.R](Code/A26/Reformatting_A26_pro.R) | `A26 Proteomics with QC variables.txt` | `pro.a26.adj` in `Data/Adjusted_Protein_A26.RData` |
| A27 | [Reformatting_A27_pro.R](Code/A27/Reformatting_A27_pro.R) | `A27_SNE_combine_protein.tsv` and `20260120_122736_Identifications.tsv` | `pro.a27.adj` in `Data/Adjusted_Protein_A27.RData` |

The P1 scripts retain legacy relative input paths. Set the input paths to your file locations before running them; save destinations follow the shared mode configuration. They also assume some functions are already attached; load `readr`, `dplyr`, and `tidyr` in a clean session. The A27 identification file supplies plate and SPQC annotations and is an additional required input; its availability must be checked separately from the two large quantification files.

A26 P1 produces 428 samples and 9,990 protein variables. A27 P1 produces 566 samples and 10,253 protein variables. Metadata columns are additional.

### M1 and P2. Preprocess and correct batches

| Animal / modality | Rmd under `Code/` | Post-ComBat workbook | Samples × analytes |
| --- | --- | --- | --- |
| A26 metabolomics | `A26/Analysis_BatchEffect_A26_met_adj.Rmd` | `Data/PostCombat_A26_Metabolite_beforeHeldout.xlsx` | 403 × 328 |
| A26 proteomics | `A26/Analysis_BatchEffect_A26_pro_adj.Rmd` | `Data/PostCombat_A26_Protein.xlsx` | 427 × 9,990 |
| A27 metabolomics | `A27/Analysis_BatchEffect_A27_met_adj.Rmd` | `Data/PostCombat_a27_Metabolite_beforeHeldout.xlsx` | 548 × 300 |
| A27 proteomics | `A27/Analysis_BatchEffect_A27_pro_adj.Rmd` | `Data/PostCombat_a27_Protein.xlsx` | 566 × 10,253 |

File-name capitalization in this table matches the repository. Dimensions exclude metadata columns. Downstream protein analyses read the `combat_corrected_unscaled` worksheet. The alternative `combat_corrected_scaled` worksheet applies plate-wise SPQC subtraction before ComBat; "scaled" here does not mean z-scoring. A27 assigns adjusted values back to their original row positions and checks sample keys against batch labels before ComBat and export, including distinct row identities for repeated SPQC labels. Run `Rscript Code/tests/test_protein_spqc_order.R` to check interleaved plates, repeated SPQC labels, metadata alignment, and missing-QC handling.

For metabolomics, the scripts apply missingness thresholds of 70% for TG analytes and 85% for other analytes, imputation, dilution adjustment, log2 transformation, SPQC adjustment, and ComBat. A26 excludes samples 244, 96, and 188; A27 excludes sample 212. QC exclusions account for differences between pre- and post-ComBat sample counts. A26 proteomics removes outlier `ID117530`.

PCA exports at this stage precede heldout selection. Each Rmd saves a reproducible pre/post matrix pair with matching sample metadata:

| Saved PCA object | Pre-ComBat | Post-ComBat |
| --- | --- | --- |
| `Data/met_a26_PCAdata.RData` | 433 × 328 | 403 × 328 |
| `Data/pro_a26_PCAdata.RData` | 427 × 9,990 | 427 × 9,990 |
| `Data/met_a27_PCAdata.RData` | 581 × 300 | 548 × 300 |
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

Outputs are `Data/Heldout_Updated_a26.xlsx` (267 selected samples) and `Data/Heldout_Updated_a27.xlsx` (527). A27 metabolomics uses all 300 retained analytes, excluding its two metadata columns.

### C4–C5. Select samples and cluster analytes

```r
source("Code/A26/CombineClustering.R")
source("Code/A27/CombineClustering.R")
```

Each script performs C4 and C5 sequentially.

| Animal | C4 output | Selected metabolomics | Selected proteomics | C5 output |
| --- | --- | --- | --- | --- |
| A26 | `Data/a26.RData` | 267 × 328 | 267 × 9,990 | `Data/A26_PChclust.RData` |
| A27 | `Data/a27.RData` | 527 × 300 | 527 × 10,253 | `Data/A27_PChclust.RData` |

C5 uses Spearman correlations, average linkage, and a 70% PC1 variance threshold. Constant protein columns are excluded before clustering. The saved cluster objects contain group assignments, scores, loadings, cluster summaries, and variance explained.

### C6. Fit spatial models and export component maps

```sh
Rscript Code/regenerate_C6_outputs.R A26 /absolute/path/to/CorrectedFigure
Rscript Code/regenerate_C6_outputs.R A27 /absolute/path/to/CorrectedFigure
```

The driver runs the spatial-analysis chunks of the corresponding `*_analysis_sCCA_laplacians.Rmd`, loading and validating the saved Delaunay graph by default (or rebuilding it when explicitly requested). It uses seed 123, 25 tuning permutations, and 10 components. Spatially matched sample counts are 233 for A26 and 508 for A27.

The driver saves spatial-model weights, and graphs only when explicitly rebuilt, to the configured data-output directory, exports paired component maps into `A26_components/` and `A27_components/`, and writes diagnostics under `Figure/C6_diagnostics/`. The A26 run also exports Table S3. It does not run the separate nonspatial CCA comparison. See [figure reproduction](README_Figures.md) for output names, presentation conventions, and supplementary-table annotations.

## Manuscript figures

Figure 2 and S2 use the four pre/post PCA objects. Figure 3 and the additional component maps use C6 and the cross-animal projection. Figure 4 and S5 use corrected exact-distance neighbor groups and the same loading-weighted protein-cluster representation. Figure 4 maps each sample’s average correlation with its neighbors; S5 displays the distribution of those same per-sample averages for each neighbor order. Figure 5 is a separate **metabolomics-only** U5/A26/A27 PCA; it does not use proteomics imputation.

The manuscript export workflow distinguishes component PNGs from final TIFF figures:

| Figure | Source and final export |
| --- | --- |
| 3 | Select component PNGs from `A26_components/`, `A27_components/`, and `A27_to_A26_components/` (10 PNGs per folder), then assemble and save `Figure3.tiff`. |
| S3 | Assemble all 10 component PNGs from `A26_components/` and save `Figure S3.tiff`. |
| S4 | Assemble all 10 component PNGs from `A27_components/` and save `Figure S4.tiff`. |
| 4 | `Code/regenerate_Figure4.R` directly exports `Figure4_corrected.tiff` at 600 dpi, as well as a PNG copy. |
| 5 | `Code/loadingComp_U5A26A27_met.Rmd` directly exports `Figure/fig5_pca.tiff` at 600 dpi. |
| S5 | `Code/regenerate_FigureS5.R` directly exports `FigureS5_corrected.tiff` at 600 dpi, as well as a PNG copy. |

`Code/sr_sCCA_figures.R`, invoked by `Code/regenerate_C6_outputs.R`, produces the A26/A27 paired component PNGs; it does not itself assemble the final TIFFs for Figures 3, S3, or S4. Run `Code/assemble_manuscript_components.R` for those final layouts. The cross-animal Rmd produces the projection component PNGs.

The A26/A27 C6 Rmds also contain separate 600-dpi TIFF exports for individual U (protein) and V (metabolite) maps (`a26_U*.tiff`, `a26_V*.tiff`, `a27_U*.tiff`, and `a27_V*.tiff`); the projection Rmd exports `proj_U*.tiff` and `proj_V*.tiff`. These are additional panel exports, not the PNG-based assembly workflow above. The C6 driver does not execute the Rmds' individual-map export blocks. The 600-dpi setting for those TIFF panels does not describe the resolution of the PNG-based final composites. See [README_Figures.md](README_Figures.md) for commands, output locations, and assembly details.

`Figure/` and exported `Supp_Table_S3*.xlsx` workbooks are ignored by Git. Generated files in an external output folder are not automatically published to this branch. Follow [README_Figures.md](README_Figures.md) for the current reproducible exports and the remaining presentation steps.
