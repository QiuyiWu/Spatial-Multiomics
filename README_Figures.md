# Manuscript figures and supplementary Table S3

Use the **QCFix** branch with Git LFS inputs downloaded and the R environment restored. Run commands from the repository root. Data inputs/outputs follow the shared reproduction/update mode described in [REPRODUCIBILITY.md](REPRODUCIBILITY.md#reproduction-mode-and-reference-updates). [README.md](README.md) documents the analysis stages; [DIAGRAM.md](DIAGRAM.md) shows their dependencies.

The commands below use the activated project library. Adding `--vanilla` bypasses the repository's `renv` startup and requires a separate explicit package setup.

## Figure index

| Figure | Content | Code |
| --- | --- | --- |
| 2 | Post-ComBat PCA, A26/A27 × metabolomics/proteomics | `Code/PCA_Figure.R` |
| 3a–b | A26 sr-sCCA components 1–2 | `Code/A26/A26_analysis_sCCA_laplacians.Rmd` and `Code/sr_sCCA_figures.R` |
| 3c–d | A27 sr-sCCA components 1–2 | `Code/A27/A27_analysis_sCCA_laplacians.Rmd` and `Code/sr_sCCA_figures.R` |
| 3e–f | A27-derived components projected onto A26 | `Code/A26A27_analysis_sCCA.Rmd` |
| 4 | Mean correlation with exact-order spatial neighbors | `Code/regenerate_Figure4.R` and the two SampleSize Rmds |
| 5 | Metabolomics PCA: blood-containing U5 versus perfused A26/A27 | `Code/loadingComp_U5A26A27_met.Rmd` |
| S2 | Pre-ComBat PCA | `Code/PCA_Figure.R` |
| S3–S4 | Additional sr-sCCA component maps | C6 Rmds and `Code/sr_sCCA_figures.R` |
| S5 | Metabolite and protein-cluster neighbor-correlation boxplots | `Code/regenerate_FigureS5.R` and the two SampleSize Rmds |

Generated figures under `Figure/` and workbooks matching `Supp_Table_S3*.xlsx` are ignored by Git. Scripts that accept an output directory can write outside the repository. Such files are not automatically committed or published.

## Figures 2 and S2: pre/post-ComBat PCA

Run the four preprocessing Rmds under `Code/A26/` and `Code/A27/` if the saved PCA objects need to be rebuilt, then run:

```sh
Rscript Code/PCA_Figure.R
```

Each preprocessing Rmd saves, to the configured data-output directory, numeric `pre_*` and `post_*` matrices plus matching `pre_*_metadata` and `post_*_metadata` tables. Metadata include `sample_id`, `plate`, and `sample_type`; row names match their matrix.

| Data file | Pre-ComBat samples × analytes | Post-ComBat samples × analytes |
| --- | --- | --- |
| `Data/met_a26_PCAdata.RData` | 433 × 328 | 403 × 328 |
| `Data/pro_a26_PCAdata.RData` | 427 × 9,990 | 427 × 9,990 |
| `Data/met_a27_PCAdata.RData` | 644 × 343 | 548 × 343 |
| `Data/pro_a27_PCAdata.RData` | 566 × 10,253 | 566 × 10,253 |

These exports precede heldout selection. The metabolomics preprocessing excludes A26 samples 244, 96, and 188 and A27 sample 212. Log2 transformation occurs before batch correction. QC exclusions differ by stage, so pre/post sample counts need not match.

Outputs:

- `Figure/post_PCA_4panel_shared_legend.tiff` — Figure 2.
- `Figure/pre_PCA_4panel_shared_legend.tiff` — Figure S2.

Both contain a shared plate legend, Metabolomic/Proteomic column headings, A26/A27 row labels, and a–d panel letters.

## C6 maps: Figures 3 and S3–S4

```sh
Rscript Code/regenerate_C6_outputs.R A26 /absolute/path/to/CorrectedFigure
Rscript Code/regenerate_C6_outputs.R A27 /absolute/path/to/CorrectedFigure
```

The driver runs the spatial-model chunks of each C6 Rmd and the corresponding section of `Code/sr_sCCA_figures.R`. It loads and validates the saved Delaunay graph by default (optional rebuilding uses `adegenet::spca()` and `adespatial`), uses seed 123 and 25 tuning permutations, and fits 10 components. It does not run the Rmds' separate nonspatial CCA comparison.

Outputs include:

- `A26_components/A26_comp1.png` through `A26_comp10.png`.
- `A27_components/A27_comp1.png` through `A27_comp10.png`.
- Current spatial-model weights named `a26_score.RData` and `a27_score.RData` in the configured data-output directory.
- When explicitly requested, rebuilt `sPCA_a26metCluster.RData` and `sPCA_a27metCluster.RData` in the configured data-output directory.
- Tuning summaries, plotting state, session information, and diagnostics in `Figure/C6_diagnostics/`.

The paired maps place metabolite scores V on the left and protein scores U on the right. PNGs use the supplied plotting code's orientation, plasma palette, 25 × 20 cm size, and 300 dpi resolution. The graph uses one coordinate per sample; all repeated sample cells are retained when rendering the flatmap.

### Cross-animal projection

Run `Code/A26A27_analysis_sCCA.Rmd` after the current A27 C6 model is available. It uses A27 model weights and cluster loadings with selected A26 data, checks the A26 graph/sample order, and projects through shared features.

The default output folder is `Figure/A27_to_A26_components`. To override it before rendering:

```r
options(proteinfix.projection_export_dir = "/absolute/path/to/output")
rmarkdown::render("Code/A26A27_analysis_sCCA.Rmd", knit_root_dir = getwd())
```

For components 1–10, it exports `proj_U*.tiff` and `proj_V*.tiff` at 600 dpi and paired `A27_to_A26_comp*.png` files at 300 dpi. Numerical results are saved as `A27_to_A26_projection.rds` in the configured data-output directory in the original fitted-model sign convention.

### Composite layout

The component export scripts produce paired maps. Assemble the final manuscript TIFFs separately using the following workflow:

| Final figure | Component PNGs used |
| --- | --- |
| `Figure3.tiff` | Selected PNGs from `A26_components/`, `A27_components/`, and `A27_to_A26_components/`; each folder contains components 1–10. |
| `Figure S3.tiff` | All 10 PNGs from `A26_components/`. |
| `Figure S4.tiff` | All 10 PNGs from `A27_components/`. |

The committed scripts do not perform this final selection/layout/TIFF assembly. Preserve the selection, layout, and export settings with the manuscript assets to reproduce the assembled figures. The component PNGs are exported at 300 dpi; the separate 600-dpi TIFF panel settings below do not establish the effective resolution of these assembled figures.

### Separate individual TIFF panel exports in the Rmds

The A26 and A27 C6 Rmds also contain export blocks for individual U (protein) and V (metabolite) maps, controlled by `export_maps`, with `fmt = "tiff"` and `dpi = 600`. Their default destinations are `Figure/A26_components/a26_U*.tiff` and `a26_V*.tiff`, and `Figure/A27_components/a27_U*.tiff` and `a27_V*.tiff`. The projection Rmd similarly exports `proj_U*.tiff` and `proj_V*.tiff` alongside its paired PNGs.

These individual TIFF panels are separate from the PNG inputs used for Figures 3, S3, and S4. `regenerate_C6_outputs.R` runs selected model chunks and the paired PNG plotting script; it does not execute the C6 Rmds' individual-map TIFF export blocks. Figures 4 and S5, by contrast, are exported directly as complete 600-dpi TIFF figures by their regeneration scripts; their PNG copies are additional outputs. Figure 5 is exported directly as TIFF by its PCA Rmd.

### Display signs

The corrected local presentation set also uses paired display-sign flips for A27 components 1 and 3 and their A27-to-A26 projections. These flips are presentation edits; the committed fitted models retain their original signs. They are not applied automatically by the export scripts. Record any sign orientation with the final manuscript assets and flip both U and V together.

## Figure 4: spatial neighbor maps

```sh
Rscript Code/regenerate_Figure4.R /absolute/path/to/output
```

The driver executes the raw-metabolite and raw-protein sections of `Code/A26/SampleSize_A26_spearman.Rmd` and `Code/A27/SampleSize_A27_spearman.Rmd`. It does not use the optional protein-cluster comparison for Figure 4.

`Code/exact_neighbor_orders.R` defines neighbors by exact shortest-path distance. Groups exclude self and all closer samples. A26 displays orders 1/2/3; A27 displays 1/3/5. Features are standardized across the spatial samples, then Spearman correlations across features are averaged over each sample's neighbors. Constant protein columns are omitted because their standardized values are undefined.

The script checks the graph groups against `igraph` shortest paths and the per-sample means against a separately computed sample-correlation matrix.

Outputs are `Figure4_corrected.tiff` (600 dpi), `Figure4_corrected.png`, numeric `A26_Figure4_results.rds` and `A27_Figure4_results.rds`, and diagnostic plots. The assembled figure has Metabolomics/Proteomics column headings, horizontal A26/A27 labels, a–d panel letters, ordinal map headings, a bottom Neighbor order label, and a shared blue–white–red `avg cor` legend on a fixed −1 to 1 scale.

## Figure S5: spatial gradient boxplots

```sh
Rscript Code/regenerate_FigureS5.R /absolute/path/to/output
```

This uses the preparation, standardization, and protein-cluster aggregation from the same two SampleSize Rmds as Figure 4. `Code/protein_cluster_contributions.R` reads current C5 protein loadings from `Data/A26_PChclust.RData` and `Data/A27_PChclust.RData`; the separately exported contribution workbooks are not used by this runner.

The left column summarizes metabolite correlations. The right summarizes correlations across loading-weighted protein clusters. Each box contains **unique undirected sample-pair** Spearman correlations for its exact neighbor order, rather than Figure 4's per-sample means. Median bars and 1.5-IQR whiskers are shown without outlier symbols; connected dots indicate means. Direct pairwise calculations check a subset of the matrix-derived correlations.

Outputs are `FigureS5_corrected.tiff` (600 dpi), `FigureS5_corrected.png`, numeric results for each animal, and diagnostics. Gray boxes, A26/A27 row labels, column headings, Average correlation, and the bottom Neighbor order label follow the manuscript layout.

## Figure 5: U5/A26/A27 metabolomics PCA

Run `Code/loadingComp_U5A26A27_met.Rmd` through its `PCA without QC samples` chunk to reproduce Figure 5. Subsequent chunks compare PCA loadings with sr-sCCA and are separate from this figure.

Inputs are the U5 `5583 FIA Q500 Data.xlsx` and `5583 UPLC Q500 Data.xlsx` workbooks and the A26/A27 raw Q500 workbooks. The Rmd applies its own filtering, half-minimum imputation before log2 transformation, dilution adjustments, QC3 normalization, and outlier-removal procedure. This is not the post-ComBat PCA used in Figure 2.

The current Figure 5 run uses **252 shared metabolites** and, after the Rmd's QC/outlier exclusions, **396 A26, 542 A27, and 68 U5 samples**. PC1 explains **40.2%** and PC2 **13.3%** of the variance. Proteomics imputation does not enter this PCA.

The Rmd's built-in export is `Figure/fig5_pca.tiff`, 89 × 73 mm at 600 dpi. Its final export uses orange A26, blue A27, and black U5. Earlier plotting expressions in the same chunk use cyan A26, red A27, and purple U5. The requested local `Figure5_corrected.tiff`/PNG presentation uses those earlier colors and a wider canvas; that presentation override is not committed to the Rmd.

## Supplementary Table S3

The **A26** C6 driver exports `Supp_Table_S3_protein_clusters_min10_UPDATED.xlsx`: one sheet per protein cluster containing at least 10 proteins, with signed loadings sorted by absolute magnitude. The current corrected output has **80 sheets and 3,834 protein rows**. The sheet count is derived from clustering, not fixed to the earlier manuscript count.

Use a new output directory or move the existing workbook before repeating the C6 export; the exporter protects existing files. The built-in export contains `protein_accession` and `loading`.

The locally annotated workbook in the corrected figure set additionally contains a third `Gene` column. Annotations were matched by accession from the supplied example, the local `PG.Genes` annotation table, and UniProt; protein groups can list multiple available gene names. It has 3,797 annotated rows and 37 marked `Not annotated`. Cluster membership, row order, and loadings were preserved.

The gene-annotation step and its mapping inputs are not part of the committed C6 exporter. Rerunning C6 therefore produces the two-column table, not the annotated three-column manuscript copy. Preserve the annotation mapping with the manuscript assets if exact reproduction of that copy is needed.
