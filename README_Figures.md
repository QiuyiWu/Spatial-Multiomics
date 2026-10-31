## 📊 Figures and Corresponding Code

This repository contains the code used to generate all figures for the manuscript:  
**“Integrated metabolomics and proteomics from voxelated cortical hemispheres of adult rhesus monkeys”**  
by *Qiuyi Wu\*, Alev M. Brigande\*, Michael W. Lutz, Pixu Shi, and Anita A. Disney*.

---

## Main Figures

- **Figure 2**: Post-ComBat PCA plots  
  `PCA_Figure.R`

- **Figure 3a–b**: Spatial correlation gradients for A26  
  `A26_analysis_sCCA_laplacians.Rmd`

- **Figure 3c–d**: Spatial correlation gradients for A27  
  `A27_analysis_sCCA_laplacians.Rmd`

- **Figure 3e–f**: Projection of A27-derived components onto A26  
  `A26A27_analysis_sCCA.Rmd`

- **Figure 4a–b**: Correlation with neighboring samples (A26)  
  `SampleSize_A26_spearman.Rmd`

- **Figure 4c–d**: Correlation with neighboring samples (A27)  
  `SampleSize_A27_spearman.Rmd`

- **Figure 5**: PCA of metabolomics with blood present (U5) vs. perfused samples (A26, A27)  
  `loadingComp_U5A26A27_met.Rmd`

---

## Supplementary Figures

- **Figure S2**: Pre-ComBat PCA plots  
  `PCA_Figure.R`

- **Figures S3–S4**: Additional spatial correlation gradient analyses  
  `A26_analysis_sCCA_laplacians.Rmd`  
  `A27_analysis_sCCA_laplacians.Rmd`

- **Figure S5**: Spatial gradient boxplots  
  `SampleSize_A26_spearman.Rmd`  
  `SampleSize_A27_spearman.Rmd`

---

## Notes

- To regenerate Figure 2 and S2, run the four `Analysis_BatchEffect_A26_met_adj.Rmd`, `Analysis_BatchEffect_A26_pro_adj.Rmd`, `Analysis_BatchEffect_A27_met_adj.Rmd`, and `Analysis_BatchEffect_A27_pro_adj.Rmd` analyses under `Code/A26` and `Code/A27`, then run `Code/PCA_Figure.R` from the repository root.
- Each analysis ends with a `save()` export to its `Data/{met,pro}_a{26,27}_PCAdata.RData` file. Each file contains numeric `pre_*` and `post_*` matrices (samples by analytes), plus matching `pre_*_metadata` and `post_*_metadata` tables with `sample_id`, `plate`, and `sample_type`. Matrix row names match metadata row names.
- Metabolomics inputs use 328 analytes for A26 and 343 for A27, with the specified outliers excluded and log2 transformation applied before correction. Pre-ComBat matrices retain the QC samples used at that stage; post-ComBat matrices follow the analyses' QC exclusions. Proteomics post matrices use ordinary ComBat without SPQC subtraction and include the corrected half-minimum imputation.
- `PCA_Figure.R` writes `Figure/post_PCA_4panel_shared_legend.tiff` (Figure 2) and `Figure/pre_PCA_4panel_shared_legend.tiff` (Figure S2), using each stage's own metadata. Both figures include Metabolomic/Proteomic column headings, A26/A27 row labels, and bold lowercase a–d panel letters.
- Spatial correlation gradients are derived from spatial CCA (sCCA) using Laplacian-based spatial weighting.
- Neighbor correlation analyses (Figure 4 and S5) quantify spatial decay of metabolomic similarity across cortical samples.



## Regenerating C6 component maps and Table S3

From the repository root, run the spatial C6 analysis and the supplied plotting code for each animal:

```sh
Rscript --vanilla Code/regenerate_C6_outputs.R A26 /path/to/CorrectedFigure
Rscript --vanilla Code/regenerate_C6_outputs.R A27 /path/to/CorrectedFigure
```

The driver evaluates the spatial-model chunks in the corresponding C6 Rmd, rebuilding the Delaunay graph from the current selected samples and coordinates. It retains seed 123, 25 tuning permutations, and 10 components. It then runs the appropriate animal section of `Code/sr_sCCA_figures.R`, preserving the supplied map orientation, labels, PNG format, 25 × 20 cm size, and 300 dpi resolution.

Outputs are `A26_components/A26_comp1.png` through `A26_comp10.png`, and the analogous A27 files. The A26 run also creates `Supp_Table_S3_protein_clusters_min10_UPDATED.xlsx`: one worksheet per A26 protein cluster containing at least 10 proteins, with signed loadings ranked by absolute magnitude. Its worksheet count is derived from the current clustering, rather than fixed to the old count of 73. Choose a new output directory if the supplementary workbook already exists.

The spatial model weights and rebuilt sPCA graphs are saved under `Data`. Intermediate plotting state, tuning summaries, and diagnostics are saved under the ignored `Figure/C6_diagnostics` folder. The driver does not run the Rmds' separate nonspatial CCA comparison.

## A27-derived components projected onto A26

`Code/A26A27_analysis_sCCA.Rmd` uses the current A27 C6 weights, A27 cluster loadings, and selected A26 data to project through the shared features. It verifies that the A26 spatial graph matches the current sample/coordinate order. The default export folder is `Figure/A27_to_A26_components`; it can be changed with `options(proteinfix.projection_export_dir = "/path/to/output")` before running the Rmd.

For components 1–10 it exports `proj_U*.tiff` and `proj_V*.tiff` (600 dpi), plus paired `A27_to_A26_comp*.png` figures (300 dpi, V on the left and U on the right). Numerical results are saved to `Data/A27_to_A26_projection.rds` in the original fitted-model sign convention. Any optional display-only sign orientation should be recorded separately from these analytical results.

### Reproduce corrected Figure 4

From the repository root, run:

```sh
Rscript --vanilla Code/regenerate_Figure4.R /path/to/output
```

This executes the raw metabolomics and proteomics Figure 4 sections of both
`SampleSize_A26_spearman.Rmd` and `SampleSize_A27_spearman.Rmd`, then assembles the
labeled TIFF (600 dpi) and PNG. The optional protein-cluster sections are not used.
Neighbor order is exact shortest-path distance on the current spatial graph:
self and all lower-order neighbors are excluded. A26 displays orders 1/2/3;
A27 displays orders 1/3/5. Features are standardized across spatial samples before
Spearman correlations across features are averaged over each sample's neighbors.
Constant protein columns are omitted because their standardized values are undefined.
Graph shells are checked against igraph shortest paths, and per-sample means are
checked against a separately computed sample correlation matrix. Result RDS files
and diagnostic plots are saved alongside the figure.
