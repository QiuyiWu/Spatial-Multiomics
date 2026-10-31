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
- `PCA_Figure.R` writes `Figure/post_PCA_4panel_shared_legend.tiff` (Figure 2) and `Figure/pre_PCA_4panel_shared_legend.tiff` (Figure S2), using each stage's own metadata. As in the existing figure layout, panel letters are added separately in Photoshop.
- Spatial correlation gradients are derived from spatial CCA (sCCA) using Laplacian-based spatial weighting.
- Neighbor correlation analyses (Figure 4 and S5) quantify spatial decay of metabolomic similarity across cortical samples.


