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
  `loadingComp_U6A26A27_met.Rmd`

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

- PCA plots (pre- and post-ComBat) are generated from the same script (`PCA_Figure.R`) with different preprocessing settings.
- Spatial correlation gradients are derived from spatial CCA (sCCA) using Laplacian-based spatial weighting.
- Neighbor correlation analyses (Figure 4 and S5) quantify spatial decay of metabolomic similarity across cortical samples.


