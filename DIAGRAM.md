# Spatial Multiomics


## A26
### Metabolomics
- M1. Data preprocessing and Batch Correction
    - Preprocess raw data `A26-Q500 Data.xlsx` and create batch corrected data (403 samples $\times$ 328 variables data frame) `PostCombat_A26_Metabolite_beforeHeldout.xlsx` (code: `Analysis_BatchEffect_A26_met_adj.Rmd`)
    - Steps:
        - Remove sample 244, 96, 188
        - For TG metabolites, remove those missing > 70%; for other metabolites, remove those missing > 85%. 
        - Impute LOD/2
        - Dilution factor
        - Remove 3 SPQC outliers
        - Log transform
        - Subtract SPQC
        - Combat 
        - PCA analysis (remove HeldOut_AnalysisIDs in Mapped_Samples sheet)

### Proteomics
- P1. Data Reformatting
    - Convert log2transformed long data `A26 Proteomics with QC variables.txt` into wide data `Adjusted_Protein_A26.RData` (428 samples $\times$ 9990 variables data frame) **[code: `Reformatting_A26_pro.R`]**
- P2. Data preprocessing and Batch Correction
    - Preprocess reformatted large protein level quantification data and create batch corrected data `PostCombat_A26_Protein.xlsx` (427 samples $\times$ 9990 variables data frame) **[code: `Analysis_BatchEffect_A26_pro_adj.Rmd`]**
    - Steps:
        - Remove outlier `ID117530`
        - Impute LOD/2
        - Dilution factor
        - Combat (we decided NOT to subtract the SPQC (plate-wise) before batch correction)
        - PCA analysis (remove HeldOut_AnalysisIDs in Mapped_Samples sheet)

### Combine Both
- C3. Heldout Preprocessing
    - Data driven way to select final mapping samples (267 samples selected) `Heldout_Updated.xlsx` **[code: `HeldoutPreprocessing.R`]**
    - Steps:
        - remove outliers (only Metabolitics data, because Proteomics outliers looks like another batch)
        - comparing the mapped samples and heldout samples and swap them if necessary
        - compute $k = (x-med)/mad$, for Proteomics and metabolitics data, each of them would have one or two k values (from original samples and heldout samples), and we unselect the one with maximum k values for the final mapping samplez. Here we use both PC1 and PC2, and compute the robust distance of each sample from the center:
        $d_i = \sqrt{(PC1_i - med(PC1))^2 + (PC2_i - med(PC2))^2}$. Then normalize this by the Median Absolute Deviation (MAD) of the distances: $k_i = \frac{d_i}{MAD(d)}$. This way, k represents how far a sample is from the robust center of the PCA space (PC1, PC2). We select samples with the smaller k values. 
- C4. Sample Selecting
    - Select mapping samples in both datasets and created data `a26.RData` (two data frames:  metabolomics -- 267 samples $\times$ 328 variables data frame; proteomics -- 267 samples $\times$ 9990 variables data frame) **[code: `CombineClustering.R`]**
- C5. Clustering Analysis
    - Variable clustering in both datasets and created selected data `A26_PChclust.RData` (two nested lists: metabolomics list and proteomics list, inside of each big list has 5 small lists with information of group, score, loading, cluster summary, and percentage of variance) **[code: `CombineClustering.R`]**
- C6. Spatial Sparse CCA analysis
    - Mapping the data then do Spatial SCCA analysis **[code: `A26_analysis_sCCA_laplacians.Rmd`]**
    - We ultimately obtained 233 matched samples across both datasets




## A27
### Metabolomics
- M1. Data preprocessing and Batch Correction
    - Preprocess raw data `A27-Q500 Data.xlsx` and create batch corrected data (548 samples $\times$ 343 variables data frame) `PostCombat_A27_Metabolite_beforeHeldout.xlsx` (code: `Analysis_BatchEffect_A27_met_adj.Rmd`)
    - Steps:
        - For TG metabolites, remove those missing > 70%; for other metabolites, remove those missing > 85%, remove outlier 212
        - Impute LOD/2
        - Dilution factor
        - Log transform
        - Subtract SPQC
        - Combat
        - PCA analysis (remove HeldOut_AnalysisIDs in Mapped_Samples sheet)

### Proteomics
- P1. Data Reformatting
    - Convert log2transformed long data `A27_SNE_combine_protein.tsv` into wide data `Adjusted_Protein_A27.RData` (566 samples $\times$ 10253 variables data frame) **[code: `Reformatting_A27_pro.R`]**
- P2. Data preprocessing and Batch Correction
    - Preprocess reformatted large protein level quantification data and create batch corrected data `PostCombat_A27_Protein.xlsx` (566 samples $\times$ 10253 variables data frame) **[code: `Analysis_BatchEffect_A27_pro_adj.Rmd`]**
    - Steps:
        - Impute LOD/2
        - Dilution factor
        - Combat (we decided NOT to subtract the SPQC (plate-wise) before batch correction)
        - PCA analysis (remove HeldOut_AnalysisIDs in Mapped_Samples sheet)

### Combine Both
- C3. Heldout Preprocessing
    - Data driven way to select final mapping samples (527 samples selected) `Heldout_Updated_a27.xlsx` **[code: `HeldoutPreprocessing.R`]**
    - Steps:
        - remove outliers (only Metabolitics data, because Proteomics outliers looks like another batch)
        - comparing the mapped samples and heldout samples and swap them if necessary
        - compute $k = (x-med)/mad$, for Proteomics and metabolitics data, each of them would have one or two k values (from original samples and heldout samples), and we unselect the one with maximum k values for the final mapping samplez. Here we use both PC1 and PC2, and compute the robust distance of each sample from the center:
        $d_i = \sqrt{(PC1_i - med(PC1))^2 + (PC2_i - med(PC2))^2}$. Then normalize this by the Median Absolute Deviation (MAD) of the distances: $k_i = \frac{d_i}{MAD(d)}$. This way, k represents how far a sample is from the robust center of the PCA space (PC1, PC2). We select samples with the smaller k values. 
- C4. Sample Selecting
    - Select mapping samples in both datasets and created data `a27.RData` (two data frames:  metabolomics -- 527 samples $\times$ 342 variables data frame; proteomics -- 527 samples $\times$ 10253 variables data frame) **[code: `CombineClustering.R`]**
- C5. Clustering Analysis
    - Variable clustering in both datasets and created selected data `A27_PChclust.RData` (two nested lists: metabolomics list and proteomics list, inside of each big list has 5 small lists with information of group, score, loading, cluster summary, and percentage of variance) **[code: `CombineClustering.R`]**
- C6. Spatial Sparse CCA analysis
    - Mapping the data then do Spatial SCCA analysis **[code: `A27_analysis_sCCA_laplacians.Rmd`]**
    - We ultimately obtained 508 matched samples across both datasets










