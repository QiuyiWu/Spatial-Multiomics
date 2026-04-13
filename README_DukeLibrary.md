# README for Duke Library Data Deposit

## Spatial-Multiomics Proteomics Input Files

This deposit accompanies the GitHub repository for the manuscript:

**Spatial-Multiomics**
GitHub repository: `https://github.com/QiuyiWu/Spatial-Multiomics`

## Overview

This Duke Library deposit contains two large protein level quantification proteomics input files that were not included in the GitHub repository because their file sizes exceed GitHub's recommended and/or allowed limits for convenient version control. These files are required only at the initial preprocessing stage of the proteomics workflow.

All downstream analysis code, processed data used in the manuscript, and paper figures are available in the GitHub repository above.

## Files included in this deposit

1. `A26 Proteomics with QC variables.txt`
2. `A27_SNE_combine_protein.tsv`


## Role of these files in the workflow

These two files are used only in the initial preprocessing scripts:

* `Reformatting_A26_pro.R`
* `Reformatting_A27_pro.R`

These scripts read the proteomics input files `A26 Proteomics with QC variables.txt` for A26 and `A27_SNE_combine_protein.tsv` for A27, and reformat them into the corresponding A26 and A27 wide-format proteomics data objects `Adjusted_Protein_A26.RData` and `Adjusted_Protein_A27.RData`, which are used for downstream preprocessing and analysis.

After this preprocessing step, the downstream analyses rely on the processed outputs rather than the original large protein level quantification files.



## Relationship to the GitHub repository

The GitHub repository contains:

* analysis scripts
* processed data objects used in downstream analysis
* figure-generation code
* figures included in the manuscript

The two large protein level quantification proteomics files are deposited separately here in Duke Library solely because of file-size constraints for GitHub hosting.

## Reproducibility notes

To fully reproduce the proteomics preprocessing workflow from protein level quantification input, users should:

1. download the two large files from this Duke Library deposit
2. place them in the expected local data directory referenced by:

   * `Reformatting_A26_pro.R`
   * `Reformatting_A27_pro.R`
3. run those preprocessing scripts before running the downstream analysis scripts in the GitHub repository

Users interested only in reproducing the manuscript analyses from the processed data may use the GitHub repository directly, without rerunning the original data preprocessing step, provided the processed files included there are sufficient for their purpose.

## Notes

These deposited files are included to support full computational reproducibility of the proteomics preprocessing pipeline associated with the manuscript.
