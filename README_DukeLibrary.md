# External proteomics inputs

This document describes inputs used by the ProteinFix branch of [Spatial-Multiomics](https://github.com/QiuyiWu/Spatial-Multiomics/tree/ProteinFix). It is not a deposit accession or a confirmation that every required file has been published.

## Quantification files

Two large long-format quantification files are described as deposited separately in Duke Library:

| Input | Used by | Role |
| --- | --- | --- |
| `A26 Proteomics with QC variables.txt` | `Code/A26/Reformatting_A26_pro.R` | A26 protein quantification, including the log2-scale `Y` values used for reformatting |
| `A27_SNE_combine_protein.tsv` | `Code/A27/Reformatting_A27_pro.R` | A27 intensity-scale `PG.Quantity` values; the script applies `log2()` |

A deposit URL or accession is not recorded here. Obtain the files from the study's data deposit or contact the authors for access details.

## Additional A27 annotation input

`Code/A27/Reformatting_A27_pro.R` also requires **`20260120_122736_Identifications.tsv`**, which provides `FileName`, `Plate.`, and `SPQC.type`. It is not tracked under `Data/`. Confirm its availability with the authors; the quantification file alone is insufficient to reproduce A27 P1.

## Running P1

1. Obtain the quantification files and the A27 annotation file.
2. Load `readr`, `dplyr`, and `tidyr` in a clean R session.
3. Adjust the scripts' legacy relative paths to the downloaded inputs.
4. Choose the reproduction/update mode described in [REPRODUCIBILITY.md](REPRODUCIBILITY.md#reproduction-mode-and-reference-updates). P1 saves `Adjusted_Protein_A26.RData` and `Adjusted_Protein_A27.RData` to that mode’s data-output directory.
5. Run the scripts, then continue with M1/P2 and C3–C6 as described in [README.md](README.md).

The saved objects contain metadata in addition to 9,990 A26 or 10,253 A27 protein variables. Both outputs are on the log2 scale before P2 imputation.

## Starting from processed inputs

Users reproducing downstream analyses can start from the processed objects already tracked in the GitHub repository. Download them with Git LFS. The external long-format inputs are needed to reproduce P1, not to read the included downstream objects.

See [README_Figures.md](README_Figures.md) for figure commands. Generated figures and supplementary-table exports are ignored by Git; their presence in a local output directory does not imply that they are included in the repository.
