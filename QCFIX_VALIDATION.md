# QCFix correction validation

Validated on October 31, 2026. This update addresses early A27 plasma-QC exclusion, Figure 5 outlier-row identity, and Figure 5 metabolite-name matching.

## Before and after

| Quantity | Previous reference | Corrected result |
| --- | --- | --- |
| Commercial plasma QC rows entering A27 SPQC/ComBat feature filtering and imputation | 63 | 0 (21 rows from each QC level removed early) |
| A27 pre-ComBat PCA matrix | 644 × 343 | 581 × 300 |
| A27 post-ComBat metabolomics | 548 × 343 | 548 × 300 |
| A27 selected samples | 527 | 527; no selected IDs changed |
| A27 metabolite clusters | 117 | 83 |
| A27 protein clusters | 6,933 | 6,933; numerical cluster results unchanged |
| Spatially mapped A27 samples | 508 | 508; saved graph validated and reused |
| Projection shared metabolites | 289 | 277 |
| Figure 5 flagged / removed outlier rows | 14 / 21 | 14 / 14 |
| Figure 5 retained samples | 1,006 | 1,013 (398 A26, 547 A27, 68 U5) |
| Figure 5 metabolites | 252 | 252 |
| Figure 5 PC1 / PC2 variance | 40.2% / 13.3% | 40.1% / 13.3% |
| PCA/sCCA loading matches, A26 | 52 literal matches | 252 normalized matches |
| PCA/sCCA loading matches, A27 | 52 literal matches | 249 normalized matches with corrected A27 features |

A27 filtering dropped 44 previously retained features and added `C16:2-OH`. The heldout workbook's distance scores changed even though its selected IDs did not. All 33 animal-specific SPQC rows remain available for A27 plate normalization and are removed before ComBat.

The corrected A27 loading comparison genuinely lacks three Figure 5 metabolites: `DG(14:1_18:1)`, `ProBetaine`, and `GCDCA`. No features or zero weights were invented to force a 252-feature match. Original labels and normalized match keys are exported.

The seven unflagged rows incorrectly removed by the old Figure 5 expression were A26 samples 212 and 308, and A27 samples 2, 8, 201, 246, and 413. Removal now follows the exact aligned row flags, so an outlier in one animal cannot remove a same-numbered sample in another.

## Updated reference data

- `Data/PostCombat_a27_Metabolite_beforeHeldout.xlsx`
- `Data/met_a27_PCAdata.RData`
- `Data/Heldout_Updated_a27.xlsx`
- `Data/a27.RData`
- `Data/A27_PChclust.RData`
- `Data/a27_score.RData`
- `Data/A27_to_A26_projection.rds`
- `Data/Figure5_PCA_results.rds` (new reproducible PCA, outlier audit, feature matching, and loading comparisons)

A26's model rerun agrees numerically with its reference. A26 Figure 4/S5 results and A27 Figure 4 protein correlations / S5 protein-cluster correlations also agree with their prior results. Raw input workbooks and proteomics preprocessing references were not changed by this correction.

## Figure package

The local `CorrectedFigure2` package contains Figures 2–5 and S2–S5, all 30 paired component PNGs, projection individual TIFF panels, numeric Figure 4/S5 results, Figure 5 loading correlations and overlap curves, and audit files. The annotated A26 supplementary Table S3 is preserved unchanged from the prior local package; the C6 exporter additionally produces its reproducible two-column workbook.

`Code/assemble_manuscript_components.R` selects components 1/2 from each of the three folders for Figure 3. S3 contains A26 components 1–10; S4 contains A27 components 1–10, in row-major order. PNG contents are preserved; TIFF rasterization at 600 dpi does not add source detail. Components and projections use the saved fitted-model signs consistently; historical display-only sign flips are not applied.

Figure 5 retains the prior local cyan A26 / red A27 / purple U5 palette and 150 × 89 mm canvas, now recorded in the notebook rather than an external override.

## Provisional Figure 5 QC3 policy

At the user's instruction, QC3 continues to participate in Figure 5 missingness filtering, half-minimum imputation, and QC3 normalization, pending PI review. QC3 is removed before the biological-sample PCA. This policy is separate from early removal of all commercial QC levels in the A27 SPQC/ComBat workflow.

## Validation and environment

- Actual preprocessing verified 63 commercial QC exclusions and finite pre/post matrices.
- Full heldout selection and both C5 modality clusterings were rerun.
- C6 validated sample/graph compatibility and fitted 10 components with seed 123 and 25 tuning permutations.
- Projection scores were checked against explicit data-by-weight multiplication.
- Figure 4 graph distances and correlation means, and S5 pairwise correlations, passed their independent checks.
- Figure 5 repeated-ID removal, name alignment, normalization, and duplicate-key rejection passed targeted tests.
- `Code/validate_QCFix.R` passed against regenerated outputs.
- The existing reproduction-mode integration tests passed: C3/C4 numerical matches, graph loading/rebuilding, output routing, and reference-file preservation.
- Final figures were visually reviewed; TIFF metadata and component-file completeness were checked.

Local execution used R 4.4.3 and the installed system library, explicitly invoked with `--vanilla`; session information is included with the local audits. `renv.lock` still records R 4.5.2. A full run under that locked R version was not performed.

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md#regenerate-the-coordinated-qcfix-corrections) for the coordinated rerun command and [README_Figures.md](README_Figures.md) for individual exports.
