# HF-signature-ML

Analysis code for the manuscript:

**A Donor-Level Single-Cell Atlas–Anchored Tissue Biomarker Signature for Heart Failure: Ensemble Machine Learning and Multi-Cohort Validation**

## Overview

This repository contains the analysis code for the study, in which the **donor, rather than the cell, serves as the statistical unit throughout**:

1. A consensus heart-failure (HF) signature (145 up / 552 down-regulated genes) is derived from **donor-level pseudobulk differential expression** in the multi-donor GSE183852 single-cell atlas (27 non-failing donors, 18 DCM patients), requiring concordant directional dysregulation in ischemic disease in the independent GSE145154 atlas.
2. Three bulk discovery cohorts (GSE16499, GSE5406, GSE57338; n = 553) are merged with ComBat under a strict **split-first, leakage-free** protocol (ComBat and z-score parameters fitted on the training set only).
3. **Twelve algorithm configurations** are benchmarked; LASSO pre-selection yields a **36-gene panel**; random forest is carried forward as the primary model. **Nested cross-validation (AUC 0.975, 95% CI 0.970–0.980)** is designated as the primary internal performance measure.
4. The frozen model is evaluated without retraining in seven independent external cohorts spanning microarray and RNA-seq platforms, with random-effects meta-analysis, leave-one-cohort-out cross-validation, post-hoc power analysis, and decision curve analysis; classifier-level external results are reported as hypothesis-generating.
5. **Decisive fingerprint-level validation in MAGNet (GSE141910; 166 DCM vs. 166 non-failing):** 35/36 panel genes replicate directionally (binomial P = 1.1 × 10⁻⁹); a label-free direction score achieves AUC 0.999 and extends to hypertrophic and peripartum cardiomyopathy.
6. Model interpretation uses kernel SHAP (gene- and pathway-level); robustness analyses include bootstrap CIs, DeLong CIs, ComBat PCA diagnostics, etiology-stratified AUCs, age/sex covariate adjustment, and LVEF regression.
7. Orthogonal analyses: spatial transcriptomic corroboration (Kuppe et al. myocardial-infarction Visium atlas) and OmniPath/CollecTRI surrogate ligand prioritization.
8. Genetic follow-up: a panel-wide cis-eQTL MR screen (nominal protective STAT3 signal) and locus-level analyses of the prior candidates PLTP/FBLN1 (GTEx heart eQTL MR, plasma pQTL MR, TWAS, colocalization, and SuSiE fine-mapping), resolving PLTP as a mechanistic rather than causal finding.

## Scripts

### Core machine-learning pipeline

| Script | Purpose |
|---|---|
| `Scripts/step20_HF_signature_ML.R` | Original (pre-split) benchmark + SHAP |
| `Scripts/step20b_HF_signature_ML_splitfirst.R` | **Primary** split-first leakage-free benchmark + SHAP |
| `Scripts/step20c_RF_SHAP_splitfirst.R` | Retrain RF primary model on split-first pipeline; RF kernel SHAP |
| `Scripts/step21_HF_signature_external_validation.R` | External validation (original pre-split model, 7 cohorts) |
| `Scripts/step21b_external_validation_splitfirst.R` | **Primary** external validation of the split-first RF model (7 cohorts) |
| `Scripts/step22_review_response_round3.R` | Bootstrap 95% CIs (internal test), DeLong 95% CIs (external), UMAP by condition/donor |
| `Scripts/step23a_ComBat_PCA_diagnostic.R` | PCA before/after ComBat, between-cohort centroid distances |
| `Scripts/step23b_nested_CV.R` | Nested cross-validation (5-fold outer × 10 repeats), LASSO re-fitted inside each outer fold |
| `Scripts/plot_fig1C_signature_scores.R` | Figure 1C: donor-level HF-up/HF-down module scores by cell type and condition |

### Robustness, meta-analysis, and interpretation

| Script | Purpose |
|---|---|
| `Scripts/step24_LOCO_meta_analysis.R` | Leave-one-cohort-out cross-validation across discovery cohorts |
| `Scripts/step24b_meta_analysis.R` | DerSimonian–Laird random-effects pooling of external cohort AUCs |
| `Scripts/step26_fig1C_stats.R` | Fig 1C module-score statistics (donor-level t tests, Cohen's d, FDR) |
| `Scripts/step27_pathway_SHAP.R` | Pathway-level SHAP aggregation over Reactome pathways |
| `Scripts/step28_subgroup.R` | Etiology-stratified (ICM/DCM) AUCs |
| `Scripts/step30_native_xgb_fullspace_rf.R` | Robustness: native-API XGBoost and full-feature-space random forest |
| `Scripts/step31_fig3_inset.R` | Meta-analysis inset figure |
| `Scripts/step32_xgb_cohort_overlap.R` | Per-cohort feature-overlap checks for XGBoost |
| `Scripts/step38_covariate_adjustment.R` | Age/sex covariate-adjusted validation |
| `Scripts/step44_dca_external.py` | Decision curve analysis in external cohorts |

### Donor-level signature derivation and cross-atlas analyses

| Script | Purpose |
|---|---|
| `Scripts/step33_crossatlas_validation.R` | GSE183852 donor-level pseudobulk differential expression (signature discovery) |
| `Scripts/step33a_build_reference_directions.R` | GSE145154 directional replication reference |
| `Scripts/step33b_extract_slim.R` | GSE183852 donor-level pseudobulk extraction |
| `Scripts/step33c_run_concordance.R` | Direction-concordance statistics between the two atlases |
| `Scripts/step34_rebuild_fig3.R` | Figure rebuild with final cohort set |
| `Scripts/step35_spatial_pltp.R` | Spatial transcriptomic validation (Kuppe et al. Visium atlas) |
| `Scripts/step36_lvef_regression.R` | LVEF regression of model scores |
| `Scripts/step37_lvef_gse183852.R` | LVEF analysis in GSE183852 |
| `Scripts/step39_nichenet_surrogate.py` | Surrogate ligand prioritization with OmniPath/CollecTRI priors |
| `Scripts/step40_niche_refine.py` | Ligand-activity refinement and sender evidence |
| `Scripts/step41_fig_S23.py` | Supplementary ligand-prioritization figure |
| `Scripts/step42_make_data_s11.py` | Supplementary Data workbook generation |
| `Scripts/step43_gse183852_ligand_replication.R` | Donor-level replication of ligand/receptor expression in GSE183852 |

### Genetic analyses

| Script | Purpose |
|---|---|
| `Scripts/step29_susie_finemap.R` | SuSiE summary-statistic fine-mapping at the PLTP locus (GWAS vs. eQTL credible sets) |

### Utilities

| Script | Purpose |
|---|---|
| `Scripts/step45_sessioninfo.R` | R session information export (reproducibility) |

## Data

All cohorts are publicly available:

- **Single-cell atlases:** GSE183852 (discovery), GSE145154 (directional replication)
- **Discovery (bulk):** GSE16499, GSE5406, GSE57338
- **External validation:** GSE116250, GSE55296, GSE42955, GSE26887, GSE46224, GSE52601, GSE48166
- **Fingerprint-level validation:** MAGNet (GSE141910)
- **Spatial:** Kuppe et al. 2022 myocardial-infarction Visium atlas (CZ CELLxGENE)
- **Genetics:** GTEx v8 / eQTL Catalogue, HERMES and GBMI HF GWAS, UKB-PPP and SCALLOP pQTL

## Requirements

- R (>= 4.4; analyses run on 4.6.0) with packages: `Seurat`, `harmony`, `limma`, `edgeR`, `glmnet`, `randomForest`, `xgboost`, `e1071`, `caret`, `pROC`, `kernlab`, `sva`, `PRROC`, `kernelshap`, `shapviz`, `susieR`, `GEOquery`, `AnnotationDbi`, `data.table`, `dplyr`, `ggplot2`, `pheatmap`, `patchwork`, `MASS`, `foreach`, `doParallel`.
- Python (>= 3.10; analyses run on 3.12) with `pandas`, `numpy`, `scipy`, `matplotlib`, `openpyxl`, `requests` for the ligand-prioritization and decision-curve scripts.

Scripts expect a working directory containing the GEO `*_symbol.csv` / `*_group.csv` matrices and the single-cell signature files under `results_scRNA/`.

## Notes

- The split-first scripts (`step20b`, `step20c`, `step21b`) are the versions reported in the manuscript; the `step20`/`step21` originals are retained for transparency.
- Two scale-incompatible cohorts (GSE46224, GSE48166) are analysed with within-cohort z-scoring and reported only as a platform-compatibility sub-analysis, not primary validation.
- Genetic results are reported as nominal, hypothesis-generating associations: colocalization and SuSiE fine-mapping did not support a shared causal variant at the PLTP locus.

## Citation

A citable Zenodo DOI will be added here upon public release of the repository.
