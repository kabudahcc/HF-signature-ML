# HF-signature-ML

Analysis code for the manuscript:

**A Single-Cell Transcriptomic Atlas–Guided Diagnostic Signature for Heart Failure: Ensemble Machine Learning, External Validation, and Genetic Association of PLTP**

## Overview

This repository contains the R scripts for the bulk-tissue machine-learning arm of the study:

1. A consensus heart-failure (HF) signature (367 up / 926 down genes) derived from the GSE145154 single-cell atlas is used as the feature space.
2. Three discovery cohorts (GSE16499, GSE5406, GSE57338) are merged with ComBat under a strict **split-first, leakage-free** protocol (ComBat and z-score parameters fitted on the training set only).
3. Twelve algorithms are benchmarked; LASSO pre-selection yields a 39-gene model; random forest is selected as the primary model (statistically tied with LDA and radial SVM).
4. The frozen model is validated in seven independent external cohorts spanning microarray and RNA-seq platforms.
5. Model interpretation uses kernel SHAP; reviewer-requested robustness analyses (bootstrap CIs, DeLong CIs, ComBat PCA diagnostics, nested cross-validation) are included.

## Scripts

| Script | Purpose |
|---|---|
| `Scripts/step20_HF_signature_ML.R` | Original (pre-split) 12-algorithm benchmark + SHAP |
| `Scripts/step20b_HF_signature_ML_splitfirst.R` | **Primary** split-first leakage-free 12-algorithm benchmark + SHAP |
| `Scripts/step20c_RF_SHAP_splitfirst.R` | Retrain RF primary model on split-first pipeline; RF kernel SHAP |
| `Scripts/step21_HF_signature_external_validation.R` | External validation (original pre-split model, 7 cohorts) |
| `Scripts/step21b_external_validation_splitfirst.R` | **Primary** external validation of the split-first RF model (7 cohorts) |
| `Scripts/step22_review_response_round3.R` | Bootstrap 95% CIs (internal test), DeLong 95% CIs (external), Figure 3 redraw, UMAP by condition/donor |
| `Scripts/step23a_ComBat_PCA_diagnostic.R` | PCA before/after ComBat, between-cohort centroid distances |
| `Scripts/step23b_nested_CV.R` | Nested cross-validation (5-fold outer × 10 repeats), LASSO re-fitted inside each outer fold |
| `Scripts/plot_fig1C_signature_scores.R` | Figure 1C: HF-up/HF-down signature scores by cell type and condition |

## Data

All cohorts are publicly available from GEO:

- **Single-cell atlas:** GSE145154
- **Discovery (bulk):** GSE16499, GSE5406, GSE57338
- **External validation:** GSE116250, GSE55296, GSE42955, GSE26887, GSE46224, GSE52601, GSE48166

## Requirements

R (>= 4.2) with packages: `glmnet`, `randomForest`, `xgboost`, `e1071`, `caret`, `pROC`, `kernlab`, `sva`, `PRROC`, `kernelshap`, `shapviz`, `Seurat`, `GEOquery`, `hugene10sttranscriptcluster.db`, `AnnotationDbi`, `ggplot2`, `dplyr`, `pheatmap`, `patchwork`, `MASS`.

Scripts expect a working directory containing the GEO `*_symbol.csv` / `*_group.csv` matrices and the single-cell signature files under `results_scRNA/step3_4_5_pathway_cross_signature/`.

## Notes

- The split-first scripts (`step20b`, `step20c`, `step21b`) are the versions reported in the manuscript; the `step20`/`step21` originals are retained for transparency.
- Two scale-incompatible cohorts (GSE46224, GSE48166) are analysed with within-cohort z-scoring and reported only as a platform-compatibility sub-analysis, not primary validation.
