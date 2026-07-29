# step22: 第三轮审稿回应补充分析
# A. 内部测试集 RF/LDA bootstrap 95% CI (2000 次重采样)
# B. 外部验证七队列 DeLong 95% CI
# C. Figure 3 重绘: 主图仅 5 个 frozen 队列; 2 个 per-cohort-z 队列单独出补充图
# D. 单细胞 UMAP 按 condition 着色 (Figure 1 补充)
# 输出: step22_review_response/

.libPaths(c("Rlib", .libPaths()))
rm(list = ls())

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(dplyr)
})

set.seed(12345)
out_dir <- "step22_review_response"
if (!dir.exists(out_dir)) dir.create(out_dir)

# A. 内部测试集 bootstrap CI ---------------------------------------------------
cat("== A. Internal test AUC bootstrap CI ==\n")
tp <- read.csv("step20b_HF_signature_ML_splitfirst/test_predictions.csv", stringsAsFactors = FALSE)
y_true <- ifelse(tp$Group == "Disease", 1, 0)
n_boot <- 2000
boot_ci <- function(prob) {
  aucs <- replicate(n_boot, {
    idx <- sample(length(y_true), replace = TRUE)
    if (length(unique(y_true[idx])) < 2) return(NA)
    as.numeric(auc(roc(y_true[idx], prob[idx], quiet = TRUE)))
  })
  aucs <- aucs[!is.na(aucs)]
  c(full = as.numeric(auc(roc(y_true, prob, quiet = TRUE))),
    lower = quantile(aucs, 0.025), upper = quantile(aucs, 0.975))
}
ci_rf  <- boot_ci(tp$RF)
ci_lda <- boot_ci(tp$LDA)
boot_df <- data.frame(
  Model = c("RF", "LDA"),
  AUC = c(ci_rf["full"], ci_lda["full"]),
  CI_lower = c(ci_rf["lower"], ci_lda["lower"]),
  CI_upper = c(ci_rf["upper"], ci_lda["upper"]))
boot_df$AUC_95CI <- sprintf("%.4f (%.4f-%.4f)", boot_df$AUC, boot_df$CI_lower, boot_df$CI_upper)
write.table(boot_df, file.path(out_dir, "internal_test_bootstrap_CI.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(boot_df)

# B. 外部验证 DeLong CI --------------------------------------------------------
cat("== B. External validation DeLong CI ==\n")
pred <- read.csv("step21b_HF_signature_external_validation_splitfirst/external_validation_predictions.csv",
                 stringsAsFactors = FALSE)
ci_rows <- lapply(split(pred, pred$Dataset), function(d) {
  y_ext <- ifelse(d$TrueLabel == "Disease", 1, 0)
  if (length(unique(y_ext)) < 2) return(NULL)
  r <- roc(y_ext, d$Prob_Disease, quiet = TRUE)
  ci <- as.numeric(ci.auc(r, method = "delong"))
  data.frame(Dataset = d$Dataset[1], n = nrow(d),
             n_Disease = sum(y_ext), n_Normal = sum(y_ext == 0),
             AUC = ci[2], CI_lower = ci[1], CI_upper = ci[3])
})
ci_df <- bind_rows(ci_rows) %>% arrange(desc(AUC))
ci_df$AUC_95CI <- sprintf("%.3f (%.3f-%.3f)", ci_df$AUC, ci_df$CI_lower, ci_df$CI_upper)
write.table(ci_df, file.path(out_dir, "external_validation_AUC_95CI.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(ci_df)

# C. Figure 3 重绘 --------------------------------------------------------------
cat("== C. Figure 3 redraw ==\n")
per_cohort_z <- c("GSE48166", "GSE46224")
pred$scaling <- ifelse(pred$Dataset %in% per_cohort_z, "per-cohort z", "frozen")
roc_df <- bind_rows(lapply(split(pred, pred$Dataset), function(d) {
  y_ext <- ifelse(d$TrueLabel == "Disease", 1, 0)
  r <- roc(y_ext, d$Prob_Disease, quiet = TRUE)
  data.frame(Dataset = paste0(d$Dataset[1], " (AUC=", round(as.numeric(auc(r)), 3), ")"),
             scaling = d$scaling[1],
             FPR = 1 - r$specificities, TPR = r$sensitivities)
}))

p_main <- ggplot(roc_df %>% filter(scaling == "frozen"),
                 aes(x = FPR, y = TPR, color = Dataset)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  theme_bw() +
  labs(title = "HF signature RF model - external validation ROC (frozen-scale cohorts)",
       x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  theme(legend.text = element_text(size = 8))
ggsave(file.path(out_dir, "HF_signature_external_validation_ROC_main.pdf"), p_main, width = 8, height = 6)

p_supp <- ggplot(roc_df %>% filter(scaling == "per-cohort z"),
                 aes(x = FPR, y = TPR, color = Dataset)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  theme_bw() +
  labs(title = "Exploratory: scale-incompatible cohorts (within-cohort z-scoring)",
       x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  theme(legend.text = element_text(size = 8))
ggsave(file.path(out_dir, "HF_signature_external_validation_ROC_scale_incompatible.pdf"),
       p_supp, width = 7, height = 6)

# D. UMAP by condition -----------------------------------------------------------
cat("== D. UMAP by condition ==\n")
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("results_scRNA/step1_QC_integration/heart_integrated.rds")
meta_cols <- colnames(obj@meta.data)
cond_col <- meta_cols[grepl("condition|group", meta_cols, ignore.case = TRUE)][1]
cat("condition column:", cond_col, "| levels:",
    paste(unique(obj@meta.data[[cond_col]]), collapse = ", "), "\n")
red <- if ("umap" %in% Reductions(obj)) "umap" else Reductions(obj)[1]
p_cond <- DimPlot(obj, reduction = red, group.by = cond_col) +
  ggtitle("UMAP by condition") + theme_bw()
ggsave(file.path(out_dir, "scatlas_umap_by_condition.pdf"), p_cond, width = 7, height = 6)

donor_col <- meta_cols[grepl("donor", meta_cols, ignore.case = TRUE)][1]
if (is.na(donor_col)) {
  sample_col <- meta_cols[grepl("sample|orig", meta_cols, ignore.case = TRUE)][1]
  obj$donor <- sub("^([A-Za-z]+-[0-9]+).*$", "\\1", obj@meta.data[[sample_col]])
  donor_col <- "donor"
}
p_donor <- DimPlot(obj, reduction = red, group.by = donor_col, label = TRUE, repel = TRUE) +
  ggtitle("UMAP by donor") + theme_bw()
ggsave(file.path(out_dir, "scatlas_umap_by_donor.pdf"), p_donor, width = 8, height = 6)

cat("\nDone. Outputs in", out_dir, "\n")
