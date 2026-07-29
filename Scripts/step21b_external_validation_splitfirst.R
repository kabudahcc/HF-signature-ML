# step21b: HF 签名 RF 模型外部验证 (split-first 修正版)
# 与 step21 的区别: 模型改为 step20b split-first 流程重训的 RF(39 基因),
# 冻存 z-score 参数直接取自 step20b best_model.rds 中的 combat_params
# (g_mean/g_sd, 仅由训练集估计, 无测试集泄露)。
# RF 与 LDA 在 split-first benchmark 中并列最优(ΔAUC = 0.0004),
# 基于 SHAP 可解释性与叙事连续性选择 RF 为主模型(见 step20b 输出)。
# 输出: step21b_HF_signature_external_validation_splitfirst/ 下的性能表与 ROC/PR/校准曲线

.libPaths(c("Rlib", .libPaths()))
rm(list = ls())

suppressPackageStartupMessages({
  library(randomForest)
  library(pROC)
  library(caret)
  library(PRROC)
  library(ggplot2)
  library(dplyr)
  library(GEOquery)
  library(hugene10sttranscriptcluster.db)
  library(AnnotationDbi)
})

set.seed(12345)

out_dir <- "step21b_HF_signature_external_validation_splitfirst"
if (!dir.exists(out_dir)) dir.create(out_dir)

# 1. 载入 split-first 模型包, 重训 RF 主模型 --------------------------------
bundle <- readRDS("step20b_HF_signature_ML_splitfirst/best_model.rds")
selected <- bundle$features
train_center <- bundle$combat_params$g_mean
train_scale  <- bundle$combat_params$g_sd
cat("Benchmark best:", bundle$method, "| RF retrained as primary model | features:", length(selected), "\n")
stopifnot(all(selected %in% names(train_center)))

train_df <- bundle$train_set
model_rf <- randomForest(Group ~ ., data = train_df, ntree = 500)

# 2. 评估函数(与 step21 相同) ------------------------------------------------
evaluate_model <- function(prob, y_true) {
  roc_obj <- roc(y_true, prob, quiet = TRUE)
  auc <- as.numeric(auc(roc_obj))
  pr_obj <- pr.curve(scores.class0 = prob[y_true == 1],
                     scores.class1 = prob[y_true == 0], curve = TRUE)
  auprc <- pr_obj$auc.integral
  best_coords <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"),
                        best.method = "youden", transpose = FALSE)
  opt_threshold <- as.numeric(best_coords$threshold[1])
  pred_05 <- ifelse(prob > 0.5, 1, 0)
  cm_05 <- confusionMatrix(factor(pred_05, levels = c(0, 1)),
                           factor(y_true, levels = c(0, 1)), positive = "1")
  cal_df <- data.frame(prob = prob, y = y_true)
  cal_summary <- tryCatch({
    breaks <- unique(quantile(prob, probs = seq(0, 1, length.out = 11)))
    breaks <- unique(c(min(prob) - 1e-8, breaks, max(prob) + 1e-8))
    cal_df$bin <- cut(cal_df$prob, breaks = breaks, include.lowest = TRUE, dig.lab = 4)
    cal_df %>%
      group_by(bin) %>%
      summarise(n = n(), mean_pred = mean(prob), obs_rate = mean(y), .groups = "drop") %>%
      filter(!is.na(bin))
  }, error = function(e) {
    data.frame(bin = "all", n = length(prob), mean_pred = mean(prob), obs_rate = mean(y_true))
  })
  list(AUC = auc, AUPRC = auprc, OptThreshold = opt_threshold,
       Sensitivity_0.5 = as.numeric(cm_05$byClass["Sensitivity"]),
       Specificity_0.5 = as.numeric(cm_05$byClass["Specificity"]),
       Accuracy_0.5 = as.numeric(cm_05$overall["Accuracy"]),
       Kappa_0.5 = as.numeric(cm_05$overall["Kappa"]),
       roc = roc_obj, pr_curve = pr_obj, calibration = cal_summary)
}

# 3. 外部验证核心函数 ----------------------------------------------------------
validate_external <- function(name, exp_mat, group_vec) {
  cat("\n==========", name, "==========\n")
  cat("Samples:", length(group_vec), "| Disease:", sum(group_vec == "Disease"),
      "| Normal:", sum(group_vec == "Normal"), "\n")
  y_ext <- ifelse(group_vec == "Disease", 1, 0)
  common_features <- intersect(selected, rownames(exp_mat))
  cat("Model features present:", length(common_features), "/", length(selected), "\n")
  if (length(common_features) < 10) { cat("Too few features, skipping\n"); return(NULL) }

  # GSE48166 / GSE46224 量纲不兼容, 改用队列内部 z-score(与 step21 相同约定)
  per_cohort_z_datasets <- c("GSE48166", "GSE46224")
  use_frozen <- !(name %in% per_cohort_z_datasets)
  cat("Scaling:", ifelse(use_frozen, "frozen(train)", "per-cohort z-score"), "\n")

  X_full <- matrix(0, nrow = ncol(exp_mat), ncol = length(selected),
                   dimnames = list(colnames(exp_mat), selected))
  if (use_frozen) {
    X_full[, common_features] <- t(exp_mat[common_features, , drop = FALSE])
    X_ext_scaled <- scale(X_full, center = train_center[selected], scale = train_scale[selected])
    X_ext_scaled[is.nan(X_ext_scaled)] <- 0
  } else {
    X_full[, common_features] <- t(scale(t(exp_mat[common_features, , drop = FALSE])))
    X_ext_scaled <- X_full
    X_ext_scaled[is.nan(X_ext_scaled)] <- 0
  }

  prob <- as.numeric(predict(model_rf, newdata = as.data.frame(X_ext_scaled), type = "prob")[, "Disease"])
  names(prob) <- colnames(exp_mat)
  ev <- evaluate_model(prob, y_ext)
  cat("AUC:", round(ev$AUC, 3), "| AUPRC:", round(ev$AUPRC, 3), "\n")
  list(Dataset = name, prob = prob, y = y_ext, eval = ev, n_features = length(common_features),
       scaling = ifelse(use_frozen, "frozen", "per_cohort_z"))
}

read_symbol_cohort <- function(ds) {
  exp <- read.csv(paste0(ds, "_symbol.csv"), header = TRUE, row.names = 1, check.names = FALSE)
  grp <- read.csv(paste0(ds, "_group.csv"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  grp$group <- trimws(as.character(grp$group))
  common <- intersect(colnames(exp), grp$geo_accession)
  list(exp_mat = as.matrix(exp[, common, drop = FALSE]),
       group_vec = setNames(grp$group[match(common, grp$geo_accession)], common))
}

results <- list()

# 4. 四个 _symbol.csv 队列 ------------------------------------------------------
for (ds in c("GSE26887", "GSE46224", "GSE52601", "GSE48166")) {
  dat <- read_symbol_cohort(ds)
  results[[ds]] <- validate_external(ds, dat$exp_mat, dat$group_vec)
}

# 5. GSE116250(RNA-seq RPKM, 本地 series matrix 离线取分组) ----------------------
exp116250 <- read.table("GSE116250_rpkm.txt.gz", header = TRUE, sep = "\t", check.names = FALSE, quote = "")
exp116250 <- exp116250 %>% filter(Common_name != "")
exp116250 <- exp116250[!duplicated(exp116250$Common_name), ]
rownames(exp116250) <- exp116250$Common_name
exp116250 <- as.matrix(exp116250[, -(1:2)])

gse <- getGEO(filename = "GSE116250_series_matrix.txt.gz", getGPL = FALSE)
pd <- pData(phenoData(gse))
disease_col <- grep("disease", names(pd), value = TRUE)[1]
pd$disease <- pd[[disease_col]]
pd <- pd %>% filter(disease %in% c("non-failing", "ischemic cardiomyopathy"))
pd$group <- ifelse(pd$disease == "non-failing", "Normal", "Disease")

title_to_geo <- data.frame(title = pd$title, geo_accession = rownames(pd), group = pd$group,
                           stringsAsFactors = FALSE)
colnames(exp116250) <- trimws(colnames(exp116250))
matched <- title_to_geo %>% filter(title %in% colnames(exp116250))
exp_val <- log2(exp116250[, matched$title, drop = FALSE] + 1)
results[["GSE116250"]] <- validate_external("GSE116250", exp_val, matched$group)

# 6. GSE55296(RNA-seq counts) ---------------------------------------------------
counts <- read.delim(gzfile("GSE55296_count_data.txt.gz"), check.names = FALSE, stringsAsFactors = FALSE)
colnames(counts)[1:2] <- c("Ensembl", "Symbol")
colnames(counts) <- make.names(colnames(counts), unique = TRUE)
counts <- counts[, colSums(!is.na(counts) & counts != "") > 0]
counts <- counts %>% filter(Symbol != "") %>% filter(!duplicated(Symbol))
rownames(counts) <- counts$Symbol
exp_counts <- as.matrix(counts[, -(1:2)])
storage.mode(exp_counts) <- "numeric"

readme <- read.delim("GSE55296_processed_data_readme.txt", check.names = FALSE, stringsAsFactors = FALSE)
colnames(readme) <- c("geo_accession", "title", "col_header")
readme$group <- ifelse(grepl("Control", readme$title, ignore.case = TRUE), "Normal",
                       ifelse(grepl("Ischemic", readme$title, ignore.case = TRUE), "Disease", NA))
readme <- readme %>% filter(!is.na(group)) %>% filter(col_header %in% colnames(exp_counts))
exp55296 <- log2(exp_counts[, readme$col_header, drop = FALSE] + 1)
results[["GSE55296"]] <- validate_external("GSE55296", exp55296, readme$group)

# 7. GSE42955(affymetrix hugene10st) --------------------------------------------
lines <- readLines(gzfile("GSE42955_series_matrix.txt.gz"))
table_start <- which(grepl("!series_matrix_table_begin", lines))
expr <- read.delim(gzfile("GSE42955_series_matrix.txt.gz"),
                   skip = table_start, comment.char = "!",
                   check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr) <- expr[, 1]
expr <- as.matrix(expr[, -1, drop = FALSE])

probe_symbols <- mapIds(hugene10sttranscriptcluster.db, keys = rownames(expr), column = "SYMBOL",
                        keytype = "PROBEID", multiVals = "first")
probe_df <- data.frame(probe = rownames(expr), symbol = probe_symbols, stringsAsFactors = FALSE) %>%
  filter(!is.na(symbol) & symbol != "")
expr_sym <- expr[probe_df$probe, , drop = FALSE]
rownames(expr_sym) <- probe_df$symbol
expr_sym <- expr_sym[order(rowMeans(expr_sym, na.rm = TRUE), decreasing = TRUE), ]
expr_sym <- expr_sym[!duplicated(rownames(expr_sym)), ]

parse_meta <- function(line) gsub('"', "", strsplit(line, "\t")[[1]][-1])
meta_lines <- lines[1:table_start]
source_vals <- parse_meta(meta_lines[grep("^!Sample_source_name_ch1", meta_lines)])
geo_vals <- parse_meta(meta_lines[grep("^!Sample_geo_accession", meta_lines)])
pd42955 <- data.frame(geo_accession = geo_vals, source = source_vals, stringsAsFactors = FALSE)
pd42955$group <- ifelse(grepl("Normal heart", pd42955$source, ignore.case = TRUE), "Normal",
                        ifelse(grepl("Ischemic cardiomyopathy", pd42955$source, ignore.case = TRUE), "Disease", NA))
pd42955 <- pd42955 %>% filter(!is.na(group)) %>% filter(geo_accession %in% colnames(expr_sym))
results[["GSE42955"]] <- validate_external("GSE42955", expr_sym[, pd42955$geo_accession, drop = FALSE],
                                           pd42955$group)

# 8. 汇总 ----------------------------------------------------------------------
results <- results[!sapply(results, is.null)]
perf_df <- bind_rows(lapply(results, function(r) {
  data.frame(Dataset = r$Dataset, n_samples = length(r$y),
             n_Disease = sum(r$y == 1), n_Normal = sum(r$y == 0),
             features_used = r$n_features, scaling = r$scaling,
             AUC = r$eval$AUC, AUPRC = r$eval$AUPRC,
             Sensitivity_0.5 = r$eval$Sensitivity_0.5,
             Specificity_0.5 = r$eval$Specificity_0.5,
             Accuracy_0.5 = r$eval$Accuracy_0.5,
             Kappa_0.5 = r$eval$Kappa_0.5,
             Threshold_youden = r$eval$OptThreshold)
}))
write.table(perf_df, file.path(out_dir, "external_validation_HF_signature_performance.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n"); print(perf_df)

# ROC
roc_df <- bind_rows(lapply(results, function(r) {
  data.frame(Dataset = paste0(r$Dataset, " (AUC=", round(r$eval$AUC, 3), ")"),
             FPR = 1 - r$eval$roc$specificities, TPR = r$eval$roc$sensitivities)
}))
p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Dataset)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  theme_bw() +
  labs(title = "HF signature RF model - external validation ROC (split-first)",
       x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  theme(legend.text = element_text(size = 8))
ggsave(file.path(out_dir, "HF_signature_external_validation_ROC.pdf"), p_roc, width = 8, height = 6)

# PR
pr_df <- bind_rows(lapply(results, function(r) {
  data.frame(Dataset = r$Dataset,
             Recall = r$eval$pr_curve$curve[, 1], Precision = r$eval$pr_curve$curve[, 2])
}))
p_pr <- ggplot(pr_df, aes(x = Recall, y = Precision, color = Dataset)) +
  geom_line(linewidth = 0.9) + theme_bw() +
  labs(title = "HF signature RF model - external validation PR curve (split-first)",
       x = "Recall", y = "Precision", color = NULL) +
  theme(legend.text = element_text(size = 8))
ggsave(file.path(out_dir, "HF_signature_external_validation_PR.pdf"), p_pr, width = 8, height = 6)

# Calibration
cal_df <- bind_rows(lapply(results, function(r) data.frame(Dataset = r$Dataset, r$eval$calibration)))
p_cal <- ggplot(cal_df, aes(x = mean_pred, y = obs_rate, color = Dataset)) +
  geom_point(aes(size = n), alpha = 0.7) + geom_line(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Dataset, ncol = 4) + theme_bw() +
  labs(title = "HF signature RF model - external validation calibration (split-first)",
       x = "Mean predicted probability", y = "Observed event rate", size = "n") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
ggsave(file.path(out_dir, "HF_signature_external_validation_calibration.pdf"), p_cal, width = 12, height = 6)

# 逐样本预测明细
pred_all <- bind_rows(lapply(results, function(r) {
  data.frame(Dataset = r$Dataset, Sample = names(r$prob),
             Prob_Disease = unname(r$prob), TrueLabel = ifelse(r$y == 1, "Disease", "Normal"))
}))
write.csv(pred_all, file.path(out_dir, "external_validation_predictions.csv"), row.names = FALSE)

cat("\nDone. Outputs in", out_dir, "\n")
