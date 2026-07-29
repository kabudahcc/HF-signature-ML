# step20b: HF 签名 -> 12 算法诊断模型 -> SHAP (split-first 修正版)
# 与 step20 的区别: 修复数据泄露 —— 先在原始矩阵上划分训练/测试集,
# ComBat 只在训练集上拟合, 其每基因x每批次的仿射校正系数冻存后应用到测试集;
# z-score 的均值/标准差同样只从训练集估计再应用到测试集。
# 输入: results_scRNA HF_signature_up/down.txt, GSE16499/5406/57338 _symbol.csv + _group.csv
# 输出: step20b_HF_signature_ML_splitfirst/ 下的 AUC 矩阵、ROC 曲线、模型文件、SHAP 图表

.libPaths(c("Rlib", .libPaths()))
rm(list = ls())

suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
  library(xgboost)
  library(e1071)
  library(caret)
  library(pROC)
  library(kernlab)
  library(sva)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(kernelshap)
  library(shapviz)
})

set.seed(12345)

out_dir <- "step20b_HF_signature_ML_splitfirst"
if (!dir.exists(out_dir)) dir.create(out_dir)

datasets <- c("GSE16499", "GSE5406", "GSE57338")
max_features <- 50        # Lasso 预选后进入 12 算法的最大基因数
train_ratio <- 0.7

# 1. 读取 HF 签名 --------------------------------------------------------
sig_up <- read.table("results_scRNA/step3_4_5_pathway_cross_signature/HF_signature_up.txt",
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)$gene
sig_down <- read.table("results_scRNA/step3_4_5_pathway_cross_signature/HF_signature_down.txt",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)$gene
sig_genes <- unique(c(sig_up, sig_down))
cat("HF signature genes:", length(sig_genes), "\n")

# 2. 读取三个 GEO 矩阵并取共有签名基因 ------------------------------------
exp_list <- list()
group_list <- list()
for (ds in datasets) {
  exp <- read.csv(paste0(ds, "_symbol.csv"), header = TRUE, row.names = 1, check.names = FALSE)
  grp <- read.csv(paste0(ds, "_group.csv"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  grp$group <- trimws(as.character(grp$group))
  common <- intersect(colnames(exp), grp$geo_accession)
  exp_list[[ds]] <- exp[, common, drop = FALSE]
  group_list[[ds]] <- setNames(grp$group[match(common, grp$geo_accession)], common)
  cat(ds, "samples:", length(common), "| Disease:", sum(group_list[[ds]] == "Disease"),
      "| Normal:", sum(group_list[[ds]] == "Normal"), "\n")
}

common_genes <- Reduce(intersect, lapply(exp_list, rownames))
feature_genes <- intersect(sig_genes, common_genes)
cat("Signature genes present in all 3 datasets:", length(feature_genes), "\n")
if (length(feature_genes) < 10) stop("共有签名基因太少，无法建模")

merged <- do.call(cbind, lapply(exp_list, function(m) m[feature_genes, , drop = FALSE]))
merged <- as.matrix(merged)   # cbind data.frame 得到的是 data.frame, 统一成矩阵
colnames(merged) <- unlist(lapply(exp_list, colnames))
sample_dataset <- setNames(unlist(lapply(datasets, function(ds) rep(ds, ncol(exp_list[[ds]])))),
                           unlist(lapply(exp_list, colnames)))
y <- factor(unlist(lapply(datasets, function(ds) unname(group_list[[ds]]))),
            levels = c("Normal", "Disease"))
stopifnot(length(y) == ncol(merged), !any(is.na(y)))

# 3. 先划分训练/测试集(在原始矩阵上, 同一种子) -----------------------------
train_idx <- createDataPartition(y, p = train_ratio, list = FALSE)[, 1]
raw_train <- merged[, train_idx, drop = FALSE]
raw_test  <- merged[, -train_idx, drop = FALSE]
y_train <- y[train_idx]
y_test  <- y[-train_idx]
batch_train <- factor(sample_dataset[colnames(raw_train)], levels = datasets)
batch_test  <- factor(sample_dataset[colnames(raw_test)], levels = datasets)
cat("Train:", ncol(raw_train), "| Test:", ncol(raw_test), "\n")
cat("Train batch counts:", paste(names(table(batch_train)), table(batch_train), collapse = ", "), "\n")

# 4. ComBat 只拟合训练集, 仿射校正系数冻存后应用到测试集 -------------------
mod_train <- model.matrix(~ y_train)
combat_train <- ComBat(dat = as.matrix(raw_train), batch = batch_train,
                       mod = mod_train, par.prior = TRUE)

# ComBat 的位置-尺度调整对每(基因, 批次)是仿射变换: adjusted = a * raw + c
# 用训练集样本按批次恢复 a, c, 再应用到测试集对应批次
combat_test <- matrix(NA_real_, nrow = nrow(raw_test), ncol = ncol(raw_test),
                      dimnames = dimnames(raw_test))
for (b in datasets) {
  tr_cols <- which(batch_train == b)
  te_cols <- which(batch_test == b)
  if (length(te_cols) == 0) next
  m_raw <- rowMeans(raw_train[, tr_cols, drop = FALSE])
  m_adj <- rowMeans(combat_train[, tr_cols, drop = FALSE])
  v_raw <- apply(raw_train[, tr_cols, drop = FALSE], 1, var)
  covv  <- sapply(seq_len(nrow(raw_train)), function(g)
                    cov(raw_train[g, tr_cols], combat_train[g, tr_cols]))
  a <- ifelse(v_raw > 0, covv / v_raw, 1)
  c0 <- m_adj - a * m_raw
  combat_test[, te_cols] <- sweep(sweep(raw_test[, te_cols, drop = FALSE], 1, a, `*`), 1, c0, `+`)
}

# 5. z-score: 均值/标准差只从训练集(ComBat 后)估计, 应用到测试集 -----------
g_mean <- rowMeans(combat_train)
g_sd   <- apply(combat_train, 1, sd)
g_sd[g_sd == 0] <- 1
X_train <- scale(t(combat_train), center = g_mean, scale = g_sd)
X_test  <- scale(t(combat_test),  center = g_mean, scale = g_sd)

train_set <- as.data.frame(X_train)
train_set$Group <- y_train
train_set$Dataset <- batch_train
test_set <- as.data.frame(X_test)
test_set$Group <- y_test
test_set$Dataset <- batch_test

x_train <- as.matrix(train_set[, feature_genes])
y_train_bin <- ifelse(train_set$Group == "Disease", 1, 0)

# 6. Lasso 预选特征(仅训练集) ---------------------------------------------
cv_lasso <- cv.glmnet(x_train, y_train_bin, family = "binomial", alpha = 1,
                      nfolds = 10, type.measure = "auc")
coef_l <- coef(cv_lasso, s = "lambda.min")
sel <- rownames(coef_l)[coef_l[, 1] != 0]
sel <- setdiff(sel, "(Intercept)")
if (length(sel) > max_features) {
  ord <- order(abs(coef_l[sel, 1]), decreasing = TRUE)
  sel <- sel[ord][seq_len(max_features)]
}
if (length(sel) < 2) {
  warning("Lasso 选中基因过少, 退而用全部签名基因")
  sel <- feature_genes
}
cat("Lasso selected", length(sel), "features\n")
write.table(data.frame(gene = sel, coef = as.numeric(coef_l[sel, 1])),
            file.path(out_dir, "model.genes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# 7. 12 算法训练与评估 ----------------------------------------------------
tr_ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                        classProbs = TRUE, summaryFunction = twoClassSummary)

fit_one <- function(method_id, train_df) {
  switch(method_id,
    "Lasso"      = glmnet(as.matrix(train_df[, sel]), ifelse(train_df$Group == "Disease", 1, 0),
                          family = "binomial", alpha = 1, lambda = cv.glmnet(as.matrix(train_df[, sel]),
                          ifelse(train_df$Group == "Disease", 1, 0), family = "binomial",
                          alpha = 1, nfolds = 5, type.measure = "auc")$lambda.min),
    "Ridge"      = glmnet(as.matrix(train_df[, sel]), ifelse(train_df$Group == "Disease", 1, 0),
                          family = "binomial", alpha = 0, lambda = cv.glmnet(as.matrix(train_df[, sel]),
                          ifelse(train_df$Group == "Disease", 1, 0), family = "binomial",
                          alpha = 0, nfolds = 5, type.measure = "auc")$lambda.min),
    "Enet"       = glmnet(as.matrix(train_df[, sel]), ifelse(train_df$Group == "Disease", 1, 0),
                          family = "binomial", alpha = 0.5, lambda = cv.glmnet(as.matrix(train_df[, sel]),
                          ifelse(train_df$Group == "Disease", 1, 0), family = "binomial",
                          alpha = 0.5, nfolds = 5, type.measure = "auc")$lambda.min),
    "Stepglm"    = step(glm(Group ~ ., data = train_df[, c(sel, "Group")], family = binomial),
                        direction = "both", trace = 0),
    "SVM_radial" = train(Group ~ ., data = train_df[, c(sel, "Group")], method = "svmRadial",
                         metric = "ROC", trControl = tr_ctrl, tuneLength = 3),
    "SVM_linear" = train(Group ~ ., data = train_df[, c(sel, "Group")], method = "svmLinear",
                         metric = "ROC", trControl = tr_ctrl, tuneLength = 3),
    "LDA"        = train(Group ~ ., data = train_df[, c(sel, "Group")], method = "lda",
                         metric = "ROC", trControl = tr_ctrl),
    "RF"         = randomForest(Group ~ ., data = train_df[, c(sel, "Group")], ntree = 500),
    "XGB_tree"   = train(Group ~ ., data = train_df[, c(sel, "Group")], method = "xgbTree",
                         metric = "ROC", trControl = tr_ctrl, tuneLength = 3, verbosity = 0),
    "XGB_linear" = train(Group ~ ., data = train_df[, c(sel, "Group")], method = "xgbLinear",
                         metric = "ROC", trControl = tr_ctrl, tuneLength = 3, verbosity = 0),
    "KNN"        = train(Group ~ ., data = train_df[, c(sel, "Group")], method = "knn",
                         metric = "ROC", trControl = tr_ctrl, tuneLength = 5),
    "NaiveBayes" = naiveBayes(Group ~ ., data = train_df[, c(sel, "Group")]),
    stop("unknown method")
  )
}

predict_prob <- function(fit, method_id, new_df) {
  x_new <- new_df[, sel, drop = FALSE]
  p <- switch(method_id,
    "Lasso"      = as.numeric(predict(fit, newx = as.matrix(x_new), type = "response")),
    "Ridge"      = as.numeric(predict(fit, newx = as.matrix(x_new), type = "response")),
    "Enet"       = as.numeric(predict(fit, newx = as.matrix(x_new), type = "response")),
    "Stepglm"    = as.numeric(predict(fit, newdata = x_new, type = "response")),
    "RF"         = as.numeric(predict(fit, newdata = x_new, type = "prob")[, "Disease"]),
    "NaiveBayes" = as.numeric(predict(fit, newdata = x_new, type = "raw")[, "Disease"]),
    as.numeric(predict(fit, newdata = x_new, type = "prob")[, "Disease"])
  )
  p
}

methods12 <- c("Lasso", "Ridge", "Enet", "Stepglm", "SVM_radial", "SVM_linear",
               "LDA", "RF", "XGB_tree", "XGB_linear", "KNN", "NaiveBayes")

fits <- list()
auc_rows <- list()
test_probs <- data.frame(row.names = rownames(test_set))
for (m in methods12) {
  cat("Training:", m, "...\n")
  fit <- tryCatch(fit_one(m, train_set), error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) next
  fits[[m]] <- fit
  prob_test <- predict_prob(fit, m, test_set)
  test_probs[[m]] <- prob_test
  roc_test <- roc(test_set$Group, prob_test, levels = c("Normal", "Disease"), direction = "<", quiet = TRUE)
  auc_rows[[m]] <- data.frame(Model = m, AUC_test = as.numeric(auc(roc_test)))
  for (ds in datasets) {
    idx <- test_set$Dataset == ds
    if (length(unique(test_set$Group[idx])) < 2 || sum(idx) < 4) next
    r <- roc(test_set$Group[idx], prob_test[idx], levels = c("Normal", "Disease"), direction = "<", quiet = TRUE)
    auc_rows[[m]][[paste0("AUC_", ds)]] <- as.numeric(auc(r))
  }
}

auc_df <- bind_rows(auc_rows) %>% arrange(desc(AUC_test))
write.table(auc_df, file.path(out_dir, "model.AUCmatrix.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
print(auc_df)

# 8. ROC 曲线 + AUC 条形图 -------------------------------------------------
roc_plot_df <- bind_rows(lapply(names(fits), function(m) {
  prob <- test_probs[[m]]
  r <- roc(test_set$Group, prob, levels = c("Normal", "Disease"), direction = "<", quiet = TRUE)
  data.frame(Model = paste0(m, " (AUC=", round(as.numeric(auc(r)), 3), ")"),
             FPR = 1 - r$specificities, TPR = r$sensitivities)
}))
p_roc <- ggplot(roc_plot_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  theme_bw() +
  labs(title = "HF signature diagnostic models (test set, split-first)", x = "1 - Specificity", y = "Sensitivity") +
  theme(legend.position = "right", legend.text = element_text(size = 7))
ggsave(file.path(out_dir, "ROC_12models_test.pdf"), p_roc, width = 8, height = 6)

p_auc <- ggplot(auc_df, aes(x = reorder(Model, AUC_test), y = AUC_test, fill = AUC_test)) +
  geom_col() + coord_flip() +
  scale_fill_gradient(low = "#92C5DE", high = "#B2182B") +
  theme_bw() + labs(x = NULL, y = "AUC (test set)", title = "12 algorithm comparison (split-first)") +
  geom_text(aes(label = round(AUC_test, 3)), hjust = -0.1, size = 3) +
  ylim(0, 1.1)
ggsave(file.path(out_dir, "AUC_barplot_12models.pdf"), p_auc, width = 7, height = 5)

# 9. 最优模型 + SHAP ------------------------------------------------------
best_method <- auc_df$Model[1]
cat("\nBest model:", best_method, "AUC =", round(auc_df$AUC_test[1], 3), "\n")
best_fit <- fits[[best_method]]
saveRDS(list(model = best_fit, method = best_method, features = sel,
             train_set = train_set[, c(sel, "Group")], lasso_cv = cv_lasso,
             combat_params = list(g_mean = g_mean, g_sd = g_sd)),
        file.path(out_dir, "best_model.rds"))

pred_fun <- function(object, newdata) {
  predict_prob(object, best_method, as.data.frame(newdata))
}
X_bg <- train_set[, sel]
shap_fit <- kernelshap(best_fit, X = X_bg, bg_X = X_bg, pred_fun = pred_fun)
sv <- shapviz(shap_fit, X_pred = X_bg, X = X_bg)

shap_importance <- sort(colMeans(abs(sv$S)), decreasing = TRUE)
write.csv(data.frame(Gene = names(shap_importance), MeanAbsSHAP = as.numeric(shap_importance)),
          file.path(out_dir, "shap_importance.csv"), row.names = FALSE)
write.csv(as.data.frame(sv$S), file.path(out_dir, "shap_values_train.csv"), row.names = TRUE)

pdf(file.path(out_dir, "SHAP_importance_barplot.pdf"), width = 6, height = 6)
print(sv_importance(sv, kind = "bar", show_numbers = TRUE) + theme_bw())
dev.off()

pdf(file.path(out_dir, "SHAP_beeswarm.pdf"), width = 7, height = 6)
print(sv_importance(sv, kind = "bee", show_numbers = TRUE) + theme_bw())
dev.off()

top_n <- min(6, length(sel))
pdf(file.path(out_dir, "SHAP_dependence_top.pdf"), width = 10, height = 7)
print(sv_dependence(sv, v = names(shap_importance)[seq_len(top_n)]) + theme_bw())
dev.off()

pdf(file.path(out_dir, "SHAP_waterfall_sample1.pdf"), width = 7, height = 6)
print(sv_waterfall(sv, row_id = 1) + theme_bw())
dev.off()

# 10. 测试集预测明细 -------------------------------------------------------
pred_out <- data.frame(Sample = rownames(test_set), Dataset = test_set$Dataset,
                       Group = test_set$Group, test_probs, check.names = FALSE)
write.csv(pred_out, file.path(out_dir, "test_predictions.csv"), row.names = FALSE)

cat("\nDone. Outputs in", out_dir, "\n")
