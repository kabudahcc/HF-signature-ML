# step23b: 嵌套交叉验证 (审稿意见: LASSO 预选与模型训练在同一训练集上,
# 需要 nested CV 以获得无偏的泛化误差估计, 并与原测试集 AUC 比较)
# 外层 5 折 x 10 次重复; 每个外层折内: z-score 参数仅由外层训练折估计,
# LASSO 重新选特征(上限 50), 再训练 RF(500 树) / LDA / 径向 SVM,
# 在外层测试折上评估 AUC。特征选择完全不接触外层测试折。

.libPaths(c("Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
  library(MASS)
  library(kernlab)
  library(caret)
  library(pROC)
  library(sva)
  library(ggplot2)
  library(dplyr)
})

set.seed(12345)
out_dir <- "step22_review_response"
if (!dir.exists(out_dir)) dir.create(out_dir)

datasets <- c("GSE16499", "GSE5406", "GSE57338")
train_ratio <- 0.7
max_features <- 50
outer_k <- 5
outer_repeats <- 10

# 1-3. 复现 step20b 至 ComBat 后的训练集 ------------------------------------
sig_up <- read.table("results_scRNA/step3_4_5_pathway_cross_signature/HF_signature_up.txt",
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)$gene
sig_down <- read.table("results_scRNA/step3_4_5_pathway_cross_signature/HF_signature_down.txt",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)$gene
sig_genes <- unique(c(sig_up, sig_down))

exp_list <- list(); group_list <- list()
for (ds in datasets) {
  exp <- read.csv(paste0(ds, "_symbol.csv"), header = TRUE, row.names = 1, check.names = FALSE)
  grp <- read.csv(paste0(ds, "_group.csv"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  grp$group <- trimws(as.character(grp$group))
  common <- intersect(colnames(exp), grp$geo_accession)
  exp_list[[ds]] <- exp[, common, drop = FALSE]
  group_list[[ds]] <- setNames(grp$group[match(common, grp$geo_accession)], common)
}
common_genes <- Reduce(intersect, lapply(exp_list, rownames))
feature_genes <- intersect(sig_genes, common_genes)
merged <- do.call(cbind, lapply(exp_list, function(m) m[feature_genes, , drop = FALSE]))
merged <- as.matrix(merged)
colnames(merged) <- unlist(lapply(exp_list, colnames))
sample_dataset <- setNames(unlist(lapply(datasets, function(ds) rep(ds, ncol(exp_list[[ds]])))),
                           unlist(lapply(exp_list, colnames)))
y <- factor(unlist(lapply(datasets, function(ds) unname(group_list[[ds]]))),
            levels = c("Normal", "Disease"))

train_idx <- createDataPartition(y, p = train_ratio, list = FALSE)[, 1]
raw_train <- merged[, train_idx, drop = FALSE]
y_train <- y[train_idx]
batch_train <- factor(sample_dataset[colnames(raw_train)], levels = datasets)
mod_train <- model.matrix(~ y_train)
combat_train <- ComBat(dat = as.matrix(raw_train), batch = batch_train,
                       mod = mod_train, par.prior = TRUE)
# combat_train: 基因 x 样本 (ComBat 校正后, 未 z-score —— z-score 在嵌套折内重新估计)
cat("Training matrix for nested CV:", nrow(combat_train), "genes x", ncol(combat_train), "samples\n")

# 4. 嵌套 CV -----------------------------------------------------------------
run_nested <- function(rep_id) {
  folds <- createFolds(y_train, k = outer_k, list = TRUE)
  rows <- list()
  for (fi in seq_along(folds)) {
    te <- folds[[fi]]; tr <- setdiff(seq_len(ncol(combat_train)), te)
    # z-score 参数仅来自外层训练折
    m <- rowMeans(combat_train[, tr, drop = FALSE])
    s <- apply(combat_train[, tr, drop = FALSE], 1, sd); s[s == 0] <- 1
    Xtr <- scale(t(combat_train[, tr, drop = FALSE]), center = m, scale = s)
    Xte <- scale(t(combat_train[, te, drop = FALSE]), center = m, scale = s)
    ytr <- y_train[tr]; yte <- y_train[te]
    ytr_bin <- ifelse(ytr == "Disease", 1, 0)
    # LASSO 仅在外层训练折内选特征
    cv_l <- cv.glmnet(Xtr, ytr_bin, family = "binomial", alpha = 1,
                      nfolds = 10, type.measure = "auc")
    co <- coef(cv_l, s = "lambda.min")
    sel <- setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
    if (length(sel) > max_features) sel <- sel[order(abs(co[sel, 1]), decreasing = TRUE)][seq_len(max_features)]
    if (length(sel) < 2) sel <- feature_genes
    dtr <- data.frame(Xtr[, sel, drop = FALSE], Group = ytr)
    dte <- data.frame(Xte[, sel, drop = FALSE])
    # RF
    fit_rf <- randomForest(Group ~ ., data = dtr, ntree = 500)
    p_rf <- predict(fit_rf, dte, type = "prob")[, "Disease"]
    # LDA
    fit_lda <- lda(Group ~ ., data = dtr)
    p_lda <- as.numeric(predict(fit_lda, dte)$posterior[, "Disease"])
    # radial SVM
    fit_svm <- ksvm(Group ~ ., data = dtr, kernel = "rbfdot", prob.model = TRUE)
    p_svm <- as.numeric(predict(fit_svm, dte, type = "probabilities")[, "Disease"])
    for (md in c("RF", "LDA", "SVM_radial")) {
      p <- switch(md, RF = p_rf, LDA = p_lda, SVM_radial = p_svm)
      r <- roc(yte, as.numeric(p), levels = c("Normal", "Disease"), direction = "<", quiet = TRUE)
      rows[[length(rows) + 1]] <- data.frame(Repeat = rep_id, Fold = fi, Model = md,
                                             AUC = as.numeric(auc(r)), N_features = length(sel))
    }
  }
  bind_rows(rows)
}

res <- bind_rows(lapply(seq_len(outer_repeats), function(r) {
  cat("Outer repeat", r, "of", outer_repeats, "\n")
  run_nested(r)
}))

write.table(res, file.path(out_dir, "nested_CV_AUC_by_fold.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 5. 汇总: 每次重复内先对 5 折取均值, 再跨重复汇总 (减少折间相关的过度自信) ----
per_rep <- res %>% group_by(Repeat, Model) %>%
  summarise(AUC = mean(AUC), N_features = mean(N_features), .groups = "drop")
summ <- per_rep %>% group_by(Model) %>%
  summarise(mean_AUC = mean(AUC),
            sd_AUC = sd(AUC),
            CI_lo = mean_AUC - 1.96 * sd_AUC / sqrt(n()),
            CI_hi = mean_AUC + 1.96 * sd_AUC / sqrt(n()),
            mean_N_features = mean(N_features), .groups = "drop")
print(summ)
write.table(summ, file.path(out_dir, "nested_CV_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 6. 图 ---------------------------------------------------------------------
p <- ggplot(per_rep, aes(x = Model, y = AUC, fill = Model)) +
  geom_boxplot(width = 0.5, alpha = 0.8, outlier.size = 0.8) +
  geom_hline(yintercept = 0.995, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.6, y = 0.995, label = "internal test AUC ≈ 0.995",
           hjust = 0, vjust = -0.6, size = 3, color = "grey40") +
  scale_fill_manual(values = c(RF = "#B2182B", LDA = "#2166AC", SVM_radial = "#4DAF4A")) +
  theme_bw(base_size = 12) +
  labs(title = "Nested cross-validation (5-fold outer × 10 repeats)",
       subtitle = "LASSO feature selection re-fitted inside every outer training fold",
       y = "AUC (mean over outer folds, per repeat)", x = NULL, fill = NULL)
ggsave(file.path(out_dir, "nested_CV_AUC_boxplot.pdf"), p, width = 6.5, height = 5)

cat("\nDone. Nested CV outputs in", out_dir, "\n")
