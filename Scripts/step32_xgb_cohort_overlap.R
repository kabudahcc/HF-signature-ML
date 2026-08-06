# step32: reviewer round-4 computational items.
#  (A) per-cohort test AUCs for the native-API XGBoost (Table 2 row)
#  (B) feature overlap between the full-987-space RF top-39 importance genes
#      and the LASSO-preselected 39 genes (equifinality check)
# Run with cwd = loco_meta_analysis/. Rebuilds the step30 pipeline (same seed).
suppressPackageStartupMessages({
  library(glmnet); library(randomForest); library(pROC); library(sva)
  library(xgboost); library(caret)
})
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
set.seed(12345)
datasets <- c("GSE16499", "GSE5406", "GSE57338")

sig_genes <- unique(c(
  read.table("data/HF_signature_up.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)$gene,
  read.table("data/HF_signature_down.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)$gene))
exp_list <- list(); group_list <- list()
for (ds in datasets) {
  exp <- read.csv(paste0("data/", ds, "_symbol.csv"), header = TRUE,
                  row.names = 1, check.names = FALSE)
  grp <- read.csv(paste0("data/", ds, "_group.csv"), header = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE)
  grp$group <- trimws(as.character(grp$group))
  common <- intersect(colnames(exp), grp$geo_accession)
  exp_list[[ds]] <- as.matrix(exp[, common, drop = FALSE])
  group_list[[ds]] <- setNames(grp$group[match(common, grp$geo_accession)], common)
}
feature_genes <- intersect(sig_genes, Reduce(intersect, lapply(exp_list, rownames)))
merged <- do.call(cbind, lapply(exp_list, function(m) m[feature_genes, , drop = FALSE]))
sample_dataset <- setNames(unlist(lapply(datasets, function(ds) rep(ds, ncol(exp_list[[ds]])))),
                           unlist(lapply(exp_list, colnames)))
y <- factor(unlist(lapply(datasets, function(ds) unname(group_list[[ds]]))),
            levels = c("Normal", "Disease"))
train_idx <- createDataPartition(y, p = 0.7, list = FALSE)[, 1]
raw_train <- merged[, train_idx, drop = FALSE]
raw_test  <- merged[, -train_idx, drop = FALSE]
y_train <- y[train_idx]; y_test <- y[-train_idx]
batch_train <- factor(sample_dataset[colnames(raw_train)], levels = datasets)
batch_test  <- factor(sample_dataset[colnames(raw_test)], levels = datasets)

combat_train <- ComBat(dat = raw_train, batch = batch_train,
                       mod = model.matrix(~ y_train), par.prior = TRUE)
combat_test <- matrix(NA_real_, nrow = nrow(raw_test), ncol = ncol(raw_test),
                      dimnames = dimnames(raw_test))
for (b in datasets) {
  tr_cols <- which(batch_train == b); te_cols <- which(batch_test == b)
  m_raw <- rowMeans(raw_train[, tr_cols, drop = FALSE])
  m_adj <- rowMeans(combat_train[, tr_cols, drop = FALSE])
  v_raw <- apply(raw_train[, tr_cols, drop = FALSE], 1, var)
  covv  <- sapply(seq_len(nrow(raw_train)), function(g) cov(raw_train[g, tr_cols], combat_train[g, tr_cols]))
  a <- ifelse(v_raw > 0, covv / v_raw, 1); c0 <- m_adj - a * m_raw
  combat_test[, te_cols] <- sweep(sweep(raw_test[, te_cols, drop = FALSE], 1, a, `*`), 1, c0, `+`)
}
g_mean <- rowMeans(combat_train); g_sd <- apply(combat_train, 1, sd); g_sd[g_sd == 0] <- 1
X_train <- scale(t(combat_train), center = g_mean, scale = g_sd)
X_test  <- scale(t(combat_test),  center = g_mean, scale = g_sd)
ytr <- ifelse(y_train == "Disease", 1, 0)

cv_lasso <- cv.glmnet(X_train, ytr, family = "binomial", alpha = 1,
                      nfolds = 10, type.measure = "auc")
coef_l <- coef(cv_lasso, s = "lambda.min")
sel <- setdiff(rownames(coef_l)[coef_l[, 1] != 0], "(Intercept)")
stopifnot(length(sel) == 39)

## (A) XGBoost native: per-cohort test AUCs -----------------------------------
dcv <- xgb.DMatrix(data = X_train[, sel], label = ytr)
xcv <- xgb.cv(params = list(objective = "binary:logistic", eval_metric = "auc",
                            max_depth = 3, eta = 0.1, subsample = 0.8,
                            colsample_bytree = 0.8),
              data = dcv, nrounds = 400, nfold = 5,
              early_stopping_rounds = 30, verbose = 0)
best_n <- xcv$best_iteration
if (is.null(best_n) || length(best_n) == 0) best_n <- which.max(xcv$evaluation_log$test_auc_mean)
xgb_final <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc",
                                     max_depth = 3, eta = 0.1, subsample = 0.8,
                                     colsample_bytree = 0.8),
                       data = dcv, nrounds = best_n, verbose = 0)
prob <- predict(xgb_final, xgb.DMatrix(X_test[, sel]))
per_cohort <- sapply(datasets, function(ds) {
  idx <- batch_test == ds
  r <- roc(factor(ifelse(y_test[idx] == "Disease", 1, 0), levels = c(0, 1)),
           prob[idx], quiet = TRUE, direction = "<")
  round(as.numeric(auc(r)), 3)
})
cat("XGBoost per-cohort test AUCs:\n"); print(per_cohort)

## (B) full-space RF top-39 importance vs LASSO 39 ----------------------------
rf_full <- randomForest(x = X_train, y = factor(ytr), ntree = 500,
                        importance = TRUE)
imp <- importance(rf_full, type = 1)[, 1]          # MeanDecreaseAccuracy
top39_full <- names(sort(imp, decreasing = TRUE))[1:39]
overlap <- intersect(top39_full, sel)
cat("full-space RF top-39 vs LASSO 39 overlap:", length(overlap), "genes\n")
cat("overlap genes:", paste(sort(overlap), collapse = ", "), "\n")

out <- data.frame(
  item = c(paste0("XGBoost_test_AUC_", datasets),
           "fullspaceRF_top39_overlap_n", "fullspaceRF_top39_overlap_genes"),
  value = c(as.character(per_cohort), as.character(length(overlap)),
            paste(sort(overlap), collapse = "; ")))
write.csv(out, "../revision_round2/step32_xgb_cohort_overlap.csv", row.names = FALSE)

# full importance table for the record
imp_df <- data.frame(gene = names(sort(imp, decreasing = TRUE)),
                     MeanDecreaseAccuracy = sort(imp, decreasing = TRUE),
                     in_LASSO39 = names(sort(imp, decreasing = TRUE)) %in% sel)
write.csv(imp_df, "../revision_round2/step32_fullspace_RF_importance.csv", row.names = FALSE)
cat("step32 done\n")
