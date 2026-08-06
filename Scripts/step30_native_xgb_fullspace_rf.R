# step30: reviewer round-3 computational items.
#  (A) native xgboost API benchmark on the 39 LASSO-preselected features
#      (replacing the failed caret xgbTree/xgbLinear interfaces), and
#  (B) random forest on the FULL 987-feature space (no LASSO pre-selection),
#      to test whether LASSO pre-selection biased the algorithm comparison.
# Same split-first protocol as step20b: identical seed/partition, ComBat fitted
# on train only with frozen affine coefficients applied to test, frozen z-score.
# Run with cwd = loco_meta_analysis/.
suppressPackageStartupMessages({
  library(glmnet); library(randomForest); library(pROC); library(sva)
  library(xgboost); library(caret)
})
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
set.seed(12345)
datasets <- c("GSE16499", "GSE5406", "GSE57338")

sig_up <- read.table("data/HF_signature_up.txt", header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)$gene
sig_down <- read.table("data/HF_signature_down.txt", header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)$gene
sig_genes <- unique(c(sig_up, sig_down))

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
common_genes <- Reduce(intersect, lapply(exp_list, rownames))
feature_genes <- intersect(sig_genes, common_genes)
cat("feature space:", length(feature_genes), "genes\n")

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
cat("Train:", ncol(raw_train), "| Test:", ncol(raw_test), "\n")

## frozen ComBat (affine coefficients recovered per gene x cohort from train)
combat_train <- ComBat(dat = raw_train, batch = batch_train,
                       mod = model.matrix(~ y_train), par.prior = TRUE)
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
g_mean <- rowMeans(combat_train)
g_sd <- apply(combat_train, 1, sd); g_sd[g_sd == 0] <- 1
X_train <- scale(t(combat_train), center = g_mean, scale = g_sd)
X_test  <- scale(t(combat_test),  center = g_mean, scale = g_sd)
ytr <- ifelse(y_train == "Disease", 1, 0)

## LASSO pre-selection (train only) -----------------------------------------
cv_lasso <- cv.glmnet(X_train, ytr, family = "binomial", alpha = 1,
                      nfolds = 10, type.measure = "auc")
coef_l <- coef(cv_lasso, s = "lambda.min")
sel <- setdiff(rownames(coef_l)[coef_l[, 1] != 0], "(Intercept)")
cat("LASSO selected", length(sel), "features\n")
stopifnot(length(sel) == 39, ncol(raw_train) == 388, ncol(raw_test) == 165)

## helpers -------------------------------------------------------------------
rep_cv_auc <- function(make_pred, folds_list) {
  # folds_list: list of test-index vectors; make_pred(tr_i, te_i) -> prob on te
  aucs <- c()
  for (te_i in folds_list) {
    tr_i <- setdiff(seq_len(nrow(X_train)), te_i)
    prob <- make_pred(tr_i, te_i)
    r <- roc(factor(ytr[te_i], levels = c(0, 1)), prob, quiet = TRUE, direction = "<")
    aucs <- c(aucs, as.numeric(auc(r)))
  }
  aucs
}
set.seed(2024)
folds <- createMultiFolds(y_train, k = 5, times = 3)  # 15 resamples, like caret repeatedcv

test_auc_ci <- function(prob) {
  r <- roc(factor(ifelse(y_test == "Disease", 1, 0), levels = c(0, 1)), prob,
           quiet = TRUE, direction = "<")
  ci <- ci.auc(r, method = "delong")
  c(AUC = as.numeric(auc(r)), lo = ci[1], hi = min(ci[3], 1))
}

results <- list()

## (A) native XGBoost on the 39 LASSO features -------------------------------
xgb_fit <- function(Xtr, ytr_v, nrounds = 150) {
  d <- xgb.DMatrix(data = Xtr, label = ytr_v)
  xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc",
                          max_depth = 3, eta = 0.1, subsample = 0.8,
                          colsample_bytree = 0.8),
            data = d, nrounds = nrounds, verbose = 0)
}
# nrounds via xgb.cv on the full training set (39 features)
dcv <- xgb.DMatrix(data = X_train[, sel], label = ytr)
xcv <- xgb.cv(params = list(objective = "binary:logistic", eval_metric = "auc",
                            max_depth = 3, eta = 0.1, subsample = 0.8,
                            colsample_bytree = 0.8),
              data = dcv, nrounds = 400, nfold = 5,
              early_stopping_rounds = 30, verbose = 0)
best_n <- xcv$best_iteration
if (is.null(best_n) || length(best_n) == 0 || is.na(best_n)) {
  best_n <- which.max(xcv$evaluation_log$test_auc_mean)
}
cat("XGBoost best nrounds:", best_n, "| cv auc:", max(xcv$evaluation_log$test_auc_mean), "\n")

cv_aucs_xgb <- rep_cv_auc(function(tr_i, te_i) {
  m <- xgb_fit(X_train[tr_i, sel], ytr[tr_i], nrounds = best_n)
  predict(m, xgb.DMatrix(X_train[te_i, sel]))
}, folds)
xgb_final <- xgb_fit(X_train[, sel], ytr, nrounds = best_n)
prob_xgb <- predict(xgb_final, xgb.DMatrix(X_test[, sel]))
ta_xgb <- test_auc_ci(prob_xgb)
results[["XGBoost_native_39feat"]] <- data.frame(
  analysis = "XGBoost (native API), 39 LASSO features",
  cv_auc_mean = mean(cv_aucs_xgb), cv_auc_sd = sd(cv_aucs_xgb),
  test_auc = ta_xgb["AUC"], test_lo = ta_xgb["lo"], test_hi = ta_xgb["hi"])

## (B) random forest on the full 987-feature space ---------------------------
# x/y interface: gene symbols like HLA-B are not valid formula names
cv_aucs_full <- rep_cv_auc(function(tr_i, te_i) {
  rf <- randomForest(x = X_train[tr_i, ], y = factor(ytr[tr_i]), ntree = 500)
  predict(rf, newdata = X_train[te_i, ], type = "prob")[, "1"]
}, folds)
rf_full <- randomForest(x = X_train, y = factor(ytr), ntree = 500)
prob_full <- predict(rf_full, newdata = X_test, type = "prob")[, "1"]
ta_full <- test_auc_ci(prob_full)
results[["RF_full_987"]] <- data.frame(
  analysis = "Random forest, full 987-feature space (no LASSO)",
  cv_auc_mean = mean(cv_aucs_full), cv_auc_sd = sd(cv_aucs_full),
  test_auc = ta_full["AUC"], test_lo = ta_full["lo"], test_hi = ta_full["hi"])

## reference: LASSO->RF on 39 features (should reproduce ~0.995) -------------
df39 <- as.data.frame(X_train[, sel]); df39$Group <- factor(ytr)
rf39 <- randomForest(Group ~ ., data = df39, ntree = 500)
prob_39 <- predict(rf39, newdata = as.data.frame(X_test[, sel]), type = "prob")[, "1"]
ta_39 <- test_auc_ci(prob_39)
results[["RF_39feat_reference"]] <- data.frame(
  analysis = "Random forest, 39 LASSO features (primary model, reference)",
  cv_auc_mean = NA_real_, cv_auc_sd = NA_real_,
  test_auc = ta_39["AUC"], test_lo = ta_39["lo"], test_hi = ta_39["hi"])

res <- do.call(rbind, results)
rownames(res) <- NULL
res <- transform(res,
  cv_auc_mean = round(cv_auc_mean, 4), cv_auc_sd = round(cv_auc_sd, 4),
  test_auc = round(test_auc, 4), test_lo = round(test_lo, 4), test_hi = round(test_hi, 4))
write.csv(res, "../revision_round2/step30_xgb_fullspace_results.csv", row.names = FALSE)
print(res)

# per-sample test predictions for the record
pred_out <- data.frame(Sample = colnames(raw_test),
                       Group = as.character(y_test),
                       Prob_XGBoost_39 = as.numeric(prob_xgb),
                       Prob_RF_full987 = as.numeric(prob_full),
                       Prob_RF_39 = as.numeric(prob_39))
write.csv(pred_out, "../revision_round2/step30_test_predictions.csv", row.names = FALSE)
cat("step30 done\n")
