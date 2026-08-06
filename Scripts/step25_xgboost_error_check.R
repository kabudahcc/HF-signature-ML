# step25: reproduce the caret XGBoost failure to document the actual error.
# Run with cwd = loco_meta_analysis/.
suppressPackageStartupMessages({
  library(glmnet); library(pROC); library(sva); library(caret)
})
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

merged <- do.call(cbind, lapply(exp_list, function(m) m[feature_genes, , drop = FALSE]))
sample_dataset <- setNames(unlist(lapply(datasets, function(ds) rep(ds, ncol(exp_list[[ds]])))),
                           unlist(lapply(exp_list, colnames)))
y <- factor(unlist(lapply(datasets, function(ds) unname(group_list[[ds]]))),
            levels = c("Normal", "Disease"))
train_idx <- createDataPartition(y, p = 0.7, list = FALSE)[, 1]
raw_train <- merged[, train_idx, drop = FALSE]
y_train <- y[train_idx]
batch_train <- factor(sample_dataset[colnames(raw_train)], levels = datasets)
mod_train <- model.matrix(~ y_train)
combat_train <- ComBat(dat = raw_train, batch = batch_train,
                       mod = mod_train, par.prior = TRUE)
g_mean <- rowMeans(combat_train)
g_sd <- apply(combat_train, 1, sd); g_sd[g_sd == 0] <- 1
X_train <- scale(t(combat_train), center = g_mean, scale = g_sd)
y_bin <- ifelse(y_train == "Disease", 1, 0)
cv_lasso <- cv.glmnet(X_train, y_bin, family = "binomial", alpha = 1,
                      nfolds = 10, type.measure = "auc")
coef_l <- coef(cv_lasso, s = "lambda.min")
sel <- setdiff(rownames(coef_l)[coef_l[, 1] != 0], "(Intercept)")
cat("LASSO selected", length(sel), "features\n")

train_df <- as.data.frame(X_train[, sel, drop = FALSE])
train_df$Group <- y_train
tr_ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                        classProbs = TRUE, summaryFunction = twoClassSummary)

for (m in c("xgbTree", "xgbLinear")) {
  cat("\n=== trying", m, "===\n")
  warns <- character(0)
  res <- tryCatch(
    withCallingHandlers({
      fit <- train(Group ~ ., data = train_df, method = m,
                   metric = "ROC", trControl = tr_ctrl, tuneLength = 3)
      cat("SUCCESS, best ROC:", max(fit$results$ROC), "\n")
    }, warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      cat("ERROR MESSAGE:\n")
      cat(conditionMessage(e), "\n")
    })
  if (length(warns)) {
    cat("WARNINGS (unique):\n")
    for (w in unique(warns)) cat("-", w, "\n")
  }
}
