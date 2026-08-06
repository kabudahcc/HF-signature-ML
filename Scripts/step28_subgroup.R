# step28: etiology subgroup analysis (ICM vs DCM vs non-failing).
# Re-runs the two well-powered LOCO folds (GSE5406, GSE57338) saving per-sample
# probabilities, then computes etiology-stratified AUCs; also stratifies the
# split-first internal test set. Sex/age metadata are not publicly available
# for these cohorts (GEO fetch blocked), so subgrouping is by etiology only.
# Run with cwd = loco_meta_analysis/.
suppressPackageStartupMessages({
  library(glmnet); library(randomForest); library(pROC); library(sva)
})
set.seed(12345)
datasets <- c("GSE16499", "GSE5406", "GSE57338")
max_features <- 50

sig_up <- read.table("data/HF_signature_up.txt", header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)$gene
sig_down <- read.table("data/HF_signature_down.txt", header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)$gene
sig_genes <- unique(c(sig_up, sig_down))

exp_list <- list(); group_list <- list(); title_list <- list()
for (ds in datasets) {
  exp <- read.csv(paste0("data/", ds, "_symbol.csv"), header = TRUE,
                  row.names = 1, check.names = FALSE)
  grp <- read.csv(paste0("data/", ds, "_group.csv"), header = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE)
  grp$group <- trimws(as.character(grp$group))
  common <- intersect(colnames(exp), grp$geo_accession)
  exp_list[[ds]] <- as.matrix(exp[, common, drop = FALSE])
  group_list[[ds]] <- setNames(grp$group[match(common, grp$geo_accession)], common)
  title_list[[ds]] <- setNames(grp$title[match(common, grp$geo_accession)], common)
}
common_genes <- Reduce(intersect, lapply(exp_list, rownames))
feature_genes <- intersect(sig_genes, common_genes)

parse_etiology <- function(ds, titles) {
  if (ds == "GSE5406") {
    ifelse(grepl("NonFailing", titles), "NonFailing",
           ifelse(grepl("Ischemic", titles), "ICM", "DCM"))
  } else if (ds == "GSE57338") {
    ifelse(grepl("Non-failing", titles), "NonFailing",
           ifelse(grepl("Ischemic", titles), "ICM", "DCM"))
  } else {
    ifelse(grepl("Normal", titles), "NonFailing", "ICM")
  }
}

loco_pred <- function(held) {
  train_ds <- setdiff(datasets, held)
  raw_train <- do.call(cbind, lapply(train_ds, function(ds) exp_list[[ds]][feature_genes, , drop = FALSE]))
  y_train <- factor(unlist(lapply(train_ds, function(ds) unname(group_list[[ds]]))),
                    levels = c("Normal", "Disease"))
  batch_train <- factor(unlist(lapply(train_ds, function(ds) rep(ds, ncol(exp_list[[ds]])))),
                        levels = train_ds)
  combat_train <- ComBat(dat = raw_train, batch = batch_train,
                         mod = model.matrix(~ y_train), par.prior = TRUE)
  g_mean <- rowMeans(combat_train)
  g_sd <- apply(combat_train, 1, sd); g_sd[g_sd == 0] <- 1
  X_train <- scale(t(combat_train), center = g_mean, scale = g_sd)
  cv_lasso <- cv.glmnet(X_train, ifelse(y_train == "Disease", 1, 0),
                        family = "binomial", alpha = 1, nfolds = 10,
                        type.measure = "auc")
  coef_l <- coef(cv_lasso, s = "lambda.min")
  sel <- setdiff(rownames(coef_l)[coef_l[, 1] != 0], "(Intercept)")
  if (length(sel) > max_features) {
    ord <- order(abs(coef_l[sel, 1]), decreasing = TRUE)
    sel <- sel[ord][seq_len(max_features)]
  }
  train_df <- as.data.frame(X_train[, sel, drop = FALSE]); train_df$Group <- y_train
  rf <- randomForest(Group ~ ., data = train_df, ntree = 500)
  raw_test <- exp_list[[held]][feature_genes, , drop = FALSE]
  X_test <- scale(t(raw_test), center = g_mean, scale = g_sd)
  prob <- predict(rf, newdata = as.data.frame(X_test[, sel, drop = FALSE]),
                  type = "prob")[, "Disease"]
  data.frame(Sample = colnames(raw_test), Dataset = held,
             Group = unname(group_list[[held]]),
             Etiology = parse_etiology(held, unname(title_list[[held]])),
             Prob = as.numeric(prob), stringsAsFactors = FALSE)
}

preds <- do.call(rbind, lapply(c("GSE5406", "GSE57338"), loco_pred))
write.csv(preds, "../revision_round2/subgroup_LOCO_predictions.csv", row.names = FALSE)

## internal test-set predictions (split-first)
tp <- read.csv("../Data/HF_signature_ML_test_predictions.csv", stringsAsFactors = FALSE)
tp$Etiology <- NA
for (ds in datasets) {
  idx <- tp$Dataset == ds
  tp$Etiology[idx] <- parse_etiology(ds, unname(title_list[[ds]][tp$Sample[idx]]))
}
tp$Prob <- tp$RF

auc_ci <- function(df, label) {
  n_d <- sum(df$Group == "Disease"); n_n <- sum(df$Group == "Normal")
  if (n_d < 2 || n_n < 2) return(NULL)
  r <- roc(factor(df$Group, levels = c("Normal", "Disease")), df$Prob,
           levels = c("Normal", "Disease"), direction = "<", quiet = TRUE)
  ci <- ci.auc(r, method = "delong")
  data.frame(Setting = label, n_Disease = n_d, n_Normal = n_n,
             AUC = round(as.numeric(auc(r)), 3),
             CI_lower = round(ci[1], 3), CI_upper = round(min(ci[3], 1), 3))
}

rows <- list()
for (ds in c("GSE5406", "GSE57338")) {
  sub <- preds[preds$Dataset == ds, ]
  for (et in c("ICM", "DCM")) {
    d2 <- sub[sub$Etiology %in% c(et, "NonFailing"), ]
    rows[[length(rows) + 1]] <- auc_ci(d2, paste0("LOCO held-out ", ds, ": ", et, " vs NonFailing"))
  }
}
for (ds in datasets) {
  sub <- tp[tp$Dataset == ds, ]
  for (et in c("ICM", "DCM")) {
    d2 <- sub[sub$Etiology %in% c(et, "NonFailing"), ]
    if (sum(d2$Etiology == et) < 2) next
    rows[[length(rows) + 1]] <- auc_ci(d2, paste0("Internal test ", ds, ": ", et, " vs NonFailing"))
  }
}
res <- do.call(rbind, rows)
write.csv(res, "../revision_round2/subgroup_etiology_AUC.csv", row.names = FALSE)
print(res)
cat("subgroup done\n")
