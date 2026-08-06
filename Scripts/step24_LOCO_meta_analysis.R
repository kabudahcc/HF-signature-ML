# step24: Leave-one-cohort-out CV (discovery cohorts) + random-effects
# meta-analysis of external validation AUCs.
# Run with cwd = loco_meta_analysis/. Pure ASCII to avoid encoding issues.
# Inputs:
#   data/GSE*_symbol.csv, data/GSE*_group.csv (copies of discovery matrices)
#   data/HF_signature_up.txt, data/HF_signature_down.txt
#   ../Data/HF_signature_external_validation_predictions.csv
# Outputs:
#   LOCO_CV_results.txt, meta_analysis_results.txt,
#   S17_forest_plot_external_AUC.png, S6_dataset_LOCO_meta.csv

suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
  library(pROC)
  library(sva)
  library(ggplot2)
})
set.seed(12345)

datasets <- c("GSE16499", "GSE5406", "GSE57338")
max_features <- 50

## ---- load signature and expression matrices -----------------------------
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
  cat(ds, "samples:", length(common),
      "| Disease:", sum(group_list[[ds]] == "Disease"),
      "| Normal:", sum(group_list[[ds]] == "Normal"), "\n")
}
common_genes <- Reduce(intersect, lapply(exp_list, rownames))
feature_genes <- intersect(sig_genes, common_genes)
cat("Signature genes present in all 3 datasets:", length(feature_genes), "\n")

## ---- Part 1: leave-one-cohort-out CV -------------------------------------
loco_rows <- list()
for (held in datasets) {
  train_ds <- setdiff(datasets, held)
  cat("\n=== LOCO fold: held-out =", held, "| train =", paste(train_ds, collapse = "+"), "===\n")

  raw_train <- do.call(cbind, lapply(train_ds, function(ds) exp_list[[ds]][feature_genes, , drop = FALSE]))
  y_train <- factor(unlist(lapply(train_ds, function(ds) unname(group_list[[ds]]))),
                    levels = c("Normal", "Disease"))
  batch_train <- factor(unlist(lapply(train_ds, function(ds) rep(ds, ncol(exp_list[[ds]])))),
                        levels = train_ds)

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
  if (length(sel) > max_features) {
    ord <- order(abs(coef_l[sel, 1]), decreasing = TRUE)
    sel <- sel[ord][seq_len(max_features)]
  }
  if (length(sel) < 2) sel <- feature_genes
  cat("LASSO selected", length(sel), "features\n")

  train_df <- as.data.frame(X_train[, sel, drop = FALSE])
  train_df$Group <- y_train
  rf <- randomForest(Group ~ ., data = train_df, ntree = 500)

  raw_test <- exp_list[[held]][feature_genes, , drop = FALSE]
  X_test <- scale(t(raw_test), center = g_mean, scale = g_sd)
  y_test <- factor(unname(group_list[[held]]), levels = c("Normal", "Disease"))
  prob <- predict(rf, newdata = as.data.frame(X_test[, sel, drop = FALSE]),
                  type = "prob")[, "Disease"]

  roc_obj <- roc(response = y_test, predictor = as.numeric(prob),
                 levels = c("Normal", "Disease"), direction = "<", quiet = TRUE)
  auc_v <- as.numeric(auc(roc_obj))
  ci_v <- ci.auc(roc_obj, method = "delong")
  cat(held, "AUC =", round(auc_v, 3), "| 95% CI", round(ci_v[1], 3), "-", round(ci_v[3], 3), "\n")
  loco_rows[[held]] <- data.frame(
    Held_out_cohort = held,
    Train_cohorts = paste(train_ds, collapse = " + "),
    n_Disease = sum(y_test == "Disease"),
    n_Normal = sum(y_test == "Normal"),
    n_LASSO_features = length(sel),
    AUC = round(auc_v, 3),
    CI_lower = round(ci_v[1], 3),
    CI_upper = round(ci_v[3], 3),
    stringsAsFactors = FALSE)
}
loco_df <- do.call(rbind, loco_rows)
write.table(loco_df, "LOCO_CV_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nLOCO mean AUC:", round(mean(loco_df$AUC), 3),
    "| range:", round(min(loco_df$AUC), 3), "-", round(max(loco_df$AUC), 3), "\n")

## ---- Part 2: meta-analysis of external AUCs ------------------------------
pred <- read.csv("../Data/HF_signature_external_validation_predictions.csv",
                 stringsAsFactors = FALSE)
compat <- c("GSE116250", "GSE55296", "GSE26887", "GSE42955", "GSE52601")

boot_se <- function(labels, scores, B = 2000) {
  idx_d <- which(labels == "Disease"); idx_n <- which(labels == "Normal")
  vals <- numeric(B)
  for (b in seq_len(B)) {
    i_d <- sample(idx_d, replace = TRUE); i_n <- sample(idx_n, replace = TRUE)
    lb <- c(labels[i_d], labels[i_n]); sc <- c(scores[i_d], scores[i_n])
    vals[b] <- as.numeric(auc(roc(lb, sc, levels = c("Normal", "Disease"),
                                  direction = "<", quiet = TRUE)))
  }
  sd(vals)
}

ma_rows <- list()
for (ds in compat) {
  sub <- pred[pred$Dataset == ds, ]
  lb <- factor(sub$TrueLabel, levels = c("Normal", "Disease"))
  r <- roc(lb, sub$Prob_Disease, levels = c("Normal", "Disease"),
           direction = "<", quiet = TRUE)
  a <- as.numeric(auc(r))
  v <- tryCatch(as.numeric(var(r, method = "delong")), error = function(e) NA_real_)
  se_src <- "DeLong"
  if (is.na(v) || v <= 1e-10) {
    se <- boot_se(sub$TrueLabel, sub$Prob_Disease)
    se_src <- "bootstrap"
  } else {
    se <- sqrt(v)
  }
  n_d <- sum(lb == "Disease"); n_n <- sum(lb == "Normal")
  # continuity correction at the boundary (AUC == 1) for logit pooling
  n_tot <- n_d + n_n
  a_c <- if (a >= 1) (a * (n_tot - 1) + 0.5) / n_tot else a
  ci_d <- tryCatch(ci.auc(r, method = "delong"), error = function(e) c(NA, a, NA))
  ma_rows[[ds]] <- data.frame(
    Cohort = ds, n_Disease = n_d, n_Normal = n_n,
    AUC = round(a, 3), SE = round(se, 4), SE_source = se_src,
    CI_lower = round(ci_d[1], 3), CI_upper = round(min(ci_d[3], 1), 3),
    AUC_logit_input = round(a_c, 4), stringsAsFactors = FALSE)
}
ma_df <- do.call(rbind, ma_rows)
print(ma_df)

dl_pool <- function(a_c, se, label) {
  yi <- log(a_c / (1 - a_c))
  vi <- (se / (a_c * (1 - a_c)))^2
  wi <- 1 / vi
  fixed <- sum(wi * yi) / sum(wi)
  Q <- sum(wi * (yi - fixed)^2)
  k <- length(yi)
  C <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- max(0, (Q - (k - 1)) / C)
  wi_r <- 1 / (vi + tau2)
  pooled <- sum(wi_r * yi) / sum(wi_r)
  se_p <- sqrt(1 / sum(wi_r))
  lo <- pooled - 1.96 * se_p; hi <- pooled + 1.96 * se_p
  I2 <- max(0, (Q - (k - 1)) / Q) * 100
  p_het <- pchisq(Q, df = k - 1, lower.tail = FALSE)
  invlogit <- function(x) 1 / (1 + exp(-x))
  data.frame(Analysis = label, k = k,
             pooled_AUC = round(invlogit(pooled), 3),
             CI_lower = round(invlogit(lo), 3),
             CI_upper = round(invlogit(hi), 3),
             Q = round(Q, 2), I2_percent = round(I2, 1),
             p_heterogeneity = signif(p_het, 3),
             tau2 = round(tau2, 4))
}

pool_all5 <- dl_pool(ma_df$AUC_logit_input, ma_df$SE,
                     "Random-effects (DL), 5 scale-compatible cohorts")
big2 <- ma_df[ma_df$Cohort %in% c("GSE116250", "GSE55296"), ]
pool_big2 <- dl_pool(big2$AUC_logit_input, big2$SE,
                     "Random-effects (DL), two largest cohorts only")
pool_df <- rbind(pool_all5, pool_big2)
print(pool_df)

sink("meta_analysis_results.txt")
cat("Per-cohort AUC and standard errors (5 scale-compatible external cohorts)\n")
cat("SE source: DeLong where finite; 2000-replicate stratified bootstrap otherwise.\n")
cat("Boundary continuity correction applied to AUC = 1 cohorts for logit pooling: (AUC*(n-1)+0.5)/n\n\n")
write.table(ma_df, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nRandom-effects meta-analysis (DerSimonian-Laird, logit scale)\n\n")
write.table(pool_df, sep = "\t", quote = FALSE, row.names = FALSE)
sink()

## ---- forest plot ---------------------------------------------------------
fp <- ma_df[, c("Cohort", "AUC", "CI_lower", "CI_upper")]
fp$Label <- sprintf("%s (n=%d/%d)", fp$Cohort,
                    ma_df$n_Disease, ma_df$n_Normal)
fp$Type <- "Cohort"
pooled_row <- data.frame(Cohort = "Pooled (random effects)",
                         AUC = pool_all5$pooled_AUC,
                         CI_lower = pool_all5$CI_lower,
                         CI_upper = pool_all5$CI_upper,
                         Label = sprintf("Pooled AUC (random effects) = %.3f (%.3f-%.3f), I2 = %.1f%%",
                                         pool_all5$pooled_AUC, pool_all5$CI_lower,
                                         pool_all5$CI_upper, pool_all5$I2_percent),
                         Type = "Pooled")
fp <- rbind(fp, pooled_row)
fp$Label <- factor(fp$Label, levels = rev(fp$Label))
p <- ggplot(fp, aes(x = AUC, y = Label)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.25, color = "grey30") +
  geom_point(aes(size = Type, shape = Type), color = "#2166ac") +
  scale_size_manual(values = c(Cohort = 2.6, Pooled = 3.6)) +
  scale_shape_manual(values = c(Cohort = 16, Pooled = 18)) +
  xlim(0.4, 1.02) +
  labs(x = "AUC (95% CI)", y = NULL,
       title = "External validation AUCs and random-effects pooled estimate") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face = "bold"))
ggsave("S17_forest_plot_external_AUC.png", p, width = 7.5, height = 4.2, dpi = 300)

## ---- S6 dataset export ----------------------------------------------------
write.csv(ma_df, "S6_dataset_meta_analysis.csv", row.names = FALSE)
write.csv(loco_df, "S6_dataset_LOCO.csv", row.names = FALSE)
cat("\nAll done.\n")
