# step24b: fixed meta-analysis of external validation AUCs.
# Cohorts with perfect separation (AUC = 1, zero DeLong and bootstrap
# variance) are excluded from inverse-variance pooling and reported
# individually. Run with cwd = loco_meta_analysis/.

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
})
set.seed(12345)

pred <- read.csv("../Data/HF_signature_external_validation_predictions.csv",
                 stringsAsFactors = FALSE)
compat <- c("GSE116250", "GSE55296", "GSE26887", "GSE42955", "GSE52601")

ma_rows <- list()
for (ds in compat) {
  sub <- pred[pred$Dataset == ds, ]
  lb <- factor(sub$TrueLabel, levels = c("Normal", "Disease"))
  r <- roc(lb, sub$Prob_Disease, levels = c("Normal", "Disease"),
           direction = "<", quiet = TRUE)
  a <- as.numeric(auc(r))
  v <- tryCatch(as.numeric(var(r, method = "delong")), error = function(e) NA_real_)
  n_d <- sum(lb == "Disease"); n_n <- sum(lb == "Normal")
  estimable <- !is.na(v) && v > 1e-10
  se <- if (estimable) sqrt(v) else NA_real_
  ci_d <- tryCatch(ci.auc(r, method = "delong"), error = function(e) c(NA, a, NA))
  ma_rows[[ds]] <- data.frame(
    Cohort = ds, n_Disease = n_d, n_Normal = n_n,
    AUC = round(a, 3),
    SE = ifelse(is.na(se), NA, round(se, 4)),
    Variance_estimable = estimable,
    CI_lower = round(ci_d[1], 3), CI_upper = round(min(ci_d[3], 1), 3),
    stringsAsFactors = FALSE)
}
ma_df <- do.call(rbind, ma_rows)
print(ma_df)

dl_pool <- function(a, se, label) {
  yi <- log(a / (1 - a))
  vi <- (se / (a * (1 - a)))^2
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
  I2 <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
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

est <- ma_df[ma_df$Variance_estimable, ]
pool_est <- dl_pool(est$AUC, est$SE,
  "Random-effects (DL), 3 cohorts with estimable variance")
big2 <- ma_df[ma_df$Cohort %in% c("GSE116250", "GSE55296"), ]
pool_big2 <- dl_pool(big2$AUC, big2$SE,
  "Random-effects (DL), two largest cohorts only")
pool_df <- rbind(pool_est, pool_big2)
print(pool_df)

sink("meta_analysis_results.txt")
cat("Per-cohort AUC and DeLong standard errors (5 scale-compatible external cohorts)\n")
cat("Cohorts with perfect separation (AUC = 1.000) have zero DeLong variance and\n")
cat("are excluded from inverse-variance pooling; they are reported individually.\n\n")
write.table(ma_df, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nRandom-effects meta-analysis (DerSimonian-Laird, logit scale)\n\n")
write.table(pool_df, sep = "\t", quote = FALSE, row.names = FALSE)
sink()

## ---- forest plot ---------------------------------------------------------
fp <- ma_df[, c("Cohort", "AUC", "CI_lower", "CI_upper", "Variance_estimable")]
fp$n_Disease <- ma_df$n_Disease; fp$n_Normal <- ma_df$n_Normal
fp$Label <- sprintf("%s (%d vs %d)%s", fp$Cohort, fp$n_Disease, fp$n_Normal,
                    ifelse(fp$Variance_estimable, "", " *"))
pooled_row <- data.frame(
  Cohort = "Pooled", AUC = pool_est$pooled_AUC,
  CI_lower = pool_est$CI_lower, CI_upper = pool_est$CI_upper,
  Variance_estimable = TRUE, n_Disease = NA, n_Normal = NA,
  Label = sprintf("Pooled (random effects): %.3f (%.3f-%.3f), I2 = %.1f%%",
                  pool_est$pooled_AUC, pool_est$CI_lower,
                  pool_est$CI_upper, pool_est$I2_percent))
fp$Type <- "Cohort"; pooled_row$Type <- "Pooled"
fp <- rbind(fp, pooled_row)
fp$Label <- factor(fp$Label, levels = rev(fp$Label))
p <- ggplot(fp, aes(x = AUC, y = Label)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper, orientation = "y"),
                width = 0.25, color = "grey30") +
  geom_point(aes(size = Type, shape = Type), color = "#2166ac") +
  scale_size_manual(values = c(Cohort = 2.6, Pooled = 3.8)) +
  scale_shape_manual(values = c(Cohort = 16, Pooled = 18)) +
  xlim(0.4, 1.02) +
  labs(x = "AUC (95% CI)", y = NULL,
       title = "External validation AUCs and pooled estimate",
       caption = "*Perfect separation (AUC = 1.000, zero variance); excluded from pooling.") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11.5, face = "bold"),
        plot.caption = element_text(size = 9, hjust = 0, color = "grey40"))
ggsave("S17_forest_plot_external_AUC.png", p, width = 7.5, height = 4.4, dpi = 300)

write.csv(ma_df, "S6_dataset_meta_analysis.csv", row.names = FALSE)
write.csv(pool_df, "S6_dataset_meta_pooled.csv", row.names = FALSE)
cat("done\n")
