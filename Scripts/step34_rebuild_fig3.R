# step34_rebuild_fig3.R
# Rebuild main Fig 3 as a proper VECTOR composite:
#   left: external-validation ROC curves (full 0-1 axes, unclipped labels)
#   inset (inside panel, lower-right): random-effects meta-analysis forest plot
# Replaces the raster-of-raster composite produced by step31_fig3_inset.R.
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

ws  <- "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
prd <- read.csv(file.path(ws, "Data", "HF_signature_external_validation_predictions.csv"),
                stringsAsFactors = FALSE)

cohorts <- c("GSE116250", "GSE26887", "GSE42955", "GSE52601", "GSE55296")
roc_pts <- list(); aucs <- numeric(length(cohorts))
for (i in seq_along(cohorts)) {
  co <- cohorts[i]
  d  <- prd[prd$Dataset == co, ]
  y  <- d$TrueLabel == "Disease"
  s  <- d$Prob_Disease
  stopifnot(any(y), any(!y))
  o  <- order(s, decreasing = TRUE)
  y <- y[o]; s <- s[o]
  P <- sum(y); N <- sum(!y)
  # full staircase (per observation; ties kept as steps, standard ROCR display)
  tp <- cumsum(y); fp <- cumsum(!y)
  tpr <- c(0, tp / P); fpr <- c(0, fp / N)
  # Mann-Whitney AUC with tie-averaged ranks (matches pROC)
  rk <- rank(s, ties.method = "average")
  auc <- (sum(rk[y]) - P * (P + 1) / 2) / (P * N)
  aucs[i] <- auc
  roc_pts[[i]] <- data.frame(FPR = fpr, TPR = tpr, Cohort = co)
}
roc_df <- do.call(rbind, roc_pts)
names(aucs) <- cohorts
print(round(aucs, 3))
stopifnot(abs(aucs["GSE116250"] - 0.857) < 0.01,
          abs(aucs["GSE55296"]  - 0.818) < 0.01,
          abs(aucs["GSE52601"]  - 0.812) < 0.01,
          aucs["GSE26887"] > 0.999, aucs["GSE42955"] > 0.999)

lab <- sprintf("%s (AUC=%s)", cohorts,
               ifelse(aucs > 0.9995, "1", sprintf("%.3f", aucs)))
lab_map <- setNames(lab, cohorts)
roc_df$Cohort <- factor(roc_df$Cohort, levels = cohorts, labels = lab)

roc <- ggplot(roc_df, aes(FPR, TPR, color = Cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = c(0.01, 0.01)) +
  labs(x = "1 \u2013 Specificity", y = "Sensitivity", color = NULL,
       title = "HF signature RF model \u2013 external validation ROC (frozen-scale cohorts)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8.5),
        legend.key.height = unit(0.45, "cm"),
        plot.title = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_blank())

# --- forest inset (same numbers as step31) ---
ma   <- read.csv(file.path(ws, "loco_meta_analysis", "S6_dataset_meta_analysis.csv"), stringsAsFactors = FALSE)
pool <- read.csv(file.path(ws, "loco_meta_analysis", "S6_dataset_meta_pooled.csv"), stringsAsFactors = FALSE)
pool3 <- pool[grepl("3 cohorts", pool$Analysis), ]
fp <- data.frame(
  Label = sprintf("%s (%d vs %d)%s", ma$Cohort, ma$n_Disease, ma$n_Normal,
                  ifelse(ma$Variance_estimable, "", " *")),
  AUC = ma$AUC, lo = ma$CI_lower, hi = ma$CI_upper, Type = "Cohort",
  stringsAsFactors = FALSE)
fp <- rbind(fp, data.frame(
  Label = sprintf("Pooled (RE): %.3f (%.3f-%.3f)", pool3$pooled_AUC, pool3$CI_lower, pool3$CI_upper),
  AUC = pool3$pooled_AUC, lo = pool3$CI_lower, hi = pool3$CI_upper, Type = "Pooled"))
fp$Label <- factor(fp$Label, levels = rev(fp$Label))

inset <- ggplot(fp, aes(AUC, Label)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey65", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, color = "grey30", linewidth = 0.35) +
  geom_point(aes(size = Type, shape = Type), color = "#2166ac") +
  scale_size_manual(values = c(Cohort = 1.4, Pooled = 2.2)) +
  scale_shape_manual(values = c(Cohort = 16, Pooled = 18)) +
  xlim(0.4, 1.02) +
  labs(x = "AUC (95% CI)", y = NULL, title = "Random-effects meta-analysis") +
  theme_bw(base_size = 7) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7.5, face = "bold"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3))

p <- roc + inset_element(inset, left = 0.42, bottom = 0.02, right = 0.985, top = 0.42)
out <- file.path(ws, "Figures", "HF_signature_external_validation_ROC.png")
ggsave(out, p, width = 8, height = 6, dpi = 300)
cat("fig3 rebuilt ->", out, "\n")
