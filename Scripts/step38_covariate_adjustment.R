# step38_covariate_adjustment.R
# Age/sex covariate-adjusted validation of the frozen 39-gene RF model.
# Two honest evaluation sets:
#   (1) GSE57338 held-out test subset (n=92; untouched by training, split-first)
#   (2) GSE116250 external cohort (n=27; fully independent)
# Analyses: adjusted logistic regression (Group ~ score + age + sex), incremental
# AUC, subgroup AUCs (sex, age-stratified), and confounding checks (score ~ age/sex).
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(pROC); library(ggplot2); library(patchwork) })

ws  <- "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
dd  <- file.path(ws, "data_external")
out <- file.path(ws, "revision_round2")

parse_series_pheno <- function(gz, fields) {
  lines <- readLines(gzfile(gz))
  get_row <- function(prefix) gsub('"', '', strsplit(lines[grep(paste0("^!", prefix, "\t"), lines)[1]], "\t")[[1]][-1])
  gsm <- get_row("Sample_geo_accession")
  ph <- data.frame(gsm = gsm, title = get_row("Sample_title"), stringsAsFactors = FALSE)
  char_idx <- grep("^!Sample_characteristics_ch1", lines)
  for (ci in char_idx) {
    vals <- gsub('"', '', strsplit(lines[ci], "\t")[[1]][-1])
    key <- tolower(sub(":.*$", "", vals[1]))
    if (key %in% fields) ph[[key]] <- sub("^[^:]*:\\s*", "", vals)
  }
  ph
}

# ---------- GSE57338 test subset ----------
ph338 <- parse_series_pheno(file.path(dd, "GSE57338_series_matrix.txt.gz"), c("age", "gender", "disease status"))
ph338$age <- as.numeric(ph338$age)
ph338$sex_m <- ifelse(tolower(ph338$gender) == "male", 1, 0)

pred <- read.csv(file.path(ws, "Data", "HF_signature_ML_test_predictions.csv"))
t338 <- merge(pred[pred$Dataset == "GSE57338", c("Sample", "Group", "RF")],
              ph338, by.x = "Sample", by.y = "gsm")
t338$y <- ifelse(t338$Group == "Disease", 1, 0)
cat("GSE57338 test subset:", nrow(t338), "| Disease:", sum(t338$y), "| Normal:", sum(1 - t338$y), "\n")

# ---------- GSE116250 external ----------
ph116 <- parse_series_pheno(file.path(dd, "GSE116250_series_matrix.txt.gz"), c("age", "sex", "disease"))
ph116$age <- as.numeric(ph116$age)
ph116$sex_m <- ifelse(tolower(ph116$sex) == "male", 1, 0)

ext <- read.csv(file.path(ws, "Data", "HF_signature_external_validation_predictions.csv"))
e116 <- merge(ext[ext$Dataset == "GSE116250", c("Sample", "Prob_Disease", "TrueLabel")],
              ph116, by.x = "Sample", by.y = "title")
e116$y <- ifelse(e116$TrueLabel == "Disease", 1, 0)
e116$RF <- e116$Prob_Disease
cat("GSE116250 external:", nrow(e116), "| Disease:", sum(e116$y), "| Normal:", sum(1 - e116$y), "\n")

# ---------- analysis suite ----------
analyze <- function(df, label) {
  cat("\n================ ", label, " (n=", nrow(df), ") ================\n", sep = "")
  res <- data.frame(cohort = label)
  # adjusted logistic regression
  fit <- glm(y ~ RF + age + sex_m, data = df, family = binomial)
  co <- summary(fit)$coefficients
  print(co)
  res$beta_score <- co["RF", "Estimate"]; res$p_score_adj <- co["RF", "Pr(>|z|)"]
  res$beta_age <- co["age", "Estimate"];   res$p_age <- co["age", "Pr(>|z|)"]
  res$beta_sex <- co["sex_m", "Estimate"]; res$p_sex <- co["sex_m", "Pr(>|z|)"]
  # AUCs: score alone / age+sex / combined
  a1 <- roc(df$y, df$RF, quiet = TRUE)
  df$prob_cov <- predict(glm(y ~ age + sex_m, data = df, family = binomial), type = "response")
  a2 <- roc(df$y, df$prob_cov, quiet = TRUE)
  a3 <- roc(df$y, predict(fit, type = "response"), quiet = TRUE)
  cat(sprintf("AUC score alone      = %.3f (%.3f-%.3f)\n", auc(a1), ci(a1)[1], ci(a1)[3]))
  cat(sprintf("AUC age+sex alone    = %.3f (%.3f-%.3f)\n", auc(a2), ci(a2)[1], ci(a2)[3]))
  cat(sprintf("AUC score+age+sex    = %.3f (%.3f-%.3f)\n", auc(a3), ci(a3)[1], ci(a3)[3]))
  res$auc_score <- auc(a1); res$auc_cov <- auc(a2); res$auc_comb <- auc(a3)
  # subgroup AUCs
  subs <- list()
  for (sx in c("male", "female")) {
    d <- df[ifelse(df$sex_m == 1, "male", "female") == sx, ]
    if (length(unique(d$y)) == 2 && nrow(d) >= 8) {
      aa <- roc(d$y, d$RF, quiet = TRUE)
      subs[[sx]] <- data.frame(cohort = label, subgroup = sx, n = nrow(d),
                               auc = as.numeric(auc(aa)), lo = ci(aa)[1], hi = ci(aa)[3])
      cat(sprintf("subgroup %-7s n=%3d  AUC=%.3f (%.3f-%.3f)\n", sx, nrow(d), auc(aa), ci(aa)[1], ci(aa)[3]))
    }
  }
  med <- median(df$age, na.rm = TRUE)
  for (ag in c("younger", "older")) {
    d <- df[(df$age <= med) == (ag == "younger"), ]
    if (length(unique(d$y)) == 2 && nrow(d) >= 8) {
      aa <- roc(d$y, d$RF, quiet = TRUE)
      subs[[ag]] <- data.frame(cohort = label, subgroup = paste0(ag, " (age", ifelse(ag=="younger","<=",">"), round(med), ")"),
                               n = nrow(d), auc = as.numeric(auc(aa)), lo = ci(aa)[1], hi = ci(aa)[3])
      cat(sprintf("subgroup %-7s n=%3d  AUC=%.3f (%.3f-%.3f)\n", ag, nrow(d), auc(aa), ci(aa)[1], ci(aa)[3]))
    }
  }
  # confounding: score ~ age / sex within each diagnosis group
  for (gg in c("Disease", "Normal")) {
    d <- df[df$Group == gg | df$TrueLabel == gg, ]
    if (nrow(d) >= 8) {
      sp <- suppressWarnings(cor.test(d$RF, d$age, method = "spearman"))
      wt <- wilcox.test(RF ~ sex_m, data = d)
      cat(sprintf("within %-7s: score~age rho=%+.3f (P=%.3g); score~sex Wilcoxon P=%.3g\n",
                  gg, sp$estimate, sp$p.value, wt$p.value))
    }
  }
  list(res = res, subs = do.call(rbind, subs))
}

r1 <- analyze(t338, "GSE57338 test subset")
r2 <- analyze(e116, "GSE116250 external")

res_all <- rbind(r1$res, r2$res)
subs_all <- rbind(r1$subs, r2$subs)
write.csv(res_all, file.path(out, "covariate_adjustment_summary.csv"), row.names = FALSE)
write.csv(subs_all, file.path(out, "covariate_subgroup_aucs.csv"), row.names = FALSE)

# ---------- figure ----------
subs_all$label <- paste0(subs_all$subgroup, "  (n=", subs_all$n, ")")
subs_all$label <- factor(subs_all$label, levels = unique(subs_all$label[order(subs_all$auc)]))
p1 <- ggplot(subs_all, aes(auc, label, color = cohort)) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25) +
  geom_vline(xintercept = 0.5, linetype = 3, color = "grey70") +
  facet_grid(cohort ~ ., scales = "free_y", space = "free_y") +
  xlim(0.4, 1.02) +
  scale_color_manual(values = c("GSE57338 test subset" = "#2166ac", "GSE116250 external" = "#b2182b")) +
  labs(x = "AUC of the 39-gene RF score (95% CI)", y = NULL,
       title = "Discrimination is stable across sex and age subgroups") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

both <- rbind(data.frame(t338[, c("RF", "age", "Group")], cohort = "GSE57338 test subset"),
              data.frame(RF = e116$RF, age = e116$age, Group = e116$TrueLabel, cohort = "GSE116250 external"))
p2 <- ggplot(both, aes(age, RF, color = Group)) +
  geom_point(size = 1.6, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  facet_wrap(~ cohort, scales = "free") +
  scale_color_manual(values = c(Disease = "#b2182b", Normal = "#7fbc41")) +
  labs(x = "Age (years)", y = "39-gene RF model score",
       title = "Model score vs age within each diagnosis group") +
  theme_bw(base_size = 11) + theme(legend.title = element_blank(), legend.position = "bottom")

ggsave(file.path(out, "S_covariate_adjustment.png"), p1 / p2, width = 8, height = 9, dpi = 300)
cat("\noutputs: covariate_adjustment_summary.csv, covariate_subgroup_aucs.csv, S_covariate_adjustment.png\n")
