# step36_lvef_regression.R
# LVEF regression: does the frozen 39-gene RF model score track per-sample ejection
# fraction in GSE52601 (the only cohort with deposited per-sample EF)?
# Scoring pipeline replicates step21b exactly (frozen train z-score parameters,
# RF retrained with seed 12345 on bundle$train_set); expression annotation
# replicates step13b (illuminaHumanv4.db, median collapse per symbol).
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(randomForest); library(illuminaHumanv4.db); library(AnnotationDbi)
  library(data.table); library(ggplot2)
})

ws  <- "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
dd  <- file.path(ws, "data_external")
out <- file.path(ws, "revision_round2")

# ---- 1. model (exact step21b reproduction) ----
bundle <- readRDS(file.path(dd, "best_model_step20b.rds"))
selected <- bundle$features
train_center <- bundle$combat_params$g_mean
train_scale  <- bundle$combat_params$g_sd
set.seed(12345)
model_rf <- randomForest(Group ~ ., data = bundle$train_set, ntree = 500)

# ---- 2. parse GSE52601 series matrix ----
lines <- readLines(gzfile(file.path(dd, "GSE52601_series_matrix.txt.gz")))
get_meta <- function(prefix) {
  idx <- grep(paste0("^!", prefix, "\t"), lines)
  if (length(idx) == 0) return(NULL)
  gsub('"', '', strsplit(lines[idx], "\t")[[1]][-1])
}
geo_ids <- get_meta("Sample_geo_accession")
titles  <- get_meta("Sample_title")

char_idx <- grep("^!Sample_characteristics_ch1", lines)
pheno <- list()
for (ci in char_idx) {
  vals <- gsub('"', '', strsplit(lines[ci], "\t")[[1]][-1])
  key <- sub(":.*$", "", vals[1])
  pheno[[key]] <- sub("^[^:]*:\\s*", "", vals)
}
meta <- data.frame(gsm = geo_ids, title = titles, stringsAsFactors = FALSE)
for (k in names(pheno)) if (length(pheno[[k]]) == nrow(meta)) meta[[gsub(" ", "_", k)]] <- pheno[[k]]
meta$EF <- suppressWarnings(as.numeric(meta$ejection_fraction))
meta$age_num <- suppressWarnings(as.numeric(meta$age))
meta$is_fetal <- grepl("fetal", meta$title)
meta$arm <- ifelse(meta$is_fetal, "fetal",
            ifelse(meta$disease == "None", "control",
            ifelse(meta$treatment == "LVAD Implantation", "pre-LVAD",
            ifelse(meta$treatment == "LVAD Explantation", "post-LVAD", "other"))))
cat("sample arms:\n"); print(table(meta$arm))

# ---- 3. expression: probe -> symbol (median collapse, as step13b) ----
tbl_start <- grep("!series_matrix_table_begin", lines) + 1
tbl_end   <- grep("!series_matrix_table_end", lines) - 1
tmp <- tempfile(); writeLines(lines[tbl_start:tbl_end], tmp)
expr <- fread(tmp, check.names = FALSE, data.table = FALSE); unlink(tmp)
rownames(expr) <- expr[[1]]; expr <- expr[, -1, drop = FALSE]
expr <- expr[, geo_ids, drop = FALSE]

mapped <- select(illuminaHumanv4.db, keys = rownames(expr), columns = "SYMBOL", keytype = "PROBEID")
mapped <- mapped[!duplicated(mapped$PROBEID), ]
expr$PROBEID <- rownames(expr)
ea <- merge(expr, mapped, by = "PROBEID", all.x = TRUE)
ea <- ea[!is.na(ea$SYMBOL) & ea$SYMBOL != "", ]
gene_expr <- sapply(split(ea[, geo_ids, drop = FALSE], ea$SYMBOL),
                    function(x) apply(x, 2, median, na.rm = TRUE))
gene_expr <- t(gene_expr)  # genes x samples
cat("genes annotated:", nrow(gene_expr), "\n")

# ---- 4. score all postnatal samples with the frozen pipeline ----
common <- intersect(selected, rownames(gene_expr))
cat("model genes found:", length(common), "/", length(selected), "\n")
X <- matrix(0, nrow = ncol(gene_expr), ncol = length(selected),
            dimnames = list(colnames(gene_expr), selected))
X[, common] <- t(gene_expr[common, , drop = FALSE])
Xs <- scale(X, center = train_center[selected], scale = train_scale[selected])
Xs[is.nan(Xs)] <- 0
prob <- as.numeric(predict(model_rf, newdata = as.data.frame(Xs), type = "prob")[, "Disease"])
names(prob) <- colnames(gene_expr)

meta$prob <- prob[meta$gsm]
post <- meta[!meta$is_fetal, ]
write.csv(post, file.path(out, "gse52601_lvef_scores.csv"), row.names = FALSE)

# sanity: reproduce the 8 samples used in the paper's external validation
old <- read.csv(file.path(ws, "Data", "HF_signature_external_validation_predictions.csv"))
old <- old[old$Dataset == "GSE52601", ]
chk <- merge(old, post[, c("gsm", "prob")], by.x = "Sample", by.y = "gsm")
cat("\nsanity check vs published scores (", nrow(chk), "samples ):\n", sep = "")
print(chk)

# ---- 5. EF regression ----
ef <- post[!is.na(post$EF), ]
cat("\n=== samples with EF:", nrow(ef), "===\n")
sp_all   <- suppressWarnings(cor.test(ef$prob, ef$EF, method = "spearman"))
failing  <- ef[ef$arm %in% c("pre-LVAD", "post-LVAD"), ]
sp_fail  <- suppressWarnings(cor.test(failing$prob, failing$EF, method = "spearman"))
pre      <- ef[ef$arm == "pre-LVAD", ]
sp_pre   <- suppressWarnings(cor.test(pre$prob, pre$EF, method = "spearman"))
cat(sprintf("Spearman all (n=%d): rho=%.3f P=%.4g\n", nrow(ef), sp_all$estimate, sp_all$p.value))
cat(sprintf("Spearman failing only (n=%d): rho=%.3f P=%.4g\n", nrow(failing), sp_fail$estimate, sp_fail$p.value))
cat(sprintf("Spearman pre-LVAD only (n=%d): rho=%.3f P=%.4g\n", nrow(pre), sp_pre$estimate, sp_pre$p.value))

# linear model with age/sex covariates (small n, descriptive)
ef$sex_m <- ifelse(tolower(ef$gender) == "male", 1, 0)
fit <- lm(EF ~ prob + age_num + sex_m, data = ef)
cat("\n=== linear model EF ~ prob + age + sex (all EF samples) ===\n")
print(summary(fit)$coefficients)

# ---- 6. paired LVAD analysis (same donor pre/post by sample id) ----
pairs_ids <- intersect(ef$sample_id[ef$arm == "pre-LVAD"], ef$sample_id[ef$arm == "post-LVAD"])
cat("\npaired LVAD donors:", paste(pairs_ids, collapse = ", "), "\n")
pre_p  <- ef[ef$arm == "pre-LVAD"  & ef$sample_id %in% pairs_ids, ]
post_p <- ef[ef$arm == "post-LVAD" & ef$sample_id %in% pairs_ids, ]
pre_p  <- pre_p[match(pairs_ids, pre_p$sample_id), ]
post_p <- post_p[match(pairs_ids, post_p$sample_id), ]
if (nrow(pre_p) >= 3) {
  wt <- wilcox.test(pre_p$prob, post_p$prob, paired = TRUE)
  cat(sprintf("paired prob pre vs post LVAD (n=%d): median %.3f vs %.3f, P=%.4g\n",
              nrow(pre_p), median(pre_p$prob), median(post_p$prob), wt$p.value))
  wtef <- wilcox.test(pre_p$EF, post_p$EF, paired = TRUE)
  cat(sprintf("paired EF pre vs post LVAD: median %.1f vs %.1f, P=%.4g\n",
              median(pre_p$EF), median(post_p$EF), wtef$p.value))
}

# ---- 7. figures ----
ef$arm <- factor(ef$arm, levels = c("control", "pre-LVAD", "post-LVAD"))
p1 <- ggplot(ef, aes(prob, EF, color = arm)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linetype = 2, aes(group = 1)) +
  scale_color_manual(values = c(control = "#7fbc41", `pre-LVAD` = "#b2182b", `post-LVAD` = "#f4a582")) +
  labs(x = "39-gene RF model score (P(Disease))", y = "Ejection fraction (%)",
       title = sprintf("Model score vs per-sample EF in GSE52601 (Spearman rho = %.2f, P = %.3g; n = %d)",
                       sp_all$estimate, sp_all$p.value, nrow(ef))) +
  theme_bw(base_size = 11) + theme(legend.title = element_blank())

pair_df <- rbind(
  data.frame(sample_id = pairs_ids, arm = "pre-LVAD",  prob = pre_p$prob,  EF = pre_p$EF),
  data.frame(sample_id = pairs_ids, arm = "post-LVAD", prob = post_p$prob, EF = post_p$EF))
p2 <- ggplot(pair_df, aes(arm, prob, group = sample_id)) +
  geom_line(color = "grey60") + geom_point(size = 2.5, color = "#b2182b") +
  labs(x = NULL, y = "39-gene RF model score (P(Disease))",
       title = sprintf("Paired pre/post-LVAD model scores (n = %d donors)", length(pairs_ids))) +
  theme_bw(base_size = 11)

library(patchwork)
ggsave(file.path(out, "S_lvef_regression.png"), p1 | p2, width = 11, height = 4.5, dpi = 300)
cat("\noutputs written: gse52601_lvef_scores.csv, S_lvef_regression.png\n")
