# step37_lvef_gse183852.R
# Donor-level LVEF regression in the Koenig atlas (GSE183852): does the 39-gene
# signature track per-donor pre-operative ejection fraction?
# Scores are computed from whole-tissue donor pseudobulks of the 39 model genes
# (slim counts extracted by step33b); EF from Koenig et al. Supplementary Table
# (MOESM3, Nuclei Patient Data; donor IDs match orig.ident 1:1).
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(Matrix); library(randomForest); library(ggplot2); library(patchwork) })

ws  <- "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
dd  <- file.path(ws, "data_external")
out <- file.path(ws, "revision_round2")

# ---- 1. donor pseudobulk of the 39 model genes ----
counts <- readRDS(file.path(out, "gse183852_counts39_RefMerge_RNA.rds"))
md <- read.csv(file.path(out, "gse183852_metadata_RefMerge.csv"))
md <- md[match(colnames(counts), md[[1]]), ]
stopifnot(identical(as.character(md[[1]]), colnames(counts)))

donors <- unique(md$orig.ident)
pb <- sapply(donors, function(d) Matrix::rowSums(counts[, md$orig.ident == d, drop = FALSE]))
lib <- tapply(md$nCount_RNA, md$orig.ident, sum)[donors]
logcpm <- log2(sweep(pb + 0.5, 2, lib / 1e6, `/`))
cat("pseudobulk:", nrow(logcpm), "genes x", ncol(logcpm), "donors\n")

# ---- 2. scores ----
bundle <- readRDS(file.path(dd, "best_model_step20b.rds"))
selected <- bundle$features
train_center <- bundle$combat_params$g_mean
train_scale  <- bundle$combat_params$g_sd
set.seed(12345)
model_rf <- randomForest(Group ~ ., data = bundle$train_set, ntree = 500)

g <- intersect(selected, rownames(logcpm))
cat("model genes present:", length(g), "/", length(selected), "\n")
X <- matrix(0, nrow = ncol(logcpm), ncol = length(selected),
            dimnames = list(colnames(logcpm), selected))
X[, g] <- t(logcpm[g, , drop = FALSE])

# (a) per-cohort z (as used for scale-incompatible cohorts in the paper)
Xz <- X
Xz[, g] <- t(scale(t(logcpm[g, , drop = FALSE])))
Xz[is.nan(Xz)] <- 0
prob_cohort <- as.numeric(predict(model_rf, as.data.frame(Xz), type = "prob")[, "Disease"])
# (b) frozen train z (descriptive; scale mismatch expected)
Xf <- scale(X, center = train_center[selected], scale = train_scale[selected])
Xf[is.nan(Xf)] <- 0
prob_frozen <- as.numeric(predict(model_rf, as.data.frame(Xf), type = "prob")[, "Disease"])
# (c) directional signature score: mean(z_gene * disease_direction)
tr <- bundle$train_set
dir_sign <- sapply(selected, function(gg) sign(mean(tr[[gg]][tr$Group == "Disease"]) -
                                               mean(tr[[gg]][tr$Group == "Normal"])))
# Xz is donors x genes; multiply each gene column by its direction, average per donor
sig_score <- rowMeans(sweep(Xz[, selected], 2, dir_sign, `*`))

sc <- data.frame(donor = colnames(logcpm), prob_cohort, prob_frozen, sig_score)

# ---- 3. merge EF ----
ef <- read.csv(file.path(dd, "koenig_donor_ef.csv"))
m <- merge(sc, ef, by = "donor")
m$group <- ifelse(m$Condition == "DCM", "DCM", "NF")
m$sex_m <- ifelse(tolower(m$Sex) == "male", 1, 0)
m <- m[!is.na(m$EF), ]
write.csv(m, file.path(out, "gse183852_donor_ef_scores.csv"), row.names = FALSE)
cat("donors with EF:", nrow(m), "| DCM:", sum(m$group == "DCM"), "| NF:", sum(m$group == "NF"), "\n")

# ---- 4. regressions ----
run <- function(df, label) {
  cat("\n=== ", label, " (n=", nrow(df), ") ===\n", sep = "")
  for (v in c("prob_cohort", "prob_frozen", "sig_score")) {
    sp <- suppressWarnings(cor.test(df[[v]], df$EF, method = "spearman"))
    cat(sprintf("%-12s Spearman rho=%+.3f P=%.4g\n", v, sp$estimate, sp$p.value))
  }
}
dcm <- m[m$group == "DCM", ]
run(dcm, "DCM donors only")
run(m, "all donors (descriptive)")
adult <- dcm[!is.na(dcm$Age) & dcm$Age >= 18, ]
run(adult, "adult DCM (sensitivity)")

cat("\n=== linear model, DCM only: EF ~ sig_score + age + sex ===\n")
dcm2 <- dcm[!is.na(dcm$Age), ]
if (nrow(dcm2) >= 8) print(summary(lm(EF ~ sig_score + Age + sex_m, data = dcm2))$coefficients)

# ---- 5. figures ----
m$group <- factor(m$group, levels = c("NF", "DCM"))
sp_d <- suppressWarnings(cor.test(dcm$sig_score, dcm$EF, method = "spearman"))
sp_a <- suppressWarnings(cor.test(m$sig_score, m$EF, method = "spearman"))
p1 <- ggplot(dcm, aes(sig_score, EF)) +
  geom_point(size = 2.5, color = "#b2182b") +
  geom_smooth(method = "lm", se = TRUE, linetype = 2, color = "grey30") +
  labs(x = "Directional 39-gene signature score (donor pseudobulk, z)",
       y = "Pre-operative EF (%)",
       title = sprintf("DCM donors (n = %d): rho = %+.2f, P = %.3g", nrow(dcm), sp_d$estimate, sp_d$p.value)) +
  theme_bw(base_size = 11)
p2 <- ggplot(m, aes(sig_score, EF, color = group)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(NF = "#7fbc41", DCM = "#b2182b")) +
  labs(x = "Directional 39-gene signature score (donor pseudobulk, z)",
       y = "Pre-operative EF (%)",
       title = sprintf("All donors (n = %d): rho = %+.2f, P = %.3g", nrow(m), sp_a$estimate, sp_a$p.value)) +
  theme_bw(base_size = 11) + theme(legend.title = element_blank())
ggsave(file.path(out, "S_lvef_gse183852.png"), p1 | p2, width = 11, height = 4.5, dpi = 300)
cat("\noutputs: gse183852_donor_ef_scores.csv, S_lvef_gse183852.png\n")
