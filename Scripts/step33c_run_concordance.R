# step33c_run_concordance.R
# Cross-atlas concordance using the slim GSE183852 extracts (no giant object needed).
# Donor-level pseudobulk log2FC (18 DCM vs 27 non-failing donors) for the 39
# signature genes x harmonised cell types, then sign-concordance vs GSE145154.
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(Matrix); library(ggplot2) })

ws  <- "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
out <- file.path(ws, "revision_round2")

counts <- readRDS(file.path(out, "gse183852_counts39_RefMerge_RNA.rds"))
meta   <- read.csv(file.path(out, "gse183852_metadata_RefMerge.csv"),
                   stringsAsFactors = FALSE, row.names = 1)
stopifnot(all(rownames(meta) == colnames(counts)))

# harmonise: our 10-category scheme
celltype_map <- c(
  "Fibroblasts"    = "Fibroblast",
  "Endothelium"    = "Endothelial",
  "Endocardium"    = "Endothelial",
  "Lymphatic"      = "Endothelial",
  "Myeloid"        = "Macrophage",
  "Cardiomyocytes" = "CM",
  "Pericytes"      = "Pericyte_SMC",
  "Smooth_Muscle"  = "Pericyte_SMC",
  "NK/T-Cells"     = "T_cell",
  "B-Cells"        = "B_cell",
  "Mast"           = "Mast"
)
meta$ct_std <- unname(celltype_map[meta$Names])
meta$group  <- ifelse(meta$condition == "Donor", "N",
                      ifelse(meta$condition == "DCM", "DCM", NA))
meta$donor  <- meta$orig.ident
keep <- !is.na(meta$ct_std) & !is.na(meta$group)
counts <- counts[, keep, drop = FALSE]; meta <- meta[keep, ]
cat("cells kept:", nrow(meta), "\n")
print(table(meta$ct_std, meta$group))

min_donors <- 3; min_cells <- 20
genes <- rownames(counts)

pb <- function(ct) {
  m <- meta[meta$ct_std == ct, ]
  if (length(unique(m$donor[m$group == "N"]))   < min_donors) return(NULL)
  if (length(unique(m$donor[m$group == "DCM"])) < min_donors) return(NULL)
  sub <- counts[, rownames(m), drop = FALSE]
  key <- paste(m$donor, m$group, sep = "||")
  pb <- sapply(unique(key), function(k) {
    cells <- rownames(m)[key == k]
    if (length(cells) < min_cells) return(rep(NA_real_, nrow(sub)))
    Matrix::rowSums(sub[, cells, drop = FALSE])
  })
  rownames(pb) <- genes
  grp <- sub(".*\\|\\|", "", colnames(pb))
  ok  <- !is.na(pb[1, ])
  pb  <- pb[, ok, drop = FALSE]; grp <- grp[ok]
  if (sum(grp == "DCM") < min_donors || sum(grp == "N") < min_donors) return(NULL)
  lib <- colSums(pb); lcpm <- log2(sweep(pb + 0.5, 2, lib / 1e6, `/`))
  res <- lapply(genes, function(g) {
    x <- lcpm[g, grp == "DCM"]; y <- lcpm[g, grp == "N"]
    tt <- tryCatch(t.test(x, y), error = function(e) NULL)
    data.frame(gene = g, cell_type = ct, comparison = "DCM_vs_N",
               ext_log2FC = mean(x) - mean(y),
               ext_p_donor = if (is.null(tt)) NA_real_ else tt$p.value,
               n_donor_disease = length(x), n_donor_normal = length(y))
  })
  do.call(rbind, res)
}

ext <- do.call(rbind, Filter(Negate(is.null),
                             lapply(sort(unique(meta$ct_std)), pb)))
ext$ext_FDR_donor <- p.adjust(ext$ext_p_donor, method = "BH")
write.csv(ext, file.path(out, "crossatlas_pseudobulk_log2FC.csv"), row.names = FALSE)

# concordance vs our reference (GSE183852 has no ICM -> DCM_vs_N only)
ref <- read.csv(file.path(out, "our_scRNA_reference_directions.csv"))
ref <- ref[ref$comparison == "DCM_vs_N", ]
m <- merge(ref, ext, by = c("gene", "cell_type", "comparison"))
m$sign_agree <- sign(m$our_log2FC) == sign(m$ext_log2FC)
m$concordance_class <- "our_not_significant"
m$concordance_class[m$our_FDR < 0.05 &  m$sign_agree] <- "concordant"
m$concordance_class[m$our_FDR < 0.05 & !m$sign_agree] <- "discordant"
write.csv(m, file.path(out, "crossatlas_direction_concordance.csv"), row.names = FALSE)

sig <- m[m$our_FDR < 0.05 & !is.na(m$sign_agree), ]
bt <- binom.test(sum(sig$sign_agree), nrow(sig), p = 0.5, alternative = "greater")
headline <- c("FCN3", "ASPN", "PLTP", "HMGN2", "RARRES1", "EDNRB", "ISLR", "BCL6", "LAPTM5")
sink(file.path(out, "crossatlas_concordance_summary.txt"))
cat("Cross-atlas direction concordance: GSE183852 (18 DCM vs 27 donor, donor-level pseudobulk)\n")
cat("Gene x cell-type pairs with our FDR<0.05:", nrow(sig), "\n")
cat("Sign-concordant:", sum(sig$sign_agree), sprintf("(%.1f%%)\n", 100 * mean(sig$sign_agree)))
cat("Binomial sign test vs 50%: p =", signif(bt$p.value, 3), "\n")
cat("Discordant pairs:\n")
print(sig[!sig$sign_agree, c("gene", "cell_type", "our_log2FC", "our_FDR", "ext_log2FC", "ext_FDR_donor")])
cat("\nHeadline genes:\n")
print(m[m$gene %in% headline, c("gene", "cell_type", "our_log2FC", "our_FDR",
                                "ext_log2FC", "ext_FDR_donor", "concordance_class")])
sink()
cat(readLines(file.path(out, "crossatlas_concordance_summary.txt")), sep = "\n")

p <- ggplot(sig, aes(our_log2FC, ext_log2FC, color = sign_agree)) +
  geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  scale_color_manual(values = c(`TRUE` = "#2166ac", `FALSE` = "#b2182b"),
                     name = "sign agreement") +
  labs(x = "GSE145154 log2FC (cell-level DE)",
       y = "GSE183852 log2FC (donor-level pseudobulk)",
       title = "Cross-atlas direction concordance (DCM vs non-failing)") +
  theme_bw(base_size = 12)
ggsave(file.path(out, "crossatlas_concordance_scatter.png"), p,
       width = 6.5, height = 5.5, dpi = 300)
cat("outputs written\n")
