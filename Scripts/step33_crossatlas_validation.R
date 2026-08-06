# step33_crossatlas_validation.R
# Cross-atlas validation of signature-gene regulatory directions.
#
# PURPOSE
#   Compare per-cell-type disease-vs-normal log2FC directions of the 39-gene
#   HF signature (reference: GSE145154, revision_round2/our_scRNA_reference_directions.csv)
#   against an INDEPENDENT external single-cell / single-nucleus heart atlas
#   (e.g. GSE183852 Koenig/Lavine 2022: 27 donor + 18 DCM;
#         Litviňuková 2020 Heart Cell Atlas: 14 healthy donors -> expanded normal pool).
#
#   Donor-level pseudobulk aggregation is used throughout, so that (for the first
#   time in this project) statistical testing is performed at the biological
#   replicate (donor) level rather than at the (pseudoreplicated) cell level.
#
# INPUT (choose ONE, edit the CONFIG block below):
#   A) ext_rds : path to a saved Seurat object (.rds) or RData/.Robj containing
#                exactly one Seurat object. Must have metadata columns for
#                donor id, condition and cell type (names configurable below).
#   B) ext_counts_csv + ext_meta_csv : a genes x cells count matrix CSV
#                (first column = gene symbol) and a cell metadata CSV with a
#                barcode column matching the count-matrix column names.
#
# OUTPUT (written to out_dir):
#   crossatlas_pseudobulk_log2FC.csv     donor-pseudobulk log2FC per gene x cell type x comparison
#   crossatlas_direction_concordance.csv merged with our reference; sign agreement flags
#   crossatlas_concordance_summary.txt   binomial sign-test + headline gene summary
#   crossatlas_concordance_scatter.png   external vs our log2FC scatter (gene x cell type points)
#
# NOTES
#   - Litviňuková 2020 contains ONLY healthy donors: with it, comparisons must be
#     run as "our DCM/ICM cells vs external normal cells" (cross-study mode,
#     ext_normal_only = TRUE), or used purely to re-estimate the normal baseline.
#   - GSE183852 contains donor + DCM (no ICM): comparisons DCM_vs_N only.

.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2)
})

## ----------------------------- CONFIG --------------------------------------
ws      <- "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
out_dir <- file.path(ws, "revision_round2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- input mode A: Seurat object ---
ext_rds <- file.path(ws, "data_external", "GSE183852_DCM_Integrated.rds")  # <- edit
# --- input mode B: counts + metadata CSV (used only if ext_rds does not exist) ---
ext_counts_csv <- file.path(ws, "data_external", "external_counts.csv")
ext_meta_csv   <- file.path(ws, "data_external", "external_metadata.csv")
ext_barcode_col <- "barcode"   # column in metadata matching count-matrix column names

# --- metadata column names in the EXTERNAL object (edit to match the atlas) ---
col_donor     <- "donor"       # e.g. "donor", "individual", "patient", "sample"
col_condition <- "condition"   # e.g. "condition", "disease", "group"
col_celltype  <- "cell_type"   # e.g. "cell_type", "cell_type_original", "annotation"

# values meaning "normal/control" and disease groups present:
val_normal    <- c("donor", "normal", "N", "NF", "control", "healthy")
val_disease   <- list(DCM = c("DCM"), ICM = c("ICM"))   # drop entries not present

# set TRUE for a normal-only atlas (Litviňuková): then disease cells come from
# OUR GSE145154 object and normals from the external atlas (cross-study mode).
ext_normal_only <- FALSE
our_rds <- "D:/早早生信分析学习课程/新R文件/results_scRNA/step1_QC_integration/heart_integrated.rds"

min_donors_per_group <- 3     # require >= this many donors per group to report a test
min_cells_per_donor  <- 20    # donor x cell type pseudobulk needs >= this many cells

# cell-type harmonisation: map external labels -> our 10 categories.
# edit after inspecting sort(unique(meta[[col_celltype]]))
celltype_map <- c(
  "Fibroblast" = "Fibroblast", "Fibroblasts" = "Fibroblast",
  "Endothelial" = "Endothelial", "Endothelial cell" = "Endothelial", "EC" = "Endothelial",
  "Macrophage" = "Macrophage", "Myeloid" = "Macrophage",
  "Cardiomyocyte" = "CM", "Cardiomyocytes" = "CM", "CM" = "CM",
  "Pericyte" = "Pericyte_SMC", "SMC" = "Pericyte_SMC", "Pericyte_SMC" = "Pericyte_SMC",
  "T cell" = "T_cell", "T_cell" = "T_cell", "NK cell" = "NK", "NK" = "NK",
  "B cell" = "B_cell", "B_cell" = "B_cell", "Mast cell" = "Mast", "Mast" = "Mast",
  "Neutrophil" = "Neutrophil"
)
## --------------------------- END CONFIG ------------------------------------

ref <- read.csv(file.path(out_dir, "our_scRNA_reference_directions.csv"),
                stringsAsFactors = FALSE)
genes <- unique(ref$gene)

## ---- load external data ----
load_external <- function() {
  if (file.exists(ext_rds)) {
    cat("Loading Seurat object:", ext_rds, "\n")
    if (grepl("\\.rds$", ext_rds, ignore.case = TRUE)) {
      obj <- readRDS(ext_rds)
    } else {
      e <- new.env(); load(ext_rds, envir = e)
      objs <- mget(ls(e), envir = e)
      obj <- objs[[which(vapply(objs, function(x) is(x, "Seurat"), logical(1)))[1]]]
    }
    counts <- GetAssayData(obj, layer = "counts")
    meta <- obj@meta.data
    return(list(counts = counts, meta = meta))
  }
  if (file.exists(ext_counts_csv) && file.exists(ext_meta_csv)) {
    cat("Loading counts CSV + metadata CSV\n")
    cm <- read.csv(ext_counts_csv, row.names = 1, check.names = FALSE)
    counts <- as(as.matrix(cm), "dgCMatrix")
    meta <- read.csv(ext_meta_csv, stringsAsFactors = FALSE)
    rownames(meta) <- meta[[ext_barcode_col]]
    meta <- meta[colnames(counts), , drop = FALSE]
    return(list(counts = counts, meta = meta))
  }
  stop("No external input found. Put the downloaded atlas in data_external/ and edit CONFIG.")
}

ext <- load_external()
meta <- ext$meta
stopifnot(all(c(col_donor, col_condition, col_celltype) %in% colnames(meta)))

meta$condition_std <- as.character(meta[[col_condition]])
meta$donor <- as.character(meta[[col_donor]])
meta$group <- NA_character_
meta$group[meta$condition_std %in% val_normal] <- "N"
for (dz in names(val_disease)) meta$group[meta$condition_std %in% val_disease[[dz]]] <- dz
meta$ct_std <- unname(celltype_map[as.character(meta[[col_celltype]])])
keep <- !is.na(meta$group) & !is.na(meta$ct_std)
cat("Cells retained after mapping:", sum(keep), "/", nrow(meta), "\n")
cat("Groups:", paste(names(table(meta$group[keep])), table(meta$group[keep]), collapse = "; "), "\n")
cat("Mapped cell types:", paste(unique(na.omit(meta$ct_std)), collapse = ", "), "\n")

counts <- ext$counts[, rownames(meta)[keep], drop = FALSE]
meta   <- meta[keep, , drop = FALSE]

## ---- cross-study mode: disease cells from our atlas, normals from external ----
if (ext_normal_only) {
  cat("Cross-study mode: external atlas supplies NORMAL cells only; disease cells from GSE145154.\n")
  our <- readRDS(our_rds)
  om <- our@meta.data
  om$donor <- om$patient; om$ct_std <- om$cell_type
  om$group <- om$condition
  oc <- GetAssayData(our, layer = "counts")
  common <- intersect(rownames(oc), rownames(counts))
  counts <- cbind(oc[common, rownames(om[om$group != "N", ]), drop = FALSE],
                  counts[common, rownames(meta[meta$group == "N", ]), drop = FALSE])
  meta <- rbind(om[om$group != "N", c("donor", "group", "ct_std")],
                meta[meta$group == "N", c("donor", "group", "ct_std")])
  meta$donor <- paste0(ifelse(meta$group == "N", "EXT_", "OUR_"), meta$donor)
} else {
  meta <- meta[, c("donor", "group", "ct_std")]
}

## ---- donor-level pseudobulk log2FC ----
genes_present <- intersect(genes, rownames(counts))
cat("Signature genes present in external matrix:", length(genes_present), "/", length(genes), "\n")

pb_log2fc <- function(ct, cmp) {
  g1 <- strsplit(cmp, "_vs_")[[1]][1]   # disease
  m  <- meta[meta$ct_std == ct & meta$group %in% c("N", g1), ]
  if (length(unique(m$donor[m$group == "N"]))  < min_donors_per_group) return(NULL)
  if (length(unique(m$donor[m$group == g1])) < min_donors_per_group) return(NULL)
  sub <- counts[genes_present, rownames(m), drop = FALSE]
  key <- paste(m$donor, m$group, sep = "||")
  donors <- unique(key)
  pb <- sapply(donors, function(k) {
    cells <- rownames(m)[key == k]
    if (length(cells) < min_cells_per_donor) return(rep(NA_real_, length(genes_present)))
    Matrix::rowSums(sub[, cells, drop = FALSE])
  })
  rownames(pb) <- genes_present
  grp <- sub(".*\\|\\|", "", colnames(pb))
  ok  <- !is.na(pb[1, ])
  pb  <- pb[, ok, drop = FALSE]; grp <- grp[ok]
  cpm <- sweep(pb + 0.5, 2, colSums(pb) / 1e6, `/`)
  lcpm <- log2(cpm)
  res <- lapply(genes_present, function(g) {
    x <- lcpm[g, grp == g1]; y <- lcpm[g, grp == "N"]
    if (length(x) < min_donors_per_group || length(y) < min_donors_per_group) return(NULL)
    tt <- tryCatch(t.test(x, y), error = function(e) NULL)
    data.frame(gene = g, cell_type = ct, comparison = cmp,
               ext_log2FC = mean(x) - mean(y),
               ext_p_donor = if (is.null(tt)) NA_real_ else tt$p.value,
               n_donor_disease = length(x), n_donor_normal = length(y),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), res))
}

comparisons <- if (ext_normal_only || "ICM" %in% meta$group) c("DCM_vs_N", "ICM_vs_N") else "DCM_vs_N"
if (!ext_normal_only) comparisons <- comparisons[sapply(comparisons, function(cmp)
  strsplit(cmp, "_vs_")[[1]][1] %in% meta$group)]

out <- list(); ii <- 0
for (ct in sort(unique(meta$ct_std))) for (cmp in comparisons) {
  r <- pb_log2fc(ct, cmp)
  if (!is.null(r)) { ii <- ii + 1; out[[ii]] <- r }
}
ext_res <- do.call(rbind, out)
ext_res$ext_FDR_donor <- p.adjust(ext_res$ext_p_donor, method = "BH")
write.csv(ext_res, file.path(out_dir, "crossatlas_pseudobulk_log2FC.csv"), row.names = FALSE)

## ---- concordance with our reference ----
m <- merge(ref, ext_res, by = c("gene", "cell_type", "comparison"))
m$sign_agree <- sign(m$our_log2FC) == sign(m$ext_log2FC)
# primary concordance set: our FDR-significant entries
m$concordance_class <- "our_not_significant"
m$concordance_class[m$our_FDR < 0.05 & m$sign_agree]  <- "concordant"
m$concordance_class[m$our_FDR < 0.05 & !m$sign_agree] <- "discordant"
write.csv(m, file.path(out_dir, "crossatlas_direction_concordance.csv"), row.names = FALSE)

sig <- m[m$our_FDR < 0.05, ]
bt <- binom.test(sum(sig$sign_agree), nrow(sig), p = 0.5, alternative = "greater")
headline <- c("FCN3", "ASPN", "PLTP", "HMGN2", "RARRES1", "EDNRB", "ISLR", "BCL6", "LAPTM5")
sink(file.path(out_dir, "crossatlas_concordance_summary.txt"))
cat("Cross-atlas direction concordance (donor-level pseudobulk)\n")
cat("Pairs tested (our FDR<0.05):", nrow(sig), "\n")
cat("Sign-concordant:", sum(sig$sign_agree),
    sprintf("(%.1f%%)", 100 * mean(sig$sign_agree)), "\n")
cat("Binomial sign test vs 50%: p =", signif(bt$p.value, 3), "\n\n")
cat("Headline genes:\n")
print(subset(m, gene %in% headline,
             select = c(gene, cell_type, comparison, our_log2FC, our_FDR,
                        ext_log2FC, ext_FDR_donor, concordance_class)))
sink()

p <- ggplot(sig, aes(our_log2FC, ext_log2FC, color = sign_agree)) +
  geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) +
  geom_point(alpha = 0.7) + geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  scale_color_manual(values = c(`TRUE` = "#2166ac", `FALSE` = "#b2182b"), name = "sign agreement") +
  labs(x = "GSE145154 log2FC (cell-level DE)", y = "External atlas log2FC (donor pseudobulk)",
       title = "Cross-atlas direction concordance") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "crossatlas_concordance_scatter.png"), p, width = 6.5, height = 5.5, dpi = 300)

cat("Done. Outputs in", out_dir, "\n")
