# step43_gse183852_ligand_replication.R
# Donor-level replication of the surrogate ligand axes in GSE183852 (18 DCM vs 27 NF).
# Loads the 12GB integrated object ONCE, extracts the ligand/receptor panel counts,
# computes per-donor pseudobulk log2FC + t tests (donor = statistical unit), BH FDR.
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

ws  <- "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
out <- file.path(ws, "revision_round2")
lg  <- file.path(out, "step43_ligand_replication.log")
say <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = lg, append = TRUE) }
say("=== start %s ===", format(Sys.time()))

panel <- c("APOE", "COL1A2", "COL14A1", "COL4A1", "COL4A2", "COL4A4",
           "COL6A1", "COL6A2", "LAMA2", "LAMA4", "LAMB1", "HAS2", "NBL1",
           "LDLR", "LRP1", "CD44", "ITGB1", "ITGA11")

# ---- load big object, extract panel counts (RNA assay) ----
e <- new.env()
say("loading Robj ...")
t0 <- Sys.time()
load(file.path(ws, "data_external", "GSE183852_DCM_Integrated.Robj"), envir = e)
say("loaded in %.1f min; objects: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")),
    paste(ls(e), collapse = ", "))

find_seurat <- function(env) {
  for (nm in ls(env)) {
    x <- env[[nm]]
    if (is(x, "Seurat")) return(list(obj = x, tag = nm))
    if (is.list(x) && !isS4(x))
      for (sub in names(x)) if (is(x[[sub]], "Seurat")) return(list(obj = x[[sub]], tag = paste0(nm, "_", sub)))
  }
  stop("no Seurat object found")
}
sr <- find_seurat(e)
say("using object %s", sr$tag)
obj <- sr$obj

mat <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"),
                error = function(er) GetAssayData(obj, assay = "RNA", slot = "counts"))
say("RNA counts %d x %d", nrow(mat), ncol(mat))
md <- obj@meta.data
say("meta cols: %s", paste(colnames(md), collapse = " | "))

g2 <- intersect(panel, rownames(mat))
say("panel genes present: %d/%d -> %s", length(g2), length(panel), paste(g2, collapse = ","))
slim <- mat[g2, , drop = FALSE]
rm(mat); gc(verbose = FALSE)

# ---- harmonise cell types & donors (same scheme as step33c) ----
celltype_map <- c(
  "Fibroblasts" = "Fibroblast", "Endothelium" = "Endothelial", "Endocardium" = "Endothelial",
  "Lymphatic" = "Endothelial", "Myeloid" = "Macrophage", "Cardiomyocytes" = "CM",
  "Pericytes" = "Pericyte_SMC", "Smooth_Muscle" = "Pericyte_SMC",
  "NK/T-Cells" = "T_cell", "B-Cells" = "B_cell", "Mast" = "Mast")
md$ct_std <- unname(celltype_map[md$Names])
md$group  <- ifelse(md$condition == "Donor", "N", ifelse(md$condition == "DCM", "DCM", NA))
md$donor  <- md$orig.ident
keep <- !is.na(md$ct_std) & !is.na(md$group)
slim <- slim[, keep, drop = FALSE]; md <- md[keep, ]
say("cells kept: %d", nrow(md))

min_donors <- 3; min_cells <- 20
run_ct <- function(ct) {
  m <- md[md$ct_std == ct, ]
  if (length(unique(m$donor[m$group == "N"])) < min_donors ||
      length(unique(m$donor[m$group == "DCM"])) < min_donors) return(NULL)
  sub <- slim[, rownames(m), drop = FALSE]
  key <- paste(m$donor, m$group, sep = "||")
  pb <- sapply(unique(key), function(k) {
    cells <- rownames(m)[key == k]
    if (length(cells) < min_cells) return(rep(NA_real_, nrow(sub)))
    Matrix::rowSums(sub[, cells, drop = FALSE])
  })
  rownames(pb) <- rownames(sub)
  grp <- sub(".*\\|\\|", "", colnames(pb))
  ok <- !is.na(pb[1, ]); pb <- pb[, ok, drop = FALSE]; grp <- grp[ok]
  if (sum(grp == "DCM") < min_donors || sum(grp == "N") < min_donors) return(NULL)
  lib <- colSums(pb); lcpm <- log2(sweep(pb + 0.5, 2, lib / 1e6, `/`))
  do.call(rbind, lapply(rownames(pb), function(g) {
    x <- lcpm[g, grp == "DCM"]; y <- lcpm[g, grp == "N"]
    tt <- tryCatch(t.test(x, y), error = function(er) NULL)
    data.frame(gene = g, cell_type = ct,
               ext_log2FC = mean(x) - mean(y),
               ext_p_donor = if (is.null(tt)) NA_real_ else tt$p.value,
               n_dcm = length(x), n_n = length(y))
  }))
}
res <- do.call(rbind, Filter(Negate(is.null), lapply(sort(unique(md$ct_std)), run_ct)))
res$ext_FDR_donor <- p.adjust(res$ext_p_donor, method = "BH")
write.csv(res, file.path(out, "gse183852_ligand_replication.csv"), row.names = FALSE)

key_pairs <- list(
  c("APOE", "Macrophage"), c("COL1A2", "Fibroblast"), c("COL14A1", "Fibroblast"),
  c("COL1A2", "Pericyte_SMC"), c("COL6A2", "Fibroblast"), c("COL4A2", "Fibroblast"),
  c("LAMB1", "Fibroblast"), c("CD44", "Fibroblast"), c("ITGB1", "Fibroblast"),
  c("ITGA11", "Fibroblast"), c("LDLR", "Fibroblast"), c("LRP1", "Fibroblast"))
sink(file.path(out, "gse183852_ligand_replication_summary.txt"))
cat("GSE183852 donor-level replication of ligand axes (18 DCM vs 27 NF)\n\n")
for (kp in key_pairs) {
  r <- res[res$gene == kp[1] & res$cell_type == kp[2], ]
  if (nrow(r) == 0) { cat(sprintf("%-8s %-14s : not evaluable\n", kp[1], kp[2])); next }
  cat(sprintf("%-8s %-14s log2FC=%+.2f  P=%.2e  FDR=%.3f  (n=%d vs %d)\n",
              r$gene, r$cell_type, r$ext_log2FC, r$ext_p_donor, r$ext_FDR_donor, r$n_dcm, r$n_n))
}
cat("\nAll rows (FDR<0.1):\n")
print(res[res$ext_FDR_donor < 0.1, ])
sink()
cat(readLines(file.path(out, "gse183852_ligand_replication_summary.txt")), sep = "\n")
say("=== done %s ===", format(Sys.time()))
