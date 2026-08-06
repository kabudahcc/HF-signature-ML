# step33a: Build our-side reference direction table for cross-atlas validation
# Extracts per-cell-type log2FC (DCM vs N, ICM vs N) for all 39 signature genes
# from the GSE145154 single-cell DE tables (step2_DEG_composition).
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))

de_dir   <- "D:/早早生信分析学习课程/新R文件/results_scRNA/step2_DEG_composition"
ws       <- "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
out_csv  <- file.path(ws, "revision_round2", "our_scRNA_reference_directions.csv")

genes <- read.delim(file.path(ws, "Data", "HF_signature_ML_model_genes.txt"),
                    stringsAsFactors = FALSE)$gene
genes <- unique(trimws(genes))
cat("Signature genes:", length(genes), "\n")

cell_types <- c("Fibroblast", "Endothelial", "Macrophage", "CM", "Pericyte_SMC",
                "T_cell", "B_cell", "NK", "Mast", "Neutrophil")
comparisons <- c("DCM_vs_N", "ICM_vs_N")

rows <- list()
k <- 0
for (ct in cell_types) {
  for (cmp in comparisons) {
    f <- file.path(de_dir, sprintf("DE_%s_%s.csv", ct, cmp))
    if (!file.exists(f)) next
    d <- read.csv(f, stringsAsFactors = FALSE)
    d <- d[d$gene %in% genes, , drop = FALSE]
    if (nrow(d) == 0) next
    for (i in seq_len(nrow(d))) {
      k <- k + 1
      rows[[k]] <- data.frame(
        gene = d$gene[i], cell_type = ct, comparison = cmp,
        our_log2FC = d$avg_log2FC[i], our_FDR = d$p_val_adj[i],
        pct_disease = d$pct.1[i], pct_normal = d$pct.2[i],
        stringsAsFactors = FALSE)
    }
  }
}
ref <- do.call(rbind, rows)
# flag genes absent from a given DE table (tested but not significant output / not detected)
ref$our_direction <- ifelse(ref$our_log2FC > 0, "up", "down")
write.csv(ref, out_csv, row.names = FALSE)
cat("Wrote", nrow(ref), "gene x cell_type x comparison rows ->", out_csv, "\n")

# quick summary of the headline genes
key <- c("FCN3", "ASPN", "PLTP", "HMGN2", "RARRES1", "EDNRB", "ISLR", "BCL6", "LAPTM5")
print(subset(ref, gene %in% key & our_FDR < 0.05,
             select = c(gene, cell_type, comparison, our_log2FC, our_FDR)))
