# step27: pathway-level SHAP aggregation for the 39 model genes.
# Map genes -> GO BP and Reactome via clusterProfiler enrichment
# (universe = 987 candidate features), then aggregate mean|SHAP|.
# Run with cwd = workspace root (Paper2_HF_signature_ML/).
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
})
set.seed(1)

genes <- read.table("Data/HF_signature_ML_model_genes.txt", header = TRUE,
                    sep = "\t", stringsAsFactors = FALSE)$gene
shap <- read.csv("Data/HF_signature_ML_shap_importance.csv", stringsAsFactors = FALSE)
stopifnot(all(genes %in% shap$Gene))
shap_map <- setNames(shap$MeanAbsSHAP, shap$Gene)
total_shap <- sum(shap_map)

sig <- unique(c(
  read.table("loco_meta_analysis/data/HF_signature_up.txt", header = TRUE,
             sep = "\t", stringsAsFactors = FALSE)$gene,
  read.table("loco_meta_analysis/data/HF_signature_down.txt", header = TRUE,
             sep = "\t", stringsAsFactors = FALSE)$gene))
# universe: 987 signature genes measured in all discovery cohorts;
# approximate universe = the 987-gene feature space is unavailable here,
# so use the full 1,263-gene consensus signature as the background.
cat("model genes:", length(genes), "| background:", length(sig), "\n")

g2e <- bitr(c(genes, sig), fromType = "SYMBOL", toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
gene_ent <- g2e$ENTREZID[match(genes, g2e$SYMBOL)]
gene_ent <- gene_ent[!is.na(gene_ent)]
bg_ent <- unique(g2e$ENTREZID[g2e$SYMBOL %in% sig])
ent2sym <- setNames(g2e$SYMBOL, g2e$ENTREZID)

collect <- function(enrich_df, source) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(NULL)
  rows <- list()
  for (i in seq_len(nrow(enrich_df))) {
    members <- strsplit(enrich_df$geneID[i], "/", fixed = TRUE)[[1]]
    mem_sym <- ent2sym[members]
    mem_sym <- mem_sym[!is.na(mem_sym) & mem_sym %in% genes]
    if (length(mem_sym) == 0) next
    s <- sum(shap_map[mem_sym])
    rows[[length(rows) + 1]] <- data.frame(
      Source = source,
      Pathway = enrich_df$Description[i],
      ID = enrich_df$ID[i],
      Padj = signif(enrich_df$p.adjust[i], 3),
      N_model_genes = length(mem_sym),
      Member_genes = paste(mem_sym, collapse = ", "),
      Sum_meanAbsSHAP = round(s, 4),
      SHAP_share_pct = round(100 * s / total_shap, 1),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

ego <- enrichGO(gene = gene_ent, universe = bg_ent, OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1,
                qvalueCutoff = 1, readable = FALSE)
erep <- enrichPathway(gene = gene_ent, universe = bg_ent,
                      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)

cat('GO BP terms:', nrow(as.data.frame(ego)), '| Reactome terms:', nrow(as.data.frame(erep)), '
')
res <- rbind(collect(as.data.frame(ego), "GO BP"),
             collect(as.data.frame(erep), "Reactome"))
res <- res[order(-res$Sum_meanAbsSHAP), ]
write.csv(res, "revision_round2/pathway_SHAP.csv", row.names = FALSE)
print(head(res[, c("Source", "Pathway", "Padj", "N_model_genes",
                   "Sum_meanAbsSHAP", "SHAP_share_pct")], 15))

top <- head(res[res$Source == 'Reactome', ], 12)
top$Pathway <- factor(top$Pathway, levels = rev(top$Pathway))
p <- ggplot(top, aes(x = Sum_meanAbsSHAP, y = Pathway, fill = Source)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(SHAP_share_pct, "%")), hjust = -0.15, size = 3) +
  scale_fill_manual(values = c("GO BP" = "#2166ac", "Reactome" = "#b2182b")) +
  labs(x = "Pathway-level SHAP (sum of mean |SHAP| of member model genes)",
       y = NULL,
       title = "Pathway-level SHAP aggregation of the 39-gene model",
       caption = "Labels: share of total model SHAP. Pathways from clusterProfiler enrichment\n(GO Biological Process and Reactome; background = consensus signature genes).") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        plot.caption = element_text(size = 8, hjust = 0, color = "grey35"),
        legend.position = "none") +
  xlim(0, max(top$Sum_meanAbsSHAP) * 1.15)
ggsave("revision_round2/S18_pathway_SHAP.png", p, width = 8.5, height = 5.5, dpi = 300)
cat("pathway SHAP done\n")
