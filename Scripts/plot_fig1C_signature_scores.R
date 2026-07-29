# Figure 1C: HF-up / HF-down consensus signature scores by cell type and condition
# Input:  results_scRNA/step3_4_5_pathway_cross_signature/HF_signature_scores_per_cell.csv
# Output: for_publication/Paper2_HF_signature_ML/Figures/scatlas_HF_signature_scores_by_celltype.pdf

.libPaths(c("Rlib", .libPaths()))
library(ggplot2)

df <- read.csv("results_scRNA/step3_4_5_pathway_cross_signature/HF_signature_scores_per_cell.csv")
df$condition <- factor(df$condition, levels = c("N", "DCM", "ICM"))
df$cell_type <- factor(df$cell_type,
  levels = c("CM", "Fibroblast", "Endothelial", "Pericyte/SMC", "Macrophage",
             "T_cell", "NK", "B_cell", "Mast", "Neutrophil"))

long <- data.frame(
  cell_type = rep(df$cell_type, 2),
  condition = rep(df$condition, 2),
  score     = c(df$HF_up_score, df$HF_down_score),
  panel     = factor(rep(c("HF-up signature (367 genes)", "HF-down signature (926 genes)"),
                         each = nrow(df)),
                     levels = c("HF-up signature (367 genes)", "HF-down signature (926 genes)")))
long <- long[!is.na(long$cell_type), ]

p <- ggplot(long, aes(x = cell_type, y = score, fill = condition)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, lwd = 0.3, width = 0.7,
               position = position_dodge(width = 0.8)) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(N = "#4DBBD5", DCM = "#E64B35", ICM = "#3C5488")) +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = "Signature score (AddModuleScore)", fill = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        strip.background = element_rect(fill = "grey90"))

ggsave("for_publication/Paper2_HF_signature_ML/Figures/scatlas_HF_signature_scores_by_celltype.pdf",
       p, width = 8, height = 8)
cat("saved\n")
