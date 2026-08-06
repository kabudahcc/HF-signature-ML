# step26: Fig 1C statistical tests (Wilcoxon rank-sum, BH FDR) + updated figure.
# Run with cwd = workspace root (Paper2_HF_signature_ML/).
suppressPackageStartupMessages(library(ggplot2))
set.seed(1)

df <- read.csv("Data/scatlas_HF_signature_scores_per_cell.csv")
df$condition <- factor(df$condition, levels = c("N", "DCM", "ICM"))
ct_levels <- c("CM", "Fibroblast", "Endothelial", "Pericyte/SMC", "Macrophage",
               "T_cell", "NK", "B_cell", "Mast", "Neutrophil")
df <- df[df$cell_type %in% ct_levels, ]
df$cell_type <- factor(df$cell_type, levels = ct_levels)

panels <- c("HF_up_score" = "HF-up signature (367 genes)",
            "HF_down_score" = "HF-down signature (926 genes)")

## ---- statistics: Wilcoxon DCM vs N and ICM vs N, per cell type per panel --
stats <- list()
for (pv in names(panels)) {
  for (ct in ct_levels) {
    sub <- df[df$cell_type == ct, ]
    for (cmp in c("DCM", "ICM")) {
      x <- sub[[pv]][sub$condition == cmp]
      y <- sub[[pv]][sub$condition == "N"]
      if (length(x) < 3 || length(y) < 3) next
      w <- suppressWarnings(wilcox.test(x, y))
      stats[[length(stats) + 1]] <- data.frame(
        Panel = panels[[pv]], Cell_type = ct, Comparison = paste(cmp, "vs. N"),
        n_case = length(x), n_normal = length(y),
        median_case = round(median(x), 3), median_normal = round(median(y), 3),
        P = signif(w$p.value, 3))
    }
  }
}
st <- do.call(rbind, stats)
st$Padj <- p.adjust(st$P, method = "BH")
st$Padj_label <- ifelse(st$Padj < 0.001, "Padj < 0.001",
                        paste0("Padj = ", signif(st$Padj, 2)))
write.csv(st, "revision_round2/fig1C_wilcoxon_stats.csv", row.names = FALSE)
print(st[, c("Panel", "Cell_type", "Comparison", "Padj")])

## ---- updated figure with significance stars ------------------------------
stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                     ifelse(p < 0.05, "*", "ns")))
st$star <- stars(st$Padj)

long <- data.frame(
  cell_type = rep(df$cell_type, 2),
  condition = rep(df$condition, 2),
  score     = c(df$HF_up_score, df$HF_down_score),
  panel     = factor(rep(panels, each = nrow(df)), levels = panels))

# annotation positions: one star row per comparison above each cell type
ann <- list()
for (i_pan in seq_along(panels)) {
  pan <- panels[[i_pan]]
  pv <- names(panels)[i_pan]
  sub_long <- long[long$panel == pan, ]
  for (i_ct in seq_along(ct_levels)) {
    ct <- ct_levels[i_ct]
    vals <- sub_long$score[sub_long$cell_type == ct]
    y_base <- max(vals, na.rm = TRUE)
    y_span <- diff(range(sub_long$score, na.rm = TRUE))
    s_dcm <- st$star[st$Panel == pan & st$Cell_type == ct & st$Comparison == "DCM vs. N"]
    s_icm <- st$star[st$Panel == pan & st$Cell_type == ct & st$Comparison == "ICM vs. N"]
    ann[[length(ann) + 1]] <- data.frame(panel = pan, cell_type = ct,
      x = i_ct + 0.13, y = y_base + 0.04 * y_span, label = s_dcm)
    ann[[length(ann) + 1]] <- data.frame(panel = pan, cell_type = ct,
      x = i_ct + 0.40, y = y_base + 0.04 * y_span, label = s_icm)
  }
}
ann <- do.call(rbind, ann)
ann$panel <- factor(ann$panel, levels = panels)
ann$cell_type <- factor(ann$cell_type, levels = ct_levels)

p <- ggplot(long, aes(x = cell_type, y = score, fill = condition)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, lwd = 0.3, width = 0.7,
               position = position_dodge(width = 0.8)) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(N = "#4DBBD5", DCM = "#E64B35", ICM = "#3C5488")) +
  geom_text(data = ann, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.2, vjust = 0) +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = "Signature score (AddModuleScore)", fill = "Condition",
       caption = "Wilcoxon rank-sum tests vs. normal within each cell type, BH FDR across all tests; *** Padj < 0.001, ** Padj < 0.01, * Padj < 0.05.\nPer-cell tests are inflated by pseudoreplication (single normal donor); see Limitations.") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        strip.background = element_rect(fill = "grey90"),
        plot.caption = element_text(size = 7.5, hjust = 0, color = "grey35"))

ggsave("Figures/scatlas_HF_signature_scores_by_celltype.png", p,
       width = 8, height = 8.6, dpi = 300)
ggsave("revision_round2/fig1C_with_stats.pdf", p, width = 8, height = 8.6)
cat("fig1C updated\n")
