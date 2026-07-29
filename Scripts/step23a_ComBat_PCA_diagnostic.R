# step23a: ComBat 前后 PCA 诊断图 (审稿意见: 验证批次效应是否被充分移除)
# 复现 step20b 的数据读取与 split-first ComBat 流程(同一种子), 然后对
# (1) ComBat 前的原始合并矩阵  (2) ComBat 后的训练集矩阵
# 分别做 PCA, 按队列和按诊断着色, 输出 4 面板图。

.libPaths(c("Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(sva)
  library(caret)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(12345)
out_dir <- "step22_review_response"
if (!dir.exists(out_dir)) dir.create(out_dir)

datasets <- c("GSE16499", "GSE5406", "GSE57338")
train_ratio <- 0.7

# 1. HF 签名 --------------------------------------------------------------
sig_up <- read.table("results_scRNA/step3_4_5_pathway_cross_signature/HF_signature_up.txt",
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)$gene
sig_down <- read.table("results_scRNA/step3_4_5_pathway_cross_signature/HF_signature_down.txt",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)$gene
sig_genes <- unique(c(sig_up, sig_down))

# 2. 读取三个 GEO 矩阵 -----------------------------------------------------
exp_list <- list(); group_list <- list()
for (ds in datasets) {
  exp <- read.csv(paste0(ds, "_symbol.csv"), header = TRUE, row.names = 1, check.names = FALSE)
  grp <- read.csv(paste0(ds, "_group.csv"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  grp$group <- trimws(as.character(grp$group))
  common <- intersect(colnames(exp), grp$geo_accession)
  exp_list[[ds]] <- exp[, common, drop = FALSE]
  group_list[[ds]] <- setNames(grp$group[match(common, grp$geo_accession)], common)
}
common_genes <- Reduce(intersect, lapply(exp_list, rownames))
feature_genes <- intersect(sig_genes, common_genes)
merged <- do.call(cbind, lapply(exp_list, function(m) m[feature_genes, , drop = FALSE]))
merged <- as.matrix(merged)
colnames(merged) <- unlist(lapply(exp_list, colnames))
sample_dataset <- setNames(unlist(lapply(datasets, function(ds) rep(ds, ncol(exp_list[[ds]])))),
                           unlist(lapply(exp_list, colnames)))
y <- factor(unlist(lapply(datasets, function(ds) unname(group_list[[ds]]))),
            levels = c("Normal", "Disease"))
y_named <- setNames(y, colnames(merged))   # unlist 的名字已被 unname 破坏, 重建命名索引

# 3. 与 step20b 相同的 split + ComBat(仅训练集拟合) ------------------------
train_idx <- createDataPartition(y, p = train_ratio, list = FALSE)[, 1]
raw_train <- merged[, train_idx, drop = FALSE]
y_train <- y[train_idx]
batch_train <- factor(sample_dataset[colnames(raw_train)], levels = datasets)
mod_train <- model.matrix(~ y_train)
combat_train <- ComBat(dat = as.matrix(raw_train), batch = batch_train,
                       mod = mod_train, par.prior = TRUE)

# 4. PCA 函数 ---------------------------------------------------------------
pca_df <- function(mat, label) {
  pc <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  ve <- round(100 * summary(pc)$importance[2, 1:2], 1)
  data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
             Dataset = sample_dataset[colnames(mat)],
             Group = y_named[colnames(mat)],
             panel = label, ve1 = ve[1], ve2 = ve[2], check.names = FALSE)
}
d_before <- pca_df(merged, "Before ComBat (all samples)")
d_after  <- pca_df(combat_train, "After ComBat (training set)")

mk_plot <- function(d, color_var, title_suffix, pal) {
  ggplot(d, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    geom_point(size = 1.6, alpha = 0.8) +
    theme_bw(base_size = 11) +
    labs(title = paste0(unique(d$panel), " — ", title_suffix),
         x = paste0("PC1 (", unique(d$ve1), "%)"),
         y = paste0("PC2 (", unique(d$ve2), "%)"),
         color = color_var) +
    scale_color_manual(values = pal)
}
pal_ds <- c(GSE16499 = "#E41A1C", GSE5406 = "#377EB8", GSE57338 = "#4DAF4A")
pal_grp <- c(Normal = "#2166AC", Disease = "#B2182B")

p1 <- mk_plot(d_before, "Dataset", "by cohort", pal_ds)
p2 <- mk_plot(d_before, "Group", "by diagnosis", pal_grp)
p3 <- mk_plot(d_after, "Dataset", "by cohort", pal_ds)
p4 <- mk_plot(d_after, "Group", "by diagnosis", pal_grp)

p_all <- (p1 | p2) / (p3 | p4)
ggsave(file.path(out_dir, "ComBat_PCA_before_after.pdf"), p_all, width = 12, height = 9)

# 5. 定量: 队列间中心距离 (PCA 前两主成分, 按队列均值的 pairwise 欧氏距离)
cent_dist <- function(d, tag) {
  cent <- aggregate(cbind(PC1, PC2) ~ Dataset, data = d, FUN = mean)
  m <- as.matrix(dist(cent[, c("PC1", "PC2")]))
  dimnames(m) <- list(cent$Dataset, cent$Dataset)
  cat(tag, " between-cohort centroid distances:\n"); print(round(m, 2))
  m
}
m_before <- cent_dist(d_before, "Before ComBat")
m_after  <- cent_dist(d_after, "After ComBat")

sink(file.path(out_dir, "ComBat_PCA_centroid_distances.txt"))
cat("Between-cohort centroid distances in PC1-PC2 space (smaller = better mixed)\n\n")
cat("Before ComBat:\n"); print(round(m_before, 3))
cat("\nAfter ComBat (training set):\n"); print(round(m_after, 3))
sink()

cat("\nDone. ComBat_PCA_before_after.pdf written to", out_dir, "\n")
