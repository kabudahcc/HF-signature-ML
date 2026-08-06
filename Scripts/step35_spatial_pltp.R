# step35_spatial_pltp.R
# Spatial validation of signature-gene directions using Kuppe et al. 2022 (Nature)
# Visium slides: 4 controls (P1/P7/P8/P17), IZ_P10 (ischemic zone, MI), FZ_P20 (fibrotic zone, MI).
# Predictions from our GSE145154 scRNA results:
#   PLTP: fibroblast UP in HF, macrophage DOWN in HF  ->
#         Fibroblast-dominant spots: MI > control; Myeloid-dominant spots: MI < control
#   ASPN/ISLR/HMGN2: fibroblast UP   RARRES1/EDNRB: fibroblast DOWN   FCN3: endothelial DOWN
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(hdf5r); library(Matrix); library(ggplot2) })

ws  <- "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
dd  <- file.path(ws, "data_external")
out <- file.path(ws, "revision_round2")

slides <- c(control_P1 = "kuppe_control_P1.h5ad", control_P7 = "kuppe_control_P7.h5ad",
            control_P8 = "kuppe_control_P8.h5ad", control_P17 = "kuppe_control_P17.h5ad",
            IZ_P10 = "kuppe_IZ_P10.h5ad", FZ_P20 = "kuppe_FZ_P20.h5ad")
genes39 <- unique(trimws(read.delim(file.path(ws, "Data", "HF_signature_ML_model_genes.txt"))$gene))
ct_cols <- c("Fibroblast", "Myeloid", "Endothelial", "Cardiomyocyte", "Pericyte", "vSMCs")

read_col <- function(node) {
  if (inherits(node, "H5Group")) {
    node[["categories"]]$read()[node[["codes"]]$read() + 1L]
  } else {
    node$read()
  }
}

read_slide <- function(fn) {
  f <- H5File$new(file.path(dd, fn), mode = "r")
  on.exit(f$close_all())
  gname <- read_col(f[["var"]][["feature_name"]])
  X <- f[["X"]]
  dat <- X[["data"]]$read(); idx <- X[["indices"]]$read(); ptr <- X[["indptr"]]$read()
  shp <- h5attributes(X)$shape
  M <- sparseMatrix(i = rep(seq_len(length(ptr) - 1), diff(ptr)), j = idx + 1L, x = dat,
                    dims = shp, index1 = TRUE)
  colnames(M) <- gname
  ob <- f[["obs"]]
  props <- sapply(ct_cols, function(ct) ob[[ct]]$read())
  colnames(props) <- ct_cols
  disease <- read_col(ob[["disease"]])
  lib <- Matrix::rowSums(M)
  spatial <- f[["obsm"]][["spatial"]]$read()
  # hdf5r returns arrays in column-major order: orient to spots x 2
  if (nrow(spatial) != length(lib) && ncol(spatial) == length(lib)) spatial <- t(spatial)
  # keep only spots with complete deconvolution proportions (max.col returns NA
  # for rows containing NA, which previously leaked NA rows into every group)
  keep <- complete.cases(props) & lib > 0
  if ("in_tissue" %in% names(ob)) keep <- keep & (ob[["in_tissue"]]$read() == 1)
  cat("  ", fn, ": kept", sum(keep), "/", length(keep),
      "spots (complete deconvolution proportions, lib > 0)\n")
  list(M = M[keep, , drop = FALSE], gname = gname,
       props = as.data.frame(props[keep, , drop = FALSE]), disease = disease[keep],
       lib = lib[keep], spatial = spatial[keep, , drop = FALSE])
}

spots_all <- list()
for (sl in names(slides)) {
  s <- read_slide(slides[sl])
  g <- intersect(genes39, s$gname)
  cat(sl, ": disease =", unique(s$disease), "| genes found:", length(g), "/", length(genes39), "\n")
  logcpm <- log2(sweep(as.matrix(s$M[, g, drop = FALSE]) + 0.5, 1, s$lib / 1e6, `/`))
  colnames(logcpm) <- g
  dom <- colnames(s$props)[max.col(s$props)]
  df <- data.frame(slide = sl, zone = sub("_.*$", "", sl), disease = s$disease[1],
                   dom = dom, s$props, logcpm, check.names = FALSE)
  spots_all[[sl]] <- df
}

all_spots <- do.call(rbind, spots_all)
all_spots$zone <- factor(all_spots$zone, levels = c("control", "IZ", "FZ"))

# ---- test 1: gene x dominant-cell-type contrasts vs control (pooled spots + donor-level) ----
targets <- list(
  PLTP    = c("Fibroblast", "Myeloid"),
  ASPN    = "Fibroblast", ISLR = "Fibroblast", HMGN2 = "Fibroblast",
  RARRES1 = "Fibroblast", EDNRB = "Fibroblast",
  FCN3    = "Endothelial"
)
our_dir <- c(PLTP_Fibroblast = "up", PLTP_Myeloid = "down", ASPN_Fibroblast = "up",
             ISLR_Fibroblast = "up", HMGN2_Fibroblast = "up", RARRES1_Fibroblast = "down",
             EDNRB_Fibroblast = "down", FCN3_Endothelial = "down")
rows <- list(); k <- 0
for (g in names(targets)) {
  if (!g %in% colnames(all_spots)) next
  for (ct in targets[[g]]) {
    sub <- all_spots[all_spots$dom == ct, ]
    y <- sub[[g]][sub$zone == "control"]
    if (length(y) < 10) next
    # per-donor means (each slide = one donor)
    dm <- tapply(sub[[g]], sub$slide, mean)
    ctrl_means <- dm[grepl("^control", names(dm))]
    for (z in c("IZ", "FZ")) {
      x <- sub[[g]][sub$zone == z]
      if (length(x) < 10) next
      z_mean <- mean(dm[grepl(paste0("^", z), names(dm))])
      wt <- suppressWarnings(wilcox.test(x, y))
      k <- k + 1
      rows[[k]] <- data.frame(
        gene = g, dominant_ct = ct, zone = z,
        n_zone = length(x), n_control = length(y),
        mean_zone = z_mean, mean_control = mean(y),
        ctrl_donor_min = min(ctrl_means), ctrl_donor_max = max(ctrl_means),
        outside_ctrl_donor_range = z_mean < min(ctrl_means) | z_mean > max(ctrl_means),
        delta = z_mean - mean(y), wilcox_p = wt$p.value)
    }
  }
}
tt <- do.call(rbind, rows)
tt$FDR <- p.adjust(tt$wilcox_p, method = "BH")
tt$expected <- our_dir[paste(tt$gene, tt$dominant_ct, sep = "_")]
tt$observed <- ifelse(tt$delta > 0, "up", "down")
tt$match <- tt$expected == tt$observed
write.csv(tt, file.path(out, "spatial_zone_contrasts.csv"), row.names = FALSE)
cat("\n=== zone contrasts (MI zone vs 4-donor control, within dominant cell type) ===\n")
print(tt[, c("gene", "dominant_ct", "zone", "n_zone", "n_control", "mean_zone",
             "ctrl_donor_min", "ctrl_donor_max", "outside_ctrl_donor_range",
             "delta", "FDR", "expected", "match")])

# ---- test 2: per-spot Spearman correlation with cell-type proportions ----
cors <- list(); k <- 0
for (z in levels(all_spots$zone)) {
  sub <- all_spots[all_spots$zone == z, ]
  for (g in c("PLTP", "ASPN", "FCN3")) {
    if (!g %in% colnames(sub)) next
    for (ct in c("Fibroblast", "Myeloid", "Endothelial")) {
      k <- k + 1
      cors[[k]] <- data.frame(zone = z, gene = g, prop = ct,
                              rho = suppressWarnings(cor(sub[[g]], sub[[ct]], method = "spearman")))
    }
  }
}
cc <- do.call(rbind, cors)
write.csv(cc, file.path(out, "spatial_proportion_correlations.csv"), row.names = FALSE)
cat("\n=== spot-level Spearman correlations (expression vs deconvolution proportion) ===\n")
print(cc)

# ---- figure: PLTP spatial maps (6 slides) + zone boxplot ----
library(patchwork)
global_mid <- median(all_spots$PLTP)
maps <- list()
for (sl in names(slides)) {
  s <- read_slide(slides[sl])
  i <- which(s$gname == "PLTP")
  pltp <- log2((s$M[, i] + 0.5) / (s$lib / 1e6))
  dd_ <- data.frame(x = s$spatial[, 1], y = s$spatial[, 2], pltp = as.vector(pltp))
  maps[[sl]] <- ggplot(dd_, aes(x, y, color = pltp)) +
    geom_point(size = 0.8) +
    scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b", midpoint = global_mid,
                          name = "PLTP log2CPM") +
    coord_fixed() + theme_void(base_size = 10) +
    labs(title = paste0(sl, " (", s$disease[1], ")")) +
    theme(plot.title = element_text(size = 10, face = "bold"), legend.position = "none")
}
box <- ggplot(all_spots[all_spots$dom %in% c("Fibroblast", "Myeloid"), ],
              aes(zone, PLTP, fill = zone)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6) +
  facet_wrap(~ dom) +
  scale_fill_manual(values = c(control = "#7fbc41", IZ = "#f4a582", FZ = "#b2182b")) +
  labs(y = "PLTP log2 CPM", x = NULL, title = "PLTP by dominant cell type and zone") +
  theme_bw(base_size = 10) + theme(legend.position = "none",
                                   plot.title = element_text(size = 10, face = "bold"))
map_grid <- wrap_plots(maps, nrow = 2) +
  plot_annotation(title = "PLTP spatial expression across Kuppe et al. Visium slides")
fig <- map_grid / box + plot_layout(heights = c(2, 1))
ggsave(file.path(out, "S_spatial_PLTP.png"), fig, width = 12, height = 9, dpi = 300)
cat("\nfigure + tables written\n")
