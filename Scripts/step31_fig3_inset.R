# step31: add a random-effects meta-analysis inset (forest-style) to main Fig 3.
# Reviewer round-3 item: external-validation ROC should display the pooled estimate.
# Run with cwd = loco_meta_analysis/.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(magick)
})
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))

ma_df <- read.csv("S6_dataset_meta_analysis.csv", stringsAsFactors = FALSE)
pool <- read.csv("S6_dataset_meta_pooled.csv", stringsAsFactors = FALSE)
pool3 <- pool[grepl("3 cohorts", pool$Analysis), ]

fp <- ma_df[, c("Cohort", "AUC", "CI_lower", "CI_upper", "Variance_estimable")]
fp$Label <- sprintf("%s (%d vs %d)%s", fp$Cohort, ma_df$n_Disease, ma_df$n_Normal,
                    ifelse(fp$Variance_estimable, "", " *"))
pooled_row <- data.frame(
  Cohort = "Pooled", AUC = pool3$pooled_AUC,
  CI_lower = pool3$CI_lower, CI_upper = pool3$CI_upper,
  Variance_estimable = TRUE,
  Label = sprintf("Pooled (RE): %.3f (%.3f-%.3f)", pool3$pooled_AUC,
                  pool3$CI_lower, pool3$CI_upper))
fp$Type <- "Cohort"; pooled_row$Type <- "Pooled"
fp <- rbind(fp, pooled_row)
fp$Label <- factor(fp$Label, levels = rev(fp$Label))

inset <- ggplot(fp, aes(x = AUC, y = Label)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey65", linewidth = 0.3) +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper, orientation = "y"),
                width = 0.2, color = "grey30", linewidth = 0.35) +
  geom_point(aes(size = Type, shape = Type), color = "#2166ac") +
  scale_size_manual(values = c(Cohort = 1.5, Pooled = 2.4)) +
  scale_shape_manual(values = c(Cohort = 16, Pooled = 18)) +
  xlim(0.4, 1.02) +
  labs(x = "AUC (95% CI)", y = NULL,
       title = "Random-effects meta-analysis") +
  theme_bw(base_size = 7) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7.5, face = "bold"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3))

roc_img <- magick::image_read("../Figures/HF_signature_external_validation_ROC.png")
roc_grob <- grid::rasterGrob(as.raster(roc_img), interpolate = TRUE)
base <- ggplot() + annotation_custom(roc_grob, -Inf, Inf, -Inf, Inf) +
  coord_fixed() + theme_void()

p <- base + inset_element(inset, left = 0.46, bottom = 0.02, right = 0.99, top = 0.42)
ggsave("../Figures/HF_signature_external_validation_ROC.png", p,
       width = 8, height = 6, dpi = 300)
cat("fig3 inset done\n")
