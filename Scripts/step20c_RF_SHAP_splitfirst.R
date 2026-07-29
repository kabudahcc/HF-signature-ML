# step20c: 在 split-first 流程上重训 RF 主模型并计算 kernel SHAP
# 背景: step20b benchmark 中 LDA 与 RF 并列最优(ΔAUC = 0.0004),
#       基于 SHAP 可解释性与外部验证一致性选择 RF 为主模型;
#       step20b 保存的 SHAP 是 LDA 的, 这里补算 RF 版本。
# 输出: step20b_HF_signature_ML_splitfirst/ 下 RF_ 前缀的 SHAP 图表与数值

.libPaths(c("Rlib", .libPaths()))
rm(list = ls())

suppressPackageStartupMessages({
  library(randomForest)
  library(kernelshap)
  library(shapviz)
  library(ggplot2)
})

set.seed(12345)

out_dir <- "step20b_HF_signature_ML_splitfirst"

bundle <- readRDS(file.path(out_dir, "best_model.rds"))
sel <- bundle$features
train_df <- bundle$train_set
cat("Retrain RF on", nrow(train_df), "samples,", length(sel), "features\n")

model_rf <- randomForest(Group ~ ., data = train_df, ntree = 500)

pred_fun <- function(object, newdata) {
  as.numeric(predict(object, newdata = as.data.frame(newdata), type = "prob")[, "Disease"])
}

X_bg <- train_df[, sel]
cat("Running kernelshap ...\n")
shap_fit <- kernelshap(model_rf, X = X_bg, bg_X = X_bg, pred_fun = pred_fun)
sv <- shapviz(shap_fit, X_pred = X_bg, X = X_bg)

shap_importance <- sort(colMeans(abs(sv$S)), decreasing = TRUE)
write.csv(data.frame(Gene = names(shap_importance), MeanAbsSHAP = as.numeric(shap_importance)),
          file.path(out_dir, "RF_shap_importance.csv"), row.names = FALSE)
write.csv(as.data.frame(sv$S), file.path(out_dir, "RF_shap_values_train.csv"), row.names = TRUE)
print(round(shap_importance, 4))

pdf(file.path(out_dir, "RF_SHAP_importance_barplot.pdf"), width = 6, height = 6)
print(sv_importance(sv, kind = "bar", show_numbers = TRUE) + theme_bw())
dev.off()

pdf(file.path(out_dir, "RF_SHAP_beeswarm.pdf"), width = 7, height = 6)
print(sv_importance(sv, kind = "bee", show_numbers = TRUE) + theme_bw())
dev.off()

top_n <- min(6, length(sel))
pdf(file.path(out_dir, "RF_SHAP_dependence_top.pdf"), width = 10, height = 7)
print(sv_dependence(sv, v = names(shap_importance)[seq_len(top_n)]) + theme_bw())
dev.off()

pdf(file.path(out_dir, "RF_SHAP_waterfall_sample1.pdf"), width = 7, height = 6)
print(sv_waterfall(sv, row_id = 1) + theme_bw())
dev.off()

# 保存 RF 模型本体, 供外部验证/投稿包使用
saveRDS(list(model = model_rf, method = "RF", features = sel,
             train_set = train_df, lasso_cv = bundle$lasso_cv,
             combat_params = bundle$combat_params),
        file.path(out_dir, "best_model_RF.rds"))

cat("\nDone. RF SHAP outputs in", out_dir, "\n")
