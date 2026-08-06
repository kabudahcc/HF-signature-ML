si <- sessionInfo()
cat(si$R.version$version.string, "\n")
pkgs <- c("Seurat","SeuratObject","harmony","glmnet","randomForest","caret",
          "kernelshap","susieR","data.table","dplyr","ggplot2","patchwork",
          "pROC","foreach","doParallel","Matrix","reticulate")
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  cat(p, ifelse(is.na(v), "NOT_INSTALLED", v), "\n")
}
