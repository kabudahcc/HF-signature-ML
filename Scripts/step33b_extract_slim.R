# step33b_extract_slim.R
# Load the giant GSE183852 integrated Seurat object ONCE, immediately extract:
#   1) full cell metadata -> CSV
#   2) counts rows for the 39 signature genes -> slim RDS (sparse, tiny)
# Progress is appended to step33b_extract.log so a timeout kill still leaves a trail.
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))

ws  <- "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
dat <- file.path(ws, "data_external")
out <- file.path(ws, "revision_round2")
lg  <- file.path(out, "step33b_extract.log")

say <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = lg, append = TRUE) }
say("=== start %s ===", format(Sys.time()))

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
say("Seurat loaded %s", format(Sys.time()))

genes <- read.delim(file.path(ws, "Data", "HF_signature_ML_model_genes.txt"))$gene
genes <- unique(trimws(genes))

e <- new.env()
say("loading RData ... %s", format(Sys.time()))
t0 <- Sys.time()
load(file.path(dat, "GSE183852_DCM_Integrated.Robj"), envir = e)
say("load done in %.1f min; objects: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")),
    paste(ls(e), collapse = ", "))

handle_seurat <- function(obj, tag) {
  say("[%s] class=%s assays=%s", tag, class(obj)[1], paste(Assays(obj), collapse = ","))
  md <- obj@meta.data
  say("[%s] cells=%d meta cols: %s", tag, nrow(md), paste(colnames(md), collapse = " | "))
  write.csv(md, file.path(out, sprintf("gse183852_metadata_%s.csv", tag)))
  say("[%s] metadata csv written", tag)
  for (a in Assays(obj)) {
    mat <- tryCatch(GetAssayData(obj, assay = a, layer = "counts"),
                    error = function(er) tryCatch(GetAssayData(obj, assay = a, slot = "counts"),
                                                  error = function(er2) NULL))
    if (is.null(mat)) { say("[%s] assay %s: no counts layer", tag, a); next }
    say("[%s] assay %s dims %d x %d class %s", tag, a, nrow(mat), ncol(mat), class(mat)[1])
    g2 <- intersect(genes, rownames(mat))
    say("[%s] assay %s: %d/39 signature genes present", tag, a, length(g2))
    if (length(g2) > 0) {
      slim <- mat[g2, , drop = FALSE]
      saveRDS(slim, file.path(out, sprintf("gse183852_counts39_%s_%s.rds", tag, a)))
      say("[%s] assay %s slim counts saved (%.1f MB)", tag, a, object.size(slim) / 1e6)
    }
    rm(mat); gc(verbose = FALSE)
  }
}

for (nm in ls(e)) {
  x <- e[[nm]]
  if (is(x, "Seurat")) {
    handle_seurat(x, nm)
  } else if (is.list(x) && !isS4(x)) {
    for (sub in names(x)) if (is(x[[sub]], "Seurat")) handle_seurat(x[[sub]], paste0(nm, "_", sub))
  } else {
    say("[%s] non-Seurat object class=%s size=%.1fMB", nm, class(x)[1], object.size(x) / 1e6)
  }
}
say("=== all done %s ===", format(Sys.time()))
