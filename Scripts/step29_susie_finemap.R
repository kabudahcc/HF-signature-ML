# step29: SuSiE fine-mapping of the PLTP locus (GTEx Heart-AA eQTL and
# HERMES HF GWAS), with 1000 Genomes EUR LD. Run with cwd = revision_round2/.
.libPaths(c("D:/早早生信分析学习课程/新R文件/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(rtracklayer); library(GenomicRanges)
  library(susieR); library(ggplot2)
})
set.seed(1)

plink <- "D:/早早生信分析学习课程/新R文件/tools/plink.exe"
kg_bfile <- "D:/早早生信分析学习课程/新R文件/MR_data/1kg_ref/EUR_chr20"

## ---- eQTL: liftover hg38 -> hg19 -----------------------------------------
eq <- fread("data/PLTP_AA_allpairs.tsv")
gr <- GRanges(seqnames = paste0("chr", eq$chromosome),
              ranges = IRanges(eq$position, width = 1))
names(gr) <- as.character(seq_along(gr))
chain <- import.chain("data/hg38ToHg19.over.chain")
lo <- unlist(liftOver(gr, chain))
eq$BP19 <- NA_integer_
hit <- as.integer(names(lo))
eq$BP19[hit] <- start(lo)
eq <- eq[!is.na(eq$BP19) & !is.na(eq$rsid) & eq$rsid != ""]
cat("eQTL variants after liftover:", nrow(eq), "\n")

## ---- GWAS ------------------------------------------------------------------
gw <- fread("data/HERMES_PLTP_region.tsv")
gw <- gw[!is.na(gw$SNP) & gw$SNP != ""]
gw$A1 <- toupper(gw$A1); gw$A2 <- toupper(gw$A2)
cat("GWAS variants in region:", nrow(gw), "\n")

## ---- harmonize by rsid -----------------------------------------------------
common <- intersect(eq$rsid, gw$SNP)
eq <- eq[match(common, eq$rsid)]
gw <- gw[match(common, gw$SNP)]
stopifnot(all(eq$rsid == gw$SNP))
# alleles: eQTL effect allele = alt; GWAS effect allele = A1
same <- eq$alt == gw$A1 & eq$ref == gw$A2
flip <- eq$ref == gw$A1 & eq$alt == gw$A2
keep <- same | flip
eq <- eq[keep]; gw <- gw[keep]; flip <- flip[keep]
gw$b[flip] <- -gw$b[flip]
gw$freq[flip] <- 1 - gw$freq[flip]
cat("harmonized variants:", nrow(gw), "\n")
z_eq <- eq$beta / eq$se
z_gw <- gw$b / gw$se

## ---- LD from 1000G EUR ------------------------------------------------------
writeLines(gw$SNP, "data/pltp_harmonized_snps.txt")
cmd <- sprintf('"%s" --bfile "%s" --extract data/pltp_harmonized_snps.txt --r square --out data/PLTP_AA_ld',
               plink, kg_bfile)
system(cmd, intern = TRUE)
bim_all <- fread("D:/早早生信分析学习课程/新R文件/MR_data/1kg_ref/EUR_chr20.bim", header = FALSE)
bim <- bim_all[bim_all$V2 %in% gw$SNP, ]
ld <- as.matrix(fread("data/PLTP_AA_ld.ld", header = FALSE))
idx <- match(bim$V2, gw$SNP)
keep2 <- !is.na(idx)
ld <- ld[keep2, keep2]
gw2 <- gw[idx[keep2]]; eq2 <- eq[idx[keep2]]
bim2 <- bim[keep2, ]
# orient LD to the harmonized effect allele (gw2$A1 == eQTL alt);
# plink --r codes genotypes by bim A1 (V5)
eff_a1 <- gw2$A1; eff_a2 <- gw2$A2
bim_a1 <- bim2$V5; bim_a2 <- bim2$V6
orient <- ifelse(bim_a1 == eff_a1 & bim_a2 == eff_a2, 1,
                 ifelse(bim_a1 == eff_a2 & bim_a2 == eff_a1, -1, NA))
ok <- !is.na(orient)
cat("variants dropped for allele mismatch with LD panel:", sum(!ok), "\n")
# drop strand-ambiguous A/T and C/G SNPs (letter matching cannot detect strand flips)
amb <- (eff_a1 %in% c("A", "T") & eff_a2 %in% c("A", "T")) |
       (eff_a1 %in% c("C", "G") & eff_a2 %in% c("C", "G"))
cat("strand-ambiguous (A/T, C/G) variants dropped:", sum(ok & amb), "\n")
ok <- ok & !amb
ld <- ld[ok, ok]; orient <- orient[ok]
gw2 <- gw2[ok]; eq2 <- eq2[ok]
D <- diag(orient)
ld <- D %*% ld %*% D
z_gw2 <- (gw2$b / gw2$se); z_eq2 <- (eq2$beta / eq2$se)
ld <- (ld + t(ld)) / 2
diag(ld) <- 1
# project to nearest positive-definite correlation matrix (reference-panel
# LD is only approximately PSD, which breaks IBSS convergence)
suppressPackageStartupMessages(library(Matrix))
ld <- as.matrix(nearPD(ld, corr = TRUE, maxit = 200)$mat)
cat("z/LD consistency diagnostic s (estimate_s_rss):\n")
cat("  GWAS:", susieR:::estimate_s_rss(z_gw2, ld, n = 964057), "\n")
cat("  eQTL:", susieR:::estimate_s_rss(z_eq2, ld, n = max(eq2$an) / 2), "\n")
cat("final variant set:", length(z_gw2), "\n")

## ---- SuSiE RSS ----------------------------------------------------------------
# reference-panel LD: disable prior-variance estimation (with only 503 LD
# samples vs n=964k GWAS the MLE can blow up); fix prior at default scale
fit_gw <- susie_rss(z = z_gw2, R = ld, n = 964057, L = 5, max_iter = 1000,
                    estimate_prior_variance = FALSE, check_prior = FALSE)
fit_eq <- susie_rss(z = z_eq2, R = ld, n = max(eq2$an) / 2, L = 5,
                    max_iter = 1000,
                    estimate_prior_variance = FALSE, check_prior = FALSE)
cat("GWAS converged:", fit_gw$converged, "| eQTL converged:", fit_eq$converged, "\n")

summ_gw <- summary(fit_gw)
summ_eq <- summary(fit_eq)
cat("\n=== GWAS credible sets ===\n"); print(summ_gw$cs)
cat("\n=== eQTL credible sets ===\n"); print(summ_eq$cs)

pip <- data.frame(rsid = gw2$SNP, BP = gw2$BP,
                  PIP_GWAS = fit_gw$pip, PIP_eQTL = fit_eq$pip,
                  z_GWAS = z_gw2, z_eQTL = z_eq2)
pip <- pip[order(-pip$PIP_GWAS), ]
write.csv(pip, "susie_PLTP_AA_PIP.csv", row.names = FALSE)
cat("\nTop GWAS PIP:\n"); print(head(pip, 8))
cat("\nTop eQTL PIP:\n")
print(head(pip[order(-pip$PIP_eQTL), c("rsid", "BP", "PIP_eQTL", "PIP_GWAS")], 8))

cs_out <- rbind(
  if (!is.null(summ_gw$cs)) data.frame(trait = "GWAS", summ_gw$cs),
  if (!is.null(summ_eq$cs)) data.frame(trait = "eQTL", summ_eq$cs))
if (!is.null(cs_out)) write.csv(cs_out, "susie_PLTP_AA_credible_sets.csv", row.names = FALSE)

# overlap of 95% credible sets (variant indices from fit$sets$cs)
gw_members <- if (!is.null(fit_gw$sets$cs)) unique(unlist(fit_gw$sets$cs)) else integer(0)
eq_members <- if (!is.null(fit_eq$sets$cs)) unique(unlist(fit_eq$sets$cs)) else integer(0)
shared <- intersect(gw_members, eq_members)
cat("\nGWAS 95% CS variants:", length(gw_members),
    "| eQTL 95% CS variants:", length(eq_members),
    "| shared:", length(shared), "\n")
if (length(gw_members)) cat("GWAS CS SNPs:", paste(gw2$SNP[gw_members], collapse = ", "), "\n")
if (length(eq_members)) cat("eQTL CS SNPs:", paste(gw2$SNP[eq_members], collapse = ", "), "\n")

sink("susie_PLTP_AA_summary.txt")
cat("SuSiE fine-mapping, PLTP locus, GTEx Heart-AA eQTL vs HERMES HF GWAS\n")
cat("variants:", length(z_gw2), "| LD: 1000G EUR chr20\n\n")
cat("GWAS 95% CS size:", length(gw_members), " SNPs:",
    paste(gw2$SNP[gw_members], collapse = ", "), "\n")
cat("eQTL 95% CS size:", length(eq_members), " SNPs:",
    paste(gw2$SNP[eq_members], collapse = ", "), "\n")
cat("shared variants in both CS:", length(shared),
    ifelse(length(shared), paste(":", paste(gw2$SNP[shared], collapse = ", ")), ""), "\n\n")
cat("Top-5 GWAS PIP:\n"); print(head(pip[, c("rsid", "BP", "PIP_GWAS", "PIP_eQTL")], 5))
cat("Top-5 eQTL PIP:\n"); print(head(pip[order(-pip$PIP_eQTL), c("rsid", "BP", "PIP_eQTL", "PIP_GWAS")], 5))
sink()

## ---- S19 figure: PIP scatter, two panels -----------------------------------
pip$panel_z <- pmin(-log10(2 * pnorm(-abs(pip$z_GWAS))), 20)
p1 <- ggplot(pip, aes(x = BP / 1e6, y = PIP_GWAS)) +
  geom_point(aes(color = PIP_eQTL > 0.1), size = 1.4, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#b2182b"),
                     labels = c("FALSE" = "eQTL PIP <= 0.1", "TRUE" = "eQTL PIP > 0.1")) +
  labs(x = NULL, y = "SuSiE PIP (HF GWAS)", color = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "top")
p2 <- ggplot(pip, aes(x = BP / 1e6, y = PIP_eQTL)) +
  geom_point(aes(color = PIP_GWAS > 0.1), size = 1.4, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#2166ac"),
                     labels = c("FALSE" = "GWAS PIP <= 0.1", "TRUE" = "GWAS PIP > 0.1")) +
  labs(x = "Position on chr20 (Mb, hg19)", y = "SuSiE PIP (Heart-AA eQTL)", color = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "top")
library(patchwork)
p12 <- p1 / p2 + plot_annotation(
  title = "SuSiE fine-mapping of the PLTP locus (GTEx Heart-AA eQTL vs HERMES HF GWAS)")
ggsave("S19_susie_PLTP_PIP.png", p12, width = 7.5, height = 5.5, dpi = 300)
cat("susie done\n")
