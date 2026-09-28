# =============================================================================
# 12e_gse125583_dd_genes_corrected.R — Section 3.6 re-run with correctly aligned covariates
#
# PRE-DECLARED (DEVIATIONS 2026-09-28, covariate-alignment correction), before running:
#   Cohort: results/intermediate/gse125583_bulk.rds rebuilt by R/12_0 (SRX-matched age and sex).
#   Age balance: point-biserial r(group, age) and Wilcoxon test.
#   Genes: the 20 DD-signature genes with their original categories (from the earlier table).
#   Models (limma-trend): M0 ~ diagnosis; M1 ~ diagnosis + age + sex;
#     M2 = M1 + BRETIGEA neuron, astrocyte and microglia scores (nMarker 50, scale = TRUE).
#     (The original M2 used marker scores whose gene sets were not recorded; BRETIGEA replaces them.)
#   Reported: per-gene log2FC and FDR under M0-M2; genome-wide DEG counts; no decision rule.
# Output: outputs/A4_GSE125583_DD_genes_corrected.csv, outputs/A4_GSE125583_DD_genes_corrected.rds
# =============================================================================
suppressPackageStartupMessages({ library(limma); library(BRETIGEA) })
g <- readRDS(here::here("results/intermediate/gse125583_bulk.rds"))
stopifnot(grepl("R/12_0", g$note))                      # must be the SRX-matched cohort
m <- g$meta

# age balance
ad <- as.numeric(m$diagnosis == "AD")
r_age <- cor(ad, m$age); w <- wilcox.test(age ~ diagnosis, data = m)
med <- tapply(m$age, m$diagnosis, median)
message(sprintf("Age: median control %.1f, AD %.1f; r(group, age) = %.3f; Wilcoxon P = %.3g",
                med["Control"], med["AD"], r_age, w$p.value))

# DD-signature genes and categories from the earlier table (gene list only; values are recomputed)
old <- read.csv(here::here("outputs/EkTabloS4_GSE125583_age_sex_adjusted.csv"), check.names = FALSE)
gcol <- names(old)[sapply(old, function(x) is.character(x) && mean(x %in% rownames(g$expr)) > 0.5)][1]
ccol <- grep("categ", names(old), ignore.case = TRUE, value = TRUE)[1]
stopifnot(!is.na(gcol), !is.na(ccol))
dd <- data.frame(gene = old[[gcol]], category = old[[ccol]], stringsAsFactors = FALSE)
message("DD genes: ", nrow(dd), " (", paste(names(table(dd$category)), table(dd$category), sep = "=", collapse = ", "), ")")

bret <- as.data.frame(BRETIGEA::brainCells(inputMat = g$expr, nMarker = 50, species = "human",
                      celltypes = c("neu", "ast", "mic"), scale = TRUE))
colnames(bret) <- paste0("bret_", colnames(bret)); m <- cbind(m, bret[m$sample_id, , drop = FALSE])
fit <- function(f) {
  des <- model.matrix(as.formula(f), data = m)
  tt <- topTable(eBayes(lmFit(g$expr, des), trend = TRUE), coef = "diagnosisAD", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt); tt }
M <- list(M0 = fit("~ diagnosis"), M1 = fit("~ diagnosis + age + sex"),
          M2 = fit("~ diagnosis + age + sex + bret_neu + bret_ast + bret_mic"))
n_deg <- sapply(M, function(t) sum(t$adj.P.Val < 0.05)); print(n_deg)
tab <- dd
for (k in names(M)) { t <- M[[k]][match(dd$gene, M[[k]]$gene), ]
  tab[[paste0("log2FC_", k)]] <- t$logFC; tab[[paste0("FDR_", k)]] <- t$adj.P.Val }
print(tab, row.names = FALSE, digits = 3)
S <- tab$category == tab$category[grepl("stress", tab$category, ignore.case = TRUE)][1]
message(sprintf("Stress genes under M2: %d/%d up, %d/%d up with FDR < 0.05",
                sum(tab$log2FC_M2[S] > 0), sum(S), sum(tab$log2FC_M2[S] > 0 & tab$FDR_M2[S] < 0.05), sum(S)))
message(sprintf("Group-score correlations: neuron %.2f, astrocyte %.2f, microglia %.2f",
                cor(ad, m$bret_neu), cor(ad, m$bret_ast), cor(ad, m$bret_mic)))
write.csv(tab, here::here("outputs/A4_GSE125583_DD_genes_corrected.csv"), row.names = FALSE)
saveRDS(list(table = tab, n_deg = n_deg, age = list(median = med, r = r_age, wilcox_p = w$p.value),
             score_group_r = c(neu = cor(ad, m$bret_neu), ast = cor(ad, m$bret_ast), mic = cor(ad, m$bret_mic)),
             run_time = Sys.time()), here::here("outputs/A4_GSE125583_DD_genes_corrected.rds"))
message("Saved outputs/A4_GSE125583_DD_genes_corrected.*")
