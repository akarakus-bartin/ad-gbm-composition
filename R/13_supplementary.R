# =============================================================================
# 13_supplementary.R — Supplementary Tables S1-S4 and Figure S1, from locked result files only
# Output: outputs/supplementary/Supplementary_Tables_S1-S4.xlsx (or CSVs), FigureS1_H1_replication.pdf/.png
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
f   <- function(p) readRDS(here::here(p))
out <- here::here("outputs/supplementary"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
b4 <- f("outputs/stage4_H1_bretigea.rds"); m4 <- f("outputs/stage4_H1_music.rds")
a1 <- f("outputs/A1_H1_replication_gse125583.rds"); a2 <- f("outputs/A2_DD_decomposition.rds")
a3 <- f("outputs/A3_GBM_limma_sensitivity.rds");     qc <- f("outputs/stage3_quadrants_corrected.rds")
gb <- f("results/intermediate/gbm_bulk.rds")

pick <- function(tt, lab) {
  d <- data.frame(gene = tt$gene, tt$logFC, tt$adj.P.Val)
  names(d)[2:3] <- paste0(c("log2FC_", "FDR_"), lab); d
}
merge_all <- function(lst) Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), lst)
shared_genes <- function(a, g) {
  x <- merge(pick(a, "a"), pick(g, "g"), by = "gene")
  x$gene[x$FDR_a < 0.05 & x$FDR_g < 0.05 & sign(x$log2FC_a) == sign(x$log2FC_g)]
}

# ---- Table S1a: shared DEGs, primary analysis (Model A), with all models --------------
sA <- shared_genes(b4$models$AD_A$tt, b4$models$GBM_A$tt)
S1a <- merge_all(list(pick(b4$models$AD_A$tt, "AD_ModelA"), pick(b4$models$AD_B$tt, "AD_ModelB"),
                      pick(m4$models$AD_C$tt, "AD_ModelC"), pick(b4$models$GBM_A$tt, "GBM_ModelA"),
                      pick(b4$models$GBM_B$tt, "GBM_ModelB"), pick(m4$models$GBM_C$tt, "GBM_ModelC")))
S1a <- S1a[S1a$gene %in% sA, ]; S1a <- S1a[order(S1a$FDR_AD_ModelA), ]
stopifnot(nrow(S1a) == 1458)

# ---- Table S1b: shared DEGs, replication cohort (GSE125583, Model A) -------------------
sR <- shared_genes(a1$models$A$tt, b4$models$GBM_A$tt)
S1b <- merge_all(list(pick(a1$models$A$tt, "GSE125583_ModelA"), pick(a1$models$B$tt, "GSE125583_ModelB"),
                      pick(b4$models$GBM_A$tt, "GBM_ModelA"), pick(b4$models$GBM_B$tt, "GBM_ModelB")))
S1b <- S1b[S1b$gene %in% sR, ]
S1b$shared_after_adjustment <- S1b$FDR_GSE125583_ModelB < 0.05 & S1b$FDR_GBM_ModelB < 0.05 &
                               sign(S1b$log2FC_GSE125583_ModelB) == sign(S1b$log2FC_GBM_ModelB)
S1b <- S1b[order(!S1b$shared_after_adjustment, S1b$FDR_GSE125583_ModelA), ]
sB <- shared_genes(a1$models$B$tt, b4$models$GBM_B$tt)
stopifnot(nrow(S1b) == a1$shared$fdr$nA, length(sB) == a1$shared$fdr$nB)
S1b$shared_after_adjustment <- S1b$gene %in% sB
n_in <- sum(S1b$shared_after_adjustment)
S1c <- merge_all(list(pick(a1$models$A$tt, 'GSE125583_ModelA'), pick(a1$models$B$tt, 'GSE125583_ModelB'),
                      pick(b4$models$GBM_A$tt, 'GBM_ModelA'), pick(b4$models$GBM_B$tt, 'GBM_ModelB')))
S1c <- S1c[S1c$gene %in% setdiff(sB, sR), ]
sig <- function(tt) tt$gene[tt$adj.P.Val < 0.05]
adA <- sig(a1$models$A$tt); adB <- sig(a1$models$B$tt); gA <- sig(b4$models$GBM_A$tt); gB <- sig(b4$models$GBM_B$tt)
message(sprintf('Replication shared DEGs: A %d, B %d; of the B genes, %d were shared under A (%.1f%% of A) and %d are newly shared',
                length(sR), length(sB), n_in, 100 * n_in / length(sR), length(sB) - n_in))
message(sprintf('GSE125583 AD DEGs: A %d, B %d, significant in both %d (%.1f%% of A)',
                length(adA), length(adB), length(intersect(adA, adB)), 100 * length(intersect(adA, adB)) / length(adA)))
message(sprintf('GBM DEGs: A %d, B %d, significant in both %d (%.1f%% of A)',
                length(gA), length(gB), length(intersect(gA, gB)), 100 * length(intersect(gA, gB)) / length(gA)))

# ---- Table S3: RRHO2 quadrant maxima, all analyses ------------------------------------
S2a <- data.frame(analysis = qc$analiz, source = "Stage 3 (corrected boundaries)", quadrant = qc$ceyrek,
                  max_neg_log10_P = qc$max, row = qc$row, column = qc$col, adjusted_P = qc$p_adj,
                  boundary_rows = qc$sinir_satir, boundary_columns = qc$sinir_sutun)
t2 <- as.data.frame(a2$table); t2$cell <- rownames(a2$table)
S2b <- rbind(data.frame(analysis = paste0("Braak II vs 0: ", t2$cell), source = "2x2 decomposition (A2)", quadrant = "DD",
                        max_neg_log10_P = t2$DD_max, top200_overlap = t2$DD_top200_overlap,
                        top200_expected = t2$DD_top200_expected, top200_P_unadjusted = t2$DD_top200_p),
             data.frame(analysis = paste0("Braak II vs 0: ", t2$cell), source = "2x2 decomposition (A2)", quadrant = "DU",
                        max_neg_log10_P = t2$DU_max, top200_overlap = NA, top200_expected = NA, top200_P_unadjusted = NA))
S2c <- do.call(rbind, lapply(names(c(plan_B6 = 1, deviation_B2 = 1)), function(k) {
  x <- a3[[k]]; lab <- if (k == "plan_B6") "Plan: Braak VI vs 0" else "Deviation: Braak II vs 0"
  data.frame(analysis = paste(lab, "(GBM limma-trend)"), source = "GBM sensitivity (A3)",
             quadrant = c("DU", "DD"), max_neg_log10_P = c(x$DU_max, x$DD_max), adjusted_P = c(x$DU_adj_p, NA),
             OR_top200 = x$OR, P_bonferroni = x$p_bonf, P_permutation = x$p_perm, decision = x$decision)
}))

# ---- Table S4: GSE125583 DD-signature and stress genes (existing exploratory table) ----
S3 <- read.csv(here::here("outputs/A4_GSE125583_DD_genes_corrected.csv"), check.names = FALSE)

# ---- Table S2: BRETIGEA astrocyte markers in the GBM cohort ---------------------------
mk  <- BRETIGEA::markers_df_human_brain
ast <- head(mk$markers[mk$cell == "ast"], 50)
sc  <- BRETIGEA::brainCells(inputMat = gb$expr, nMarker = 50, species = "human",
                            celltypes = c("ast", "end", "mic", "neu", "oli", "opc"), scale = TRUE)
ast <- ast[ast %in% rownames(gb$expr)]
tum <- colnames(gb$expr) %in% gb$meta$sample_id[gb$meta$diagnosis == "GBM"]
S4 <- data.frame(marker = ast, rank_in_BRETIGEA = match(ast, mk$markers[mk$cell == "ast"]),
                 r_with_astrocyte_score_all = sapply(ast, function(g) cor(gb$expr[g, ], sc[, "ast"])),
                 r_with_astrocyte_score_tumours = sapply(ast, function(g) cor(gb$expr[g, tum], sc[tum, "ast"])),
                 r_with_astrocyte_score_controls = sapply(ast, function(g) cor(gb$expr[g, !tum], sc[!tum, "ast"])))
message("Check against text (AQP4 -0.22, GFAP -0.27, ETNPPL +0.87):")
print(S4[S4$marker %in% c("AQP4", "GFAP", "ETNPPL"), ], row.names = FALSE, digits = 3)

# ---- write tables ----------------------------------------------------------------------
sheets <- list(S1a_shared_DEGs_primary = S1a, S1b_shared_DEGs_replication = S1b, S1c_newly_shared_ModelB = S1c,
               S2_BRETIGEA_astro_markers = S4,
               S3a_RRHO2_quadrants = S2a, S3b_DD_decomposition = S2b, S3c_GBM_limma = S2c,
               S4_GSE125583_DD_genes = S3)
if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(sheets, file.path(out, "Supplementary_Tables_S1-S4.xlsx"), rowNames = FALSE)
  message("Saved Supplementary_Tables_S1-S4.xlsx")
} else if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(sheets, file.path(out, "Supplementary_Tables_S1-S4.xlsx")); message("Saved xlsx (writexl)")
} else {
  for (n in names(sheets)) write.csv(sheets[[n]], file.path(out, paste0(n, ".csv")), row.names = FALSE)
  message("openxlsx/writexl not installed: saved one CSV per table")
}

# ---- Figure S1: H1 replication ----------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", trim = TRUE)
nd  <- function(tt) sum(tt$adj.P.Val < 0.05)
d <- data.frame(
  cohort = factor(rep(c("Primary AD meta-cohort (n = 51)", "Replication, GSE125583 (n = 195)"), each = 4),
                  levels = c("Primary AD meta-cohort (n = 51)", "Replication, GSE125583 (n = 195)")),
  what  = factor(rep(rep(c("AD DEGs", "Shared DEGs"), each = 2), 2), levels = c("AD DEGs", "Shared DEGs")),
  model = rep(c("A: unadjusted", "B: + BRETIGEA"), 4),
  n = c(nd(b4$models$AD_A$tt), nd(b4$models$AD_B$tt), m4$shared$ALL_fdr$nA, m4$shared$ALL_fdr$nB,
        a1$deg_A, a1$deg_B, a1$shared$fdr$nA, a1$shared$fdr$nB))
d$ret <- ave(d$n, d$cohort, d$what, FUN = function(v) v[2] / v[1])
p <- ggplot(d, aes(model, n, fill = model)) + geom_col(width = 0.65) +
  geom_text(aes(label = fmt(n)), vjust = -0.3, size = 2.8) +
  facet_grid(what ~ cohort, scales = "free_y") +
  geom_text(data = unique(d[, c("cohort", "what", "ret")]), inherit.aes = FALSE,
            aes(x = 1.5, y = Inf, label = sprintf("B/A: %.1f%%", 100 * ret)), vjust = 1.5, size = 2.8) +
  scale_fill_manual(values = c("#999999", "#0072B2"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "Genes at FDR < 0.05") +
  theme_classic(base_size = 9) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(out, "FigureS1_H1_replication.pdf"), p, width = 6.5, height = 4.6, device = cairo_pdf)
ggsave(file.path(out, "FigureS1_H1_replication.png"), p, width = 6.5, height = 4.6, dpi = 300)
message("Saved FigureS1_H1_replication; all files in ", out)
