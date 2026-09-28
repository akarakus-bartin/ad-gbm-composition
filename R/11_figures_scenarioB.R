# =============================================================================
# 11_figures_scenarioB.R — Figures 2-5 (Scenario B), generated only from locked result files.
# Figure 1 (study design) needs no data and is produced separately.
# Output: outputs/figures_scenarioB/FigureN_*.pdf and .png
# =============================================================================
for (p in c("ggplot2", "patchwork", "RRHO2", "here")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p, " — install.packages('", p, "')")
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

out <- here::here("outputs/figures_scenarioB")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
f <- function(p) readRDS(here::here(p))

b4   <- f("outputs/stage4_H1_bretigea.rds")
m4   <- f("outputs/stage4_H1_music.rds")
s1   <- f("outputs/stage1_plan_conformant.rds")
st   <- f("outputs/sensitivity_analysis_full_state.rds")
nc   <- f("outputs/stage4_negative_control.rds")
ad   <- f("results/intermediate/ad_bulk.rds")
gb   <- f("results/intermediate/gbm_bulk.rds")
leng <- f("results/intermediate/leng_deg_primary.rds")
neft <- f("results/intermediate/neftel_deg_primary.rds")
PARAMS_step <- 100   # plan 6.3: RRHO2 step size

# Okabe-Ito colour-blind-safe palette
col <- c(A = "#999999", B = "#0072B2", C = "#E69F00", DD = "#D55E00", DU = "#0072B2",
         orijinal = "#999999", donor = "#000000")
th <- theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom", plot.tag = element_text(face = "bold", size = 11),
        plot.title = element_text(size = 9, face = "bold"))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  message("Saved: ", name)
}
ndeg <- function(m) sum(m$tt$adj.P.Val < 0.05)
fmt  <- function(x) format(x, big.mark = ",", trim = TRUE)
dec  <- function(x, d = 2) formatC(x, format = "f", digits = d)

# -----------------------------------------------------------------------------
# FIGURE 2 — H1
# -----------------------------------------------------------------------------
d2a <- data.frame(
  kohort = rep(c("AD (n = 51)", "GBM (n = 338)"), each = 3),
  model  = rep(c("A", "B", "C"), 2),
  n      = c(ndeg(b4$models$AD_A), ndeg(b4$models$AD_B), ndeg(m4$models$AD_C),
             ndeg(b4$models$GBM_A), ndeg(b4$models$GBM_B), ndeg(m4$models$GBM_C)))
p2a <- ggplot(d2a, aes(model, n, fill = model)) +
  geom_col(width = 0.7) + geom_text(aes(label = fmt(n)), vjust = -0.3, size = 2.7) +
  facet_wrap(~ kohort, scales = "free_y") +
  scale_fill_manual(values = col[c("A", "B", "C")], guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Model (A: unadjusted, B: BRETIGEA, C: MuSiC)", y = "DEGs (FDR < 0.05)",
       title = "DEGs per cohort") + th

keys <- c(ALL_fdr = "All controls", ALL_lfc = "All controls\n+ log2FC threshold",
          BA9_fdr = "BA9 only", BA9_lfc = "BA9 only\n+ log2FC threshold",
          HIP_fdr = "Hippocampus only", HIP_lfc = "Hippocampus only\n+ log2FC threshold")
d2b <- do.call(rbind, lapply(names(keys), function(k) {
  s <- m4$shared[[k]]
  data.frame(analiz = keys[[k]], model = c("A", "B", "C"), n = c(s$nA, s$nB, s$nC))
}))
d2b$analiz <- factor(d2b$analiz, levels = unname(keys))
p2b <- ggplot(d2b, aes(analiz, n, fill = model)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = fmt(n)), position = position_dodge(0.8), vjust = -0.3, size = 2.3) +
  scale_fill_manual(values = col[c("A", "B", "C")], name = "Model") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Shared DEGs", title = "Shared DEGs (primary analysis: 'All controls')") +
  th + theme(axis.text.x = element_text(size = 7))

grp_ad  <- as.character(ad$meta$diagnosis[match(rownames(m4$proportions$ad), colnames(ad$expr))])
grp_gbm <- as.character(gb$meta$diagnosis[match(rownames(m4$proportions$gbm), gb$meta$sample_id)])
relab <- c(AD = "AD", Control = "Control", GBM = "GBM tumour")
grp_ad <- relab[grp_ad]; grp_gbm <- relab[grp_gbm]
prop_long <- function(p, g, lab) {
  a <- aggregate(p, list(grup = g), mean)
  do.call(rbind, lapply(2:ncol(a), function(j)
    data.frame(kohort = lab, grup = a$grup, sinif = sub("^mus_", "", colnames(a)[j]), oran = a[[j]])))
}
d2c <- rbind(prop_long(m4$proportions$ad, grp_ad, "AD cohort"),
             prop_long(m4$proportions$gbm, grp_gbm, "GBM cohort"))
d2c$sinif <- factor(d2c$sinif, levels = c("neu", "ast", "mic", "oli", "opc", "end"),
                    labels = c("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Endothelial"))
p2c <- ggplot(d2c, aes(grup, oran, fill = sinif)) +
  geom_col(width = 0.7) + facet_wrap(~ kohort, scales = "free_x") +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"),
                    name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  labs(x = NULL, y = "Mean MuSiC proportion", title = "MuSiC proportions (see text)") +
  th + guides(fill = guide_legend(nrow = 2))

save_fig((p2a | p2c) / p2b + plot_annotation(tag_levels = "A"), "Figure2_H1", 7.2, 6.4)

# -----------------------------------------------------------------------------
# FIGURE 3 — H2: RRHO2 maps (recomputed from the same inputs)
# -----------------------------------------------------------------------------
univ <- intersect(leng$de_table$gene, neft$de_table$gene)
if (length(univ) != 3857) warning("Shared universe has ", length(univ), " genes; Stage 3 used 3,857")
rrho_mat <- function(ad_tab) {
  a <- ad_tab[ad_tab$gene %in% univ, ]; g <- neft$de_table[neft$de_table$gene %in% univ, ]
  cg <- intersect(a$gene, g$gene); a <- a[match(cg, a$gene), ]; g <- g[match(cg, g$gene), ]
  rr <- RRHO2::RRHO2_initialize(
    list1 = data.frame(gene = a$gene, value = sign(a$logFC) * -log10(a$PValue)),
    list2 = data.frame(gene = g$gene, value = sign(g$logFC) * -log10(g$PValue)),
    stepsize = PARAMS_step, labels = c("AD", "GBM"), method = "hyper", log10.ind = TRUE)
  rr$hypermat
}
H <- list("Original implementation\n(pseudoreplicated, no age, Braak II)" = rrho_mat(leng$de_table),
          "Deviation: Braak II vs 0\n(donor level, age)"                  = rrho_mat(s1$primary_B2$tt),
          "Plan: Braak VI vs 0\n(donor level, age)"                       = rrho_mat(s1$primary_B6$tt))
stopifnot(abs(max(H[[1]], na.rm = TRUE) - 6.54) < 0.01)   # reproduction of the original result

# True quadrant boundaries: empty (NA) bands in the RRHO2 map at the sign change of each list
bounds <- function(m) {
  br <- which(rowMeans(is.na(m)) > 0.5); bc <- which(colMeans(is.na(m)) > 0.5)
  stopifnot(length(br) > 0, length(bc) > 0)
  list(r_up = 1:(min(br) - 1), r_dn = (max(br) + 1):nrow(m),
       c_up = 1:(min(bc) - 1), c_dn = (max(bc) + 1):ncol(m), br = range(br), bc = range(bc))
}
qmax <- function(m, r, c) { v <- max(m[r, c], na.rm = TRUE)
  pos <- which(m == v & !is.na(m), arr.ind = TRUE); pos <- pos[pos[, 1] %in% r & pos[, 2] %in% c, , drop = FALSE][1, ]
  c(max = v, row = pos[[1]], col = pos[[2]]) }
quad_tab <- do.call(rbind, lapply(names(H), function(k) {
  m <- H[[k]]; b <- bounds(m)
  q <- rbind(UU = qmax(m, b$r_up, b$c_up), UD = qmax(m, b$r_up, b$c_dn),
             DU = qmax(m, b$r_dn, b$c_up), DD = qmax(m, b$r_dn, b$c_dn))
  data.frame(analiz = gsub("\n", " ", k), ceyrek = rownames(q), q, row.names = NULL,
             p_adj = pmin(1, 10^(-q[, "max"]) * length(m)),
             sinir_satir = paste(b$br, collapse = "-"), sinir_sutun = paste(b$bc, collapse = "-"))
}))
quad_tab$Test1_DU <- ifelse(quad_tab$ceyrek == "DU", ifelse(quad_tab$p_adj < 0.001, "PASS", "FAIL"), "")
print(quad_tab, row.names = FALSE)
saveRDS(quad_tab, here::here("outputs/stage3_quadrants_corrected.rds"))
message("Corrected quadrant table: outputs/stage3_quadrants_corrected.rds")
npix <- length(H[[1]]); thr <- 3 + log10(npix)             # adjusted P < 0.001 threshold on the -log10(raw P) scale
d3 <- do.call(rbind, lapply(names(H), function(k) {
  m <- H[[k]]; e <- expand.grid(i = seq_len(nrow(m)), j = seq_len(ncol(m)))
  cbind(analiz = k, e, v = m[cbind(e$i, e$j)])
}))
d3$analiz <- factor(d3$analiz, levels = names(H))
lines <- do.call(rbind, lapply(names(H), function(k) { b <- bounds(H[[k]])
  data.frame(analiz = k, h = mean(b$br), v = mean(b$bc)) }))
lab <- do.call(rbind, lapply(names(H), function(k) { b <- bounds(H[[k]]); m <- H[[k]]
  data.frame(analiz = k, x = c(mean(b$c_up), mean(b$c_dn), mean(b$c_up), mean(b$c_dn)),
             y = c(mean(b$r_up), mean(b$r_up), mean(b$r_dn), mean(b$r_dn)),
             t = c("UU", "UD", "DU\n(H2)", "DD")) }))
lines$analiz <- factor(lines$analiz, levels = names(H)); lab$analiz <- factor(lab$analiz, levels = names(H))
p3 <- ggplot(d3, aes(j, i, fill = v)) + geom_raster() +
  geom_hline(data = lines, aes(yintercept = h), linetype = 2, linewidth = 0.3, inherit.aes = FALSE) +
  geom_vline(data = lines, aes(xintercept = v), linetype = 2, linewidth = 0.3, inherit.aes = FALSE) +
  geom_text(data = lab, aes(x, y, label = t), inherit.aes = FALSE, size = 2.6, colour = "white", fontface = "bold") +
  facet_wrap(~ analiz, nrow = 1) + scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(name = expression(-log[10](P)), limits = c(0, max(d3$v, na.rm = TRUE)), na.value = "grey85") +
  labs(x = "GBM ranking (up → down)", y = "AD ranking (up → down)",
       caption = paste0("Dashed lines: sign change (quadrant boundary). Pre-specified threshold (DU quadrant): -log10(P) > ", dec(thr))) +
  coord_fixed() + th + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                             legend.position = "right")
save_fig(p3, "Figure3_H2_RRHO2", 7.2, 3.0)

# -----------------------------------------------------------------------------
# FIGURE 4 — Exploratory: sensitivity of the DD signal to corrections
# -----------------------------------------------------------------------------
pr <- st$pseudoreplication
d4a <- data.frame(yontem = factor(c("Original\n(no donor)", "voom +\ndonor block", "Donor-\nsummed"),
                                  levels = c("Original\n(no donor)", "voom +\ndonor block", "Donor-\nsummed")),
                  n = as.numeric(pr$n_sig[c("original", "voom_dupcor", "donor_sum")]))
p4a <- ggplot(d4a, aes(yontem, n)) + geom_col(width = 0.6, fill = "#555555") +
  geom_text(aes(label = n), vjust = -0.3, size = 2.7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Genes at FDR < 0.05", title = paste0("Pseudoreplication\n(within-donor correlation ", dec(pr$icc_donor), ")")) + th

o  <- leng$de_table; tb <- pr$ttB
cg <- intersect(o$gene, tb$gene)
d4b <- data.frame(orig = o$logFC[match(cg, o$gene)], donor = tb$logFC[match(cg, tb$gene)],
                  dd = cg %in% unlist(st$three_categories))
p4b <- ggplot(d4b, aes(orig, donor)) +
  geom_point(data = subset(d4b, !dd), size = 0.4, alpha = 0.3, colour = "#999999") +
  geom_point(data = subset(d4b, dd), size = 1.2, colour = col[["DD"]]) +
  geom_abline(linetype = 2, linewidth = 0.3) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 2.7,
           label = paste0("r = ", dec(cor(d4b$orig, d4b$donor), 3))) +
  labs(x = "log2FC, original", y = "log2FC, donor-summed", title = "Effect estimates are preserved\n(orange: original DD genes)") + th

tp <- st$threshold_profile
tp$analiz <- factor(ifelse(tp$analysis == "original", "original", "donor level"), levels = c("original", "donor level"))
p4c <- ggplot(tp, aes(n, -log10(p_BH), colour = quadrant, linetype = analiz)) +
  geom_line() + geom_point(size = 1.4) +
  geom_hline(yintercept = -log10(0.05), linetype = 3, linewidth = 0.3) +
  scale_colour_manual(values = col[c("DD", "DU")], name = "Quadrant") +
  scale_linetype_manual(values = c(original = 2, `donor level` = 1), name = NULL) +
  scale_x_continuous(breaks = c(100, 200, 300, 500)) +
  labs(x = "Top N genes", y = expression(-log[10](P[BH])), title = "Pre-specified\nthreshold profile") + th

d4d <- subset(quad_tab, ceyrek == "DD")
d4d$analiz <- factor(c("Original", "Deviation\n(Braak II)", "Plan\n(Braak VI)"), levels = c("Original", "Deviation\n(Braak II)", "Plan\n(Braak VI)"))
d4d$v <- d4d$max
p4d <- ggplot(d4d, aes(analiz, v)) + geom_col(width = 0.6, fill = col[["DD"]]) +
  geom_text(aes(label = dec(v)), vjust = -0.3, size = 2.7) +
  geom_hline(yintercept = thr, linetype = 2, linewidth = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = expression("Maximum " * -log[10](P)),
       title = "RRHO2 DD quadrant\n(dashed: adjusted P = 0.001)") + th

save_fig((p4a | p4b) / (p4c | p4d) + plot_annotation(tag_levels = "A"), "Figure4_DD_sensitivity", 7.2, 6.4)

# -----------------------------------------------------------------------------
# FIGURE 5 — H1 negative control
# -----------------------------------------------------------------------------
it <- nc$iterations; nA <- ndeg(b4$models$AD_A)
p5a <- ggplot(it, aes(arm1_ad)) + geom_histogram(bins = 40, fill = "#999999") +
  geom_vline(xintercept = nA, colour = col[["B"]], linewidth = 0.9) +
  geom_vline(xintercept = 0.8 * nA, linetype = 2) +
  labs(x = "AD DEGs", y = "Iterations", title = paste0("Arm 1: noise covariates\n(retention ", dec(nc$retain_arm1), ")")) + th
p5b <- ggplot(it, aes(arm2_ad)) + geom_histogram(bins = 40, fill = "#999999") +
  geom_vline(xintercept = nA, colour = col[["B"]], linewidth = 0.9) +
  labs(x = "AD DEGs", y = "Iterations", title = paste0("Arm 2: random-gene PCs\n(retention ", dec(nc$retain_arm2), ")")) + th
p5c <- ggplot(it, aes(arm2_compR2)) + geom_histogram(bins = 40, fill = "#999999") +
  geom_vline(xintercept = 0.25, linetype = 2) +
  labs(x = expression("Mean " * R^2 * " of BRETIGEA scores"), y = "Iterations",
       title = paste0("Arm 2: composition captured\n(median ", dec(nc$median_compR2), ")")) + th
save_fig((p5a | p5b | p5c) + plot_annotation(tag_levels = "A",
         caption = sprintf("Blue line: AD DEGs under Model A (%s). Dashed lines: pre-specified thresholds.", fmt(nA))),
         "Figure5_negative_control", 7.2, 3.0)

message("All figures: ", out)
