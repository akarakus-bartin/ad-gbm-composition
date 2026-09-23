# Methods — Draft Outline

**Manuscript:** Cross-disease transcriptomic convergence of Alzheimer's disease and glioblastoma
**Author:** Ahmet Karakuş, Bartın Üniversitesi
**Language:** Turkish (draft), edilgen ağırlıklı akademik üslup
**Started:** 2026-09-23 afternoon
**Status:** Outline — draft to follow

---

## 2.1 Data sources (Veri kaynakları)

### 2.1.1 Discovery cohorts (single-cell)
- **Leng et al. 2021** — GEO GSE147528 (snRNA-seq)
  - Post-mortem human brain: entorhinal cortex + superior frontal cortex
  - n=42 donors across Braak stages 0-VI
  - Cell type: RORB+ excitatory neurons (Braak stage ≥IV vulnerable)
- **Neftel et al. 2019** — GEO GSE131928 (snRNA-seq + smart-seq2)
  - Adult and pediatric GBM tumors
  - n=28 samples
  - Cell state: NPC-like, OPC-like, AC-like, MES-like

### 2.1.2 Validation cohorts (bulk)
- **GSE125583** — Bulk RNA-seq, fusiform gyrus
  - n=289 (Advanced AD, All AD, Control)
  - Braak stage stratified
- **TCGA-GBM** — recount3 pipeline
  - n=157 primary GBM tumors
- **GTEx cortex** — recount3 pipeline
  - n=510 non-tumor cortical brain samples

---

## 2.2 Stage 1: AD-vulnerable neuron signature (Leng)

- **Cell filter:** RORB+ excitatory neurons (Braak stage ≥IV subset)
- **Pseudo-bulk aggregation:** Sum counts per donor per cell type
- **DE analysis:** edgeR quasi-likelihood F-test
- **Contrast:** braak_groupBraak2 - braak_groupBraak0 (advanced vs early)
- **Threshold:** FDR<0.05 (no logFC cutoff)
- **Output:** 486 DEGs (104 up, 382 down)

---

## 2.3 Stage 2: GBM neural-mimicry signature (Neftel)

- **Cell filter:** Malignant cells only (CopyKAT / inferCNV verified)
- **Pseudo-bulk aggregation:** Sum counts per sample per cell state
- **DE analysis:** edgeR quasi-likelihood F-test
- **Contrast:** NeuralLineage (NPC + OPC-like) vs Other (AC + MES-like)
- **Threshold:** FDR<0.05 (no logFC cutoff)
- **Output:** 1,576 DEGs (765 up, 811 down)

---

## 2.4 Stage 3: Cross-disease convergence testing

### 2.4.1 Shared gene universe
- 3,857 genes intersected between two signatures
- Used as background for all downstream tests

### 2.4.2 Pre-registered decision rule
- H2 acceptance criteria: At least 2 of 3 tests must pass pre-specified thresholds

### 2.4.3 Three independent tests
- **Test 1 — RRHO2:** Rank-rank hypergeometric overlap (Cahill 2018)
  - Signed -log10(P) ranks
  - Bonferroni threshold at max signal
- **Test 2 — Hypergeometric:** Top-N overlap enrichment
  - N=200 top genes per signature
  - Fisher's exact test
- **Test 3 — Permutation:** 10,000 iterations
  - Preserves marginal directional distributions
  - Empirical null OR distribution

### 2.4.4 Post-hoc 4-quadrant analysis
- After H2 rejection, top-N ranked overlap decomposed into 4 quadrants
- UU (both UP), UD (AD UP × GBM DOWN), DU (AD DOWN × GBM UP) [H2a], DD (both DOWN) [novel]
- Hypergeometric enrichment per quadrant
- Bonferroni correction across quadrants

### 2.4.5 Pathway enrichment (clusterProfiler)
- enrichGO (GO BP) + enrichKEGG
- DD signature (20 genes) and DU signature (18 genes)
- Explicit filter (see DEVIATIONS.md v3):
  - Count ≥ 3, p_raw < 0.01 (DD) or < 0.05 (DU)
  - Exclude tissue-irrelevant terms
  - Exclude generic parent terms

---

## 2.5 Stage 4: Independent cohort validation

### 2.5.1 AD validation (GSE125583)
- Log CPM normalization
- Braak-stratified expression: Advanced AD (Braak V-VI) vs Control
- Wilcoxon signed-rank test for dose-response
- Direction concordance analysis (17/20 more extreme in Advanced)

### 2.5.2 GBM validation (TCGA + GTEx via recount3)
- recount3 pipeline for cross-cohort standardization
- limma-voom differential expression
- TCGA-GBM vs GTEx cortex logFC computed

### 2.5.3 Cross-disease convergence quantification
- Spearman rank correlation between AD and GBM logFC (20 DD genes)
- Pearson correlation as secondary
- Three-category gene classification: neuronal-loss shared, stress shared, tumor-divergent

---

## 2.6 Pre-registration + DEVIATIONS transparency

- Analysis plan pre-registered via Bartın University institutional email archive
- Timestamp preceded first data access to Stage 3+
- All post-hoc modifications documented in DEVIATIONS.md
- ~600+ lines of methodology decisions and justifications

---

## 2.7 Software + Reproducibility

- R 4.6.0 on macOS ARM64
- Key packages: edgeR, limma, RRHO2, clusterProfiler, recount3
- Full analysis code + git history: [repo URL — pending submission]
- 5 commits on 2026-09-23 alone (dev disipilini)
