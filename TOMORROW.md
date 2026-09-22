# TOMORROW.md — 2026-09-23 Plan (manuscript writing başlangıcı)

**Son commit:** <bugünkü ikinci commit>
**Bugün (2026-09-22) tamamlananlar:**

### Sabah (Stage 3 v2 + AD validation)
- Stage 3 v2 kod hijyeni + pathway enrichment + RRHO2 heatmap
- AD validation: GSE125583 (n=289), dose-response Wilcoxon p=1.3e-5

### Akşam (GBM validation)
- TCGA-GBM (n=157) vs GTEx cortex (n=510)
- **Spearman ρ = 0.713, p = 4.2e-4** — güçlü cross-disease concordance
- **3 kategori discovery:**
  - Neuronal loss shared (6 gen)
  - Stress shared (6 gen)
  - Divergent (8 gen — tumor-specific)

---

## Manuscript için elimizdeki yapı

### Stage 1 (Leng AD snRNA-seq)
- 486 FDR<0.05 vulnerable neuron signature

### Stage 2 (Neftel GBM snRNA-seq)
- 1,576 FDR<0.05 neural-mimicry signature

### Stage 3 (Cross-disease convergence)
- H2 formally REJECTED (1/3 pre-specified test)
- DU (H2a partial): 18 genes, neurogenesis-enriched
- **DD (novel): 20 genes, HIF-1/AMPK/glycolysis KEGG**

### Stage 4 AD validation
- Direction check: 13/20 DOWN (65%)
- **Dose-response: 17/20 more extreme in Braak V-VI (Wilcoxon p=1.3e-5)**
- İki alt-grup: neuronal-loss (14) + glial-stress (6)

### Stage 4 GBM validation
- **Cross-disease Spearman ρ=0.713 (p=4.2e-4)**
- **3 kategori:** Neuronal loss (6) + Stress shared (6) + Divergent (8)

---

## Yarın için öncelikli görevler

### Öncelik 1 — Manuscript writing başlangıcı (2-3 saat)

**Öncelik 1a: Introduction draft (~1 saat)**
- AD-GBM cross-disease context (neuronal-glial biology overlap)
- Cellular reprogramming hypothesis in both diseases
- Existing literature: Neftel 2019, Leng 2021, Mathys 2019
- Gap: cross-disease transcriptomic convergence tested

**Öncelik 1b: Methods draft (~1-1.5 saat)**
- Stage 1-4 methodology özet
- Statistical tests, thresholds
- Reproducibility statements (git repo, pre-registration)

### Öncelik 2 — Ana figure hazırlığı (1-2 saat)

**Figure 5 — Cross-disease scatter plot:**
- x-axis: AD logFC (advanced vs control)
- y-axis: GBM logFC (tumor vs GTEx)
- Points: 20 DD genes, colored by 3 categories
- Diagonal + quadrants
- Spearman ρ annotation

### Öncelik 3 (opsiyonel) — Category-specific enrichment
- 3 kategori için GO/KEGG
- Tümor-specific divergent genlerin pathway'i (Warburg + proliferation confirmation)

---

## Manuscript hedef dergiler (revised)

Şu anki elimizdeki veri seti ile:

**Primary target:**
- **Briefings in Bioinformatics (IF ~9.5)** — methodology + discovery + validation hikayesi
- Kabul olasılığı: %40-55 (Spearman 0.713 güçlü kanıt)

**Secondary target:**
- **Genes & Diseases (IF ~7.1)** — biology-focused reviewers
- **NPJ Genomic Medicine (IF ~6.7)** — Nature portfolio, open access
- Kabul olasılığı: %55-70

**Safety net:**
- **BMC Medical Genomics (IF ~2.8)** — güvenli, kabul olasılığı yüksek
- **Neuroinformatics (IF ~3.9)** — methodology değer görür

---

## Enerji ve tempo

Bugün ~8 saat aktif çalışıldı (sabah 5.5 + akşam 2.5)
Toplam 4 gün: 20+ commit, 3 major aşama, 2 validation cohort

Yarın hafif başla — manuscript writing farklı bir zihniyet gerektirir.
3-4 saat aktif yazım idealdir.

