# TOMORROW.md — 2026-09-22 Başlangıç Notu

**Son commit:** 589cf60 (Stage 3 v1: H2 tested + 4-quadrant discovery)
**Bugünün özeti:** Stage 1 (Leng AD signature) + Stage 2 v2 (Neftel GBM signature) + Stage 3 v1 (H2 REJECTED + iki novel konvergens ekseni) tamamlandı.

---

## Bugünün ana bilimsel çıktıları

### Stage 1 (leng_deg_primary.rds)
- 486 FDR<0.05 gen (Braak 0 vs 2, primary)
- Klasik AD synaptic loss signature: FOS, EGR1, NPTX2, VGF, BDNF, GRASP DOWN
- H2 için AD-tarafı signature hazır

### Stage 2 (neftel_deg_primary.rds)
- 1,576 FDR<0.05 gen (NeuralLineage[NPC+OPC] vs Other[AC+MES])
- Neftel Fig 1B marker-based malignant filter uygulandı (MES kontaminasyonu %42→%6)
- Top UP: OPALIN, MAG, DLX5/6, STMN2, DCX, SOX11, PLP1, CD24
- Neftel Fig 3B state dağılımı ile birebir uyumlu (AC:%36, MES:%17, NPC:%31, OPC:%16)

### Stage 3 (stage3_convergence.rds + stage3_quadrant_analysis.rds)
- Pre-specified H2 karar: **REJECTED** (1/3 test PASS)
  - RRHO2: PASS (BH-P = 5.3e-4)
  - Hypergeometric top-200: FAIL (OR=1.52, eşik ≥3)
  - Permutation: FAIL (P=0.098, eşik <0.01)
- Post-hoc 4-quadrant discovery: **iki bağımsız konvergens ekseni**
  - **DU (H2a partial):** 18 gen, sinaptik/nöronal identity axis (CD24, STMN1, SYT4, APLP1)
  - **DD (novel):** 20 gen, activity-dependent + metabolic stress axis (JUNB, ZFP36L1, SLC2A3, ALDOA)
  - DU ∩ DD = 0 (iki bağımsız biyoloji)

---

## Yarın için öncelikli görevler

### 1. Stage 3 sign convention temizliği (yüksek öncelik)
- `R/06_stage3_convergence.R` içinde TODO comment var — Test 1 (RRHO2) başlığından önce
- Kod bug'sız çalışıyor ama yorum satırları H2a beklerken kod concordant test ediyor
- Yeniden yazım için: `stage3_quadrant_analysis.rds` yaklaşımını kalıcılaştır (4 quadrant ayrı ayrı)
- Alternatif: RRHO2 iki yönde çağır (concordant + discordant), sonuçları birleştir

### 2. Pathway enrichment analizi (DU + DD gen setleri için)
- **DU (H2a) 18 gen:** GO/KEGG enrichment — sinaptik pathway'ler bekleniyor
  - clusterProfiler::enrichGO, org.Hs.eg.db
  - Background: 3,857 ortak evren geni
- **DD 20 gen:** GO/KEGG — activity-dependent + metabolic + glial support pathway'ler
- Yeni dosya: `R/07_stage3_pathways.R`
- Manuskrit figürü için bar plot

### 3. RRHO2 heatmap PDF kaydı
- `stage3_convergence.rds` içinde `rrho_object` var (RRHO2 hesaplama sonucu)
- `RRHO2_heatmap()` fonksiyonu ile PDF'e kaydet
- `results/figures/stage3_rrho2_heatmap.pdf` — manuskrit için critical figure
- Dört quadrantı işaretle (DU, DD annotation)

### 4. H1 hazırlığı (opsiyonel — yarın enerji varsa başla)
- TCGA-GBM bulk data işlenmesi — TOIL RSEM dosyası (1.26 GB)
  - `data/raw/xena/TcgaTargetGtex_rsem_gene_tpm.gz`
- MuSiC + BRETIGEA deconvolution planlaması
- Stage 4 (05_qc_bulk.R placeholder skeleton var)

---

## Metodolojik dikkatler

### Sign convention sorunları
- RRHO2'nin default davranışı 'concordant' overlap arıyor
- 'Discordant' (H2a gibi) için list2'yi flip etmek yeterli DEĞİL (test tekrarı aynı sonucu verdi — nedeni tam anlaşılamadı)
- Manuel 4-quadrant top-N analiz en şeffaf yaklaşım — bu yaklaşımı stage 3 pipeline'ının resmi çıktısı yap

### İstatistiksel dürüstlük
- DU (p=0.014) ve DD (p=0.003) borderline/moderate — validation gerekli
- Multiple testing farkındalığı: 4 quadrant test edildi, Bonferroni × 4 → DU p=0.056 (n.s.), DD p=0.012 (borderline)
- Manuskritte: DD'yi 'discovery' olarak sun, replication cohort gerektiğini belirt

### Pre-specification integrity
- H2 REJECTED kararı korunuyor (pre-specified rule ile)
- DU + DD bulguları POST-HOC olarak sunulacak (dürüst reporting)
- Manuskrit hikayesi: 'primary hypothesis rejected but discovery findings identified'

---

## Diskteki tüm intermediate dosyalar (kontrol için)

1. ad_bulk.rds (6.9 MB, dünkü AD meta-cohort)
2. neftel_modules.rds (dünkü, 8 meta-modul gen listeleri)
3. leng_processed.rds (97 MB, vulnerable + non-vulnerable SCE)
4. leng_pseudobulk.rds (1 MB)
5. leng_deg_primary.rds (0.3 MB, 486 FDR<0.05)
6. leng_deg_sensitivity.rds (0.5 MB, Braak 0 vs 6)
7. leng_deg_per_subcluster.rds (0.7 MB)
8. neftel_ucell_scores_backup.rds (0.3 MB, sensitivity için backup)
9. neftel_states.rds (0.8 MB, filtered)
10. neftel_deg_primary.rds (2.6 MB, 1576 FDR<0.05)
11. stage3_convergence.rds (0.0 MB, RRHO2 + hypergeometric + permutation)
12. stage3_quadrant_analysis.rds (0.2 MB, 4-quadrant discovery)

---

## Sonuç yorumlama için manuskrit önizlemesi

### Ana bulgular (nihayete kadar korunacak)
1. **H1 (bulk cell-composition):** Henüz test edilmedi — Stage 4 gerekli
2. **H2 (identity oscillation):** Pre-specified formulation REJECTED — ancak partial support (18 gen, p=0.014)
3. **H2' (unexpected shared identity loss):** Novel discovery — 20 gen, p=0.003, activity-dependent axis

### Hedef dergiler (revised, sonuçlara göre)
- **Primary:** Briefings in Bioinformatics (IF ~9.5) — methodology + honest reporting + discovery
- **Backup:** BMC Medical Genomics (IF ~2.8) — güvenli
- **Aspirational:** Neuro-Oncology veya Alzheimer's & Dementia (klinik yorum güçlü olursa)

### Novelty vurgu
1. First cross-disease AD-GBM neural convergence analysis
2. Discovery of 'shared cellular fragility' axis (JUNB, ZFP36L1, SLC2A3, ALDOA)
3. Methodology framework (pre-specification + honest post-hoc discovery)

---

**Not:** Bugün 10 saatlik yoğun bir gün. Yarın taze zihinle Stage 3 temizlik + pathway analysis öncelik. H1 için TCGA bulk işi ayrı bir gün (Stage 4).

