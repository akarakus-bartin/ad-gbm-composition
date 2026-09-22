# TOMORROW.md — 2026-09-23 Plan (Stage 4 devam: GBM validation)

**Son commit:** <bugünkü commit>
**Bugün (2026-09-22) tamamlananlar:**
- Stage 3 v2 kod hijyeni + pathway enrichment + RRHO2 heatmap ✓
- Synapse Certified User (yedek, kullanılmadı) ✓
- recount3 pipeline kurulu ✓
- Stage 4 AD validation TAMAMLANDI (GSE125583 / SRP181886, n=289)
- **Dose-response bulgusu: 17/20 gen daha ekstrem in Advanced AD, Wilcoxon p=1.3e-5**
- İki alt-grup keşfi: Neuronal-loss (14) + Glial-stress (6)

---

## Yarın için ana görev — Stage 4 GBM validation

### Ana hipotez
DD signature GBM'de de reproducible mi? Beklenen:
- Group A (nöronal) → GBM'de DOWN (Neftel Stage 2 ile uyumlu)
- Group B (glial-stress) → GBM'de UP? (reactive glia-like ama tumor context)

### Öncelikli veri kaynakları

1. **recount3 içindeki GBM cohortları (öncelik)** — açık, standardize, ideal
   - TCGA-GBM recount3'te var: `available_projects()` kontrol
   - CGGA: recount3 içinde olmayabilir, alternatif indirme gerek

2. **GSE108474 (REMBRANDT)** — GEO açık, ~500 sample, brain cancer + normal
   - Yedek seçenek

### Metodoloji
AD tarafındaki analiz aynen tekrarlanacak:
- 20 DD gen'i çek (aynı ENSEMBL ID'ler)
- CPM normalize
- GBM tumor vs normal brain karşılaştırması
- Direction check + dose-response (WHO grade ile stratify)
- Group A / Group B ayrı analiz

### Beklenen süreçler
- Veri indirme + preprocessing: ~1 saat
- Direction analiz: ~30 dk
- Dose-response (WHO grade II/III/IV): ~30 dk
- Kayıt + grafik: ~30 dk
- **Toplam: 2-3 saat**

---

## AD validation bulgusu özeti (referans)

### Sample özeti
- GSE125583 / SRP181886, fusiform gyrus, RNA-seq (recount3)
- 219 AD + 70 control
- Braak stratification: I-III (early), V-VI (advanced)

### Ana istatistikler
- Direction check: 13/20 DOWN (65%), binomial p=0.13 (nominal)
- Advanced AD: 14/20 DOWN (70%), binomial p=0.058 (borderline)
- **Dose-response: 17/20 gen daha ekstrem in Advanced AD**
  - Binomial p = 0.0013
  - Wilcoxon signed-rank p = 1.3e-05
- Group A (neuronal, 14 gen): +17.6% magnitude
- Group B (glial-stress, 6 gen): +25.4% magnitude

### İki alt-grup
- **Group A (neuronal-loss, 14 gen DOWN):** NRN1, PNMA2, OLFM1, IDS, ALDOA, HRH1, BAIAP2, ATP1B2, UBC, PPP2CB, PEA15, SLC2A3, CANX
- **Group B (glial-stress, 6 gen UP):** JUNB, CEBPD, ZFP36L1, PER1, FLNA, PFKFB3

---

## Manuskript güncelleme — güç kazanımı

### Önceki tahmin (Stage 3 sonrası)
IF 5-10 hedef

### Şimdiki tahmin (Stage 4 AD validation sonrası)
IF 8-12 hedef — çünkü:
- Discovery (Stage 3) + Validation (Stage 4 AD) mevcut
- Wilcoxon p=1.3e-5 hakem gözünde güçlü kanıt
- Cell-composition confounding indirect kanıt (H1 dolaylı destek)
- İki alt-grup keşfi manuscript'i daha derin yapıyor

GBM validation eklendiğinde IF potansiyeli daha da artabilir.

---

## Enerji notu

Bugün ~5.5 saat aktif çalışıldı — yoğun ama verimli.
Yarın 3-4 saat yeterli (GBM validation daha rutin, AD paradigm'i tekrarla).

