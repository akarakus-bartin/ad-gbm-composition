Initial commit hash: 3fcf8d1fa2e4fab5dde875ebe54b8a7ef6a1d7ff 2026-09-20 18:14:57 +0300

# DEVIATIONS from the pre-specified Analysis Plan

This file records any deviation from `docs/Analysis_Plan_v1.docx`.

Every deviation entry must include:
- **Date** the deviation was introduced
- **Section of the plan** affected
- **Reason** for the deviation
- **What was changed** (exactly)
- **Impact** on the primary analysis (if any)MIND (bMIND)
- **Git commit hash** where the change was made

Deviations are additive: never edit or delete existing entries; append new ones.

---

## Deviation log

## 2026-09-20 — MIND (bMIND) paketi MuSiC ile değiştirildi

**Section of the plan:** 6.4 (Stage 4, deconvolution method 2) and Section 5.2 methods list
**Commit:** <Adım 1'de aldığınız hash>
**Reason:** MIND paketi ve bağımlılığı BisqueRNA, R 4.6.0 için CRAN, Bioconductor ve GitHub üzerinden erişilebilir değildir. Her iki paket de bakımsız görünmektedir (CRAN'dan kaldırılmışlar).
**Change:** İkinci deconvolution yöntemi olarak bMIND (MIND paketi) yerine MuSiC (Wang et al., 2019, Nat Commun; xuranw/MuSiC) kullanıldı. Her ikisi de referans-tabanlı deconvolution yaklaşımıdır ve bilimsel amaç (BRETIGEA marker-tabanlı yaklaşımına bağımsız bir alternatif sağlamak) korunmuştur.
**Impact on primary analysis:** H1 karar kuralı aynı kalmıştır (iki bağımsız deconvolution yöntemi arasında ortak shared-DEG oranı karşılaştırması). Yöntemler farklı olduğu için sensitivity analizinde CIBERSORTx ile üçüncü bir triangülasyon da yapılacaktır (Section 8.4.a).
**Justification:** Paket erişilebilirliği zorunlu bir kısıttır. MuSiC, referans-tabanlı bulk deconvolution literatüründe en yaygın kullanılan yöntemlerden biridir ve aktif olarak bakılmaktadır (son güncelleme 2024).
---

## 2026-09-20 — GSE48350 diagnosis kodlaması ve Braak-tabanlı AD tanımı

**Section of the plan:** 4.1 (GSE48350 inclusion criteria) and 5.1 (bulk preprocessing)
**Commit:** <Adım 1'de aldığınız hash>
**Reason:** GSE48350'nin metadata'sı doğrudan bir diagnosis alanı içermiyor. Manuel inceleme sonrası veri yapısı şu şekilde tespit edildi:
- 80 örnek Braak stage annotasyonu taşıyor ve title alanında "_AD_" içeriyor (hepsi AD hastası, Braak I-VI dağılımında)
- 140 örnek "C" (Caucasian) koduyla, Braak yok (sağlıklı kontroller, karışık yaş dağılımı 20-99)
- 33 örnek "AA" (African American) koduyla, 30-48 yaş aralığında (analiz için uygun değil — early-onset AD yaş aralığı dışında)

**Change:** İki metodolojik karar alındı:
1. Diagnosis, hazır bir alandan okunmak yerine title/characteristics kombinasyonundan türetildi. Karar mantığı `R/02_qc_bulk.R` içinde `process_gse48350()` fonksiyonundaki yorumlarda ve kodda açıkça belgelendi.
2. AD grubu ileri Braak stage'lerine (V, V-VI, VI; n≈42) daraltıldı. Bu, planın Stage 1 snRNA-seq analizinde kullanılan Braak 0 vs Braak VI karşılaştırmasıyla metodolojik tutarlılık sağlar. Erken/orta Braak stage'leri (I-IV, n≈35) sensitivity analizinde kullanılabilir.

**Impact on primary analysis:** Final örnek sayısı beklenenden az (region + yaş + Braak filtrelemesi sonrası ~60-80 örnek yerine önceki analizde 75). Bu, GSE48350'nin tek başına gücünü azaltır ama GSE36980 meta-cohort'a eklendiğinde toplam n büyür. Analiz planının pre-specified karar kuralları değişmez.

**Justification:** GSE48350'nin karma yapısı (üç ayrı çalışmayı içeren kompozit veri seti) analiz planı yazılırken tam olarak öngörülemedi. Braak-tabanlı AD tanımı Alzheimer araştırma literatüründe standart bir yaklaşımdır (Braak & Braak 1991, *Acta Neuropathologica*).
<!--
Template for new entries:

## YYYY-MM-DD — <short title>

**Section of the plan:** e.g. 6.1 (Stage 1 procedure)
**Commit:** <git commit hash>
**Reason:** <why this deviation was necessary — data quality issue, tool bug, etc.>
**Change:** <exactly what was changed — file path, parameter value, procedure step>
**Impact on primary analysis:** <none | H1 decision affected | H2 decision affected | sensitivity only>
**Justification:** <scientific reasoning>
-->
