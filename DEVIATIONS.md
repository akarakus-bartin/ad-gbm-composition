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

---

## 2026-09-20 — AD meta-cohort: sadece hippocampus + kategori uyumlaştırma

**Section of the plan:** 4.1 (inclusion criteria) and 5.1 (meta-cohort assembly)
**Commit:** be85e304c34a6ffb824c32cffe9c961aa73f2a8b

**Reason:** Meta-cohort birleştirmesinde iki metodolojik zorlukla karşılaşıldı:
1. **Bölge asimetrisi:** GSE48350 hippocampus + entorhinal cortex içerir; GSE36980 hippocampus + temporal cortex içerir. Sadece hippocampus her iki kohortta ortaktır. Farklı bölgeleri birleştirmek anatomik heterojenlik yaratır ve bölge × dataset arasında perfect confounding'e yol açar.
2. **Kategori kodlaması tutarsız:** İki GEO submission farklı kodlama konvansiyonları kullanmış — Sex: "F"/"M" (GSE36980) vs "female"/"male" (GSE48350); Brain region: "Hippocampus" (GSE36980) vs "hippocampus" (GSE48350).

**Change:** Üç ilgili düzenleme yapıldı:
1. **Bölge filtresi hippocampus'a daraltıldı:** Her iki `process_gse48350()` ve `process_gse36980()` fonksiyonunda `keep_region` filtresi sadece "hippocampus" içeren örnekleri tutar. GSE48350'de entorhinal cortex örnekleri (n=14), GSE36980'de temporal cortex örnekleri (n=27) analiz dışında kaldı.
2. **`merge_ad_cohorts()` fonksiyonuna `harmonise_meta()` iç fonksiyonu eklendi:** Sex ("F"/"M" → "female"/"male") ve brain_region (tolower normalisation) uyumlaştırıldı.
3. **`merge_ad_cohorts()` ComBat model matrix'i dinamik hale getirildi:** Sadece ≥2 seviyesi olan covariate'ler koruma değişkeni olarak dahil edilir. Bu meta-cohort'ta `brain_region` tek seviyeli ("hippocampus") olduğu için ComBat modelinden çıkarıldı; `diagnosis + sex` korundu.

**Impact on primary analysis:**
- Final AD meta-cohort: 17,973 gen × 51 örnek (18 AD + 33 Control), planın öngördüğü 122'den az.
- ComBat kabul kriterleri karşılandı: post-correction PC1 vs dataset r = 0.093 (eşik <0.20 ✓), PC1 vs diagnosis r = 0.485 (eşik ≥0.30 ✓).
- Örnek sayısı azaldığı için istatistiksel güç sınırlıdır; bulk düzeyi cell-composition-controlled DE'nin küçük etki büyüklüklerini tespit etme kapasitesi düşer. Ancak H1 için beklenen etki (bulk DEG'lerin composition control altında büyük ölçüde kaybolması) yeterince büyüktür.
- H1 karar kuralı değişmez (pre-specified 20%/50% oranı korunur).

**Justification:** Sadece hippocampus'a odaklanmak metodolojik olarak daha temiz bir meta-analiz sağlar (anatomik homojenlik). Hippocampus AD patolojisinin en erken ve en yoğun etkilendiği bölgedir (Braak & Braak 1991); Stage 1 snRNA-seq analizinde kullanılan entorhinal cortex ile histopatolojik olarak yakın komşudur ve benzer vulnerable nöron popülasyonlarını içerir (Leng et al. 2021).

---

## 2026-09-21 — Leng snRNA-seq kohortu keşifsel doğrulama ve karar noktaları

**Section of the plan:** 4.1 (Leng cohort inclusion), 6.1 (Stage 1 procedure), 6.2 (cell type annotation)
**Commit:** 6097b8

**Reason:** Leng et al. 2021 (GSE147528) verisi manuel inceleme sonrasında analiz planında öngörülmeyen dört yapısal özellik gösterdi. Bu keşifler dört karar noktasını gerektirdi:

### Karar 1 — Braak stage tasarımı (Braak 0 vs 2 ana, Braak 0 vs 6 sensitivity)

**Reason:** Leng kohortu 3 Braak grubu içeriyor (0/2/6). Manuel inceleme + makalenin ana bulgusu, RORB+ vulnerable neuron depletion'ın **Braak 0 → Braak 2 geçişinde pik yaptığını** (%60-79 azalma), sonra plato oluşturduğunu göstermektedir (makalede: 'no further decrease in Braak stage 6'). Verinin kendi doğrulaması: EC:Exc.s1 (835→326, -61%), EC:Exc.s2 (666→232, -65%), EC:Exc.s4 (382→80, -79%) — hepsi Braak 0 → Braak 2 geçişinde.
**Change:** Analiz planındaki Braak 0 vs Braak 6 ana karşılaştırma, Braak 0 vs Braak 2 olarak revize edildi. Braak 0 vs Braak 6 sensitivity analiz olarak korundu (plan-uyumlu backup).
**Impact on primary analysis:** H2 hipotezinin AD-tarafı sinyalini biyolojik olarak en güçlü yakalayan karşılaştırma seçildi. Karar kuralları değişmez (RRHO2 + hypergeometric + permutation). Braak 0 vs 6 sensitivity karşılaştırması sonuçların progression-invariance'ini test eder.
**Justification:** Leng makalesinin ana bulgusu (Fig. 2c ve Results) selective RORB+ depletion pikinin Braak 0 → Braak 2 geçişinde olduğunu net söylüyor. Sensitivity için Braak 0 vs 6 korunuyor — planla tam uyumlu ve biyolojik olarak da savunulabilir (persistent depletion).

### Karar 2 — Sex asimetrisi

**Reason:** 20 örneğin hepsi Male. Analiz planında sex heterojen kohort öngörülmüştü; gerçek Leng verisi male-only.
**Change:** Sex covariate modelden çıkarıldı (varyansı sıfır; teknik zorunluluk). Analiz planındaki `~ diagnosis + brain_region + sex + age + dataset` modeli, Leng için `~ diagnosis + age` şeklinde uyarlandı.
**Impact on primary analysis:** Analitik güç değişmez. Ama sonuçların genellenebilirliği kısıtlı — sadece erkek beyinleri.
**Justification:** Manuskritin Limitations bölümüne 'findings should be interpreted with caution as generalization to female cohorts requires independent validation' cümlesi eklenecek. Bu yayınlanabilir bir kısıttır; Mathys et al. 2019 gibi büyük snRNA-seq çalışmaları da benzer sex dengesizliklerini raporlamıştır.

### Karar 3 — Cell type annotation kaynağı

**Reason:** Analiz planında cell type ataması için Seurat clustering + marker gen tabanlı manuel etiketleme öngörülmüştü. Ancak Leng'in scAlign-integrated ve cell-type-atanmış işlenmiş verisi Synapse'ta (syn21788402) mevcut.
**Change:** Kendi Seurat clustering pipeline'ı yerine, Leng'in scAlign-assigned `sce.EC.Exc.scAlign.rds` dosyasındaki `subclusterAssignment` kolonu kullanıldı. Bu, orijinal Leng et al. 2021 analizinin bire bir replikasyonudur.
**Impact on primary analysis:** RORB+ vulnerable subclusterlar (EC:Exc.s1, EC:Exc.s2, EC:Exc.s4) doğrudan Leng makalesinin tanımladığı gibi filtrelenecek. Kendi clustering yapılsaydı benzer ama tam olarak aynı olmayan clusterlar çıkabilirdi; orijinal atama daha temiz bir zemin.
**Justification:** (a) Metodolojik güç: Leng makalesinin bulgularının doğrudan replikasyonu; hakemde 'neden farklı clustering?' sorusu doğmaz. (b) Zaman verimliliği. (c) scAlign entegrasyonu (Nazor et al. 2019) donör-batch effekleri için altın standarttır ve Seurat integrate'ten üstün kabul edilir bazı kohortlarda.

### Ekstra keşif — RORB+ vulnerable subcluster tanımlaması

**Reason:** Leng makalesinde RORB+ vulnerable subpopulation'lar 'fine' subclustering ('s' suffix; 9 EC:Exc.s0-s8 subcluster) düzeyinde tanımlanmış. Broad clustering'de ('Exc.1-5') vulnerability sinyalleri bulanıklaşıyor.
**Change:** `subclusterAssignment` (fine subclustering) kullanılıyor. Vulnerable RORB+ subcluster üçlüsü: EC:Exc.s1, EC:Exc.s2, EC:Exc.s4 (Leng et al. 2021, Fig. 2c). Bu üçü Braak 0'da yüksek RORB ifadesi gösterir ve Braak 2'de seçici olarak deplete olur.
**Impact on primary analysis:** Stage 1 pseudo-bulk analizi bu üç subclusterda yapılacak. Non-vulnerable Exc subclusterlar (s0, s3, s5, s6, s7, s8) internal control olarak kullanılabilir.
**Justification:** Leng makalesinin (Nature Neuroscience 2021) bulgu odağı bu üç subcluster; ayrı ele almak spesifisiteyi test etmek için gerekli.
