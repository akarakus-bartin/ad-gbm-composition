# DEVIATIONS from the pre-specified Analysis Plan

This file records any deviation from `docs/Analysis_Plan_v1.docx`.

Every deviation entry must include:
- **Date** the deviation was introduced
- **Section of the plan** affected
- **Reason** for the deviation
- **What was changed** (exactly)
- **Impact** on the primary analysis (if any)
- **Git commit hash** where the change was made

Deviations are additive: never edit or delete existing entries; append new ones.

---

## Deviation log

## 2026-XX-XX — MIND (bMIND) paketi MuSiC ile değiştirildi

**Section of the plan:** 6.4 (Stage 4, deconvolution method 2) and Section 5.2 methods list
**Commit:** <bu değişikliğin git commit hash'ini buraya yazın>
**Reason:** MIND paketi ve bağımlılığı BisqueRNA, R 4.6.0 için CRAN, Bioconductor ve GitHub üzerinden erişilebilir değildir. Her iki paket de bakımsız görünmektedir (CRAN'dan kaldırılmışlar).
**Change:** İkinci deconvolution yöntemi olarak bMIND (MIND paketi) yerine MuSiC (Wang et al., 2019, Nat Commun; xuranw/MuSiC) kullanıldı. Her ikisi de referans-tabanlı deconvolution yaklaşımıdır ve bilimsel amaç (BRETIGEA marker-tabanlı yaklaşımına bağımsız bir alternatif sağlamak) korunmuştur.
**Impact on primary analysis:** H1 karar kuralı aynı kalmıştır (iki bağımsız deconvolution yöntemi arasında ortak shared-DEG oranı karşılaştırması). Yöntemler farklı olduğu için sensitivity analizinde CIBERSORTx ile üçüncü bir triangülasyon da yapılacaktır (Section 8.4.a).
**Justification:** Paket erişilebilirliği zorunlu bir kısıttır. MuSiC, referans-tabanlı bulk deconvolution literatüründe en yaygın kullanılan yöntemlerden biridir ve aktif olarak bakılmaktadır (son güncelleme 2024).

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
