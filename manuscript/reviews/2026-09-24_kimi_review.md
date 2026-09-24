# External Review — Kimi (received 2026-09-24)

**Received after:** Full manuscript draft (commit 47dcaaa)
**Reviewer:** Kimi (independent AI critique)
**Verdict:** 'Major revision sonrası kabul edilebilir düzeyde. Analitik çekirdek sağlam.'
**Overall tone:** More detailed and more positive than ChatGPT review

---

## Points confirming ChatGPT critique

1. 'İlk sistematik kanıt' claim overstated — MUST revise
2. H2 rejection ↔ DD discovery statistical tension
3. Post-hoc DD framing as exploratory
4. Braak dose-response statistics need reformulation
5. HIF/AMPK claim premature (nominal p only)

---

## NEW critical concerns from Kimi

### 1. TNFRSF12A study (Scientific Reports 2025, PubMed 40595248)

- AD-GBM Mendelian randomization + scRNA-seq + bulk + cell-cell communication + in vitro
- 413 shared tumor suppressor genes identified
- INVERSE regulation pattern reported
- DIRECTLY LIMITS our 'first evidence' claim
- CRITICAL bonus: published in Scientific Reports → editorial fit evidence!

### 2. Mitochondrial signatures study (Frontiers in Aging Neuroscience 2025)

- AD-GBM cross-disease markers + multi-omic
- Second independent 2025 study

### 3. Braak coding inconsistency — CONCRETE ERROR in Methods

- Methods 2.2: 'Braak 0 / Braak 1 / Braak 2' groups (three groups)
- Neftel Braak literal I-VI (six stages)
- Validation: 'Braak 0-II / Tüm AH / Braak V-VI'
- Reader cannot understand what contrasts are being compared
- MUST clarify with literal Braak numbers throughout

### 4. GTEx age confounder

- TCGA-GBM avg age ~60
- GTEx age distribution different
- HIF-1/hypoxia genes could be artefactually affected
- Sensitivity analysis with age covariate needed

### 5. Neftel cohort pediatric+adult mixture

- Pediatric GBM transcriptome differs from adult
- DD signature dependence on this mix unknown
- Sensitivity subset analysis recommended

### 6. Permutation test null model unclear

- What labels shuffled? What marginal distribution preserved?
- For RRHO2 max or top-200 OR?
- MUST specify precisely — determines p=0.31 reliability

### 7. Cell fraction vs expression suppression distinction

- In AD, is DD gene downregulation due to vulnerable neuron loss OR
  suppression in surviving neurons?
- Leng data allows cell fraction analysis
- This affects biological interpretation of 'neuronal identity loss' axis

### 8. Wilcoxon on ranks loses effect size info

- 17/20 |logFC| ranking test uses ranks only
- Effect size magnitudes not captured
- Report effect sizes with signed test

### 9. 3,857 gene wording in Figure 2A

- 3,857 is shared universe (Stage 3), not Stage 1 DE test universe
- Results 3.1 conflates the two
- MUST clarify

---

## Kimi's SUGGESTED journal targets (revised strategy)

### Primary targets

**1. Scientific Reports (Nature Portfolio)**
- Editorial fit CONFIRMED: 2025 TNFRSF12A paper same topic
- APC-based OA (~2000-2500 USD APC — cost concern)
- ÜAK-safe: standard 'Article' label
- Reasonable acceptance rate for major-revision manuscripts

**2. Frontiers in Aging Neuroscience**
- 2025 mitochondrial signatures paper published there
- Full OA (APC-based, ~2,950 USD)
- CAVEAT: same-year competing paper — differentiation must be sharp

**3. BMC Bioinformatics or BMC Genomics**
- Method contribution positioning (pre-registration + multi-test protocol)
- Full OA (~2,690 USD APC)

### Secondary targets (higher IF, more selective)

- Alzheimer's & Dementia: DADM — needs more AD-specific analyses (APOE, second cohort)
- Neuro-Oncology Advances (Oxford)
- npj Systems Biology and Applications

### Third-tier realistic

- Journal of Alzheimer's Disease

### AVOID (Kimi's warning)

- Genome Biology, Nucleic Acids Research — analytical novelty insufficient
- Public data reanalysis alone won't pass their threshold

---

## Kimi's overall assessment quotes

On methodological discipline:
> 'Metodolojik disiplin üst düzeyde. Ön kayıt, 2/3 karar kuralı, DEVIATIONS belgesi — 
> alanda (özellikle TEK YAZARLI, halka açık veri çalışmalarında) çok nadir.'

On pseudo-bulk approach:
> 'Pseudo-bulk + edgeR-QLF tasarımı doğru tercih. Tek hücre analizinde istatistiksel
> olarak savunulabilir altın standart yaklaşım.'

On 3-category framework:
> 'Basit bir ortak gen listesi sunmak yerine yorum katmanı ekliyor;
> bu, makalenin en özgün entelektüel katkısı olmaya aday.'

On validation architecture:
> 'Keşif → doğrulama geçişi, keşif kohortlarının küçük örneklemini telafi ediyor.'

On H2 rejection:
> 'H2 hipotezinin formal olarak reddedilip sonuçta bu reddin ötesinde yeni bir yakınsama
> ekseninin çıkarılması, çalışmanın bilimsel dürüstlük açısından en güçlü yanı.'

---

## Combined ChatGPT + Kimi action plan — priority ordered

### MUST DO before submission

1. Rewrite 'first evidence' claims (Abstract, Intro para 3, para 5, Disc 4.1, Disc 4.4)
2. Add references to Liu 2013 + TNFRSF12A 2025 + Mitochondrial signatures 2025
3. Reframe DD as exploratory finding, not confirmatory validation
4. Reformulate Braak dose-response (ordered variable + covariates, effect sizes)
5. Clarify Braak coding notation (literal stage numbers throughout)
6. Add pathway multi-test correction (q-values in Figure 4 supplementary)
7. Specify permutation null model precisely
8. Fix '3,857 test edilen gen' wording in Results 3.1
9. Soften HIF/AMPK claim + revise title if needed
10. Add cell fraction sensitivity analysis (Leng data)
11. Add GTEx age-covariate sensitivity analysis

### STRONGLY RECOMMENDED

12. Neftel pediatric vs adult sensitivity subset
13. Individual gene logFC + CI table (supplementary)
14. Complete [VERIFY] references in Discussion
15. Abstract 360 → 300 words

### OPTIONAL but value-adding

16. Top-100/300 robustness in 4-quadrant
17. DU quadrant (n=18) validation OR justification

---

## Journal decision update

Based on Kimi's discovery that Scientific Reports published 2025 TNFRSF12A paper
on same topic:

- Editorial fit STRONGLY evidenced
- Competitive positioning: our niche is cell-type + pre-registration + 3-axes framework
- APC cost concern (~2000-2500 USD) — check Bartın R&P agreements
- ÜAK-safe article label

REVISED PRIMARY TARGET CANDIDATE: Scientific Reports (subject to APC cost verification)
