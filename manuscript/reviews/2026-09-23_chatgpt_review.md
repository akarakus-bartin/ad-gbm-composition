# External Review — ChatGPT (received 2026-09-23 late evening)

**Received after:** Full manuscript draft completion (commit 47dcaaa)
**Reviewer:** ChatGPT (independent AI critique)
**Verdict:** Major revision required, but publishable potential

---

## Critical concerns (major revision — structural)

### 1. Internal inconsistency: '20-gene shared down' claim vs 3-category finding

- Discovery: 20 genes selected as DD (both directions in cell-type contrasts)
- Bulk validation: 6 down + 6 up + 8 mixed direction
- Title 'shared stress/hypoxia axis' NOT directly supported by validation
- The three axes framework partially rescues this, but title/abstract/discussion
  still frame as 'shared down regulation confirmed' — must be rewritten

### 2. Discovery and validation contrasts measure different biology

- Leng discovery: within-neuron Braak progression effect
- Neftel discovery: between-malignant-state contrast
- GSE125583 validation: Advanced AD vs Control (mixed neuron loss + glia)
- TCGA-GTEx validation: tumor vs normal cortex (bulk tissue, not cell state)
- Cross-contrast comparison amplifies interpretation gap

### 3. Post-hoc DD signature framed as confirmatory

- DD selected AFTER RRHO2 heatmap inspection
- Bonferroni correction across 4 quadrants doesn't convert exploratory to confirmatory
- Framing must be exploratory throughout
- DU BH-P=0.014 vs 'Bonferroni failed' terminology needs table clarification

### 4. Braak dose-response statistics flawed

- 'Advanced AD' ⊂ 'All AD' → paired comparisons are NOT independent
- 17/20 more extreme |logFC| ≠ monotonic dose-response
- Correct approach: ordered Braak (0-6) modeling + covariates (age, sex, tech)
- Report each gene direction + uncertainty separately

### 5. 'First systematic evidence' claim NOT supportable

- AD-GBM transcriptomic comparisons published 2013 and 2017
- 2025 Scientific Reports paper: single-cell AD-GBM + omics
- 2025 additional paper: single-cell + experimental validation
- Our literature review inadequate — 'ilk kanıt' claim must be withdrawn
- Reframe defensible contribution: specific RORB+ vs NPC/OPC-like comparison

---

## Moderate concerns (methodological rewording)

### 6. Three convergence tests NOT statistically independent

- Fisher (top-200) and permutation test same overlap statistic
- RRHO2 uses same signed rankings
- '2/3 independent tests' overstates independence
- Rewrite as 'multiple summaries of same data'
- Pre-registered decision rule can remain, but scope softened

### 7. TCGA-GTEx interpretation overreach

- recount3 reduces alignment/count batch, doesn't eliminate sample-source effects
- Bulk tumor vs normal cortex ≠ malignant cell-state contrast (discovery)
- Spearman ρ=0.713 cannot be called 'confirmation of shared mechanism in malignant cells'
- Suggest: second GBM cohort + malignant-cell-level reanalysis

### 8. HIF-1 / AMPK emphasis premature

- Raw p=0.041, 0.049 nominally significant
- Multi-test correction (FDR/q-values) not shown for pathway results
- Cannot claim 'shared hypoxia mechanism' from two nominal p-values
- Title 'reveals a shared stress/hypoxia axis' must be softened until verified

---

## Minor concerns (clarification)

### 9. Pre-registration timeline inconsistency

- Methods 2.6: 'plan registered before any data observation'
- Methods elsewhere: 'before Stage 3 and Stage 4 access'
- Must reconcile: what was already visible when plan drafted?

### 10. BH-P vs Bonferroni naming inconsistency

- RRHO2 result: '43x43 Bonferroni' vs 'BH-P' — must align terminology
- 4-quadrant analysis: BH-P values must match reported method

---

## Reviewer verdict

- Research question: Interesting and publishable
- Analytic design: Good foundation, contrasts need biological refinement
- Findings reliability: 20-gene validation + pathway interpretation need re-examination
- Originality: Moderate; specific cell-state comparisons distinctive
- Submission maturity: MAJOR REVISION REQUIRED

---

## Response strategy (planned tomorrow with fresh mind)

### To be addressed BEFORE submission:

1. REWRITE Title — remove 'reveals a shared stress/hypoxia axis' → more descriptive
2. REWRITE Abstract — 3-category finding as primary, not 'shared down validation'
3. REWRITE Introduction para 5 — remove 'first systematic evidence'
4. UPDATE Introduction — literature review including 2013, 2017, 2025 papers
5. REWRITE Results 3.2-3.3 — post-hoc DD framing throughout, exploratory language
6. REWRITE Results 3.3 — Braak analysis reformulated (ordered variable, covariates)
7. REWRITE Discussion 4.1 — main findings honest positioning
8. REWRITE Discussion 4.2.2 — HIF/AMPK claim softened, multi-test acknowledged
9. REWRITE Discussion 4.4 — 'first systematic' language removed throughout
10. UPDATE Methods 2.4 — three tests NOT independent statistically
11. UPDATE Methods 2.5 — TCGA-GTEx caveat added
12. UPDATE Methods 2.6 — pre-registration timeline reconciled
13. UPDATE DEVIATIONS.md — this critique's revision commitments

### To be considered:

- Second GBM cohort for validation (CGGA?)
- Ordered Braak modeling in GSE125583 with covariates
- Pathway multi-test correction (q-values) added to Figure 4 supplementary
- Literature review paragraph in Introduction acknowledging prior AD-GBM work

---

## Bilim insanı olarak dürüst self-assessment

This critique is largely correct. Our manuscript oversold the
validation strength and understated the discovery-validation contrast
mismatch. The 3-category framework rescues part of this but the
title/abstract/conclusion still overstated shared-signal claims.

Response: Accept critique, plan systematic revision, address each
point with methodological rigor. This is exactly the kind of external
review a manuscript needs BEFORE submission — much better to fix
now than face desk-reject or reviewer rejection later.
