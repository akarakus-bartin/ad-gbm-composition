# AD–GBM transcriptomic overlap under cell-composition adjustment

Code, pre-specified analysis plan, deviation register and locked results for:

> Karakuş A. *Alzheimer's disease–glioblastoma transcriptomic overlap is markedly reduced after adjustment for estimated cell composition: a pre-specified re-evaluation.* Submitted to *Neurobiology of Aging*.

Author: Ahmet Karakuş, Bartın University (akarakus@bartin.edu.tr)

## What is here

| Path | Content |
|---|---|
| `docs/Analysis_Plan_v1.docx` | Pre-specified analysis plan, first committed before any analysis (commit `3fcf8d1`, 20 September 2026) and never modified since |
| `DEVIATIONS.md` | Dated register of every departure from the plan, including errors found in an audit of the original implementation |
| `config/params.yml` | Pre-specified parameters and thresholds |
| `R/` | Analysis scripts (see the map below) |
| `outputs/` | Locked result files (`.rds`), figures, supplementary tables |
| `manuscript/facts_scenarioB.md` | Every number in the manuscript, generated from the locked results by `R/10_facts_scenarioB.R` |

## Script map

Scripts are numbered in the order they were written and run. The original implementation is kept, with its errors, so that the audit can be checked.

**Data preparation**

- `00_setup.R` packages, parameters, helpers
- `01_download_data.R` public data from GEO and UCSC Xena
- `02_qc_bulk.R` bulk QC, normalisation, ComBat (AD meta-cohort)
- `03_qc_singlecell.R` snRNA-seq and scRNA-seq QC

**Original implementation (superseded; see DEVIATIONS.md)**

- `04_stage1_ad_signature.R`, `05_stage2_gbm_signature.R`, `06_stage3_convergence.R`, `07_stage4_bulk_reanalysis.R`, `08_sensitivity.R`, `09_figures.R`
- An audit against the plan found that Stage 1 treated subclusters from the same donor as independent samples and omitted the planned age covariate, and that H1 had not been tested. These scripts are retained for transparency only.

**Plan-conformant re-analysis (used for all hypothesis decisions)**

- `04b_stage1_plan_conformant.R` Stage 1, donor-level, with age
- `06b_stage3_plan_conformant.R` Stage 3 (H2), three pre-specified tests
- `07a_build_gbm_bulk_plan.R` GBM bulk cohort as specified in the plan
- `07b_stage4_H1_bretigea.R` H1 with BRETIGEA composition covariates
- `07c_allen_reference_subset.R` MuSiC reference selection (committed before the expression matrix was read)
- `07d_stage4_H1_music.R` H1 with MuSiC composition covariates
- `07e_stage4_negative_control.R` negative control for H1

**Analyses added after external review (pre-declared before running)**

- `12_build_gse125583_bulk.R` saves the GSE125583 cohort (see limitation below)
- `12a_H1_replication_gse125583.R` replication of H1 in an independent AD cohort
- `12b_DD_decomposition.R` donor × age decomposition of the original convergence signal
- `12c_GBM_limma_sensitivity.R` GBM limma-trend sensitivity analysis

**Outputs**

- `10_facts_scenarioB.R` facts sheet; `11_figures_scenarioB.R` Figures 2–5; `13_supplementary.R` Tables S1–S4 and Figure S1; `14_build_submission.R` submission package

For every hypothesis test, the commit that fixed the analysis decisions precedes the commit that contains the results. The commit history should not be rewritten.

## Pre-specified decision rules

**H1 (Stage 4).** Supported if the number of shared DEGs under composition adjustment is at most 20% of the unadjusted number with both deconvolution methods; rejected if at least 50% with either method; inconclusive otherwise. *Implemented deviation:* MuSiC replaced the planned bMIND (package availability).

**H2 (Stage 3).** Test 1, RRHO2 at the DU position, adjusted P < 0.001 (*implemented deviation:* P multiplied by the number of pixels instead of BH, which is more conservative); Test 2, hypergeometric overlap of the top 200 genes, odds ratio ≥ 3 and Bonferroni P < 0.001; Test 3, permutation (1,000), empirical P < 0.01. Supported if all three pass, rejected if two or more fail, inconclusive if exactly one fails.

Outcome: H1 supported; H2 rejected in every analysis variant.

## Data

All input data are public; raw data are not redistributed here.

| Data | Accession / source | Use |
|---|---|---|
| Leng et al. (2021) snRNA-seq | GEO GSE147528 | Stage 1 |
| Neftel et al. (2019) scRNA-seq | GEO GSE131928 (Smart-seq2) | Stage 2 |
| AD microarrays | GEO GSE48350, GSE36980 | H1, primary AD cohort |
| TCGA-GBM and GTEx brain | UCSC Xena, TOIL recompute | H1, GBM cohort |
| Human M1 snRNA-seq (Bakken et al., 2021) | Allen Brain Map | MuSiC reference |
| Fusiform cortex RNA-seq | GEO GSE125583, via recount3 | H1 replication (exploratory) |

`results/intermediate/` contains the five processed cohort objects read by the figure, facts-sheet and supplementary scripts (`ad_bulk.rds`, `gbm_bulk.rds`, `gse125583_bulk.rds`, `leng_deg_primary.rds`, `neftel_deg_primary.rds`; about 49 MB). Larger upstream objects, such as `leng_processed.rds` (93 MB), are not included and are regenerated from the public data by the data-preparation scripts.

## Software

R 4.6.0 with edgeR 4.10.3, limma 3.68.2, sva 3.60.0, BRETIGEA 1.0.4, MuSiC 1.0.0 and RRHO2 1.0; random seed 20260101. A `Dockerfile` from the planning stage is included but has not been re-tested against these versions.

## Known limitations of this repository

- The download and log-CPM steps for GSE125583 are not yet scripted; `12_build_gse125583_bulk.R` saves the cohort from objects created interactively. The saved cohort is included in `results/intermediate/`.
- The analysis plan was committed before analysis, but the repository was made public only at submission. Commit timestamps were generated locally and are not independently verified.

## License

Code: see `LICENSE`. Please cite the article if you use this material.

