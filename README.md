# AD-GBM Cell-Composition-Controlled Re-evaluation

Analytical code accompanying the pre-specified analysis plan for the study:
**"Cell-composition-controlled re-evaluation of transcriptomic convergence between Alzheimer's disease and glioblastoma multiforme: a hypothesis-driven bulk and single-cell integrative analysis"**

Author: Ahmet Karakuş (Bartın University)

---

## What this repository is

This repository implements the four-stage analytical pipeline pre-specified in `docs/Analysis_Plan_v1.docx` (**this file is the reference; do not modify the code in ways that deviate from it without logging the deviation in `DEVIATIONS.md`**).

## Repository structure

```
.
├── README.md                    # this file
├── DEVIATIONS.md                # log of any deviation from the analysis plan
├── LICENSE
├── .gitignore
├── Dockerfile                   # reproducible R environment
├── renv.lock                    # exact R package versions (created on first run)
├── config/
│   └── params.yml               # ALL parameters (edit here, not in scripts)
├── R/
│   ├── 00_setup.R               # loads packages, params, helpers
│   ├── 01_download_data.R       # data acquisition (GEO + UCSC Xena)
│   ├── 02_qc_bulk.R             # bulk QC + normalisation + ComBat
│   ├── 03_qc_singlecell.R       # snRNA-seq + scRNA-seq QC
│   ├── 04_stage1_ad_signature.R # Stage 1: AD-vulnerable neuron signature
│   ├── 05_stage2_gbm_signature.R# Stage 2: GBM neural-mimicry signature
│   ├── 06_stage3_convergence.R  # Stage 3: H2 tests
│   ├── 07_stage4_bulk_reanalysis.R # Stage 4: H1 test
│   ├── 08_sensitivity.R         # all sensitivity analyses (Section 8)
│   └── 09_figures.R             # publication figures
├── data/
│   └── raw/                     # raw downloads (gitignored)
├── results/
│   ├── intermediate/            # analysis outputs (gitignored, Zenodo-deposited)
│   └── logs/                    # per-stage logs and session info
├── manuscript/
│   └── figures/                 # publication-ready figures
└── docs/
    └── Analysis_Plan_v1.docx    # the pre-specified analysis plan
```

## Prerequisites

- **R** ≥ 4.3 (tested with 4.4)
- **Bioconductor** ≥ 3.18
- ~ 50 GB disk space for raw data downloads
- 32 GB RAM recommended for bulk analyses; 64 GB recommended for single-cell

## Setup

### Option 1: Docker (recommended for reproducibility)

```bash
docker build -t adgbm-analysis .
docker run -it --rm -v $(pwd):/workspace adgbm-analysis
```

Inside the container, R is preconfigured with all required packages.

### Option 2: Local R with renv

```r
install.packages("renv")
renv::restore()   # installs exact versions from renv.lock
```

If `renv.lock` does not exist yet, initialise it:

```r
renv::init()
# Then install packages listed in R/00_setup.R
```

## Running the pipeline

Scripts are numbered in execution order. Run them one at a time and inspect outputs before proceeding:

```r
source("R/01_download_data.R")   # Weeks 1-2 (once, may take hours)
source("R/02_qc_bulk.R")         # ~30 min
source("R/03_qc_singlecell.R")   # ~1-2 hours (single-cell QC is slow)
source("R/04_stage1_ad_signature.R")   # ~10 min
source("R/05_stage2_gbm_signature.R")  # ~15 min
source("R/06_stage3_convergence.R")    # ~5 min
source("R/07_stage4_bulk_reanalysis.R")# ~30 min
source("R/08_sensitivity.R")           # varies; several analyses are TODO stubs
source("R/09_figures.R")               # ~5 min
```

Each stage saves intermediate objects to `results/intermediate/` as `.rds` files and produces a per-stage log at `results/logs/<stage>.log`.

## Pre-specified decision rules

The two primary hypotheses have decision thresholds defined in `config/params.yml`. **Do not modify these values without logging the change in `DEVIATIONS.md` with justification.**

### H1 (Stage 4)

- Supported if shared DEG count under cell-composition control ≤ 20% of uncorrected count in BOTH BRETIGEA and bMIND models
- Rejected if ≥ 50% in either method
- Inconclusive otherwise

### H2 (Stage 3)

Three tests must be evaluated:

- Test 1 (RRHO2): BH-adjusted P < 0.001 at concordant maximum
- Test 2 (Hypergeometric top-200): odds ratio ≥ 3 AND Bonferroni P < 0.001
- Test 3 (Permutation, 1000 replicates): empirical P < 0.01

Composite: SUPPORTED if all 3 pass; REJECTED if ≥ 2 fail; INCONCLUSIVE if exactly 1 fails.

## What is a "deviation"

Any of the following requires an entry in `DEVIATIONS.md`:

- Changing any value in `config/params.yml`
- Skipping a step or filter defined in the plan
- Adding an analysis not described in the plan (unless labelled "exploratory" in the manuscript)
- Using different software versions than those recorded in `renv.lock`

## Data sources and expected file structures

| Source | Accession | Where to find | Notes |
|---|---|---|---|
| Leng snRNA-seq | GSE147528 | GEO or Broad Single Cell Portal (SCP1198) | Get from SCP for clean cell annotations |
| Neftel scRNA-seq | GSE131928 | GEO or Broad Single Cell Portal (SCP503) | Get from SCP for pre-computed state assignments |
| AD microarray 1 | GSE48350 | GEO | Raw CEL files |
| AD microarray 2 | GSE36980 | GEO | Raw CEL files |
| TCGA-GBM | via UCSC Xena TOIL | https://xenabrowser.net | Manual download recommended (large files) |
| GTEx brain | via UCSC Xena TOIL | https://xenabrowser.net | Filter to hippocampus + frontal cortex BA9 |
| Allen M1 reference | Allen Brain Atlas | https://portal.brain-map.org/atlases-and-data/rnaseq | Required for bMIND |

## Known TODOs in the code

The code contains explicit `TODO:` markers where local judgment or dataset-specific adaptation is required:

- Metadata parsing for GSE48350 and GSE36980 (column names vary by GEO submission)
- Neftel cellular state column name (differs between GEO and SCP versions)
- bMIND Allen reference preparation (multi-step, dataset-specific)
- CIBERSORTx triangulation in sensitivity analysis (external service)

These are expected and normal for real-data analysis; they are marked so they are not overlooked.

## Timestamping

Per the analysis plan (Section 10.5), the timestamp of the first commit of this repository containing `docs/Analysis_Plan_v1.docx` serves as the pre-specification timestamp. Do not squash or rewrite commit history.

## Contact

Ahmet Karakuş — akarakus@bartin.edu.tr

## License

See `LICENSE`.
