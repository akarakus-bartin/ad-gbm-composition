# Yarın için 3 iş — 2026-09-21

## 1. Leng metodolojik kararlar (30 dk, kahveden sonra)
- Braak 2 grubu: sadece 0 vs 6 (plan-uyumlu) veya + 0 vs 2 (sensitivity)?
- Sex asimetrisi (hepsi male) → Limitations bölümüne not
- Cell type annotation: Broad SCP metadata var mı? Yoksa Seurat clustering + marker gen
- DEVIATIONS.md'ye Leng kohort keşif kaydı

## 2. Neftel Smart-seq2 pipeline (2-3 saat)
- TPM matrisini yükle (GSM3828672)
- Adult filter (5,745 hücre, 21 tümör)
- UCell ile 8 modül skorla
- Neftel Fig 2 quadrant kuralları ile state ata
- Cycling hücreleri ayır (G1/S veya G2/M > 1.0)
- results/intermediate/neftel_states.rds
- Git commit

## 3. GBM bulk TOIL processing (1-2 saat)
- Xena phenotype yükle
- TCGA-GBM tümörler (~170) + GTEx hippocampus/frontal cortex (~50) filtrele
- TPM matrisini yükle, örnekleri süz
- results/intermediate/gbm_bulk.rds
- Git commit
