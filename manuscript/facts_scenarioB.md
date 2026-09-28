# Olgu sayfası — Senaryo B (R/10_facts_scenarioB.R ile üretildi; elle düzenlemeyin)
Üretim: 2026-09-28 08:42 | git HEAD: f2ad237

## H1 (birincil; Stage 4)
- AD kohortu: n = 51 (AD 18, Control 33); bölge: hippocampus 51; veri seti: GSE36980 17, GSE48350 34
- GBM kohortu: n = 338 (Control 186, GBM 152); kontrol bölgesi: Brain - Frontal Cortex (Ba9) 102, Brain - Hippocampus 84
- DEG FDR<0,05 — AD: A 4449, B 0, C 0 | GBM: A 14811, B 7787, C 14853
- Paylaşılan DEG: A 1458, B 0, C 0; oran_B 0.000, oran_C 0.000 -> **SUPPORTED**
- Tanı VIF: AD_B 2.2, GBM_B 13.4

## H2 (ikincil; Stage 3) — plan karar kuralı
| Analiz | Test 1 DU | Test 1 DD | Test 2 | Test 3 | Karar |
|---|---|---|---|---|---|
| Orijinal (psödo-replike, yaşsız, Braak II) | DU 4.21 (adj p 0.114) | DD 6.54 (adj p 0.000531) | OR 1.52, Bonf p 0.278 | emp p 0.098 | **REJECTED** |
| Plan: Braak VI vs 0, donör, yaş | DU 1.28 (adj p 1) | DD 2.06 (adj p 1) | OR 1.40, Bonf p 0.458 | emp p 0.157 | **REJECTED** |
| Sapma: Braak II vs 0, donör, yaş | DU 3.37 (adj p 0.796) | DD 1.67 (adj p 1) | OR 1.40, Bonf p 0.458 | emp p 0.157 | **REJECTED** |

## Stage 1
- Braak VI vs 0: donör 3 vs 3, artık df 3, FDR<0,05 = 0, imza = 0 gen
- Braak II vs 0: donör 3 vs 4, artık df 4, FDR<0,05 = 1, imza = 0 gen
- Psödo-replikasyon: donör ICC 0.346; FDR<0,05 orijinal 486, voom+dupCor 27, donör-toplam 0; logFC r = 0.981

## Negatif kontrol (keşifsel)
- Kol 1 (gürültü): AD korunma 0.851 -> power loss does NOT explain collapse
- Kol 2 (rastgele gen PC): AD korunma 0.000, medyan kompozisyon R² 0.501 -> NON-DISCRIMINATING (random genes capture composition)

## GSE125583 (keşifsel)
- Yaş-grup r = 0.071; stres genleri yaş+cinsiyet+kompozisyon ayarlı: 6/6 yukarı, 5/6 FDR<0,05
