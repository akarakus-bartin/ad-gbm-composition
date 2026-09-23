# Introduction — Draft (v2, refined)

**Manuscript:** Cross-disease transcriptomic convergence of Alzheimer's disease and glioblastoma
**Author:** Ahmet Karakuş, Bartın Üniversitesi
**Language:** Turkish (draft), will translate to English later
**Started:** 2026-09-23
**Status:** COMPLETED — 5/5 paragraphs, v2 refinements applied

---

## ¶1 — Broad context

Alzheimer hastalığı (AH) ve glioblastoma multiforme (GBM), yüzeyde birbirinden çok farklı iki hastalık olarak görünmektedir. İlki, yavaş ilerleyen nörodejenerasyon ve nöron kaybı ile karakterize; ikincisi ise agresif, hızla büyüyen bir beyin tümörüdür. Buna karşın, her iki hastalık da serebral korteksin nöroepitelyal biyolojisini derinden etkiler ve prognozları hâlâ dramatik biçimde olumsuzdur: AH şu an dünyada en yaygın demans biçimi olarak yaklaşık 55 milyon kişiyi etkilerken¹, GBM tanısından sonraki median sağkalım 15 ayı geçmez². Farklı klinik yörüngelere karşın, son yıllarda ortaya çıkan tek hücre transkriptomik çalışmaları, her iki hastalıkta da olgun kortikal nöron kimliğinin bir 'kimlik krizi' yaşadığı yönünde çarpıcı ipuçları sunmuştur³,⁴. Bu paralel gözlem, sistematik olarak sınanmamış bir soruyu ortaya çıkarmaktadır: iki hastalık, kortikal hücre kimliği bozukluğunda ortak transkripsiyonel bir eksen mi paylaşmaktadır?

---

## ¶2 — Existing knowledge

AH nörodejenerasyonunun tek hücre analizi, hastalığın belirli nöron alt tiplerini seçici olarak vurduğunu ortaya koymuştur. Leng ve arkadaşları³ post-mortem entorhinal ve superior frontal korteks örneklerinden elde ettikleri snRNA-seq verisinde, üst-katman RORB+ eksitatör nöronların hastalığın erken evrelerinden itibaren savunmasız olduğunu göstermiştir; bu bulgu Mathys ve arkadaşları⁵ ile Green ve arkadaşları⁶ tarafından bağımsız kohortlarda tekrarlanmıştır. Savunmasız RORB+ nöronlar, NPTX2, VGF, EGR1 ve FOS gibi sinaptik plastisiteye özgü aktivite-bağımlı gen programlarını kaybederken, aynı korteks bölgesindeki diğer nöronlar bu programı korumaktadır — bu seçici savunmasızlık, AH patolojisinin mikroskobik dinamikleri açısından merkezi bir bulgu haline gelmiştir. GBM tarafında ise Neftel ve arkadaşları⁴ 20 pediatrik ve erişkin tümör örneği üzerinde yaptıkları tek hücre analizinde, tek bir tümörün bile içinde dört farklı hücresel durum barındırdığını göstermiştir: nöral-öncü benzeri (NPC-like), oligodendrosit-öncü benzeri (OPC-like), astrositik (AC-like) ve mezenkimal (MES-like). Bu dört durum, tümör-içi plastisite ile birbirine dönüşebilmekte ve normal beyin gelişimi sırasında görülen olgunlaşma programlarını yeniden aktive etmektedir. GBM'in bu 'nöral taklit' bileşeni özellikle NPC-like ve OPC-like durumlarda belirgin olup, malign hücrelere sağladığı proliferasyon ve tedaviye direnç avantajları giderek daha fazla belgelenmektedir⁷. Böylece, AH ve GBM tam olarak zıt yönlerde ilerlemektedir: birinde olgun nöronal kimliğin kaybı, diğerinde immatür öncü kimliğinin geri kazanılması — ancak ortaya çıkan iki farklı sonuç arasında paylaşılan bir transkripsiyonel altyapı olup olmadığı henüz sistematik olarak incelenmemiştir.

---

## ¶3 — Knowledge gap

AH ve kanser arasındaki ilişki uzun süredir dikkat çekmektedir, ancak literatürde bu ilişki büyük ölçüde epidemiyolojik boyutta ele alınmıştır: AH hastalarında çeşitli kanserlerin (özellikle akciğer, meme ve kolon) görülme sıklığının beklenenden düşük olduğu birden çok kohort çalışmasında raporlanmıştır⁸. Bu ters ilişkinin moleküler kökenlerini araştıran çeşitli aday-gen çalışmaları p53, PIN1 ve Wnt yolağı gibi ortak düzenleyicileri öne sürmüş⁹, fakat bu yaklaşımlar hem odaklandıkları gen setleriyle sınırlı kalmış hem de genellikle bulk RNA-seq verisine dayanmıştır. Bulk analizin temel sınırlaması ise, farklı hastalıklarda gözlenen 'ortak' gen ekspresyon değişikliklerinin, gerçek biyolojik yakınsamadan mı yoksa doku düzeyindeki hücre kompozisyonu farklılıklarından mı kaynaklandığını ayırt edememesidir — AH'de nöron kaybı ve reaktif glia artışı, GBM'de tümör hücresi + tümör mikroçevresi karışımı, bulk sinyalleri tipik olarak birbirine karıştırır. Bu sorun, son yıllarda tek hücre transkriptomiği ve büyük ölçekli standardize edilmiş kohort veri kaynakları¹⁰ ile aşılabilir bir zemine kavuşmuştur; buna karşın, AH ve GBM'in tek hücre imzaları arasındaki yakınsamayı önceden belirlenmiş hipotezler altında ve hücre-tipi spesifik çözünürlükte sistematik olarak sınayan bir çalışma bilgimiz dahilinde bulunmamaktadır. Böyle bir analiz, hem paylaşılan transkripsiyonel eksenleri hem de hastalığa özgü ayrışmaları görünür kılabilir ve iki hastalık arasındaki 'zıt-yön ortak biyoloji' hipotezini formal olarak test etme imkânı sunar.

---

## ¶4 — This study

Bu çalışmada, AH ve GBM'in tek hücre transkripsiyonel imzaları arasındaki yakınsama, dört aşamalı bir çerçevede sistematik olarak sınanmıştır. İlk aşamada Leng ve arkadaşlarının³ AH snRNA-seq verisinden Braak V-VI (ileri patoloji) ile Braak 0-II (kontrol/erken) arasında farklı ifade edilen genler çıkarılarak RORB+ eksitatör nöron imzası (486 gen, FDR<0.05) elde edilmiştir. İkinci aşamada Neftel ve arkadaşlarının⁴ GBM snRNA-seq verisinden malign hücre filtresi uygulanarak nöral-taklit durumlar (NPC-like + OPC-like) ile diğer durumlar (AC-like + MES-like) arasındaki karşılaştırma yapılmış ve 1.576 genlik nöral-taklit imzası (FDR<0.05) oluşturulmuştur. Üçüncü aşamada bu iki imza arasındaki çapraz-hastalık yakınsama, üç bağımsız istatistiksel test — RRHO2 rank-rank hipergeometrik örtüşme, hipergeometrik zenginleşme testi ve 10.000 permütasyon — ile önceden belirlenmiş bir karar kuralı çerçevesinde (H2 kabulü için üç testten en az ikisi geçmeli) değerlendirilmiştir. Dördüncü aşamada, üçüncü aşamada ortaya çıkan çekirdek imzalar, tamamen bağımsız iki kohortta doğrulanmıştır: GSE125583 (n=289, fuziform korteks, Braak evresi stratifiye) AH tarafı için¹¹ ve TCGA-GBM (n=157 birincil tümör) ile GTEx korteks (n=510) karşılaştırması GBM tarafı için, hepsi recount3 boru hattı ile standardize edilmiş şekilde¹⁰. Analiz planı, herhangi bir sonuç görülmeden önce Bartın Üniversitesi kurumsal e-posta arşivi üzerinden zaman-damgalı olarak kaydedilmiş; sonrasında ortaya çıkan tüm sapmalar bir DEVIATIONS belgesinde şeffaf biçimde raporlanmıştır, bu sayede post-hoc rasyonelleştirme riski en aza indirilmiştir.

---

## ¶5 — Main findings preview

Elde edilen bulgular, çapraz-hastalık yakınsama hipotezinin ilk formülasyonunun (H2) formal olarak reddedildiğini göstermiştir; buna karşın veriler, önceden öngörülmeyen fakat biyolojik olarak yorumlanabilir bir yakınsama eksenini ortaya koymuştur. Nöral-taklit yönlü klasik yakınsama (AH'de erişkin kimliği kaybı → GBM'de immatürite kazanımı, DU çeyreği) yalnızca sınırda destek almış (18 gen, nörojenez zenginleşmesi p<0.001) ve önceden belirlenmiş Bonferroni eşiğinin altında kalmıştır. Ancak, üç testin ortak sinyal noktası her iki hastalıkta da azalış gösteren bir çekirdek gen kümesi üzerinde yoğunlaşmıştır: DD imzası olarak adlandırılan bu 20 genlik küme (BAIAP2, ALDOA, NRN1, PNMA2, OLFM1, IDS ve diğerleri), RRHO2 heatmap'inde en yoğun örtüşme bölgesini oluşturmuş ve HIF-1, AMPK ve fruktoz-mannoz metabolizması gibi metabolik-stres yolaklarında anlamlı zenginleşme sergilemiştir. Bağımsız GSE125583 kohortunda DD imzasının Braak evresi ile doz-yanıt ilişkisi (Wilcoxon işaretli-sıra p=1.3×10⁻⁵) gösterdiği, TCGA-GBM ile GTEx korteks karşılaştırmasında ise iki hastalığın logFC değerleri arasında çok güçlü sıralama uyumu (Spearman ρ=0.713, p=4.2×10⁻⁴) bulunduğu tespit edilmiştir. Genler kategoriye ayrıldığında üç ayrı biyolojik pattern ortaya çıkmıştır: her iki hastalıkta paylaşılan nöronal kimlik kaybı (6 gen), her ikisinde paylaşılan stres yanıtı aktivasyonu (6 gen: JUNB, CEBPD, ZFP36L1, PER1, FLNA, PFKFB3), ve tümöre özgü ayrışan metabolik-proliferatif regülasyon (8 gen). Bu üç-kategorili çerçeve, AH ve GBM arasındaki paylaşılan biyolojiye yeni bir kavramsal çerçeve — hücresel kırılganlık ekseni — katmakta ve bu yaklaşım için ilk sistematik kanıtları sunmaktadır.

---

## Revision notes (v2, 2026-09-23 afternoon)

- ¶1: Son iki cümle yumuşatıldı; RORB+ / NPC-like detayları ¶2'ye bırakıldı (tekrar önlendi)
- ¶5 açılış: 'cross-disease convergence' → 'çapraz-hastalık yakınsama' (¶3 ile tutarlılık)
- ¶5 son cümle: 'ilk kez formal olarak / şu ana kadar gözlenmemiş' iddialar yumuşatıldı
- Dosya sonu NA + gereksiz tekrar satırları temizlendi
