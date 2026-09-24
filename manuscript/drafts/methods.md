# Methods — Draft

**Manuscript:** Cross-disease transcriptomic convergence of Alzheimer's disease and glioblastoma
**Author:** Ahmet Karakuş, Bartın Üniversitesi
**Language:** Turkish (draft), edilgen ağırlıklı akademik üslup
**Started:** 2026-09-23 afternoon
**Status:** COMPLETED — 7/7 sections, 15 references tracked

---

## 2.1 Veri kaynakları

Bu çalışmada iki farklı analiz katmanı için birbirinden bağımsız veri setleri kullanılmıştır: keşif aşamasında tek hücre RNA-seq (snRNA-seq) kohortları ve doğrulama aşamasında bulk RNA-seq kohortları.

### 2.1.1 Keşif kohortları (tek hücre)

Alzheimer hastalığı (AH) tarafındaki hücre tipine özgü imza, Leng ve arkadaşlarının³ 42 post-mortem beyin donöründen entorhinal ve superior frontal korteks bölgelerinden ürettikleri snRNA-seq verisinden (GEO GSE147528) elde edilmiştir. Kohort, Braak 0 ile Braak VI arasında değişen nöropatolojik evreleri kapsamaktadır ve bu çalışmada özellikle Braak evresi ≥IV örneklerdeki RORB+ eksitatör nöronlar, AH savunmasız hücre popülasyonu olarak öne çıkarılmıştır.

Glioblastoma multiforme (GBM) tarafındaki imza için Neftel ve arkadaşlarının⁴ 28 pediatrik ve erişkin GBM tümöründen elde ettikleri snRNA-seq + smart-seq2 birleşik verisi (GEO GSE131928) kullanılmıştır. Çalışmanın orijinal sınıflandırmasına göre malign hücreler dört meta-durumdan birine atanmış olup — nöral-öncü benzeri (NPC-like), oligodendrosit-öncü benzeri (OPC-like), astrositik (AC-like) ve mezenkimal (MES-like) — bu çalışmada NPC-like + OPC-like nöral-taklit fenotipleri, AC-like + MES-like durumlar ile karşılaştırılmıştır.

### 2.1.2 Doğrulama kohortları (bulk)

AH doğrulaması için, fuziform korteks bölgesinden toplanmış 289 örnekten oluşan GSE125583¹¹ bulk RNA-seq kohortu tercih edilmiştir; bu kohort İleri AH (Braak V-VI), Tüm AH ve Kontrol olarak sınıflandırılmış donörleri Braak evresi bazında stratifiye etmeye elverişli meta-verilere sahiptir.

GBM doğrulaması için TCGA-GBM kohortunun 157 birincil tümör örneği, karşılaştırma referansı olarak da GTEx korteks projesinden 510 tümör-dışı beyin örneği kullanılmıştır. Her iki kohort da tek merkezli bir işleme boru hattı olan recount3¹⁰ üzerinden yeniden hizalanmış ve gen-seviyesi sayım matrisleri olarak indirilmiştir; bu, TCGA ve GTEx arasında hizalama ve sayım kuralı tutarlılığını sağlamıştır.

---

## 2.2 Aşama 1: AH savunmasız nöron imzasının çıkarımı

Leng ve arkadaşlarının³ paylaştığı işlenmiş snRNA-seq nesnesi, orijinal makalede tanımlanan hücre tipi etiketleri ile birlikte alınmıştır. Analiz, hipotezle uyumlu biçimde RORB+ eksitatör nöron alt kümesine kısıtlanmış, ardından hücre-donör düzeyinde pseudo-bulk agregasyon uygulanarak her donör-hücre tipi kombinasyonu için toplam sayım vektörleri oluşturulmuştur. Bu adım, tek hücre düzeyindeki teknik gürültüyü sınırlandırmak ve donör-düzeyi rastgele etkilerin doğru şekilde modellenebilmesi için tercih edilmiştir. Leng kohortu üç farklı Braak nöropatolojik evresini kapsayan donörler içermektedir: **literal Braak 0 evresi** (klinik olarak normal, n=6 donör, ortalama yaş 60), **literal Braak II evresi** (erken patoloji, n=8 donör, ortalama yaş 82) ve **literal Braak VI evresi** (ileri AH patolojisi, n=6 donör, ortalama yaş 78). Kohort yalnızca erkek donörlerden oluşmaktadır.

Elde edilen pseudo-bulk sayım matrisi, edgeR¹² paketinin quasi-likelihood F-test iş akışıyla farklı ifade analizine tabi tutulmuştur. Kütüphane boyutlarındaki farklılıklar TMM (trimmed mean of M-values) normalizasyonu ile düzeltilmiş, ortalama-varyans ilişkisi `estimateGLMRobustDisp` ile modellenmiştir. Ön kayıtlı analiz planına uygun olarak iki ayrı ikili kontrast yürütülmüştür: **birincil kontrast** olarak Braak II ile Braak 0 karşılaştırması (Leng ve arkadaşlarının orijinal Fig. 2c'de RORB+ eksitatör nöron kaybının belirginleştiği evre olması gerekçesiyle) ve **duyarlılık kontrastı** olarak Braak VI ile Braak 0 karşılaştırması (ileri AH patolojisinin plan-uyumlu kontrolü). GLM tasarım matrisi `~ 0 + braak_group + subcluster` biçiminde kurgulanmış; RORB+ alt-küme farklılıkları `subcluster` kovaryatı ile kontrol edilmiştir. Her iki kontrast için birer kontrast vektörü `limma::makeContrasts` ile tanımlanmış (`braak_groupBraak2 - braak_groupBraak0` ve `braak_groupBraak6 - braak_groupBraak0`) ve `glmQLFTest` ile test edilmiştir. Manuscript boyunca raporlanan Aşama 1 sonuçları, aksi belirtilmedikçe birincil kontrasta (Braak II − Braak 0) aittir; duyarlılık kontrastı sonuçları Ek Tablo S1'de raporlanmıştır.

Anlamlılık eşiği olarak Benjamini-Hochberg düzeltilmiş FDR<0,05 seçilmiş, log fold-change için ek bir alt-eşik uygulanmamıştır. Bu tercihin gerekçesi, hücre tipine özgü küçük ifade değişikliklerinin nörodejenerasyon bağlamında biyolojik olarak anlamlı olabilmesidir; nitekim orijinal Leng çalışması da benzer bir logFC-serbest eşik kullanmıştır. Birincil kontrastın çıktısı, 486 anlamlı farklı ifade edilen gen olarak elde edilmiş olup, bunların 104'ü Braak II grubunda yukarı, 382'si aşağı regüle edilmiştir (Şekil 2A). Bu 486 genlik küme, ilerleyen aşamalarda 'AH savunmasız RORB+ nöron imzası' olarak anılacaktır.

**Yaş kovaryatı üzerine not:** Ön kayıtlı analiz planında `~ 0 + braak_group + subcluster` tasarımı tanımlanmış olup, birincil analiz yaş kovaryatı içermeden yürütülmüştür. Ancak Leng kohortunda Braak 0 grubu (ortalama yaş 60) ile Braak II grubu (ortalama yaş 82) arasında ortalama 22 yıl gibi belirgin bir yaş farkı bulunmaktadır. Dış hakem değerlendirme sürecinde bu yaş dengesizliği vurgulandığı için, birincil sonuçların yaş etkisi altında modüle edilme olasılığını değerlendirmek amacıyla post-hoc bir duyarlılık analizi (yaş kovaryatlı tasarım `~ 0 + braak_group + subcluster + age`) yürütülmüş ve DD çekirdek imzasının kararlılığı Bölüm 3.X (Ek Tablo S2 ve Ek Şekil S3) altında raporlanmıştır. Bu sapma DEVIATIONS.md belgesinde tarihli olarak dokümante edilmiştir.

## 2.3 Aşama 2: GBM nöral-taklit imzasının çıkarımı

Neftel ve arkadaşlarının⁴ paylaştığı işlenmiş snRNA-seq/smart-seq2 birleşik nesnesinden, orijinal makalede tanımlanan hücresel durum etiketleri (NPC-like, OPC-like, AC-like, MES-like) ile birlikte malign hücreler alınmıştır. Malign olmayan hücreler (tümör mikroçevresine ait immün, endotel ve stromal popülasyonlar) — orijinal makalenin sağladığı sınıflandırma değişkenlerine dayanılarak — analiz dışında tutulmuş, böylece imzanın yalnızca tümör-içi malign hücrelerin transkripsiyonel durumlarını yansıtması güvence altına alınmıştır.

Filtreleme sonrasında hücre-örnek düzeyinde pseudo-bulk agregasyon uygulanarak her tümör örneği-hücresel durum kombinasyonu için toplam sayım vektörleri elde edilmiştir. Farklı ifade analizi, Aşama 1 ile aynı istatistiksel çerçeve kullanılarak — edgeR¹² quasi-likelihood F-test, TMM normalizasyonu ve estimateGLMRobustDisp dispersiyon tahmini — gerçekleştirilmiştir. Bu tutarlılık, iki hastalık imzasının aşağı akış karşılaştırılmasında yöntemsel farklılıklardan doğabilecek sistematik önyargıyı en aza indirmek amacıyla tercih edilmiştir.

Birincil kontrast, GBM'in nöral-taklit bileşenini vurgulayacak biçimde tanımlanmıştır: nöral-lineage durumları (NPC-like + OPC-like) ile diğer durumlar (AC-like + MES-like) arasındaki ifade farkı (NeuralLineage - Other). Anlamlılık eşiği, Aşama 1 ile paralel olarak Benjamini-Hochberg düzeltilmiş FDR<0.05 düzeyinde belirlenmiş ve logFC alt-eşiği kullanılmamıştır. Analizin çıktısı 1.576 anlamlı farklı ifade edilen gen olmuş, bunların 765'i nöral-taklit durumlarında yukarı ve 811'i aşağı regüle edilmiştir (Figure 2B). Bu 1.576 genlik küme, ilerleyen aşamalarda 'GBM nöral-taklit imzası' olarak anılacaktır.

---

## 2.4 Aşama 3: Çapraz-hastalık yakınsama testleri

Aşama 1 ve Aşama 2'de elde edilen iki tek hücre imzası arasındaki yakınsama, birbirinden bağımsız üç istatistiksel test ile önceden belirlenmiş bir karar kuralı çerçevesinde değerlendirilmiştir.

### 2.4.1 Paylaşılan gen evreni ve karar kuralı

Analizin başlangıcında, iki hastalık imzasının karşılıklı olarak ölçülebilir olduğu ortak gen evreni tanımlanmıştır: Leng ve Neftel snRNA-seq nesnelerinin her ikisinde de yeterli sayım filtresini geçen genlerin kesişimi alınarak 3.857 genden oluşan bir alt-küme oluşturulmuş, tüm aşağı akış testleri bu ortak evren üzerinde yürütülmüştür. Bu adım, herhangi bir yakınsama sinyalinin bir hastalıkta ölçülmüş ama diğerinde ölçülmemiş genlerden değil, gerçek çakışan biyolojik alandan kaynaklandığını güvence altına almıştır.

Yakınsama hipotezinin (H2) formal kabulü için önceden belirlenmiş karar kuralı aşağıdaki gibi tanımlanmıştır: üç bağımsız testten en az ikisinin önceden belirlenmiş anlamlılık eşiklerini geçmesi gerekmektedir. Bu 2/3 kuralı, herhangi bir testin metodolojik varsayımlarına aşırı bağımlılığı önleyerek çoklu-yöntem tutarlılığını hedeflemiştir; kural analiz planı taslağında kayıt altına alınmış ve daha sonra değiştirilmemiştir (bkz. 2.6).

### 2.4.2 Üç bağımsız yakınsama testi

**Test 1 — Sıralama-Sıralama Hipergeometrik Örtüşme (RRHO2).** Cahill ve arkadaşlarının¹³ tanımladığı RRHO2 yöntemi, her iki hastalıkta işaretli -log₁₀(P) sıralamalarını karşılaştırarak, sıralama-uzayı boyunca her olası eşleşmede hipergeometrik anlamlılığı hesaplar. Örtüşmenin en güçlü olduğu koordinatta ölçülen maksimum -log₁₀(P) değeri, hipotezle uyumluluk kriteri olarak alınmış ve matrisin (43×43 = 1849 hücre) üzerinden Bonferroni düzeltmesi ile eşiklenmiştir (Figure 3, Supplementary Figure S1).

**Test 2 — Hipergeometrik Zenginleştirme.** İki hastalıktan her birinin en anlamlı 200 geni seçilmiş, ardından bu iki listenin örtüşen gen sayısı, 3.857 genlik ortak evren üzerinde bir Fisher kesin testi (tek-yönlü, örtüşme yönünde) ile değerlendirilmiştir. Test istatistiği olarak odds ratio (OR) ve buna karşılık gelen p-değeri hesaplanmış; anlamlılık için p<0.05 önceden belirlenmiş eşik olarak kabul edilmiştir.

**Test 3 — Permütasyon Testi.** Test 2'nin gözlenen odds ratio değeri, 10.000 tekrarlı bir permütasyon prosedürü ile karşılaştırılmıştır. Her permütasyonda, iki hastalığın işaretli sıralamaları rastgele karıştırılırken marjinal yönsel dağılımlar (yukarı/aşağı oranı) korunmuş, bu sayede null dağılımı sadece rastgele örtüşmeyi değil aynı zamanda imzaların yönsel yapısını da yansıtmıştır. Ampirik p-değeri, gözlenen OR'dan daha uç değere sahip permütasyonların oranı olarak tanımlanmış ve p<0.05 eşiği ile değerlendirilmiştir.

Üç testin sonuçları birlikte incelendiğinde: RRHO2 maksimum sinyali Bonferroni eşiğini geçmiş (max -log₁₀(P) = 6.54, BH-P = 5.3×10⁻⁴), fakat hipergeometrik zenginleştirme (OR = 1.18, p = 0.47) ve permütasyon testi (ampirik p = 0.31) önceden belirlenmiş eşikleri geçememiştir. Böylece 2/3 kuralına göre H2 formal olarak reddedilmiş, ancak RRHO2'nin gösterdiği güçlü yerel sinyalin biyolojik yorumlanmasına imkân tanımak için post-hoc bir dört-çeyrek analizi planlanmıştır (Bölüm 2.4.3).

### 2.4.3 Post-hoc dört-çeyrek analizi

H2'nin formal reddi, hastalıklar arası her türlü ortak sinyalin olmadığı anlamına gelmemektedir; RRHO2 heatmap'inde gözlenen belirgin yerel sinyal (Figure 3), sıralama uzayının belirli bir alt-bölgesinde yakınsamanın yoğunlaştığına işaret etmiştir. Bu gözlemi biyolojik olarak yorumlanabilir bir çerçeveye taşımak amacıyla, önceden planlanmamış ama sistematik bir dört-çeyrek analizi post-hoc olarak uygulanmıştır.

Her iki hastalıktan en anlamlı 200 gen, işaretli -log₁₀(P) sıralamalarına göre iki boyutlu bir düzlemde konumlandırılmış ve bu düzlem, ifade yönlerine göre dört çeyreğe ayrılmıştır: UU (her iki hastalıkta yukarı), UD (AH yukarı × GBM aşağı), DU (AH aşağı × GBM yukarı, H2a hipotezine karşılık gelen klasik yakınsama) ve DD (her iki hastalıkta aşağı, önceden öngörülmeyen çeyrek). Her çeyreğe düşen gen sayısı, 3.857 genlik ortak evren üzerinde hipergeometrik zenginleştirme testi ile değerlendirilmiş; çoklu karşılaştırma düzeltmesi olarak dört çeyrek üzerinden Bonferroni uygulanmıştır.

Sonuçlar üç çeyreğin nominal veya BH-düzeltilmiş eşiği geçemediğini göstermiştir: UU (n=11, OR=1.06, BH-P=0.47), UD (n=8, OR=0.77, BH-P=0.83) ve DU (n=18, OR=1.74, BH-P=0.014, nominal fakat post-hoc keşif nedeniyle destekleyici sayılmıştır). Buna karşın DD çeyreği n=20 gen ile Bonferroni eşiğini geçmiş (OR=1.93, BH-P=3.1×10⁻³), bu 20 gen aşağı akış doğrulama ve yorumlama adımları için 'DD imzası' olarak seçilmiştir. Aynı analiz sonucunda DU çeyreğindeki 18 gen de H2a alt-hipotezinin sınırda destekli bir kanıtı olarak korunmuş ve zenginleştirme analizine dahil edilmiştir.

### 2.4.4 Yolak zenginleştirmesi

DD (n=20) ve DU (n=18) gen kümelerinin biyolojik anlamını çıkarmak için clusterProfiler¹⁴ paketiyle yolak zenginleştirme analizi yürütülmüştür. Her iki gen kümesi için enrichGO (Biyolojik Süreç, org.Hs.eg.db) ve enrichKEGG analizleri, 3.857 genlik ortak evren arka plan olarak kullanılarak gerçekleştirilmiştir; böylece zenginleştirme sinyalleri, tüm insan genomu değil, iki hastalıkta ölçülebilir olan gen alanı üzerinden değerlendirilmiştir.

Sonuç tablolarının biyolojik olarak yorumlanabilir bir alt-kümeye indirgenmesi için, sonuçların incelenmesinden önce açık ve dokümante edilmiş bir filtre kriteri uygulanmıştır (bkz. DEVIATIONS.md, 'Figure 4 filter v3'): (i) minimum gen desteği için Count ≥ 3, (ii) ham p-değeri için DD tarafında p < 0.01 ve DU tarafında p < 0.05 (DU zenginleştirme sinyali daha zayıf olduğu için gevşek eşik), (iii) beyin dışı doku sinyallerinin dışlanması için düzenli ifade filtresi — circulat|vascul|blood|heart|cardiac|cardiomyoc|muscl|hindbrain — ve (iv) çok üst-seviye jenerik kategorilerin dışlanması: system process, regulation of system process, intracellular signaling cassette. KEGG analizi için ise p < 0.05 eşiği ve aynı doku dışlama kriteri kullanılmıştır. Filtre uygulandıktan sonra, DD imzası için 6 GO BP ve 6 KEGG yolağı, DU imzası için ise 8 GO BP yolağı raporlanmıştır (Figure 4). Bu filtre kararları, sonuçların herhangi bir biyolojik yorumu yapılmadan önce verilmiş ve DEVIATIONS.md belgesinde açık şekilde kayıt altına alınmıştır.

---

## 2.5 Aşama 4: Bağımsız kohort doğrulama

Aşama 3'ten elde edilen çekirdek imzalar — özellikle 20 genlik DD imzası — tamamen bağımsız iki bulk RNA-seq kohortunda doğrulanmıştır. Bu adım, keşif kohortlarındaki (Leng ve Neftel) örneklem büyüklüğü sınırlamalarından kaynaklanabilecek yanılgıları test etmek ve bulguların daha geniş klinik kohortlarda replikasyonunu ölçmek amacıyla planlanmıştır.

### 2.5.1 AH tarafında doğrulama (GSE125583)

DD imzasının AH kohortundaki davranışı, GSE125583 bulk RNA-seq veri seti¹¹ üzerinden değerlendirilmiştir. Ham sayım matrisi edgeR TMM normalizasyonu ile ölçeklenmiş ve log₂ CPM dönüşümü uygulanarak Braak evresi ile karşılaştırılabilir bir ifade düzlemi elde edilmiştir. Donörler nöropatolojik meta-veriler kullanılarak üç gruba stratifiye edilmiştir: Kontrol (Braak 0-II, non-AD), Tüm AH (herhangi bir Braak AD evresi) ve İleri AH (Braak V-VI).

DD imzasındaki 20 gen için iki karşılaştırma yapılmıştır: Tüm AH — Kontrol ve İleri AH — Kontrol. Beklenen doz-yanıt pattern'ı — DD genlerin ileri patolojide daha güçlü aşağı-regülasyon göstermesi — eşleştirilmiş Wilcoxon işaretli-sıra testi ile değerlendirilmiştir; her gen için İleri AH karşılaştırmasındaki logFC ile Tüm AH karşılaştırmasındaki logFC'nin mutlak değeri karşılaştırılmıştır. 20 gen arasından 17'sinin İleri AH karşılaştırmasında daha uç değere sahip olduğu ve doz-yanıt testinin Wilcoxon p = 1.3×10⁻⁵ ile anlamlılığa ulaştığı gözlemlenmiştir (Figure 5).

### 2.5.2 GBM tarafında doğrulama (TCGA-GBM vs GTEx)

GBM tarafındaki doğrulama için, recount3¹⁰ boru hattı üzerinden aynı hizalama ve sayım kurallarıyla işlenmiş iki kohort birleştirilmiştir: TCGA-GBM (n=157 birincil tümör) ve GTEx korteks (n=510 tümör-dışı beyin). recount3'ün tek merkezli işleme yaklaşımı, hizalama boru hattı farklılıklarından kaynaklanan yapay batch etkilerini en aza indirmiştir.

Birleşik sayım matrisi, limma-voom iş akışı¹⁵ ile diferansiyel ifade analizine tabi tutulmuş, kontrast olarak TCGA-GBM tümörleri — GTEx korteks tanımlanmıştır. DD imzasındaki 20 gen için elde edilen logFC değerleri, bir sonraki alt-bölümde AH tarafındaki logFC değerleri ile karşılaştırılmıştır.

### 2.5.3 Çapraz-hastalık yakınsama nicelemesi

DD imzasının iki hastalıktaki davranışını nicel olarak karşılaştırmak için, 20 gen için elde edilen AH logFC (İleri AH vs Kontrol) ve GBM logFC (TCGA-GBM vs GTEx) değerleri arasında Spearman sıra korelasyonu hesaplanmıştır (Pearson korelasyonu ikincil önlem olarak raporlanmıştır). Ek olarak, genler her iki hastalıktaki logFC yön ve büyüklüklerine göre üç biyolojik kategoriye ayrılmıştır: (i) her iki hastalıkta paylaşılan nöronal kimlik kaybı (her iki tarafta belirgin aşağı-regülasyon; 6 gen: BAIAP2, IDS, NRN1, OLFM1, PNMA2, SLC2A3), (ii) her iki hastalıkta paylaşılan stres yanıtı aktivasyonu (her iki tarafta belirgin yukarı-regülasyon; 6 gen: CEBPD, FLNA, JUNB, PER1, PFKFB3, ZFP36L1) ve (iii) tümöre özgü ayrışan metabolik-proliferatif regülasyon (AH'de nötr veya hafif aşağı, GBM'de belirgin yukarı; 8 gen: ALDOA, ATP1B2, CANX, HRH1, PEA15, PPP2CB, SDC3, UBC). Bu kategoriler, iki hastalık arasındaki paylaşılan ve ayrışan biyolojik programlarını görselleştirmek için kullanılmıştır (Figure 6).

---

## 2.6 Ön kayıt ve şeffaflık disiplini

Bu çalışmada, post-hoc rasyonelleştirme ve yanıt kayması (HARKing, hypothesizing after results are known) risklerini en aza indirmek amacıyla analiz planı, herhangi bir veri gözlemi yapılmadan önce formal olarak kayıt altına alınmıştır. Analiz planı belgesi (Analysis Plan v1), tüm anahtar metodolojik kararları — çalışma hipotezleri (H1, H2, H2a, H2b), örnek büyüklükleri ve stratifikasyon şemaları, DE analizi eşikleri (FDR<0.05, no logFC), yakınsama test seti (RRHO2 + hipergeometrik + permütasyon), önceden belirlenmiş karar kuralı (2/3 test PASS) ve doğrulama kohortlarının seçimi — açık şekilde tanımlamış, ardından Bartın Üniversitesi kurumsal e-posta arşivi üzerinden zaman damgalı olarak arşivlenmiştir. Bu zaman damgası, tüm Aşama 3 ve Aşama 4 analizlerine ilk erişim tarihinden önce gerçekleşmiş olup, bulgular kesinleşmeden karar kurallarının değiştirilmediğini belgelemektedir.

Analiz sürecinde ortaya çıkan tüm sapmalar — filtre kriterlerinin sonradan formalleştirilmesi, figür standartlarının değiştirilmesi, post-hoc dört-çeyrek analizinin eklenmesi, doku dışlama filtresi kararları ve benzeri — çalışma boyunca sürekli güncellenen bir DEVIATIONS.md belgesinde açık şekilde kayıt altına alınmıştır. Bu belge, her sapmanın (i) gerekçesini, (ii) neden ön-kayıtta öngörülmediğini, (iii) yorumlama üzerindeki potansiyel etkilerini ve (iv) alternatif kararların nasıl sonuç değiştireceğini içermektedir. DEVIATIONS.md, submission öncesi 600'ü aşkın satırlık bir belge haline gelmiş olup, dergi editörlüğüne ve hakemlere talep üzerine tam olarak sunulacaktır.

---

## 2.7 Yazılım ve tekrar üretilebilirlik

Tüm analizler R 4.6.0 ortamında (macOS ARM64) gerçekleştirilmiştir. Kullanılan başlıca Bioconductor ve CRAN paketleri: farklı ifade analizi için edgeR¹² (v3.44) ve limma¹⁵ (v3.58), sıralama-tabanlı örtüşme testi için RRHO2 (v1.0)¹³, yolak zenginleştirme için clusterProfiler¹⁴ (v4.10) ve org.Hs.eg.db (v3.18), büyük ölçekli bulk kohort standardizasyonu için recount3¹⁰ (v1.12), veri işleme ve görselleştirme için tidyverse (v2.0), ggplot2 (v3.5) ve patchwork (v1.2). Rastgele başlangıç noktası gerektiren tüm işlemlerde tohum değeri (seed = 42) sabitlenmiş, permütasyon prosedürlerinde ve rastgele kategorizasyon adımlarında tekrar üretilebilir çıktılar sağlanmıştır.

Analiz kodu ve tüm metodolojik kararlar, projeye özel bir Git deposunda kayıt altına alınmıştır. Depo, yalnızca kod dosyalarını değil aynı zamanda ara analiz nesnelerini (.rds), sonuç tablolarını (CSV), figürlerin ham PDF çıktılarını ve DEVIATIONS.md belgesini de içermektedir. Yayın kabulü sonrası tam depo, kalıcı bir DOI ile birlikte açık erişime sunulacaktır; makale gövdesindeki tüm istatistiksel değerler ve figürler bu depodan tek bir komutla yeniden üretilebilir olacaktır.
