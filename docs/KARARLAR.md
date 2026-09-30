# Karar Günlüğü

Projede başvuru dosyasından sapan ya da yöntemsel önemi olan her karar burada
gerekçesiyle kayıtlıdır. Amaç: jüri/hakem sorularına izlenebilir yanıt vermek.

Biçim: **K-NNN — Başlık** · tarih · ilgili WP

---

## K-001 — ICGC-PACA-CA yerine CPTAC-PDAC ana keşif kohortu
*2026-09-30 · WP1*

**Durum:** ICGC Veri Portalı (dcc.icgc.org) Haziran 2024'te kapatıldı. PACA-CA
verisine artık yalnızca DACO onaylı kontrollü erişimle (SFTP) ulaşılabiliyor.

**Karar:** PACA-CA kapsam dışı. Ana çoklu-omik keşif kohortu olarak CPTAC-PDAC
(Cao ve ark., Cell 2021) kullanılır.

**Gerekçe:** CPTAC-PDAC, başvuruda hedeflenen dört omik katmanını (mRNA, miRNA,
protein, metilasyon) *aynı hastalarda* ve *tümör + normal doku* için birlikte
sunar. MOFA+ eşleşmiş örneklerden en çok fayda sağladığı için bu, PACA-CA'dan
metodolojik olarak daha güçlüdür. Kohort sayısı dört bağımsız kaynağın altına
düşmez (CPTAC, TCGA, GEO doku, GEO serum).

## K-002 — Kan (serum) doğrulama kohortlarının eklenmesi
*2026-09-30 · WP5, WP6*

**Sorun:** Keşif verilerinin tamamı doku kökenlidir; "kan tabanlı panel"
iddiası yalnızca sekresyon tahminine dayanırsa zayıf kalır.

**Karar:** Serum miRNA kohortları GSE85589 ve GSE59856 bağımsız kan doğrulaması
olarak eklenir. Paneldeki miRNA adayları doğrudan hasta serumunda sınanır.

**Gerekçe:** Dokuda bulunan sinyalin kanda da ölçülebildiğini gösteren doğrudan
kanıt, sıvı biyopsi iddiasının en güçlü dayanağıdır.

## K-003 — Veri sızıntısını önleme ilkesi
*2026-09-30 · WP3, WP6*

**Karar:** Özellik seçimi, normalizasyon parametreleri, batch düzeltmesi ve
hiperparametre optimizasyonu yalnızca eğitim katmanında öğrenilir (iç içe çapraz
doğrulama). Doğrulama kohortları model/panel sabitlendikten sonra **tek sefer**
değerlendirilir.

**Gerekçe:** Tüm veride seçilen özelliklerle yapılan çapraz doğrulama AUC'yi
sistematik olarak şişirir. Bu, biyobelirteç literatüründeki en yaygın hatadır.

## K-004 — Hasta düzeyinde veri bölme (GroupKFold)
*2026-09-30 · WP3, WP6*

**Sorun:** CPTAC-PDAC'ta normal doku örnekleri tümör örnekleriyle *aynı
hastalardan* alınmıştır (67 hastada eşleşmiş normal). GSE15471'de de hasta başına
tümör + normal + teknik tekrar bulunur.

**Karar:** Tüm çapraz doğrulama bölmeleri hasta kimliğine göre gruplanır; bir
hastanın hiçbir örneği aynı anda hem eğitim hem test katmanında bulunamaz. Teknik
tekrarlar (GSE15471 `_rep`) önce ortalanır.

**Gerekçe:** Aksi halde model hastaya özgü (germline, bireysel) sinyali öğrenir
ve test performansı gerçekçi olmayan biçimde yükselir.

## K-005 — GEO kohortlarına özgü ön işleme
*2026-09-30 · WP1*

| Kohort | Bulgu | İşlem |
|---|---|---|
| GSE15471 | 39 tümör / 39 normal, 3 hastada teknik tekrar | Tekrarlar hasta+doku tipi bazında ortalanır |
| GSE62452 | Evre bilgisi var (PDAC: 49 erken / 19 geç); 40 örneğin başlığı hasta değil numune no. (ör. `9727`/`9728T`) | Erken evre alt-analizi; eşlenemeyen örnekler ayrı hasta sayılır (yalnızca dış doğrulamada kullanıldığından sızıntı riski yok) |
| GSE16515 | Değerler log dönüşümsüz (maks. 65 045); bir örnekte hasta no. boş | log2(x) uygulanır; boş kimliğe GSM no. atanır |
| GSE85589 | Serum; PDAC 88, sağlıklı 19, diğer GİS kanserleri; **CA19-9 durumu** | CA19-9 ile doğrudan kıyas (K-006) |
| GSE59856 | Serum; PDAC 100, sağlıklı 150, benign pankreas/safra 21, diğer kanserler 300; evre + operabilite | Özgüllük ve benign ayrımı analizi |

## K-006 — CA19-9 ile doğrudan karşılaştırma
*2026-09-30 · WP6*

**Karar:** GSE85589'da örnek bazında mevcut CA19-9 durumu (>37 / <37 U/mL)
kullanılarak panelin performansı CA19-9 ile aynı hastalarda kıyaslanır; ayrıca
"CA19-9 negatif" PDAC hastalarında panelin yakalama oranı raporlanır.

**Gerekçe:** Başvurunun ihtiyaç gerekçesi CA19-9'un yetersizliğidir; bu iddia
kendi verimizle sınanmalıdır. CA19-9 negatif hastaları yakalamak, klinik
açıdan en anlamlı katkıdır.

## K-007 — Doku bileşimi karıştırıcısı
*2026-09-30 · WP3, WP5*

**Sorun:** Tümör-komşu normal pankreas dokusu asiner hücrelerden zengindir;
tümör/normal ayrımında öne çıkan genlerin bir kısmı kanser biyolojisi değil doku
bileşimi farkını yansıtabilir (ör. asiner enzim genlerinin tümörde "azalması").
CPTAC her örnek için neoplastik, asiner, stromal vb. fraksiyonları sunar.

**Karar (K-007):** Aday belirteçler, doku fraksiyonlarıyla korelasyon açısından
denetlenir; kanda ölçüm için tümörde *artan* ve salgılanan adaylar önceliklendirilir.

## K-008 — TCGA-PAAD: birincil tanı seçimi ve histoloji filtresi
*2026-09-30 · WP1*

**Sorun:** GDC bazı vakalarda birden fazla tanı kaydı tutar (ör. sonradan
gelişen metastaz). İlk kaydın alınması 35 vakada evrenin, bazı vakalarda da
histolojinin yanlış okunmasına yol açıyordu (ör. TCGA-2J-AAB9: ilk kayıt metastaz).

**Karar:** `diagnosis_is_primary_disease = True` olan kayıt kullanılır
(`gdc._primary_diagnosis`, birim testli). Yalnızca duktal adenokarsinom
histolojileri tutulur (Infiltrating duct carcinoma NOS, Adenocarcinoma NOS,
Mucinous adenocarcinoma); nöroendokrin karsinom (8), indiferansiye ve karışık
alt tipler (2) dışlanır.

**Sonuç:** 168 PDAC tümörü + 4 normal; evresi bilinmeyen vaka 35 → 2.
PDAC evre dağılımı: IA 5, IB 10, IIA 28, IIB 116, III 3, IV 4 (erken evre = 159).
TCGA neredeyse tamamen erken evre olduğundan erken/geç ayrımı için değil,
**erken evre tümör sinyalinin bağımsız doğrulaması** için kullanılır.

## K-009 — LinkedOmics CPTAC miRNA verisinde teknik parti etkisi
*2026-09-30 · WP1*

**Bulgu:** LinkedOmics CPTAC-PDAC miRNA tablosunda örnekler PC1'de (varyansın
%51'i) birbirinden tamamen kopuk iki kümeye ayrılıyordu. Kümeler tespit
derinliğinde farklıydı (~1250 vs ~480 miRNA/örnek; medyan log2 −0,05 vs 4,46)
ve tümör/normal, merkez (C3L/C3N) veya ülke ile ilişkili değildi. Ham veriden
parti etiketi, lojistik regresyonla CV-AUC = 1,000 ile tahmin edilebiliyordu.

**Karar:** ComBat ile düzeltmek yerine veri kaynağı değiştirildi (K-010). GDC'nin
harmonize işlenmiş CPTAC miRNA verisinde bu iki kümeli yapı bulunmamaktadır;
dolayısıyla etki LinkedOmics işleme sürecine özgüdür.

## K-010 — CPTAC normal miRNA verisi: kopya tespiti ve GDC'ye geçiş
*2026-09-30 · WP1*

**Bulgu:** LinkedOmics `microRNA_TPM_log2_Normal.cct` dosyasındaki 21 normal
örneğin **21'i de aynı hastanın tümör verisinin birebir kopyasıdır**
(Pearson r = 1,0000; `np.allclose` = True). Bilinen belirteçlerde (miR-21,
miR-216a, miR-217, miR-375) tümör/normal farkı hiç görülmüyordu. Aynı test
protein (0/67) ve metilasyon (0/32) için kopya göstermedi; mRNA'da asiner genler
normalde ~64 kat yüksek (biyolojik olarak beklenen).

**Karar:** CPTAC miRNA'sı GDC'deki açık erişimli harmonize isoform dosyalarından
(CPTAC-3, Pancreas: 146 T / 38 N) yeniden türetildi. Tümörler çalışmanın 140
hastasıyla sınırlandı. Aynı işlem TCGA-PAAD için de uygulandı.

**Doğrulama (GDC verisi, 140 T / 38 N):** miR-21-5p ↑ (p = 4×10⁻¹⁶),
miR-196a-5p ↑ (p = 1×10⁻¹³), miR-210-3p ↑ (p = 2×10⁻¹⁴), miR-375-3p ↓
(p = 1×10⁻⁷), miR-216b-5p ↓ — literatürle uyumlu. Kopya: 0/38.

**Ek kazanım:** Olgun miRNA düzeyinde MIMAT erişim numarası anahtar olarak
kullanıldı; doku (CPTAC, TCGA) ve serum (GSE85589, GSE59856) kohortları miRBase
ad değişikliklerinden bağımsız olarak eşleştirilebilir.

## K-011 — MOFA+ tasarımı: denetimsiz keşif, CV içinde yeniden eğitim
*2026-09-30 · WP2, WP3*

**Karar:**
- MOFA+ yalnızca CPTAC-PDAC'ta (eşleşmiş çoklu-omik tek kohort) eğitilir; tek
  grup (tümör + normal birlikte), 4 görünüm, en fazla 15 faktör, ARD + spike-and-slab
  ağırlık önselleri, %0,5'ten az varyans açıklayan faktörler budanır.
  Yakınsama: 'medium' ölçütü, en fazla 1000 iterasyon ('slow' + 3000 iterasyon bu
  veride >65 dk sürüyordu; 'medium' 716. iterasyonda yakınsadı).
- Özellik ön seçimi yalnızca varyansa göredir (mRNA 3000, protein 2000,
  metilasyon 3000, miRNA 500); sınıf etiketi kullanılmaz. Metilasyon beta
  değerleri M değerine dönüştürülür (Gauss olabilirliği için).
- Ölçülmemiş görünümler eksik veri olarak modele verilir (MOFA+ bunu
  olasılıksal olarak ele alır); örnek atılmaz.
- Faktörlerin klinik anlamı (tümör/normal, erken/geç evre, doku bileşimi)
  model kurulduktan sonra, BH-FDR düzeltmeli testlerle değerlendirilir.
- WP3'te MOFA faktörleri sınıflandırıcı girdisi olarak kullanılırsa MOFA+
  çapraz doğrulamanın her eğitim katmanında yeniden eğitilir; test örnekleri
  eğitilmiş ağırlıklarla faktör uzayına yansıtılır.

**Gerekçe:** Denetimsiz faktör analizi etiket sızıntısı yaratmaz ve biyolojik
yorumu kolaylaştırır. Etiketle ilişkilendirilmiş faktör/özellik *seçimi* ise
denetimli bir adımdır ve yalnızca eğitim verisinde yapılmalıdır (K-003).

## K-012 — MOFA+ bulguları: baskın tümör/normal ekseni doku bileşimidir
*2026-09-30 · WP2 → WP3, WP5*

**Model:** 13 aktif faktör; açıklanan toplam varyans mRNA %48, protein %56,
miRNA %59, metilasyon %70. Kendi R² hesabımız MOFA+'ın iç hesabıyla birebir
örtüşür (en büyük fark 0,000 puan).

**Bulgular:**
1. **F1** tümör/normali neredeyse kusursuz ayırır (AUC 0,988; yalnızca 4 omiği
   tam örneklerde 0,983; her omik alt kümesinde >0,98 → eksik görünüm deseninden
   kaynaklanan yapay bir sinyal değildir). Ancak en yüksek ağırlıklı proteinler
   asiner sindirim enzimleridir (PLA2G1B, CPA2, PNLIP, AMY2A/B, PRSS3, SYCN).
2. **F2** aynı asiner programın mRNA/miRNA yüzüdür (CPA1, CELA3A/B, CTRB1/2;
   miR-216a/b, miR-217, miR-802) ve tümörlerde asiner fraksiyonla güçlü
   korelasyon gösterir (Spearman ρ = 0,70).
3. **Hiçbir faktör erken ve geç evreyi ayırmaz** (en iyi AUC 0,59; tüm q > 0,2).
4. Doku bileşimiyle ilişkili diğer faktörler: F5 ↔ neoplastik hücre oranı
   (ρ = 0,45), F9 ↔ stromal oran (ρ = −0,41).

**Yorum ve karar:**
- Dokuda tümörü normalden ayıran en büyük sinyal, tümörün asiner (enzim üreten)
  dokuyu ortadan kaldırmasıdır. Bu, K-007'de öngörülen karıştırıcıyı gerçek
  veride doğrular. Asiner enzim kaybı **kan testi için uygun değildir**:
  kanda "azalan" doku sinyali güvenilir ölçülemez; üstelik pankreatitte bu
  enzimler kanda artar ve yanlış pozitife yol açar.
- **WP3/WP5'e yansıması:** Asiner program genleri/miRNA'ları (F1/F2 yüksek
  ağırlıklı özellikler ve asiner fraksiyonla |ρ| > 0,5 korelasyon gösteren
  özellikler) aday panelden dışlanır ya da ayrıca işaretlenir; tümörde **artan**
  ve salgılanan adaylar önceliklendirilir.
- Evre ekseni bulunmaması, erken evre tümörlerin tümör sinyalini zaten tam
  taşıdığını düşündürür: PDAC-normal ayrımı için geliştirilen bir panelin erken
  evrede de çalışması beklenir. Bu, WP6'da erken/geç evre duyarlılığı ayrı
  raporlanarak sınanacaktır.

## K-013 — WP3 modelleme protokolü ve doku performansının doğru yorumu
*2026-09-30 · WP3*

**Protokol:** Omik başına tümör/normal sınıflandırması; 5 katlı × 5 tekrarlı
hasta-gruplu dış CV; her dış eğitim katmanında 3 katlı grup-CV + Optuna (15
deneme, TPE). Ön filtre (en fazla 200 özellik) pipeline'ın parçasıdır ve her iç
katmanda yeniden öğrenilir. Modeller: L1-lojistik (referans), XGBoost, LightGBM.
Füzyon: kapılı dikkatli MIL (torba = örnek, parça = ölçülmüş omik görünümü;
görünüm başına ön filtre 30–50 özellik; 3 tohumlu topluluk; eğitim hastalarının
%20'si ile erken durdurma). Aynı dış katmanlar tüm modellerde kullanılır.

**İki senaryo:** `free` (tüm özellikler, iki yön) ve `blood` (yalnızca tümörde
artan + asiner program dışlaması: tümörlerde asiner fraksiyonla ρ > 0,5).
Kan yönelimli kısıtlar dokuda anlamlı AUC kaybına yol açmadı (tüm omiklerde
fark güven aralıkları içinde).

**Hesaplama notu:** Tek değişkenli filtre istatistikleri aynı eğitim alt kümesi
için önbelleklenir (Optuna denemeleri aynı iç katmanı yeniden kullanır);
sonuç değişmez (birim testli), süre ~25 kat kısalır.

**Yorum uyarısı (sunum için bağlayıcı):** Doku düzeyinde tümör/normal AUC'leri
(0,96–0,99) *kolay* bir problemi ölçer ve kan testi performansı değildir.
Panelin gerçek sınavı bağımsız serum kohortlarıdır (GSE85589, GSE59856; WP6).
Raporda doku ve serum sonuçları ayrı başlıklarda verilir.

**Aday seçimi:** Kararlılık seçimi (blood senaryosu, 75 model = 25 dış katman ×
3 model); seçim sıklığı ile sıralama. Tek değişkenli AUC/log2FC betimleyicidir.

## K-014 — Attention-MIL'de ölçüm deseni kısayolu ve çoklu-omik füzyonun katkısı
*2026-09-30 · WP3*

**Bulgu:** CPTAC'ta ölçülen omik sayısı etiketle tamamen örtüşür: yalnızca 1–2
omiği ölçülmüş 71 örneğin **tamamı normaldir** (PDAC: 0). Tüm örneklerle
eğitilen MIL genel AUC 0,986 verdi, ancak 4 omiği tam 117 örnekte AUC 0,935'e
düştü (tek omikli modellerin tümünün altında). Model "az görünüm = normal"
kısayolunu öğrenmiş olabilir; genel AUC bu nedenle **raporlanmaz**.

**Karar:** MIL yalnızca 4 omiği tam örneklerde (100 T / 17 N) eğitilir ve
değerlendirilir (`--complete-only`). Sonuç: AUC 0,973 (%95 GA 0,918–1,000;
tekrarlar arası 0,952 ± 0,033), erken evre 0,975.

**Karşılaştırma (aynı 117 örnek, blood senaryosu):** mRNA-LightGBM 0,992;
miRNA-LightGBM 0,987; protein-L1 0,974; metilasyon-XGBoost 0,964; MIL 0,973.

**Yorum:** Dokuda çoklu-omik füzyon en iyi tek omiğin üzerine anlamlı katkı
sağlamaz; tümör/normal sinyali her omikte zaten güçlüdür. Bu, panelin tek
ölçüm platformunda (serum miRNA → RT-qPCR veya plazma proteini → ELISA)
kurulabileceğini destekler; bu da klinik maliyet ve uygulanabilirlik açısından
avantajdır. MIL'in katkısı yorumlanabilirliktir (örnek başına omik dikkati).
Genel ilke: eksik görünümlerin deseni etiketle ilişkiliyse füzyon modelleri
yalnızca tam örneklerde ya da desenden bağımsız bir tasarımla değerlendirilmelidir.

## K-015 — WP4 yorumlanabilirlik: dış katman SHAP ve yolak analizi
*2026-09-30 · WP4*

**SHAP:** Her omikte XGBoost (blood senaryosu) WP3'ün 1. tekrarındaki 5 dış
katmanda yeniden eğitildi; SHAP değerleri yalnızca test örneklerinde, XGBoost'un
yerleşik kesin TreeSHAP'ı ile (`pred_contribs`) hesaplandı. En etkili özellikler
kararlılık seçimiyle büyük ölçüde örtüşür (miRNA: miR-135b, miR-21, miR-203a,
miR-196a; protein: SERPINB5, SULF1, SDCBP2, SAMD9; metilasyon: CCNA1, HOXD10,
TFPI2). Not: en etkili özelliklerin tamamının "yüksek değer → PDAC" yönünde
olması **bulgu değildir**; 'up' yön filtresinin doğal sonucudur.

**GSEA (tüm genler, Welch t sıralaması):** Tümörde EMT (protein NES 3,27),
interferon yanıtları, koagülasyon (2,80), KRAS sinyali (2,23), E2F/G2-M
proliferasyon programları (mRNA) aktif; oksidatif fosforilasyon (protein) ve
pankreas beta hücre programı baskılanmış — PDAC biyolojisiyle tutarlı.

**ORA (arka plan = test edilen tüm genler):** Protein adayları (sel_freq ≥ 0,5,
n = 39) p53 sinyal yolağı (FDR 0,002) ve ECM-reseptör etkileşiminde zenginleşir.

**Sınırlılıklar:**
- Myc hedefleri ve oksidatif fosforilasyon mRNA'da pozitif, proteinde negatif
  NES gösterir; asiner doku (yüksek protein sentezi/mitokondri) bileşim farkının
  protein düzeyine yansıması olası açıklamadır (K-007, K-012).
- miRNA hedef zenginleştirmesi (miRTarBase 2017; ≥2 aday miRNA hedefi) keşif
  amaçlıdır: veritabanı kanser çalışmalarına yanlı olduğundan "Pathways in
  cancer" gibi genel terimlerin zenginleşmesi beklenen bir yanlılıktır.
- Gen seti kütüphaneleri erişim tarihiyle birlikte data/raw/genesets/ altında
  saklanır.

## K-016 — miRNA panelinin serum etiketlerine bakılmadan kilitlenmesi
*2026-09-30 · WP5*

**İlke:** Panel üyeleri ve skorlama kuralı yalnızca doku verisinden ve
*etiketsiz* serum ölçülebilirlik kontrolünden türetilir; serum kohortlarının
PDAC/kontrol etiketleri panel kilitlenene kadar hiçbir adımda kullanılmaz.
Kilit dosyası (`results/panel/mirna_panel_v1.json`) doğrulamadan önce ayrı bir
commit ile git'e işlenir; commit zaman damgası ön-kayıt kanıtıdır.

**Filtre hunisi:** 40 doku adayı (sel_freq ≥ 0,5) → asiner güvenliği (ρ < 0,3)
40 → kan hücresi/hemoliz kökenli miRNA dışlaması 40 → serumda ölçülebilirlik
(GSE59856 örneklerinin ≥ %50'sinde 3D-Gene tespit sınırı 2,612 üzerinde +
GSE85589 platformunda bulunma) **15** → fazlalık giderme (dokuda |ρ| > 0,7;
miR-221-3p ve miR-210-5p → miR-210-3p, miR-1247-5p → miR-1247-3p) 12 →
panel boyutu 5 (doku grup-CV AUC eğrisi: 1→0,923, 3→0,970, 5→0,978 tepe; 10→0,974).

**Kilitli panel v1:** miR-21-5p, miR-210-3p, miR-27a-3p, miR-935, miR-24-3p
(hepsi tümörde artan yönde).

**Skorlama kuralı (kilitli):** Panel miRNA'larının kohort-içi z-skor ortalaması
(z-skor parametreleri etiketsiz, tüm örneklerden). Serumda ağırlık öğrenilmez.

**Gözlem:** Dokuda en güçlü adayların (miR-135b, miR-196a, miR-203a) serumda
çoğunlukla tespit sınırı altında kalması, "dokuda iyi = kanda ölçülebilir"
varsayımının geçersiz olabileceğini gösterir; ölçülebilirlik filtresi bu
nedenle zorunludur.

**Sınırlılık:** Skorlama kuralının kohort-içi z-skor kullanması, skorun
kohortun vaka karışımına (PDAC/kontrol oranı) bağlı olmasına yol açar; klinik
uygulamada sabit referans aralıkları (ör. sağlıklı referans popülasyonu
dağılımı) gerekecektir.

## K-017 — Fazlalık giderme hatası ve panellerin v2 olarak yeniden kilitlenmesi; protein aday paneli
*2026-09-30 · WP5*

**Hata:** Panel fazlalık giderme adımında Spearman korelasyonu tümör + normal
örnekler birlikte kullanılarak hesaplanmıştı. Bu durumda sınıf farkı, tümörde
artan tüm özellikleri yapay olarak yüksek korelasyonlu gösterir. Etkisi protein
panelinde belirgin oldu: 9 Düzey-1 (salgılanan) adayın 7'si, kana salgılanan
LGALS3BP dahil, "fazlalık" gerekçesiyle elenmişti.

**Düzeltme:** Fazlalık korelasyonu yalnızca tümörler içinde hesaplanır (iki
belirtecin sınıf farkının ötesinde aynı bilgiyi taşıyıp taşımadığı).

**Zamanlama ve şeffaflık:** Hata, serum etiketleri hiçbir analizde
kullanılmadan önce fark edildi. miRNA paneli v1 (commit fca9e9b) ve v2 **aynı
üyelere** sahiptir (miR-21-5p, miR-210-3p, miR-27a-3p, miR-935, miR-24-3p);
düzeltilmiş kuralla v2 olarak yeniden kilitlendi. v1 dosyası geçmiş kaydı
olarak depoda tutulur. WP6 doğrulaması **v2** dosyasını kullanır.

**Protein aday paneli v2 (8):** SULF1, LGALS3BP, NMI, LAMC2, THBS1, CD55
(Düzey 1: salgılanan + kanda ölçülmüş), GLA, SERPINB5 (Düzey 2: kanda ölçülmüş).
Doku grup-CV AUC 0,992 (8 üye). Salgılanma/plazma kanıtı: HPA (secretome
location, kanda MS konsantrasyonu) + UniProt (sinyal peptidi, 'Secreted');
ham yanıtlar data/raw/annotations/ altında tarih damgalı.

**Başvurudan sapma:** SignalP 6.0 ve SecretomeP yerine UniProt sinyal peptidi
anotasyonu ve HPA salgılanma sınıflandırması kullanıldı (programatik,
tekrarlanabilir, lisans gerektirmez; HPA sınıfı klasik dışı salgılanmayı da kapsar).

**Sınırlılıklar:**
- Protein paneli kanda **doğrulanmamıştır** (plazma proteomik kohortu yok);
  "aday panel" olarak raporlanır.
- THBS1 trombosit α-granüllerinde boldur ve pıhtılaşma sırasında serbestlenir;
  ölçüm **plazmada** yapılmalıdır (serum değerleri trombosit aktivasyonunu yansıtır).
- Hiçbir adayda HPA immünoanaliz (IM) kan ölçümü bulunmadı; ELISA
  uygulanabilirliği ticari kit/antikor taramasıyla ayrıca teyit edilmelidir.

## K-018 — Serum doğrulaması için ön-kayıtlı analiz planı
*2026-09-30 · WP6*

Analiz planı (`scripts/wp6_01_serum_validation.py` başlığı) sonuçlar
görülmeden önce ayrı bir commit ile kaydedilmiştir. Birincil sonlanım: her iki
serum kohortunda kilitli panel v2 ile PDAC vs sağlıklı AUC. İkincil sonlanımlar
(erken/ameliyat edilebilir/geç evre, benign, diğer kanserler, %90/%95 özgüllükte
duyarlılık, CA19-9 kıyası, rastgele panel kontrolü, DCA) önceden
listelenmiştir. Plan dışı analizler raporda "keşif amaçlı" olarak ayrıca etiketlenir.

## K-019 — Ön-kayıtlı serum doğrulamasının sonucu: kilitli miRNA paneli başarısız
*2026-09-30 · WP6*

**Sonuç (değiştirilmeden raporlanır):**
- Birincil sonlanım, PDAC vs sağlıklı AUC: GSE85589 **0,536** (%95 GA
  0,405–0,666); GSE59856 **0,396** (0,324–0,470). %95 özgüllükte duyarlılık
  %16 ve %7.
- Rastgele 5'li panel kontrolü (serumda ölçülebilir 880 miRNA havuzu):
  ampirik p = 0,43 ve 0,71 → panel şanstan iyi değil.
- Erken evre PDAC vs sağlıklı (GSE59856, n = 19): 0,433. Ameliyat edilebilir: 0,427.
- CA19-9 (GSE85589): duyarlılık %74, özgüllük %100; panel aynı özgüllükte %16.
  CA19-9 negatif 23 PDAC'ın %17'si panelle %95 özgüllük eşiğinde pozitif.
- Tek miRNA düzeyinde miR-24-3p GSE59856'da ters yönde (AUC 0,22).

**Yorum:** Dokuda tümörde artan ve serumda ölçülebilen miRNA'lar serumda
PDAC'ı ayırt etmedi. Doku sinyalinin kana doğrudan aktarılabileceği varsayımı
bu veride desteklenmedi. Olası nedenler (sonraki keşif analizleriyle sınanacak):
(1) dolaşımdaki miRNA havuzuna kan hücreleri ve diğer dokuların katkısının tümör
katkısını baskılaması; (2) serum miRNA düzeylerinin tümör dışı (sistemik,
inflamatuvar) süreçlerden etkilenmesi; (3) platform ve ön-analitik farklılıklar.

**Durum:** Bu sonuç projenin temel varsayımlarından birini çürütür ve
raporda birincil bulgu olarak sunulacaktır. Bundan sonraki analizler
**keşif amaçlıdır** ve ön-kayıtlı sonucu değiştirmez.
