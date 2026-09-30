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
