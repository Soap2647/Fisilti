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
