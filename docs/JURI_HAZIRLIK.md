# Fısıltı — Jüri Hazırlık Rehberi

Bu dosya projeyi sıfırdan anlatır. Sıra: önce hikâye, sonra kavramlar, sonra yaptıklarımız, sonra rakamlar, en sonda jürinin sorabileceği sorular ve cevapları. Sunumdan önce en az iki kez baştan sona oku; "Ezberle" diye işaretli kutuları yüksek sesle tekrar et.

---

## 1. Otuz saniyede proje (ezberle)

> "Pankreas kanseri sessiz ilerler; çoğu hasta ameliyat edilemeyecek evrede tanı alır ve beş yıllık sağkalım %10'un altındadır. Fısıltı'da, dünyanın açık paylaştığı yedi hasta veri setini yapay zekâ ile inceleyip, kanda ölçülebilecek ve erken evreyi yakalayabilecek 8 proteinlik bir aday panel bulduk. Bu panel, hiç görmediği üç bağımsız hastane verisinde tümörü normalden AUC 0,87 ile 0,95 arasında ayırdı; erken evrede de 0,86. Kanda çalışmayan miRNA panelimizi saklamadık, nedenini gösterdik. Sıradaki adım, 8 proteini gerçek hasta plazmasında ELISA ile ölçmek."

## 2. Tek cümlelik özet

Açık verilerden, yapay zekâ ile, erken pankreas kanseri için kanda ölçülebilir bir biyobelirteç adayı bulduk ve bunu ön-kayıtlı, dürüst bir şekilde test ettik.

---

## 3. Problem: neden pankreas kanseri?

- **Pankreas** karnın arkasında, midenin arkasında duran bir organ. İki işi var: sindirim enzimi üretmek (asiner hücreler) ve insülin üretmek (Langerhans adacıkları).
- **PDAC** (pankreas duktal adenokarsinomu): pankreas kanserlerinin yaklaşık %90'ı. Pankreas kanallarını döşeyen hücrelerden çıkar.
- **Neden ölümcül?** Belirti vermez; belirti çıktığında (sarılık, kilo kaybı, sırt ağrısı) çoğunlukla yayılmıştır. Hastaların yalnızca küçük bir kısmı ameliyat edilebilir evrede yakalanır.
- **Beş yıllık sağkalım %10'un altında** (kaynak: Klein AP, Nat Rev Gastroenterol Hepatol 2021). Erken evrede yakalananlarda bu oran belirgin şekilde daha yüksek.
- **Tek kan belirteci CA19-9.** Sorunları:
  - Erken evrede sık yükselmez (duyarlılık düşük).
  - İyi huylu durumlarda da yükselir (sarılık, pankreatit) → yanlış alarm.
  - İnsanların yaklaşık %5–10'u (Lewis antijeni negatif) CA19-9'u hiç üretemez.
  - Bu yüzden tarama testi olarak önerilmez; teşhis sonrası takipte kullanılır.

**İhtiyaç:** Basit bir kan testiyle, özellikle yüksek riskli kişilerde (yeni başlayan diyabet, kronik pankreatit, ailede pankreas kanseri), kanseri erken evrede yakalamak.

---

## 4. Bilmen gereken kavramlar (sade dille)

| Kavram | Sade anlamı |
|---|---|
| **Biyobelirteç** | Ölçülebilen ve bir hastalığa işaret eden molekül. Ör. şeker hastalığında kan şekeri. |
| **Panel** | Tek molekül yerine birkaç molekülün birlikte değerlendirilmesi. |
| **Omik** | Hücredeki bir molekül türünün toplu ölçümü. |
| **mRNA (transkriptom)** | Genden çıkan "tarif kopyası"; hangi genin ne kadar çalıştığını gösterir. |
| **Protein (proteom)** | Tariflerden üretilen asıl iş gören moleküller. Kanda ELISA ile ölçülen şey. |
| **miRNA** | Küçük RNA'lar; diğer genlerin üretimini kısar. Kanda kararlıdır, RT-qPCR ile ölçülür. |
| **Metilasyon** | DNA üzerindeki "aç/kapa düğmeleri". |
| **Çoklu-omik** | Bu katmanları birlikte kullanmak. |
| **Kohort** | Bir çalışmada toplanmış hasta grubu / veri seti. |
| **ELISA** | Hastane laboratuvarında kandaki bir proteini ölçen standart, ucuz test. |
| **RT-qPCR** | Kandaki RNA miktarını ölçen standart test (COVID testinin yöntemi). |
| **AUC** | Model bir hasta ile bir sağlıklıyı ne sıklıkla doğru sıralıyor: **0,5 = yazı tura, 1,0 = kusursuz.** 0,8 iyi, 0,9 çok iyi kabul edilir. |
| **Duyarlılık** | Hastaların yüzde kaçını yakalıyor. |
| **Özgüllük** | Sağlamların yüzde kaçına "temiz" diyor (yanlış alarm vermiyor). |
| **Çapraz doğrulama** | Veriyi parçalara bölüp modeli bir kısmında eğitip diğerinde test etmek; ezber yapmasını önler. |
| **İç içe çapraz doğrulama** | Ayar yapma (hiperparametre) ile test etmeyi ayrı katmanlarda yapmak; test verisi ayar için hiç kullanılmaz. |
| **Hasta gruplu bölme** | Aynı hastanın tümör ve normal örneği hep aynı tarafta kalır; model "bu hastayı tanıyorum" diye kopya çekemez. |
| **Sızıntı (leakage)** | Test bilgisinin farkında olmadan eğitime karışması; sahte yüksek başarı üretir. Projenin en çok korunduğu tuzak. |
| **Ön-kayıt** | Sonucu görmeden panel ve analiz planını kilitleyip tarih damgasıyla kaydetmek. "Cevapları sınavdan önce zarfa mühürlemek." |
| **Bağımsız doğrulama** | Modeli, eğitimde hiç görmediği başka bir hastanenin/çalışmanın verisinde test etmek. |

---

## 5. Veriler (yedi kohort)

Hepsi kamuya açık ve anonim; bu yüzden etik kurul onayı gerekmedi.

| Kohort | Rol | Örnek | Omik |
|---|---|---|---|
| **CPTAC-PDAC** | Keşif (modeller burada eğitildi) | 140 tümör · 92 normal | 4 omik, aynı hastalarda |
| **TCGA-PAAD** | Doku doğrulama | 175 tümör · 10 normal | 4 omik |
| **GSE15471** | Bağımsız doku doğrulama | 36 tümör · 36 normal | mRNA |
| **GSE62452** | Bağımsız doku doğrulama, evre bilgili | 69 tümör · 61 normal | mRNA |
| **GSE16515** | Bağımsız doku doğrulama | 36 tümör · 16 normal | mRNA |
| **GSE85589** | Serum (kan) | 88 PDAC · 19 sağlıklı · 125 diğer | miRNA |
| **GSE59856** | Serum (kan) | 100 PDAC · 150 sağlıklı · 321 diğer | miRNA |

- **CPTAC**: ABD Ulusal Kanser Enstitüsü'nün proteomik projesi. Seçme sebebimiz: dört omik aynı hastalarda ölçülmüş, protein verisi var.
- **TCGA**: Kanser Genom Atlası.
- **GEO**: ABD'nin açık gen ifadesi veri arşivi; GSE ile başlayan numaralar oradaki çalışmalar.
- **ICGC** başvuruda geçiyordu ama portalı 2024'te kapandı; yerine CPTAC kullanıldı (karar günlüğünde kayıtlı).

---

## 6. Ne yaptık? Beş adım

**1. Veri toplama ve temizlik.** Yedi kohort indirildi, her dosyanın SHA-256 özeti kaydedildi (aynı veriyi herkes tekrar indirebilir). Altı ciddi veri tuzağı bulundu ve düzeltildi (bölüm 7).

**2. MOFA+ (çoklu-omik faktör analizi).** Dört omik katmanını etiket kullanmadan birleştiren bir yöntem; verideki ana "eksenleri" (13 faktör) bulur. En önemli bulgusu bir uyarıydı (bölüm 8).

**3. Yapay zekâ modelleri.**
- **L1-lojistik regresyon:** basit, açıklanabilir; gereksiz özellikleri sıfırlar.
- **XGBoost ve LightGBM:** karar ağaçlarına dayalı güçlü modeller.
- **Attention-MIL:** dikkat mekanizmalı derin öğrenme; her hastanın dört omiğini bir "torba" olarak görür ve hangi omiğe ne kadar güveneceğini öğrenir.
- Ayarlar **Optuna** ile yapıldı; hepsi **5×5 iç içe, hasta gruplu çapraz doğrulama** ile değerlendirildi (toplam yaklaşık 600 model).

**4. Yorum.**
- **SHAP:** her molekülün karara katkısını gösterir; yalnızca test örneklerinde hesaplandı.
- **GSEA / yolak analizi:** modelin seçtiği genler hangi biyolojik süreçlere ait? (KRAS, EMT, p53 çıktı.)

**5. Panel ve doğrulama.**
- Adaylar şu filtrelerden geçti: tümörde artmalı, kana salgılanmalı (Human Protein Atlas + UniProt kanıtı), kanda ölçülebilir olmalı, kan hücresi kökenli olmamalı (hemoliz riski), birbirinin kopyası olmamalı.
- Paneller **kilitlendi**, analiz planı **sonuçtan önce** kaydedildi, sonra bağımsız kohortlarda test edildi.

---

## 7. Altı veri tuzağı (jüri bunu sever)

1. **CPTAC "normal" miRNA dosyaları tümörün birebir kopyasıydı** (21/21 çift, korelasyon r = 1,0000). Fark edilmeseydi model "kopyayı" öğrenirdi. GDC'deki ham veriden doğru normaller (38 gerçek normal) alındı.
2. **TCGA'da metastaz kaydı birincil tanının önündeydi.** 35 hastanın evresi yanlış okunuyordu; düzeltildi, eksik evre 35'ten 2'ye indi.
3. **Aynı hastanın tümörü ve normali** eğitimde ve testte ayrı düşebilirdi → hasta bazlı bölme.
4. **Eksik omik deseni etiketi ele veriyordu.** Normal örneklerde bazı omikler hep eksikti; Attention-MIL bunu öğrenip 0,986'ya çıktı. Sadece dört omiği tam olan hastalarla adil sonuç: 0,973.
5. **Panel seçiminde yanlış korelasyon.** Gereksizlik kontrolü tümör+normal birlikte yapılınca salgılanan iyi adaylar elenmişti. Tümör içinde yapılacak şekilde düzeltildi ve bu, **kan verisine bakmadan** yapıldı.
6. **Ölçek ve tekrar hataları:** GSE16515 log dönüşümsüz yüklenmişti, GSE15471'de teknik tekrarlar vardı; düzeltildi.

Hepsi `docs/KARARLAR.md` dosyasında (22 karar, K-001 … K-022) gerekçesiyle yazılı.

---

## 8. MOFA+: "en güçlü sinyal bir tuzaktı"

- Tümörü normalden en iyi ayıran eksen (AUC 0,99) aslında **kanser değil, pankreasın enzim üreten dokusunun (asiner doku) kaybıydı.** Tümör büyüdükçe normal doku azalıyor; model bunu "kanser işareti" sanabilir.
- Sorun: kanda "azalan doku" ölçülemez; ayrıca pankreatitte de asiner doku kaybolur → yanlış alarm.
- Karar: bu genler adaylardan çıkarıldı; yalnızca **tümörde artan** belirteçler seçildi.
- Mesaj: "Yüksek AUC her zaman doğru biyoloji demek değil; biz neyi öğrendiğine baktık."

---

## 9. Sonuçlar (rakamlar)

### 9a. Doku içinde yapay zekâ (CPTAC, iç içe çapraz doğrulama)
| Omik | En iyi AUC | Model |
|---|---|---|
| Protein | 0,99 | L1-lojistik |
| mRNA | 0,99 | LightGBM |
| Metilasyon | 0,99 | XGBoost |
| miRNA | 0,98 | XGBoost |
- Dört omiği birleştiren Attention-MIL: 0,973 → tek omikten iyi değil. Anlamı: **tek, ucuz bir test türü yeterli olabilir.**
- Erken evre performansı geç evreyle aynı.
- **Dikkat:** doku kolay bir problem; bu rakamlar kan testi performansı değil. Bunu kendin söyle, jüri söylemeden.

### 9b. Model bilinen biyolojiyi kendisi buldu
- **KRAS** sinyali (PDAC'ların %90'ından fazlasında KRAS mutasyonu var), **epitel-mezenkimal geçiş** (yoğun tümör stroması), **p53 yolağı** (FDR 0,002).
- Literatürde bilinen belirteçler önceden söylenmeden bulundu: **miR-21, miR-196a, SERPINB5, TFPI2 metilasyonu.**

### 9c. Paneller
- **miRNA paneli (5):** miR-21-5p, miR-210-3p, miR-27a-3p, miR-935, miR-24-3p. Eleme: 40 → 15 → 14 → 5.
- **Protein paneli (8):** SULF1, LGALS3BP, NMI, LAMC2, THBS1, CD55, GLA, SERPINB5. Eleme: 39 → 31 → 8.

### 9d. Ön-kayıt (30 Eylül 2026)
18:20 miRNA paneli kilitlendi → 18:26 hata düzeltmesi → 19:25 serum analiz planı kaydedildi → 19:28 serum sonucu → 19:36 doku doğrulama planı kaydedildi → sonra doku sonucu. Kanıt: GitHub commit geçmişi (değiştirilemez zaman damgaları).

### 9e. Gerçeğin anı 1: miRNA paneli kanda çalışmadı
- AUC **0,54** (GSE85589) ve **0,40** (GSE59856); 0,5 yazı tura.
- Rastgele 5 miRNA'dan iyi değil (p = 0,43 ve 0,71).
- Aynı serum verisinde **CA19-9:** duyarlılık %74, özgüllük %100.

### 9f. Neden çalışmadı?
- İki serum kohortu kendi içinde "mükemmel" (AUC ≈ 1,00) ama **birbirini doğrulamıyor**: birinde eğitilen model diğerinde 0,40 / 0,48; iki kohortun sinyalleri arasındaki korelasyon ρ = 0,02 (sıfır).
- Dokuda artan miRNA kanda artmıyor (ρ = 0,05).
- Yorum: kandaki miRNA'nın büyük kısmı kan hücrelerinden gelir; tümörün katkısı çok küçük. Ayrıca kohortlar farklı platform ve laboratuvarlarda ölçülmüş; "mükemmel" tek kohort sonuçları büyük ihtimalle teknik farklılıktan.
- Mesaj: **"Tek kohortta yüksek AUC'ye güvenmeyin."** Bu, alan için de bir uyarı.

### 9g. Gerçeğin anı 2: protein paneli bağımsız dokuda doğrulandı
- Hiç görmediği üç kohortta AUC **0,92** (GSE15471) · **0,87** (GSE62452) · **0,95** (GSE16515).
- Erken evre (GSE62452): **0,86.**
- 24/24 gen-kohort çiftinde değişim beklenen yönde.
- Rastgele 8 genlik panellerden üstün (**p ≤ 0,004**; demoda 2000 rastgele panelin hiçbiri yetişemedi).
- Sabit eşikle (skor > 0, etiket görülmeden) isabet: GSE16515 %90 · GSE15471 %83 · GSE62452 %79.
- **Sınır:** bunlar mRNA ölçümü ve doku verisi; protein panelinin kanda ölçümü yapılmadı.

---

## 10. Sekiz protein (kısa kimlik kartları)

| Protein | Ne iş yapar | Neden aday | Dış doğrulama AUC |
|---|---|---|---|
| **SULF1** | Hücre dışı şekerleri düzenleyen salgılanan enzim | Tümör stromasında yüksek | 0,84–0,97 |
| **LGALS3BP** | Salgılanan glikoprotein (Mac-2 bağlayıcı) | Kana salgılanır | 0,84–0,89 |
| **NMI** | İnterferon yanıtında rol alan protein | Kanda tespit edilmiş | 0,75–0,88 |
| **LAMC2** | Laminin γ2; invazyon cephesinde | Literatürde serumda PDAC için doğrulanmış | 0,91–0,95 |
| **THBS1** | Trombospondin-1; damar oluşumu, stroma | Bir şirket panelinde de var | 0,59–0,83 |
| **CD55** | Kompleman sistemini düzenler | Kana salgılanır | 0,83–0,93 |
| **GLA** | Lizozomal α-galaktozidaz | Plazmada ölçülmüş | 0,64–0,77 |
| **SERPINB5** | Maspin; PDAC'ta belirgin artar | Bilinen PDAC belirteci | 0,82–0,94 |

---

## 11. Literatürde neredeyiz? (orijinallik sorusu için)

- **Kan testiyle erken PDAC yakalama, dünyada çok çalışılan bir alan.** İlk biz değiliz; bunu açıkça söyle.
- **Immunovia** (şirket, 2025): TIMP1, ICAM1, THBS1, CTSD + CA19-9 paneli; 128 erken evre hasta ve 465 yüksek riskli kişide, %98 özgüllükte %85 duyarlılık (Cancers 2025). **THBS1 bizim panelimizde de var.**
- **LAMC2:** 2013'te doku proteomiğiyle bulunup 425 kişilik serumda doğrulanmış (Kosanam ve ark., Mol Cell Proteomics 2013); CA19-9 + CA125 + LAMC2 paneli (Chan ve ark., Clin Cancer Res 2014).
- **Bizim özgün katkımız:**
  1. Bu **8'li kombinasyon** literatürde yayımlanmış görünmüyor.
  2. Yöntemimiz, **başkalarının kanda doğruladığı iki proteini (LAMC2, THBS1) kendi kendine buldu** → doğru yere bakıyoruz; diğer 6'sı yeni aday.
  3. **"Kanda ölçülebilirlik" baştan bir tasarım kuralı.**
  4. Ön-kayıt, hasta bazlı doğrulama, olumsuz sonucu raporlama; iki açık serum kohortunun birbirini doğrulamadığını göstermek.
  5. Tamamen açık kaynak ve tekrarlanabilir.

> **Ezberle:** "Bu alanda ilk biz değiliz; şirketler bile var. Bizim farkımız: açık verilerden, baştan kanda ölçülebilirlik şartıyla yeni bir 8'li kombinasyon bulduk. Yöntemimizin doğru çalıştığını da şuradan biliyoruz: başkalarının kanda doğruladığı iki proteini kendiliğinden seçti."

---

## 12. Zayıf yönler (önce sen söyle)

1. **Kanda henüz doğrulama yok.** Protein paneli yalnızca doku verisinde test edildi. Kanda ölçülen açık protein verisi yok.
2. **miRNA paneli kanda başarısız oldu.**
3. **Dış doğrulama mRNA ile yapıldı**, protein ile değil. mRNA ve protein seviyesi her zaman birebir örtüşmez.
4. **Kontrol grubu çoğunlukla normal doku.** Asıl zor ayrım kanser ve pankreatit. **LGALS3BP** literatürde kronik pankreatitte PDAC'tan daha yüksek bulunmuş (Maekawa ve ark., World J Gastroenterol 2016); bu yüzden pilotta pankreatit grubu şart.
5. **Örnek sayıları orta düzeyde;** bazı kohortlarda erken evre bilgisi yok.
6. **Retrospektif veri;** prospektif (ileriye dönük) çalışma gerekecek.

Bunları saklamamak projeyi güçlü gösterir. Jüri zayıf yanı senden önce bulursa güven kaybedersin; sen söylersen olgunluk olarak görür.

---

## 13. Yol haritası ve bütçe

1. **Pilot (planlanan):** PDAC, kronik pankreatit ve sağlıklı gönüllülerden plazma; 8 proteinin ticari ELISA kitleriyle ölçümü; CA19-9 ile karşılaştırma. Etik kurul onayı ve bir hastaneyle iş birliği gerekiyor (hedef: Kastamonu Eğitim ve Araştırma Hastanesi / Tıp Fakültesi). Hastaneyle henüz anlaşma yapılmadı; "planlanan" de.
2. **Çok merkezli, ileriye dönük çalışma:** erken evre odaklı, ön-kayıtlı protokol.
3. **Fon:** TÜBİTAK 1001 / 2209 / 2232 tipi destekler; aynı çerçeve akciğer, kolorektal, mide kanserine uyarlanabilir.
4. **Fikri mülkiyet:** üniversite Teknoloji Transfer Ofisi ile patent ve lisans değerlendirmesi.
- **Şu ana kadarki maliyet: 0 TL** (açık veri, açık kaynak yazılım, kişisel bilgisayar).
- **Ödül (50.000 TL) ne için?** Pilot için ELISA kitleri ve laboratuvar sarf malzemesinin bir kısmı. Rakam sorulursa "kesin bütçe hastane ve laboratuvarla birlikte çıkarılacak; sekiz kit ve sarf ilk kalem" de, uydurma rakam verme.

**Hedef kullanıcılar:** onkoloji ve gastroenteroloji klinikleri; yüksek riskli bireylerin taranması (yeni başlayan diyabet, kronik pankreatit, ailede pankreas kanseri); tanı kiti (IVD) üreticileri.

---

## 14. Teknik altyapı (yazılım sorulursa)

- Python 3.12; pandas, scikit-learn, XGBoost, LightGBM, Optuna, PyTorch (Attention-MIL), mofapy2 (MOFA+), SHAP, gseapy.
- 34 birim testi; her veri dosyası SHA-256 ile kayıtlı; rastgelelik sabit tohumla kontrollü (aynı sonuçlar tekrar üretilebilir).
- Kod: github.com/Soap2647/Fisilti. **Sunumdan önce repoyu herkese açık yap**, posterdeki QR oraya gidiyor.
- Demo: `sunum/demo/fisilti_demo.html` (internetsiz çalışır).
- Kapanış filmi: Three.js ile kodla üretilmiş 3B animasyon; anatomi modeli BodyParts3D (açık lisanslı insan anatomisi veritabanı); sesler kodla sentezlendi.

---

## 15. Jürinin sorabileceği sorular ve cevaplar

### A. Genel / anlam
**S: Projeyi bir cümlede anlatır mısın?**
C: Açık verilerden yapay zekâ ile, erken pankreas kanserini kanda yakalayabilecek 8 proteinlik bir aday panel bulduk ve bağımsız hastane verilerinde doğruladık; sıradaki adım gerçek kanda ölçmek.

**S: Neden pankreas kanseri?**
C: En ölümcül kanserlerden biri, sessiz ilerliyor, erken yakalamanın etkisi çok büyük ve elimizde iyi bir tarama testi yok. Kazanç en büyük yerde çalışmak istedim.

**S: Bu bir kan testi mi? Kanseri teşhis edebiliyor musunuz?**
C: Hayır, henüz değil. Bir **aday panel** bulduk ve dokuda doğruladık. Kan testine dönüşmesi için gerçek hasta plazmasında ölçülmesi gerekiyor; yol haritamızın ilk adımı bu.

**S: Bunu neden sen yaptın, bölümün bilgisayar mühendisliği?**
C: Bu problemin darboğazı veri analizi ve yapay zekâ. Biyolojik kararlarda literatüre ve açık veri tabanlarına dayandım; laboratuvar kısmı için hekim ve biyolog iş birliği planlıyorum.

### B. Veri
**S: Hasta verisini nereden buldunuz? Etik kurul?**
C: Tamamı kamuya açık ve anonimleştirilmiş araştırma verileri: CPTAC, TCGA ve GEO. Kimlik bilgisi içermiyor; bu yüzden etik kurul onayı gerekmedi. Kanda pilot çalışma için etik kurul alınacak.

**S: Neden bu kohortlar?**
C: CPTAC'ta dört omik aynı hastalarda ölçülmüş ve protein verisi var; keşif için ideal. Doğrulama için farklı ülkelerden, farklı laboratuvarlardan bağımsız setler seçtim ki sonuç tek bir laboratuvara özgü olmasın.

**S: Veri kalitesinden emin misiniz?**
C: Değildik, o yüzden kontrol ettik ve altı ciddi sorun bulduk. Örneğin bir kaynakta "normal" diye etiketli miRNA dosyaları tümörün birebir kopyasıydı. Hepsini düzelttik ve karar günlüğüne yazdık.

**S: Örnek sayısı az değil mi?**
C: Keşif setinde 232 örnek var; doğrulamada üç bağımsız kohortta toplam 254 örnek. Bir biyobelirteç çalışması için makul, ama klinik karar için yeterli değil; bu yüzden çok merkezli bir çalışma planlıyoruz.

### C. Yöntem / yapay zekâ
**S: Hangi yapay zekâ modellerini kullandınız, neden birden fazla?**
C: L1-lojistik, XGBoost, LightGBM ve Attention-MIL. Farklı varsayımları olan modeller aynı belirteçleri seçiyorsa sonuç bir modelin tuhaflığı değildir. Basit model (L1-lojistik) proteinde en iyi çıktı; bu da açıklanabilirlik için iyi.

**S: AUC 0,99 çok yüksek, ezber (overfitting) olmasın?**
C: Üç önlemimiz var. Birincisi, iç içe çapraz doğrulama: ayar ve test ayrı katmanlarda. İkincisi, hasta bazlı bölme. Üçüncüsü ve en önemlisi, modelin hiç görmediği üç bağımsız kohortta test ettik ve 0,87–0,95 aldık. Ayrıca dokuda tümör-normal ayrımının kolay bir problem olduğunu biliyoruz; o 0,99'u kan testi başarısı olarak sunmuyoruz.

**S: Çapraz doğrulama nedir?**
C: Veriyi beş parçaya bölüp dördünde eğitip beşincide test etmek, bunu her parça için tekrarlamak. Biz bunu beş kez farklı karıştırmayla tekrarladık (5×5).

**S: MOFA+ nedir, neden kullandınız?**
C: Dört omiği etiket kullanmadan birleştiren bir faktör analizi; verideki ana eksenleri bulur. En güçlü eksenin kanser değil normal doku kaybı olduğunu MOFA+ sayesinde fark ettik ve o tuzaktan kaçtık.

**S: Attention-MIL nedir?**
C: Dikkat mekanizmalı bir derin öğrenme modeli. Her hastanın dört omiğini bir torba gibi görür ve hangi omiğe ne kadar güveneceğini kendisi öğrenir. İlk denemede eksik veri desenini kopya olarak kullandığını yakaladık; adil versiyonda 0,973.

**S: SHAP nedir?**
C: Modelin kararını hangi moleküllerin ne kadar etkilediğini gösteren bir açıklama yöntemi. Sadece test örneklerinde hesapladık ki açıklama da sızıntı içermesin.

**S: Neden derin öğrenme ağırlıklı değil?**
C: Örnek sayısı yüzlerle sınırlı; bu ölçekte basit modeller daha kararlı. Nitekim dört omiği birleştiren derin model tek omikten daha iyi değildi.

### D. Sonuçlar
**S: miRNA paneli neden başarısız oldu?**
C: Kandaki miRNA'nın çoğu kan hücrelerinden gelir; tümörün izi çok küçük kalır. Ayrıca iki açık serum veri seti birbirini bile doğrulamıyor: birinde eğitilen model diğerinde yazı turadan farksız. Yani sorun büyük ölçüde açık serum verilerinin tutarsızlığı.

**S: Başarısız sonucu neden sunuyorsun?**
C: Çünkü sonucu görmeden planı kilitledik; sonuç kötü çıktı diye gizlemek bilimsel olarak yanlış olurdu. Ayrıca bu olumsuz sonuç alana bir uyarı: tek veri setinde çok yüksek görünen kan testlerine temkinli yaklaşılmalı.

**S: Protein paneli kanda da çalışır mı?**
C: Bilmiyoruz; dürüst cevap bu. Lehte işaretler: sekizi de kana salgılanma kanıtı olan proteinler, ikisi (LAMC2, THBS1) başka ekiplerin kan çalışmalarında işe yaramış. Cevabı ancak gerçek plazma ölçümü verecek.

**S: Dış doğrulamayı mRNA ile yaptınız, protein ile değil. Bu sorun değil mi?**
C: Bir sınırlama, kabul ediyorum. Bağımsız kohortlarda protein ölçümü yoktu. mRNA ile protein her zaman birebir örtüşmez; bu yüzden asıl test plazma proteini olacak.

**S: Erken evrede de çalışıyor mu?**
C: Evre bilgisi olan kohortta, yalnızca erken evre tümörlerle AUC 0,86. Yani sinyal erken evrede de var, biraz zayıflıyor.

**S: CA19-9 ile karşılaştırdınız mı?**
C: Protein panelini kanda ölçemediğimiz için doğrudan karşılaştıramadık. Serum veri setinde CA19-9'un duyarlılığı %74, özgüllüğü %100 idi. Pilot çalışmada panel ile CA19-9'u aynı hastalarda kıyaslayacağız; muhtemelen en iyi sonuç ikisinin birlikte kullanımından gelecek.

**S: Pankreatit ile kanseri ayırabiliyor mu?**
C: Henüz test etmedik; en önemli açık soru bu. Hatta panelimizdeki LGALS3BP'nin literatürde pankreatitte yüksek bulunduğunu biliyoruz. Bu yüzden pilotta pankreatit grubu mutlaka olacak; o protein orada elenirse panel 7'ye iner.

**S: p ≤ 0,004 ne demek?**
C: Rastgele seçilmiş 8 genlik panellerin bizim panel kadar iyi çıkma olasılığı binde 4'ten az. Yani sonuç şans eseri değil.

### E. Orijinallik
**S: Bu zaten yapılmadı mı?**
C: Alan çok çalışılıyor, ilk biz değiliz; bir şirketin 4 proteinli paneli bile var. Bizim katkımız: yeni bir 8'li kombinasyon, baştan kanda ölçülebilirlik kuralı ve ön-kayıtlı, dürüst bir doğrulama. Yöntemimizin doğruluğunu da başkalarının kanda doğruladığı iki proteini kendi kendine bulmasıyla gösteriyoruz.

**S: Literatür taraması yaptınız mı?**
C: Evet; LAMC2'nin 2013'te serumda doğrulandığını, THBS1'in Immunovia panelinde olduğunu, LGALS3BP'nin pankreatitte yüksek olduğunu biliyoruz. Bunları hem destek hem risk olarak değerlendirdik.

### F. Uygulama, ticarileşme, etki
**S: Bu ürüne dönüşürse kim kullanır?**
C: Önce yüksek riskli gruplar: 50 yaş üstü yeni diyabet tanısı alanlar, kronik pankreatit hastaları, ailesinde pankreas kanseri olanlar. Gastroenteroloji ve onkoloji klinikleri; ileride tanı kiti üreticileri.

**S: Neden herkese değil de yüksek riskliye?**
C: Pankreas kanseri genel nüfusta nadir; herkesi tararsanız mükemmel bir test bile çok sayıda yanlış alarm verir. Risk yüksek gruplarda test çok daha anlamlı.

**S: Maliyeti ne olur?**
C: ELISA hastanelerde zaten kullanılan, ucuz ve yaygın bir yöntem; özel cihaz gerekmez. Kesin maliyet pilot çalışmayla netleşecek.

**S: Ödülü alırsan ne yapacaksın?**
C: Pilot çalışma için ELISA kitleri ve laboratuvar sarfı. Ayrıca bu sonuçlarla TÜBİTAK gibi daha büyük bir fon başvurusu hazırlamak.

**S: Patent?**
C: Üniversitenin Teknoloji Transfer Ofisi ile değerlendireceğiz. Kod açık kaynak; asıl değer kanda doğrulanmış panel ve ölçüm protokolünde olacak.

**S: Bir sonraki adım ne kadar sürer?**
C: Etik kurul ve hasta toplama dahil bir pilot çalışma tipik olarak 6–12 ay. Hastaneyle iş birliği görüşmesi ilk iş.

### G. Zor / sıkıştırıcı sorular
**S: Kanda çalışmayan bir şeyin neresi başarı?**
C: İki sonucumuz var. miRNA kanda çalışmadı ve bunu gösterdik. Protein paneli bağımsız dokuda güçlü çalıştı ve kanda test edilmeye hazır. Bir AR-GE projesinde hangi yolun çıkmaz, hangisinin umut verici olduğunu net olarak ayırmak başarıdır.

**S: Neden bu kadar az şeyi kanıtlayabildiniz?**
C: Bir öğrenci projesinde, sıfır bütçeyle, açık veriyle gidilebilecek en uzak noktaya gittik: aday bulundu, bağımsız dokuda doğrulandı, ön-kayıtla test edildi. Kalan adım laboratuvar ister.

**S: Yapay zekâya sadece kod mu yazdırdın?**
C: Kodlamada yapay zekâ araçlarından yararlandım; bu artık standart bir mühendislik pratiği. Ama hangi verinin kullanılacağı, hangi tuzakların kontrol edileceği, neyin raporlanacağı gibi kararlar projenin özü ve hepsi gerekçeleriyle karar günlüğünde. Kodun tamamı açık; isteyen tekrar çalıştırabilir. *(Bu soruya dürüst cevap ver; yapay zekâ desteğini inkâr etme.)*

**S: Sonuçları güzel göstermek için seçim yaptınız mı?**
C: Hayır; tam tersi için önlem aldık. Paneller ve analiz planları sonuçtan önce kilitlendi ve GitHub'da zaman damgasıyla duruyor. Kötü çıkan miRNA sonucunu da aynen sunduk.

**S: Demoda neden bazı hastaları yanlış biliyor?**
C: Çünkü gerçek veri ve sabit bir eşik kullanıyoruz; bu kohortta isabet %79–90. Mükemmel değil, o yüzden kanda doğrulama yapacağız.

---

## 16. Sunum günü kontrol listesi

- [ ] GitHub reposunu herkese açık yap (posterdeki QR için).
- [ ] Sunumun yapılacağı bilgisayarda kapanış filminin **sesli** oynadığını dene; yedek MP4'ü USB'de taşı.
- [ ] Demoyu o bilgisayarda aç (`fisilti_demo.html` ve `data.js` aynı klasörde). Sesi kontrol et.
- [ ] Demoda: **E** (GSE16515) → 2–3 kez Boşluk → **3** → Boşluk (yarış). Toplam 1–2 dakika.
- [ ] Posterin basılı halini ve PDF'ini yanında taşı.
- [ ] Su, saat; 8–10 dakikayı aşma.
- [ ] Bilmediğin bir soru gelirse: "Bunu kesin bilmiyorum, uydurmak istemem; not alıp araştırırım." Bu cevap hiçbir zaman puan kaybettirmez.

## 17. Kapanış cümlesi (ezberle)

> "Fısıltıyı duymak için önce gürültüyü ayıklamak gerekir. Dokuda kanserin gerçek izini bulduk ve bağımsız hastanelerde doğruladık; kanda çalışmayanı saklamadık ve nedenini gösterdik. Sekiz proteinlik adayımız şimdi gerçek hastaların kanında test edilmeye hazır."
