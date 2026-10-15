# Fısıltı — Sunum Çalışma Kağıdı

Bu kağıdın amacı ezber değil, **mantığı kavramak**. Önce hikâyeyi anla (Bölüm 1–9), sonra rakamları otur (Bölüm 10), sonra soruları kendi cümlelerinle cevapla (Bölüm 11–14). Sonda kendini test et (Bölüm 15).

> **Altın kural:** Bilmediğin bir şey sorulursa uydurma. "Bunu kesin bilmiyorum, uydurmak istemem; not alıp araştırırım" de. Bu cevap puan kaybettirmez; uydurduğun bir rakam yakalanırsa kaybettirir.

---

## 0. Tek sayfalık özet (sunumdan 5 dakika önce bunu oku)

**Projenin özü:**
> "Pankreas kanseri sessiz ilerler; çoğu hasta ameliyat edilemeyecek evrede tanı alır. Açık paylaşılmış yedi hasta veri setini yapay zekâ ile inceleyip, kanda ölçülebilecek 8 proteinlik bir aday panel bulduk. Panel, hiç görmediği üç bağımsız hastane verisinde tümörü AUC 0,87–0,95 ile ayırdı; erken evrede 0,86. Kanda çalışmayan miRNA panelimizi saklamadık. Sıradaki adım, 8 proteini gerçek hasta kanında ölçmek."

**Mantık zinciri** (kaybolursan buna dön):

**Geç fark ediliyor** → **erken yakalamak lazım** → **en kolay yol kan testi** → **tümör kana iz bırakır** → **binlerce molekül var, hangisi?** → **açık veri + yapay zekâ ile aday seçtik** → **sadece kanda ölçülebilenleri tuttuk** → **kopya çekmesin diye sıkı test ettik** → **miRNA çalışmadı, saklamadık** → **8 protein bağımsız hastanelerde çalıştı** → **sıradaki adım gerçek kanda ölçmek**

**Beş altın rakam:**

| Rakam | Anlamı |
|---|---|
| **0,87 – 0,95** | 8 protein, hiç görmediği 3 hastane verisinde (AUC) |
| **0,86** | Sadece erken evre tümörlerde (AUC) |
| **0 / 2000** | 2000 rastgele panelin hiçbiri bizimki kadar iyi değil |
| **0,54 ve 0,40** | miRNA paneli kanda (yazı tura) |
| **0 TL** | Bugüne kadarki maliyet |

---

## 1. Dert ne?

- Pankreas karnın en derininde, midenin arkasında. Sindirim enzimi ve insülin üretir.
- Orada çıkan kanser (**PDAC**, pankreas kanserlerinin yaklaşık %90'ı) **ses çıkarmaz**: ağrı yok, şişlik yok.
- Belirti çıktığında (sarılık, kilo kaybı, sırt ağrısı) genelde yayılmıştır; çoğu hasta ameliyat edilemez.
- **Beş yıllık sağkalım %10'un altında** (Klein AP, Nat Rev Gastroenterol Hepatol 2021).

**Mantık:** Bu kanseri öldürücü yapan, hastalık kadar **geç fark edilmesi**. Erken yakalanırsa ameliyatla kurtarma şansı var.

**Neden kan testi?** Herkesi MR'a sokamayız: pahalı, yavaş, sırası uzun. Kan tahlili ucuz, hızlı, her yerde var.

**Mevcut kan testi neden yetmiyor?** Tek belirteç **CA19-9**:

- Erken evrede çoğu zaman yükselmez.
- Sarılık ve pankreatit gibi iyi huylu durumlarda da yükselir → yanlış alarm.
- İnsanların yaklaşık %5–10'u bu molekülü hiç üretemez.
- Bu yüzden tarama testi olarak önerilmez, teşhisten sonra takipte kullanılır.

👉 **Amacımız:** Erken pankreas kanserini kandan yakalayabilecek yeni bir test adayı bulmak.

---

## 2. Fikir: tümör kana iz bırakır

- Hırsız eve girince ayak izi bırakır. Tümör de kana bazı moleküller bırakır → **biyobelirteç**.
- Projenin adı buradan: tümör kana bir **fısıltı** bırakır, biz onu duymaya çalışıyoruz.
- **Sorun:** Vücutta binlerce molekül var. Hangileri tümörün izi? Laboratuvarda tek tek denemek yıllar sürer, çok pahalıdır.

---

## 3. Çözüm: açık veri + yapay zekâ

- Dünyada pankreas kanseri hastalarından alınmış dokular ölçülmüş ve **herkese açık, anonim** olarak paylaşılmış. Bazıları tümör, bazıları sağlıklı doku.
- Yapay zekâya sorduğumuz soru: **"Bu iki grubu ayırırken en çok hangi moleküllere bakıyorsun?"**
- Benzetme: Dedektif binlerce şüpheli arasından en güçlü delili olanları seçiyor. Laboratuvara girmeden aday listesini daraltıyoruz.
- Hücreye **dört kamerayla** baktık:

| Kamera | Sade anlamı |
|---|---|
| **mRNA** | Hangi gen ne kadar çalışıyor (genin "tarif kopyası") |
| **Protein** | Tariflerle üretilen, iş gören moleküller. Kanda ELISA ile ölçülen şey |
| **miRNA** | Diğer genleri kısan küçük RNA'lar. Kanda kararlı, RT-qPCR ile ölçülür |
| **Metilasyon** | DNA üzerindeki aç-kapa düğmeleri |

Birden fazla katmana birlikte bakmaya **çoklu-omik** denir.

### Kullandığımız yedi veri seti (kohort)

| Kohort | Rol | Örnek |
|---|---|---|
| **CPTAC-PDAC** | Keşif (modeller burada eğitildi), 4 omik aynı hastalarda | 140 tümör · 92 normal |
| **TCGA-PAAD** | Doku doğrulama | 175 tümör · 10 normal |
| **GSE15471** | Bağımsız doku doğrulama | 36 tümör · 36 normal |
| **GSE62452** | Bağımsız doku doğrulama (evre bilgisi var) | 69 tümör · 61 normal |
| **GSE16515** | Bağımsız doku doğrulama | 36 tümör · 16 normal |
| **GSE85589** | Kan (serum) miRNA | 88 PDAC · 19 sağlıklı · 125 diğer |
| **GSE59856** | Kan (serum) miRNA | 100 PDAC · 150 sağlıklı · 321 diğer |

- **CPTAC**: ABD Ulusal Kanser Enstitüsü'nün protein projesi. **TCGA**: Kanser Genom Atlası. **GEO**: ABD'nin açık gen verisi arşivi (GSE numaraları).
- Veriler anonim ve açık olduğu için **etik kurul gerekmedi**.
- Başvuruda geçen **ICGC** portalı 2024'te kapandı; yerine CPTAC kullanıldı (karar günlüğünde kayıtlı).

---

## 4. Yöntemin beş adımı

1. **Veri:** 7 kohort indirildi, her dosyanın parmak izi (SHA-256) kaydedildi; herkes aynı veriyi indirip aynı sonucu alabilir.
2. **MOFA+:** Dört omiği, etiket kullanmadan birleştiren bir yöntem. Verideki ana "eksenleri" (13 faktör) bulur.
3. **Yapay zekâ modelleri:**
   - **L1-lojistik regresyon:** basit, açıklanabilir; gereksiz molekülleri sıfırlar.
   - **XGBoost, LightGBM:** karar ağacı tabanlı güçlü modeller.
   - **Attention-MIL:** dikkat mekanizmalı derin öğrenme; bir hastanın dört omiğine birlikte bakar ve hangisine ne kadar güveneceğini öğrenir.
   - Ayarlar **Optuna** ile yapıldı. Toplam yaklaşık 600 model eğitildi.
   - **Neden birden fazla model?** Farklı varsayımlarla çalışan modeller aynı molekülleri seçiyorsa, sonuç tek bir modelin tuhaflığı değildir.
4. **Yorum:** **SHAP** (hangi molekül karara ne kadar katkı verdi) ve **yolak analizi** (seçilen genler hangi biyolojik süreçlere ait).
5. **Panel ve doğrulama:** Adaylar filtrelendi, panel kilitlendi, bağımsız verilerde test edildi.

---

## 5. En akıllı kural: kanda ölçülebilir olmalı

Tümörde harika bir işaret, kana hiç çıkmıyorsa kan testinde işe yaramaz. Üç filtre uyguladık:

1. **Tümörde artmalı.** Azalan bir şeyi kanda yakalamak zor.
2. **Kana salgılanmalı.** Kanıtı Human Protein Atlas ve UniProt veri tabanlarından aldık.
3. **Hastanede ölçülebilmeli.** ELISA veya RT-qPCR ile.

Ek olarak kan hücrelerinden gelen molekülleri çıkardık; tüp sallanınca hücreler patlarsa (hemoliz) yanıltırlar. Birbirinin kopyası olanları da attık.

**Sonuç:**

- **Protein paneli (8):** SULF1, LGALS3BP, NMI, LAMC2, THBS1, CD55, GLA, SERPINB5. Eleme: 39 → 31 → 8.
- **miRNA paneli (5):** miR-21-5p, miR-210-3p, miR-27a-3p, miR-935, miR-24-3p. Eleme: 40 → 15 → 14 → 5.

---

## 6. Kendimizi kandırmamak (projenin kalbi)

Yapay zekânın en büyük tehlikesi **kopya çekmek**. Soruları önceden gören öğrenci sınavdan 100 alır ama bir şey öğrenmemiştir. Bilimdeki adı **sızıntı (leakage)**.

**Üç önlem:**

**a) Sınav soruları ayrı tutuldu (çapraz doğrulama).** Veri 5 parçaya bölündü: 4'ünde eğitildi, 1'inde test edildi, her parça sırayla test oldu. Bu iş 5 kez farklı karıştırmayla tekrarlandı (5×5). Ayarlar ayrı bir iç katmanda yapıldı (**iç içe**), yani test verisi ayar için hiç kullanılmadı. Aynı hastanın tümörü ve normali hep aynı tarafta kaldı (**hasta bazlı bölme**).

**b) Başka okulun sınavına sokuldu (bağımsız doğrulama).** Model bir veri setinde eğitildi, sonra hiç görmediği, başka ülkelerdeki başka hastanelerin verisinde test edildi. **Gerçek başarı buradaki sonuçtur.**

**c) Cevaplar zarfa mühürlendi (ön-kayıt).** Panel ve test planı, **sonuçlar görülmeden** kilitlendi ve GitHub'a tarihiyle kaydedildi. "Sonuç kötü çıktı, biraz değiştireyim" demek imkânsız hale geldi.

**30 Eylül 2026'daki zaman çizelgesi:** 18:20 miRNA paneli kilitlendi → 18:26 hata düzeltmesi → 19:25 kan analiz planı kaydedildi → 19:28 kan sonucu → 19:36 doku doğrulama planı kaydedildi → sonra doku sonucu.

### Yakaladığımız altı gizli veri hatası

1. Bir veri setinde "sağlıklı" diye etiketli miRNA dosyaları **tümörün birebir kopyasıydı** (21/21, korelasyon 1,0000). Fark etmeseydik model kopyayı öğrenirdi. Ham veriden doğru dosyalar alındı.
2. TCGA'da metastaz kaydı birincil tanının önündeydi; 35 hastanın evresi yanlış okunuyordu. Düzeltilince eksik evre 35'ten 2'ye indi.
3. Aynı hastanın iki örneği eğitime ve teste ayrı düşebilirdi → hasta bazlı bölme.
4. Normal örneklerde bazı omikler hep eksikti; model **eksiklik desenini kopya olarak** kullanıp 0,986'ya çıktı. Adil sürümde 0,973.
5. Panel seçiminde yanlış bir benzerlik hesabı, iyi salgılanan adayları elemişti. **Kan verisine bakmadan** düzeltildi.
6. İki veri setinde ölçek ve tekrar hataları vardı; düzeltildi.

Hepsi gerekçeleriyle `docs/KARARLAR.md` dosyasında yazılı (22 karar).

### MOFA+ uyarısı: en güçlü sinyal bir tuzaktı
Tümörü normalden en iyi ayıran eksen (AUC 0,99) aslında **kanser değil, normal pankreas dokusunun kaybıydı**: tümör büyüdükçe enzim üreten doku azalıyor. Bu kanda ölçülemez ve pankreatitte de olur, yani yanlış alarm üretirdi. Bu yüzden sadece **tümörde artan** moleküller seçildi.

**Mesaj:** Yüksek skor her zaman doğru biyoloji demek değil; biz modelin neyi öğrendiğine baktık.

---

## 7. Sonuçlar: iki haber

### Önce ölçüyü bil

- **AUC:** Modelin önüne bir hasta bir de sağlıklı örnek koy; hangisinin hasta olduğunu ne sıklıkla doğru biliyor? **0,5 = yazı tura, 1,0 = kusursuz**, 0,8 iyi, 0,9 çok iyi.
- **Duyarlılık:** Hastaların yüzde kaçını yakalıyor.
- **Özgüllük:** Sağlıklıların yüzde kaçına "temiz" diyor, yani yanlış alarm vermiyor.

### Dokuda yapay zekâ (keşif verisi içinde)

| Katman | AUC |
|---|---|
| Protein | 0,99 |
| mRNA | 0,99 |
| Metilasyon | 0,99 |
| miRNA | 0,98 |

- Dört katmanı birleştirmek tek katmandan daha iyi değil (0,973) → **tek, ucuz bir test türü yeterli olabilir.**
- **Dikkat:** Dokuda tümörü normalden ayırmak kolay bir problem. Bu 0,99'u kan testi başarısı diye sunma; jüri söylemeden sen söyle.

### Model bilinen biyolojiyi kendi kendine buldu

- **KRAS** sinyali: PDAC'ların %90'ından fazlasında bozuk olan gen.
- **Epitel-mezenkimal geçiş:** tümörün sertleşen bağ dokusu.
- **p53 yolağı.**
- Literatürde bilinen belirteçler, biz söylemeden bulundu: miR-21, miR-196a, SERPINB5, TFPI2.
- **Mantık:** Model saçma bir şey öğrenseydi bunları bulamazdı.

### 😞 Kötü haber: miRNA paneli kanda çalışmadı

- Gerçek kan verisinde AUC **0,54 ve 0,40**, yani yazı tura.
- Rastgele seçilmiş 5 miRNA'dan bile iyi değil.
- **Neden?**
  - Kandaki miRNA'ların çoğu tümörden değil kan hücrelerinden geliyor; tümörün sesi kalabalıkta kayboluyor.
  - İki açık kan veri seti kendi içinde "mükemmel" (AUC ≈ 1,00) ama **birbirini tutmuyor**: birinde öğrenilen diğerinde 0,40 / 0,48. Yani o "mükemmel" sonuçlar büyük ihtimalle laboratuvar farkından.
  - Dokuda artan miRNA kanda artmıyor.
- **Mantık:** Neyin çalışmadığını bilmek de bilgidir. Bu aynı zamanda alana bir uyarı: tek bir veri setindeki çok yüksek kan testi sonuçlarına güvenmeyin.

### 😊 İyi haber: 8'li protein paneli tuttu

| Bağımsız kohort | AUC | Sabit eşikle isabet |
|---|---|---|
| GSE15471 | **0,92** | %83 |
| GSE62452 | **0,87** | %79 |
| GSE16515 | **0,95** | %90 |

- Sadece erken evre tümörlerde: **0,86**.
- 24 gen-kohort ölçümünün 24'ünde de değişim beklenen yönde.
- **2000 rastgele 8'li panelin hiçbiri** bu kadar iyi olamadı (p ≤ 0,004) → sonuç şans eseri değil.

### Dürüst sınır

- Protein paneli **dokuda** doğrulandı, **gerçek hasta kanında henüz ölçülmedi.** Bu proteinlerin kanda ölçüldüğü açık bir veri yok.
- Bağımsız doğrulama mRNA ölçümüyle yapıldı; mRNA ile protein her zaman birebir örtüşmez.

> **Benzetme:** "Şüphelileri bulduk, deliller güçlü. Ama mahkemeye, yani gerçek hasta kanına henüz çıkmadılar."

---

## 8. Sekiz proteinin kimlik kartı

| Protein | Ne iş yapar | Neden aday | AUC aralığı |
|---|---|---|---|
| **SULF1** | Hücre dışındaki şekerleri düzenleyen salgılanan enzim | Tümörün bağ dokusunda yüksek | 0,84–0,97 |
| **LGALS3BP** | Salgılanan bir glikoprotein | Kana salgılanır | 0,84–0,89 |
| **NMI** | Bağışıklık (interferon) yanıtında rol alır | Kanda tespit edilmiş | 0,75–0,88 |
| **LAMC2** | Tümörün yayıldığı uçta bulunan laminin | Başka ekipler serumda doğrulamış | 0,91–0,95 |
| **THBS1** | Damar oluşumu ve bağ dokusu | Bir şirketin panelinde de var | 0,59–0,83 |
| **CD55** | Bağışıklık sisteminin kompleman kolunu düzenler | Kana salgılanır | 0,83–0,93 |
| **GLA** | Hücre içi bir enzim | Plazmada ölçülmüş | 0,64–0,77 |
| **SERPINB5** | Maspin; PDAC'ta belirgin artar | Bilinen PDAC belirteci | 0,82–0,94 |

---

## 9. Sıradaki adım, bütçe ve etki

**Yol haritası:**

1. **Pilot (planlanan):** Hastaneden kan örnekleri (PDAC, **pankreatit** ve sağlıklı kişiler), 8 proteinin ELISA ile ölçümü, CA19-9 ile karşılaştırma. Etik kurul ve hastane iş birliği gerekiyor; hastaneyle henüz anlaşma yok, "planlanan" de.
2. **Çok merkezli, ileriye dönük çalışma:** erken evre odaklı, ön-kayıtlı.
3. **Fon:** TÜBİTAK destekleri. Aynı yöntem akciğer, kolorektal ve mide kanserine de uyarlanabilir.
4. **Fikri mülkiyet:** Üniversitenin Teknoloji Transfer Ofisi ile değerlendirme.

**Bütçe:**

- Bugüne kadarki maliyet **0 TL**: açık veri, açık kaynak yazılım, kişisel bilgisayar.
- **Ödül ne için?** Pilot için ELISA kitleri ve laboratuvar sarfı. Kesin rakam sorulursa uydurma: "Kesin bütçe hastane ve laboratuvarla birlikte çıkacak, ilk kalem sekiz kit ve sarf malzemesi."

**Kim kullanır?** Önce yüksek riskli gruplar:

- 50 yaş üstü, yeni diyabet tanısı alanlar (pankreas kanseri bazen ilk olarak diyabetle kendini gösterir)
- Kronik pankreatit hastaları
- Ailesinde pankreas kanseri olanlar ve BRCA gibi kalıtsal gen değişikliği taşıyanlar
- Pankreas kisti (IPMN) takibindeki hastalar

**Neden herkese değil?** Pankreas kanseri genel nüfusta nadir. Herkesi tararsan en iyi test bile çok fazla yanlış alarm verir. Risk yüksek grupta pozitif sonuç çok daha anlamlıdır.

---

## 10. Literatürde neredeyiz? (orijinallik)

- **Kanla erken pankreas kanseri yakalama, çok çalışılan bir alan.** İlk biz değiliz; bunu açıkça söyle.
- **Immunovia** (şirket, 2025): 4 protein + CA19-9 paneli; 128 erken evre hasta ve 465 yüksek riskli kişide, sağlamların %98'ine doğru "temiz" derken hastaların %85'ini yakalamış. **THBS1 bizim panelimizde de var.**
- **LAMC2:** 2013'te doku analiziyle bulunup 425 kişilik serumda doğrulanmış; 2014'te CA19-9 + CA125 + LAMC2 paneli yayımlanmış.
- **LGALS3BP:** Bir 2016 çalışmasında kronik pankreatitte kanserden daha yüksek bulunmuş → **zayıf noktamız**, pilotta test edilecek.

**Bizim özgün katkımız:**

1. **Bu 8'li kombinasyon** literatürde yayımlanmış görünmüyor.
2. Yöntemimiz, başkalarının kanda doğruladığı **iki proteini (LAMC2, THBS1) kendi kendine buldu.** Bu, doğru yere baktığımızı gösteriyor; diğer 6'sı yeni aday.
3. **Kanda ölçülebilirlik baştan bir kural.**
4. Ön-kayıt, hasta bazlı test, başarısız sonucu raporlama ve iki açık kan veri setinin birbirini tutmadığını göstermek.
5. Tamamen açık kaynak; herkes tekrar çalıştırabilir.

> **Söyleyiş:** "Bu alanda ilk biz değiliz; şirketler bile var. Bizim farkımız: açık verilerden, baştan kanda ölçülebilirlik şartıyla yeni bir 8'li kombinasyon bulduk. Yöntemimizin doğru çalıştığını da şuradan biliyoruz: başkalarının kanda doğruladığı iki proteini kendiliğinden seçti."

---

## 11. Doktor sorarsa: "Bunu nerede, nasıl kullanacağız?"

> "Herkese uygulanan bir tarama testi olarak düşünmüyoruz. Riski yüksek hastalarda, görüntülemeden önce kullanılan bir **kan ön-elemesi** olarak düşünüyoruz. 50 yaşından sonra yeni diyabet tanısı alan, kronik pankreatiti olan ya da ailesinde pankreas kanseri bulunan hastada rutin kan tahliliyle birlikte panel bakılır. Skor yüksekse hasta öncelikli olarak endoskopik ultrason ya da MR'a yönlendirilir, düşükse rutin takibe devam edilir. Görüntülemenin ve biyopsinin yerini almıyor; **kime ve ne zaman görüntüleme yapılacağına karar vermeye yardım ediyor.** CA19-9'un yerine değil, onunla birlikte kullanılır."

**Hasta yolculuğu:** Yüksek riskli hasta → kan tahlili + panel + CA19-9 → skor yüksekse öncelikli EUS/MR, düşükse 6–12 ayda tekrar → lezyon varsa biyopsi → erken evrede cerrahi.

**Hekimin gözünde avantajları:**

- ELISA her hastanede var, aynı tüpten bakılır, özel cihaz gerekmez.
- Ucuz ve tekrarlanabilir.
- Pahalı ve sırası uzun görüntülemelerde "kimin sırası önce gelmeli" sorusuna cevap verir.
- CA19-9'u üretemeyen hastalarda ek bilgi verebilir.

**Üstüne giderse:**

- *Pozitif çıkana ne diyeceksiniz?* → "Kanser tanısı koymuyoruz; 'görüntülemeyi öne alın' diyoruz. Tanıyı görüntüleme ve biyopsi koyar."
- *Kanseri kaçırırsa?* → "Bu hastalar zaten düzenli takipte; panel takibin yerine geçmiyor, takibe ek bilgi veriyor."
- *Yanlış pozitif hastayı korkutmaz mı?* → "Sonucu tanı değil öncelik olarak veriyoruz ve sadece zaten takip gereken gruplarda kullanıyoruz."
- *Ne zaman kullanılabilir?* → "Kanda pilot, çok merkezli doğrulama ve onay süreçleri; birkaç yıllık bir yol. Bugün sunduğumuz şey bu yolun ilk ve en ucuz adımı."

**Doktor karşısında asla "kanseri teşhis ediyoruz" deme.** De ki: "Kime öncelikli görüntüleme yapılacağına karar vermeye yardım eden bir kan testi adayı."

---

## 12. Zayıf yönler (önce sen söyle)

Jüri zayıf yanı senden önce bulursa güven kaybedersin. Sen söylersen olgunluk olarak görür.

1. **Kanda henüz doğrulama yok.**
2. **miRNA paneli kanda başarısız oldu.**
3. **Bağımsız doğrulama mRNA ile yapıldı**, protein ile değil.
4. **Kontrol grubu çoğunlukla normal doku.** Asıl zor ayrım kanser ve pankreatit. LGALS3BP pankreatitte yüksek olabilir; elenirse panel 7'ye iner.
5. **Örnek sayıları orta düzeyde**; bazı kohortlarda evre bilgisi yok.
6. **Geriye dönük veri**; ileriye dönük çalışma gerekecek.

---

## 13. Jüri soruları: kısa cevaplar

### Genel

- **Bir cümlede proje?** → Açık verilerden yapay zekâ ile, erken pankreas kanserini kanda yakalayabilecek 8 proteinlik bir aday panel bulduk ve bağımsız hastane verilerinde doğruladık; sıradaki adım gerçek kanda ölçmek.
- **Neden pankreas?** → En ölümcül kanserlerden biri, sessiz ilerliyor, erken yakalamanın etkisi çok büyük ve iyi bir tarama testi yok.
- **Bu bir kan testi mi?** → Henüz değil; aday panel. Kanda ölçüm yol haritamızın ilk adımı.
- **Bilgisayar mühendisi olarak neden bu konu?** → Bu problemin darboğazı veri analizi ve yapay zekâ. Biyolojik kararlarda literatüre dayandım, laboratuvar kısmı için hekim ve biyolog iş birliği planlıyorum.

### Veri

- **Veriyi nereden buldunuz, etik kurul?** → CPTAC, TCGA, GEO: açık ve anonim veriler; kimlik bilgisi yok, etik kurul gerekmedi. Kan pilotu için alınacak.
- **Veri kalitesinden emin misiniz?** → Değildik, kontrol ettik ve altı ciddi hata bulduk. Örneğin "sağlıklı" etiketli dosyalar tümörün kopyasıydı.
- **Örnek sayısı az değil mi?** → Keşifte 232 örnek, bağımsız doğrulamada 254. Biyobelirteç keşfi için makul, klinik karar için yetmez; bu yüzden çok merkezli çalışma planlıyoruz.

### Yöntem

- **0,99 çok yüksek, ezber olmasın?** → İç içe ve hasta bazlı doğrulama yaptık; en önemlisi hiç görmediği üç hastanede test ettik ve 0,87–0,95 aldık. Ayrıca dokuda ayrımın kolay olduğunu biliyoruz; 0,99'u kan testi başarısı diye sunmuyoruz.
- **Neden birden fazla model?** → Farklı modeller aynı molekülleri seçiyorsa sonuç bir modelin tuhaflığı değildir. En iyi protein sonucunu en basit model verdi; bu açıklanabilirlik için iyi.
- **Neden derin öğrenme ağırlıklı değil?** → Örnek sayısı yüzlerle sınırlı; bu ölçekte basit modeller daha kararlı. Derin model tek omikten daha iyi çıkmadı.
- **MOFA+ ne işe yaradı?** → En güçlü eksenin kanser değil normal doku kaybı olduğunu gösterdi; o tuzağa düşmedik.
- **SHAP nedir?** → Kararı hangi molekülün ne kadar etkilediğini gösterir; sadece test örneklerinde hesapladık.
- **p ≤ 0,004 ne demek?** → Rastgele bir 8'li panelin bizimki kadar iyi çıkma ihtimali binde 4'ten az.

### Sonuç

- **miRNA neden çalışmadı?** → Kandaki miRNA'nın çoğu kan hücrelerinden gelir. Ayrıca iki açık kan veri seti birbirini bile tutmuyor; sorun büyük ölçüde verinin tutarsızlığı.
- **Başarısız sonucu neden gösteriyorsun?** → Planı sonuçtan önce kilitledik; kötü çıkanı gizlemek bilimsel olarak yanlış olurdu. Ayrıca alana bir uyarı.
- **Protein kanda da çalışır mı?** → Bilmiyoruz, dürüst cevap bu. Lehte işaretler: sekizi de kana salgılanıyor, ikisi başka ekiplerin kan çalışmalarında işe yaramış.
- **Erken evrede de çalışıyor mu?** → Evre bilgisi olan kohortta sadece erken evre tümörlerle 0,86.
- **CA19-9 ile kıyasladınız mı?** → Paneli kanda ölçemediğimiz için doğrudan kıyaslayamadık; pilotta aynı hastalarda kıyaslayacağız. Büyük ihtimalle en iyi sonuç ikisi birlikte.
- **Pankreatitten ayırır mı?** → Henüz test etmedik; en önemli açık soru bu, pilotta pankreatit grubu olacak.

### Orijinallik ve etki

- **Bu zaten yapılmadı mı?** → Bölüm 10'daki söyleyiş.
- **Kim kullanır?** → Yüksek riskli gruplar; gastroenteroloji ve onkoloji klinikleri; ileride tanı kiti üreticileri.
- **Maliyeti?** → ELISA ucuz ve yaygın; kesin maliyet pilotla netleşecek.
- **Ödülü alırsan?** → Pilot için ELISA kitleri ve sarf; sonra TÜBİTAK gibi büyük bir fon başvurusu.
- **Patent?** → Üniversitenin Teknoloji Transfer Ofisi ile değerlendireceğiz. Kod açık; asıl değer kanda doğrulanmış panelde ve ölçüm protokolünde olacak.

### Zor sorular

- **Kanda çalışmayan bir şeyin neresi başarı?** → İki sonuç var: miRNA'nın çıkmaz olduğunu gösterdik, protein panelinin umut verici olduğunu gösterdik. AR-GE'de doğru yolu yanlıştan ayırmak başarıdır.
- **Neden bu kadar az şey kanıtlandı?** → Öğrenci projesi, sıfır bütçe, açık veri: gidilebilecek en uzak noktaya gittik. Kalan adım laboratuvar ister.
- **Yapay zekâya kod mu yazdırdın?** → Kodlamada yapay zekâ araçlarından yararlandım; bu artık standart bir mühendislik pratiği. Ama hangi verinin kullanılacağı, hangi tuzağın kontrol edileceği, neyin raporlanacağı kararları projenin özü ve hepsi gerekçeleriyle karar günlüğünde. Kod açık, isteyen tekrar çalıştırabilir. *(Dürüst cevap ver, inkâr etme.)*
- **Sonuçları güzel göstermek için seçim yaptınız mı?** → Tam tersi için önlem aldık: planlar sonuçtan önce kilitlendi, GitHub'da tarihiyle duruyor. Kötü çıkan sonucu da aynen sunduk.
- **Demoda neden bazı hastaları yanlış biliyor?** → Gerçek veri ve sabit bir eşik kullanıyoruz; isabet %79–90. Mükemmel değil, o yüzden kanda doğrulama yapacağız.

---

## 14. Demo ve film nasıl sunulur

**Demo** (`sunum/demo/fisilti_demo.html`, internetsiz çalışır):

1. **E** tuşu → GSE16515 kohortu (isabet %90).
2. **Boşluk** 2–3 kez → tek hasta testleri. Söyle: "Bu hasta örnekleri modelin hiç görmediği bir hastaneden."
3. **3** → **Boşluk** → rastgele panellerle yarış. Söyle: "2000 rastgele panelin hiçbiri bizimki kadar iyi değil; sonuç şans eseri değil."
4. Yanlış bir tahmin çıkarsa: "Mükemmel değil; bu yüzden kanda doğrulama yapacağız."
- Toplam 1–2 dakika. En etkileyici an yarış.

**Kapanış filmi:** 74 saniye, sesli, son slaytta kendiliğinden başlar. Geçiş cümlesi: "Son olarak bu yolculuğu bir dakikada göstermek istiyorum."

---

## 15. Kendini test et (sesli cevapla)

Bu soruları kağıda bakmadan, kendi cümlelerinle cevaplayabiliyorsan hazırsın:

1. Bu kanser neden bu kadar ölümcül?
2. Neden kan testi, neden CA19-9 yetmiyor?
3. Yapay zekâ bu işte tam olarak ne yaptı?
4. Binlerce molekül yerine neden sadece 8 protein?
5. Modelin kopya çekmediğini nasıl biliyoruz? (Üç önlem)
6. Ön-kayıt nedir, neden önemli?
7. MOFA+ hangi tuzağı gösterdi?
8. miRNA neden kanda çalışmadı?
9. Protein panelinin en güçlü üç kanıtı ne?
10. Bu şu an bir kan testi mi? Değilse neden?
11. Doktor "nerede kullanacağız?" derse ne dersin?
12. "Bu zaten yapılmadı mı?" sorusuna cevabın ne?
13. Projenin en büyük zayıf yanı ne?
14. Ödülü alırsan ne yapacaksın?

---

## 16. Sunum günü kontrol listesi

- ☐ GitHub reposunu herkese açık yap (posterdeki QR için).
- ☐ Sunum bilgisayarında kapanış filmini **sesli** dene; yedek MP4'ü USB'de taşı.
- ☐ Demoyu o bilgisayarda aç (`fisilti_demo.html` ve `data.js` aynı klasörde); sesi kontrol et.
- ☐ Poster basılı hali ve PDF'i yanında.
- ☐ Süre 8–10 dakika; aşma.
- ☐ Su iç, yavaş konuş, jüriye bak. Rakamlar senin elinde.

## Kapanış cümlesi

> "Fısıltıyı duymak için önce gürültüyü ayıklamak gerekir. Dokuda kanserin gerçek izini bulduk ve bağımsız hastanelerde doğruladık; kanda çalışmayanı saklamadık ve nedenini gösterdik. Sekiz proteinlik adayımız şimdi gerçek hastaların kanında test edilmeye hazır."
