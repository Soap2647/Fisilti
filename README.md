# Fısıltı

**Erken evre pankreas kanseri (PDAC Evre I/II) için çoklu-omik, yapay zekâ destekli,
kan tabanlı biyobelirteç paneli keşfi.**

> Pankreas kanseri, vücudun yıllarca fark edilmeyen, neredeyse duyulmayan bir sinyal
> göndermesiyle ilerler. Fısıltı, kandaki bu küçük sinyali yapay zekâ ile yükseltip
> zamanında duyulabilecek bir uyarıya dönüştürmeyi hedefler.

Kastamonu Üniversitesi AR-GE Proje Pazarı ve Tanıtım Günleri kapsamında
geliştirilmektedir. Yürütücü: Mert Yasin Yılmaz.

---

## Veri kaynakları

| Kohort | Rol | Doku | Omik katmanları | Örnek |
|---|---|---|---|---|
| CPTAC-PDAC | Keşif (MOFA+, model) | Tümör + komşu normal | mRNA, protein, miRNA*, metilasyon | 140 T / 92 N (83 erken evre) |
| TCGA-PAAD | Dış doğrulama (erken evre) | Tümör (+10 normal) | mRNA, miRNA, metilasyon, RPPA | 175 T (164 erken evre) |
| GSE15471 | Dış doğrulama | Tümör + eşleşmiş normal | mRNA (U133 Plus 2) | 36 T / 36 N |
| GSE62452 | Dış doğrulama, evre alt-analizi | Tümör + normal | mRNA (HuGene 1.0 ST) | 69 T / 61 N |
| GSE16515 | Dış doğrulama | Tümör + normal | mRNA (U133 Plus 2) | 36 T / 16 N |
| GSE85589 | **Kan doğrulaması**, CA19-9 kıyası | Serum | miRNA | 88 PDAC / 19 sağlıklı / 125 diğer |
| GSE59856 | **Kan doğrulaması**, özgüllük | Serum | miRNA | 100 PDAC / 150 sağlıklı / 21 benign / 300 diğer kanser |

\* CPTAC miRNA, LinkedOmics yerine GDC'den alınır: LinkedOmics normal miRNA tablosunun
tümör verisinin kopyası olduğu tespit edildi (bkz. K-010).

Tüm veriler kamuya açık ve anonimdir; yeni etik kurul onayı gerekmez.
Başvuru dosyasından sapmalar ve yöntemsel kararlar: [docs/KARARLAR.md](docs/KARARLAR.md).

## İş paketleri

| WP | İçerik | Durum |
|---|---|---|
| WP1 | Veri toplama ve küratasyon | ✅ |
| WP2 | MOFA+ çoklu-omik entegrasyon | ✅ |
| WP3 | Attention-MIL + XGBoost/LightGBM, özellik seçimi | ⏳ |
| WP4 | SHAP, pathway zenginleştirme | ⏳ |
| WP5 | Sekresyon/plazma filtresi, panel tasarımı | ⏳ |
| WP6 | Kohortlar arası doğrulama, CA19-9 kıyası, DCA | ⏳ |
| WP7 | Paketleme (Docker, Nextflow), rapor, poster | ⏳ |

## Kurulum

Gereksinim: Python 3.11+ (NVIDIA GPU isteğe bağlı).

```bash
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r env/requirements.txt
.venv\Scripts\python.exe -m pip install torch --index-url https://download.pytorch.org/whl/cu124
.venv\Scripts\python.exe -m pip install -e . --no-deps
```

## Çalıştırma (WP1)

```bash
.venv\Scripts\python.exe scripts/wp1_01_download_cptac.py      # LinkedOmics
.venv\Scripts\python.exe scripts/wp1_02_download_geo.py        # GEO doku + serum
.venv\Scripts\python.exe scripts/wp1_03_download_tcga.py       # GDC (~3,4 GB)
.venv\Scripts\python.exe scripts/wp1_04_mirna_mature_gdc.py    # GDC olgun miRNA
.venv\Scripts\python.exe scripts/wp1_05_preprocess_geo.py
.venv\Scripts\python.exe scripts/wp1_06_preprocess_cptac.py
.venv\Scripts\python.exe scripts/wp1_07_preprocess_tcga.py
.venv\Scripts\python.exe scripts/wp1_08_qc.py                  # özet tablo + PCA
```

İndirilen her dosya kaynak URL'si ve SHA-256 özetiyle `data/raw/MANIFEST.tsv`
dosyasına kaydedilir.

## Tekrarlanabilirlik ilkeleri

- Tüm rastgelelik tek bir tohumdan (`config/config.yaml → project.seed`) türetilir.
- Özellik seçimi ve hiperparametre ayarı yalnızca eğitim katmanlarında yapılır
  (iç içe çapraz doğrulama); dış kohortlar panel sabitlendikten sonra bir kez
  değerlendirilir.
- Veri bölmeleri hasta düzeyindedir; aynı hastanın tümör ve normal örnekleri
  eğitim ve test katmanlarına dağılmaz.

## Proje yapısı

```
config/      merkezi yapılandırma
src/fisilti/ yeniden kullanılabilir Python paketi
scripts/     WP bazlı çalıştırılabilir adımlar (wpN_MM_*.py)
data/        raw → interim → processed (git dışı)
results/     şekiller ve tablolar
docs/        karar günlüğü, rapor taslakları
tests/       birim testleri (pytest)
```

## Lisans

MIT
