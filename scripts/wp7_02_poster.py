"""WP7-02: Kastamonu AR-GE poster şablonunu (70×100 cm) proje sonuçlarıyla doldurur.

Girdi : sunum/poster/sablon.pptx (kirmizi-kastamonu-poster-sablon-1.ppt'nin .pptx kopyası)
Çıktı : sunum/poster/Fisilti_Poster.pptx + poster grafikleri (sunum/poster/img/)
Şablon kuralları: Times New Roman, gövde 28 pt, iki yana yaslı, tek satır aralığı; başlık 40 pt kalın büyük harf.
Tüm sayılar results/ tablolarından okunur. Bilinmeyen bilgiler PLACEHOLDER sözlüğünde, kırmızı köşeli parantezle.

Çalıştırma:  .venv\\Scripts\\python.exe scripts/wp7_02_poster.py
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import qrcode
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Pt
from sklearn.metrics import roc_curve

from fisilti.config import PROJECT_ROOT, path

POSTER = PROJECT_ROOT / "sunum" / "poster"
IMG = POSTER / "img"
IMG.mkdir(parents=True, exist_ok=True)
FONT = "Times New Roman"
RED = RGBColor(0xC0, 0x00, 0x00)
INK = RGBColor(0x1A, 0x1A, 0x1A)
GREY = RGBColor(0x55, 0x55, 0x55)
NAVY = RGBColor(0x1F, 0x3A, 0x5F)
PH = RGBColor(0xD0, 0x00, 0x00)  # yer tutucu rengi

# Kullanıcıdan teyit bekleyen bilgiler (None -> kırmızı yer tutucu basılır)
INFO = {
    "bolum": None,          # ör. "Mühendislik ve Mimarlık Fakültesi, Bilgisayar Mühendisliği"
    "danisman": None,       # ör. "Dr. Öğr. Üyesi Ad SOYAD"
    "tematik": None,        # ör. "Sağlık / Biyoinformatik ve Yapay Zekâ"
    "kategori": "Öğrenci",
    "proje_kodu": None,
    "eposta": None,
}

plt.rcParams.update({"font.family": "Times New Roman", "font.size": 13, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.facecolor": "white", "axes.facecolor": "white"})
MRED, MNAVY, MGREY = "#C00000", "#1F3A5F", "#8A8A8A"


# ----------------------------------------------------------------------------- grafikler
def fig_pipeline(file):
    steps = [("1", "Veri", "7 kohort · 4 omik"), ("2", "MOFA+", "çoklu-omik faktörler"),
             ("3", "Yapay zekâ", "600 model · iç içe CV"), ("4", "Yorum", "SHAP · GSEA"),
             ("5", "Doğrulama", "kilitli panel · ön-kayıt")]
    fig, ax = plt.subplots(figsize=(11.7, 2.35), dpi=300)
    ax.set_xlim(0, 100); ax.set_ylim(0, 20); ax.axis("off")
    w, gap = 17.2, 3.5
    for k, (n, t1, t2) in enumerate(steps):
        x = k * (w + gap)
        ax.add_patch(FancyBboxPatch((x, 1.5), w, 17, boxstyle="round,pad=0,rounding_size=1.8",
                                    fc="#FBEAEA" if k != 4 else "#C00000", ec=MRED, lw=1.4))
        c = "white" if k == 4 else "#1A1A1A"
        ax.text(x + w / 2, 14.2, n, ha="center", va="center", fontsize=17, fontweight="bold", color="white" if k == 4 else MRED)
        ax.text(x + w / 2, 9.6, t1, ha="center", va="center", fontsize=15, fontweight="bold", color=c)
        ax.text(x + w / 2, 5.0, t2, ha="center", va="center", fontsize=11.5, color=c)
        if k < 4:
            ax.add_patch(FancyArrowPatch((x + w + 0.4, 10), (x + w + gap - 0.4, 10), arrowstyle="-|>", mutation_scale=14, color=MRED, lw=1.4))
    fig.savefig(file, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def fig_external_roc(file):
    lock = json.loads((PROJECT_ROOT / "results/panel/protein_panel_v2.json").read_text(encoding="utf-8"))
    genes = [m["gene"] for m in lock["members"]]
    res = pd.read_csv(path("tables") / "wp6_tissue_external.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.1), dpi=300)
    for ax, gse in zip(axes, ["GSE15471", "GSE62452", "GSE16515"]):
        e = pd.read_parquet(path("processed") / "geo" / f"{gse}_expr.parquet")
        l = pd.read_csv(path("processed") / "geo" / f"{gse}_labels.tsv", sep="\t", index_col=0).loc[e.columns]
        z = e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1).replace(0, np.nan), axis=0)
        s = z.loc[genes].mean(axis=0)
        r = res[(res.cohort == gse) & (res.comparison == "T1_PDAC_vs_normal")].iloc[0]
        m = l.group.isin(["PDAC", "normal"])
        fpr, tpr, _ = roc_curve(l.group[m] == "PDAC", s[m])
        ax.plot(fpr, tpr, c=MRED, lw=2.6, label=f"Tüm PDAC · AUC {r.auc:.2f}".replace(".", ","))
        if l.early.notna().any():
            me = (l.group == "normal") | ((l.group == "PDAC") & (l.early == 1))
            fpr, tpr, _ = roc_curve(l.group[me] == "PDAC", s[me])
            re_ = res[(res.cohort == gse) & (res.comparison == "T2_erken_PDAC_vs_normal")].iloc[0]
            ax.plot(fpr, tpr, c=MNAVY, lw=2, ls="--", label=f"Erken evre · AUC {re_.auc:.2f}".replace(".", ","))
        ax.plot([0, 1], [0, 1], c=MGREY, lw=1, ls=":")
        ax.set_title(f"{gse} · {int(r.n_pos)} tümör / {int(r.n_neg)} normal", fontsize=13)
        ax.set_xlabel("1 − özgüllük"); ax.set_ylabel("Duyarlılık")
        ax.legend(loc="lower right", frameon=False, fontsize=11.5)
    fig.tight_layout()
    fig.savefig(file)
    plt.close(fig)


def make_qr(url, file):
    qr = qrcode.QRCode(border=1, box_size=12)
    qr.add_data(url); qr.make(fit=True)
    qr.make_image(fill_color="#1A1A1A", back_color="white").save(file)


# ----------------------------------------------------------------------------- pptx yardımcıları
def write(shape, paras, top=None, height=None, box=None):
    """paras: [dict(runs=[(metin, {bold, italic, color})], size, align, space_after)]"""
    tf = shape.text_frame
    tf.word_wrap = True
    tf.auto_size = None
    tf.vertical_anchor = MSO_ANCHOR.TOP
    for p in list(tf.paragraphs)[1:]:
        p._p.getparent().remove(p._p)
    first = True
    for spec in paras:
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        for r in list(p.runs):
            r._r.getparent().remove(r._r)
        p.alignment = spec.get("align", PP_ALIGN.JUSTIFY)
        p.space_after = Pt(spec.get("space_after", 4))
        p.line_spacing = 1.0
        for text, st in spec["runs"]:
            run = p.add_run()
            run.text = text
            f = run.font
            f.name = FONT; f.size = Pt(spec.get("size", 28))
            f.bold = st.get("bold", False); f.italic = st.get("italic", False)
            f.color.rgb = st.get("color", INK)
    if box is not None:  # beyaz kutunun içine, 24 pt iç boşlukla
        shape.left = Pt(box[0] + 24); shape.width = Pt(box[1] - 48)
    if top is not None:
        shape.top = Pt(top)
    if height is not None:
        shape.height = Pt(height)


def P(*runs, size=28, align=PP_ALIGN.JUSTIFY, space_after=4):
    return {"runs": [(r, {}) if isinstance(r, str) else r for r in runs], "size": size, "align": align, "space_after": space_after}


def B(t, **kw):
    return (t, {"bold": True, **kw})


def ph(key, label):
    v = INFO[key]
    return (v, {}) if v else (f"[{label}]", {"color": PH, "bold": True})


def heading(t):
    return P(B(t), size=32, align=PP_ALIGN.LEFT, space_after=8)


def add_table(slide, rows, left, top, width, col_w, size=20, head_fill=RED, row_h=None):
    nr, nc = len(rows), len(rows[0])
    shp = slide.shapes.add_table(nr, nc, Pt(left), Pt(top), Pt(width), Pt(row_h * nr if row_h else 30 * nr))
    tbl = shp.table
    for j, w in enumerate(col_w):
        tbl.columns[j].width = Pt(width * w)
    for i, row in enumerate(rows):
        if row_h:
            tbl.rows[i].height = Pt(row_h)
        for j, val in enumerate(row):
            cell = tbl.cell(i, j)
            cell.margin_left = cell.margin_right = Pt(6); cell.margin_top = cell.margin_bottom = Pt(3)
            tf = cell.text_frame; tf.word_wrap = True
            p = tf.paragraphs[0]; p.alignment = PP_ALIGN.LEFT
            text, color, bold = (val if isinstance(val, tuple) else (val, None, False))
            r = p.add_run(); r.text = text
            r.font.name = FONT; r.font.size = Pt(size)
            r.font.bold = (i == 0) or bold
            r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF) if i == 0 else (color or INK)
            cell.fill.solid()
            cell.fill.fore_color.rgb = head_fill if i == 0 else (RGBColor(0xF7, 0xEE, 0xEE) if i % 2 == 0 else RGBColor(0xFF, 0xFF, 0xFF))
    return shp


def pic(slide, file, left, top, width=None, height=None):
    return slide.shapes.add_picture(str(file), Pt(left), Pt(top), Pt(width) if width else None, Pt(height) if height else None)


# ----------------------------------------------------------------------------- içerik
def main() -> None:
    fig_pipeline(IMG / "is_akisi.png")
    fig_external_roc(IMG / "dis_dogrulama_roc.png")
    make_qr("https://github.com/Soap2647/Fisilti", IMG / "qr_github.png")

    ext = pd.read_csv(path("tables") / "wp6_tissue_external.tsv", sep="\t")
    single = pd.read_csv(path("tables") / "wp6_tissue_external_single_genes.tsv", sep="\t", index_col=0)
    ann = pd.read_csv(path("tables") / "wp5_protein_annotation.tsv", sep="\t", index_col=0)
    serum = pd.read_csv(path("tables") / "wp6_serum_primary_secondary.tsv", sep="\t")
    ca = pd.read_csv(path("tables") / "wp6_serum_ca19_9.tsv", sep="\t", index_col=0)["value"]
    tab = pd.read_csv(path("tables") / "wp3_tabular_summary.tsv", sep="\t")
    f2 = lambda x: f"{x:.2f}".replace(".", ",")
    t1 = ext[ext.comparison == "T1_PDAC_vs_normal"].set_index("cohort").auc
    early = ext[ext.comparison == "T2_erken_PDAC_vs_normal"].auc.iloc[0]
    p1 = serum[serum.comparison == "P1_PDAC_vs_saglikli"].set_index("cohort").auc
    best = tab[tab.config == "blood"].groupby("view").auc.max()

    prs = Presentation(POSTER / "sablon.pptx")
    s = prs.slides[0]
    sh = list(s.shapes)

    # --- Başlık alanı
    write(sh[6], [P(B("FISILTI: ERKEN EVRE PANKREAS KANSERİ İÇİN KAN TABANLI BİYOBELİRTEÇ KEŞFİ"),
                    size=40, align=PP_ALIGN.CENTER, space_after=0)], top=336, height=96)
    sh[6].left = Pt(440); sh[6].width = Pt(1140)
    write(sh[7], [P(B("Mert Yasin YILMAZ"), (" · Kastamonu Üniversitesi, ", {}), ph("bolum", "Fakülte / Bölüm"), size=26, align=PP_ALIGN.CENTER, space_after=0),
                  P(("Danışman: ", {}), ph("danisman", "Unvan Ad SOYAD"), size=24, align=PP_ALIGN.CENTER, space_after=0)],
          top=440, height=80)
    sh[7].left = Pt(440); sh[7].width = Pt(1140)
    write(sh[14], [P(B("Tematik Alan-Alt Alan"), size=24, align=PP_ALIGN.CENTER, space_after=2),
                   P(ph("tematik", "Alan / Alt alan"), size=22, align=PP_ALIGN.CENTER, space_after=0)])
    write(sh[19], [P(B("Kategori / Katılımcı Grubu:"), size=24, align=PP_ALIGN.CENTER, space_after=2),
                   P(ph("kategori", "Öğrenci / Akademisyen"), size=22, align=PP_ALIGN.CENTER, space_after=0)])
    write(sh[13], [P(B("Proje Kodu"), size=24, align=PP_ALIGN.CENTER, space_after=2),
                   P(ph("proje_kodu", "Kod"), size=22, align=PP_ALIGN.CENTER, space_after=0)])

    # --- Özet (kutu: x 64, y 626, w 934, h 781)
    ozet = (
        "Pankreas duktal adenokarsinomu (PDAC) beş yıllık sağkalımı %10'un altında olan ve çoğunlukla ameliyat "
        "edilemeyecek evrede tanı alan bir kanserdir; klinikteki tek kan belirteci CA19-9'un duyarlılığı ve "
        "özgüllüğü sınırlıdır. Fısıltı, kamuya açık yedi kohortun dört omik katmanını (mRNA, protein, miRNA, "
        "metilasyon) MOFA+ ile bütünleştirdi; L1-lojistik, XGBoost, LightGBM ve dikkat mekanizmalı çoklu-örnek "
        "öğrenimi (Attention-MIL) modellerini hasta bazlı, iç içe çapraz doğrulamayla eğitti. Paneller kan verisinin "
        "etiketlerine bakılmadan kilitlendi ve analiz planları sonuçlar görülmeden kaydedildi (ön-kayıt). "
        f"Sekiz proteinlik aday panel, modelin hiç görmediği üç bağımsız doku kohortunda AUC {f2(t1.min())}–{f2(t1.max())} "
        f"(erken evre {f2(early)}) verdi. Beş miRNA'lık serum paneli ise bağımsız serum kohortlarında doğrulanamadı "
        f"(AUC {f2(p1['GSE85589'])} ve {f2(p1['GSE59856'])}); iki kamuya açık serum kohortunun birbirini de "
        "doğrulamadığı gösterildi. Sonuçlar, plazmada ELISA ile doğrulanmaya hazır bir aday panel sunmakta ve kan "
        "tabanlı biyobelirteç araştırmalarında bağımsız doğrulamanın önemini ortaya koymaktadır."
    )
    write(sh[8], [heading("Özet"), P(ozet),
                  P(B("Anahtar kelimeler: "), ("Pankreas kanseri; Çoklu-omik analiz; Yapay zekâ; Sıvı biyopsi; Ön-kayıt", {"italic": True}),
                    size=28, align=PP_ALIGN.LEFT)], top=652, height=730, box=(64, 934))

    # --- Materyal ve Yöntem (kutu: x 1043, y 626, w 879, h 781)
    yontem = [
        heading("Projenin Materyali ve Yöntemi"),
        P("• ", B("Veri: "), "CPTAC-PDAC (keşif; 140 tümör, 92 normal, 4 omik), TCGA-PAAD, 3 GEO doku ve 2 GEO serum kohortu; tümü açık ve anonim.", space_after=3),
        P("• ", B("Entegrasyon: "), "MOFA+ ile 13 çoklu-omik faktör (denetimsiz).", space_after=3),
        P("• ", B("Modelleme: "), "5×5 iç içe, hasta gruplu çapraz doğrulama; Optuna ile ayar; kararlılık seçimi.", space_after=3),
        P("• ", B("Yorum: "), "test örneklerinde SHAP; GSEA (Hallmark, KEGG, Reactome).", space_after=3),
        P("• ", B("Panel: "), "tümörde artan, salgılanan (HPA, UniProt), kanda ölçülebilir belirteçler; kilitleme ve ön-kayıt.", space_after=3),
    ]
    write(sh[11], yontem, top=652, height=400, box=(1043, 879))
    pic(s, IMG / "is_akisi.png", 1062, 1060, width=840)
    rows = [["Omik (doku, CPTAC)", "En iyi AUC", "Model"],
            ["Protein", f2(best["protein"]), "L1-lojistik"], ["mRNA", f2(best["mrna"]), "LightGBM"],
            ["Metilasyon", f2(best["methylation"]), "XGBoost"], ["miRNA", f2(best["mirna"]), "XGBoost"]]
    add_table(s, rows, 1062, 1236, 840, [0.44, 0.24, 0.32], size=20, row_h=30)

    # --- Amaç (kutu: x 61, y 1432, w 934, h 661)
    write(sh[10], [heading("AMAÇ"),
                   P(B("Ana amaç: "), "erken evre PDAC'ı (Evre I–II) kan örneğinden yakalayabilecek, ELISA / RT-qPCR ile "
                     "ölçülebilir 5–10 belirteçlik bir panel keşfetmek ve bağımsız kohortlarda sınamak.")],
          top=1458, height=190, box=(61, 934))
    goals = [["Başvurudaki hedef", "Sonuç", "Durum"],
             ["Çoklu kohort veri entegrasyonu", "7 kohort, 4 omik; 6 veri tuzağı yakalandı", ("Karşılandı", RGBColor(0x1B, 0x6B, 0x2F), True)],
             ["MOFA+ ile ortak faktörler", "13 faktör; baskın eksen asiner doku kaybı", ("Karşılandı", RGBColor(0x1B, 0x6B, 0x2F), True)],
             ["Yapay zekâ ile aday sıralama", "600 model, hasta bazlı iç içe doğrulama", ("Karşılandı", RGBColor(0x1B, 0x6B, 0x2F), True)],
             ["Yorumlanabilirlik (SHAP, yolaklar)", "KRAS, EMT, p53; bilinen belirteçler bulundu", ("Karşılandı", RGBColor(0x1B, 0x6B, 0x2F), True)],
             ["5–10 belirteçlik panel", "5 miRNA + 8 protein, kilitli", ("Karşılandı", RGBColor(0x1B, 0x6B, 0x2F), True)],
             ["Kohortlar arası doğrulama", f"Doku: AUC {f2(t1.min())}–{f2(t1.max())}; serum miRNA: doğrulanamadı", ("Kısmen", RGBColor(0xB0, 0x5A, 0x00), True)]]
    add_table(s, goals, 84, 1655, 890, [0.36, 0.47, 0.17], size=20, row_h=46)

    # --- Uygulama Alanı ve Hedef Grubu (kutu: x 1046, y 1428, w 879, h 661)
    write(sh[18], [heading("Projenin Uygulama Alanı ve Hedef Grubu"),
                   P("• ", B("Hedef grup: "), "onkoloji ve gastroenteroloji klinikleri, yüksek riskli bireylerin taranması "
                     "(yeni başlayan diyabet, kronik pankreatit, aile öyküsü), IVD tanı kiti üreticileri.", space_after=3),
                   P("• ", B("Sıradaki adım: "), "Kastamonu Eğitim ve Araştırma Hastanesi ile PDAC, pankreatit ve sağlıklı "
                     "gönüllülerde plazma ELISA pilotu; CA19-9 ile doğrudan kıyas.", space_after=3)],
          top=1454, height=300, box=(1046, 879))
    prot = [["Protein", "Salgılanma kanıtı", "Dış doğrulama AUC"]]
    for g in ["SULF1", "LGALS3BP", "NMI", "LAMC2", "THBS1", "CD55", "GLA", "SERPINB5"]:
        secl = str(ann.at[g, "hpa_secretome"]) if pd.notna(ann.at[g, "hpa_secretome"]) else ""
        ev = "kana salgılanır" if "blood" in secl else ("salgılanır" if (ann.at[g, "signal_peptide"] or "Secreted" in secl) else "plazmada ölçülmüş")
        v = single.loc[g]
        prot.append([g, ev, f"{f2(v.min())}–{f2(v.max())}"])
    add_table(s, prot, 1066, 1762, 840, [0.28, 0.42, 0.30], size=20, row_h=33)

    # --- Özgün Değer ve Yaygın Etki (kutu: x 63, y 2116, w 934, h 610)
    write(sh[16], [heading("Projenin Özgün Değeri ve Yaygın Etkisi"),
                   P("• ", B("Ön-kayıt: "), "panel ve analiz planı sonuç görülmeden, zaman damgalı olarak kilitlendi.", space_after=2),
                   P("• ", B("Bağımsız doğrulama: "), f"protein paneli 3 kohortta AUC {f2(t1.min())}–{f2(t1.max())}; rastgele panellerden üstün (p ≤ 0,004).", space_after=2),
                   P("• ", B("Alana uyarı: "), "iki kamuya açık serum miRNA kohortu birbirini doğrulamadı (ρ = 0,02).", space_after=2)],
          top=2141, height=270, box=(63, 934))
    pic(s, IMG / "dis_dogrulama_roc.png", 92, 2432, width=876)

    # --- Referanslar (kutu: x 1043, y 2116, w 879, h 610)
    refs = [
        "Argelaguet R ve ark. MOFA+. Genome Biol 2020;21:111.",
        "Ilse M, Tomczak JM, Welling M. Attention-based deep MIL. ICML 2018.",
        "Lundberg SM, Lee SI. SHAP. NeurIPS 2017.",
        "Cao L ve ark. Proteogenomic characterization of PDAC. Cell 2021;184:5031.",
        "TCGA Research Network. Pancreatic adenocarcinoma. Cancer Cell 2017;32:185.",
        "Badea L 2008 (GSE15471); Yang S 2016 (GSE62452); Pei H 2009 (GSE16515).",
        "GEO serum kohortları GSE85589 ve GSE59856.",
        "Klein AP. Pancreatic cancer epidemiology. Nat Rev Gastroenterol Hepatol 2021.",
        "T.C. Sağlık Bakanlığı. Ulusal Kanser Kontrol Programı 2023–2028.",
    ]
    write(sh[12], [heading("Referanslar ve Bilgilendirmeler")] +
          [P(f"{k}. {r}", size=19, align=PP_ALIGN.LEFT, space_after=1) for k, r in enumerate(refs, 1)],
          top=2141, height=400, box=(1043, 879))
    pic(s, IMG / "qr_github.png", 1066, 2556, width=150)
    info = s.shapes.add_textbox(Pt(1232), Pt(2556), Pt(670), Pt(150))
    write(info, [P(B("Kod ve veri: "), "github.com/Soap2647/Fisilti", size=20, align=PP_ALIGN.LEFT, space_after=2),
                 P(B("İletişim: "), ph("eposta", "e-posta"), size=20, align=PP_ALIGN.LEFT, space_after=2),
                 P("Kod, 22 yöntem kararı ve tüm sonuçlar açık kaynaktır. Veriler kamuya açık; etik kurul onayı gerekmemiştir. Kavram kanıtı bütçesi: 0 TL.", size=18, align=PP_ALIGN.LEFT)])

    out = POSTER / "Fisilti_Poster.pptx"
    prs.save(out)
    missing = [k for k, v in INFO.items() if not v]
    print(f"[ok] {out.name} · eksik bilgi: {missing or 'yok'}")


if __name__ == "__main__":
    main()
