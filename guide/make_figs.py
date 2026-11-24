"""Diagrams for the HER-OMI study guide. Hand-built SVG so labels sit on the drawing.
Colour carries meaning throughout: blue = HER2 / ERBB2, orange = AR / luminal, grey = neutral or background."""
import pandas as pd

OUT = "guide/figs"
BLUE, BLUE_L = "#2563b8", "#dbe7f7"
ORANGE, ORANGE_L = "#c8551b", "#fbe3d4"
GREY, GREY_L = "#6b7280", "#eef0f3"
INK = "#1f2328"
FONT = "Inter, DejaVu Sans, sans-serif"


def svg(w, h, body):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
            f'font-family="{FONT}" font-size="13" fill="{INK}">'
            '<defs><marker id="a" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{GREY}"/></marker>'
            '<marker id="ab" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{BLUE}"/></marker>'
            '<marker id="ao" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{ORANGE}"/></marker></defs>'
            f'{body}</svg>')


def box(x, y, w, h, fill, stroke, r=8, sw=1.5, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>'


def text(x, y, s, size=13, weight=400, anchor="middle", fill=INK, italic=False):
    st = ' font-style="italic"' if italic else ""
    lines = s.split("\n")
    out = ""
    for i, ln in enumerate(lines):
        out += (f'<text x="{x}" y="{y + i * size * 1.3:.1f}" font-size="{size}" font-weight="{weight}" '
                f'text-anchor="{anchor}" fill="{fill}"{st}>{ln}</text>')
    return out


def arrow(x1, y1, x2, y2, col="a", color=GREY, sw=1.8, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{sw}" marker-end="url(#{col})"{d}/>'


def save(name, s):
    open(f"{OUT}/{name}.svg", "w").write(s)


# ---------------------------------------------------------------- 1. the big picture
def fig_bigpicture():
    b = ""
    cols = [
        ("DNA", "the recipe book", "The ERBB2 gene sits on\nchromosome 17. Extra copies\n(amplification) make HER2+\ncancer.", "copy number and\nbreakpoints near ERBB2\n(step 11)"),
        ("RNA", "a working copy of one recipe", "The cell copies ERBB2 into\nmessenger RNA. RNA levels\nof thousands of genes can\nbe read at once.", "the 144 gene score and\nERBB2 mRNA\n(steps 2 to 7, 13 to 15)"),
        ("Protein", "the dish that is served", "Ribosomes read the RNA and\nbuild the HER2 receptor,\nwhich sits on the cell\nsurface.", "HER2 and active HER2\n(pY1248) protein\n(step 8)"),
    ]
    x0, w, gap = 20, 190, 40
    for i, (t, sub, body, test) in enumerate(cols):
        x = x0 + i * (w + gap)
        b += box(x, 20, w, 190, BLUE_L if i != 1 else "#e8f0fb", BLUE)
        b += text(x + w / 2, 46, t, 18, 700, fill=BLUE)
        b += text(x + w / 2, 66, sub, 12, 400, italic=True)
        b += text(x + w / 2, 96, body, 12)
        b += box(x, 225, w, 70, "white", GREY, dash="4 3")
        b += text(x + w / 2, 244, "what we checked here", 11, 600, fill=GREY)
        b += text(x + w / 2, 261, test, 11)
        if i < 2:
            b += arrow(x + w + 4, 115, x + w + gap - 4, 115)
    # drug
    x = x0 + 3 * (w + gap)
    b += arrow(x - gap + 4, 115, x - 4, 115)
    b += box(x, 20, 150, 190, ORANGE_L, ORANGE)
    b += text(x + 75, 46, "The drug", 18, 700, fill=ORANGE)
    b += text(x + 75, 66, "trastuzumab deruxtecan", 11, 400, italic=True)
    b += text(x + 75, 96, "An antibody that\ngrabs HER2 protein\non the surface.\nNo protein, nothing\nto grab.", 12)
    b += text(x + 75, 250, "The claim we test:\nHER2-like RNA means\nthis drug could work.", 11, 600, fill=ORANGE)
    save("fig1_bigpicture", svg(x + 170, 310, b))


# ---------------------------------------------------------------- 2. breast cancer family tree
def fig_tree():
    b = ""
    b += box(210, 10, 340, 46, GREY_L, GREY)
    b += text(380, 31, "Breast cancer", 15, 700)
    b += text(380, 48, "pathologist stains for three proteins: ER, PR, HER2", 11)
    kids = [(30, "Hormone receptor positive", "ER or PR present\nabout 70% of patients\nhormone pills work", GREY_L, GREY),
            (275, "HER2 positive", "HER2 protein in excess\nabout 15%\nanti-HER2 drugs work", BLUE_L, BLUE),
            (520, "Triple negative (TNBC)", "none of the three\nabout 15%\nno targeted drug, chemo", "#fff", INK)]
    for x, t, s, f, c in kids:
        b += f'<path d="M380,56 L380,72 L{x + 105},72 L{x + 105},88" fill="none" stroke="{GREY}" stroke-width="1.5"/>'
        b += box(x, 88, 210, 80, f, c, sw=2 if "TNBC" in t else 1.5)
        b += text(x + 105, 108, t, 13, 700, fill=c)
        b += text(x + 105, 126, s, 11)
    # TNBC children
    sub = [(330, "Basal-like", "most TNBC; looks like\nthe outer (basal) duct cells", GREY_L, GREY),
           (530, "LAR (luminal androgen\nreceptor)", "driven by the AR hormone\nreceptor; looks luminal", ORANGE_L, ORANGE)]
    for x, t, s, f, c in sub:
        b += f'<path d="M625,168 L625,186 L{x + 95},186 L{x + 95},200" fill="none" stroke="{GREY}" stroke-width="1.5"/>'
        b += box(x, 200, 190, 82, f, c)
        b += text(x + 95, 219, t, 12, 700, fill=c)
        b += text(x + 95, 250 if "\n" in t else 238, s, 11)
    b += box(250, 300, 470, 56, "white", BLUE, dash="5 4")
    b += text(485, 322, "IBCGA (Indian cohort) proposed a new TNBC group: \"HER2-like\"", 12, 700, fill=BLUE)
    b += text(485, 340, "Our question: is it really HER2, or is it the orange LAR box under another name?", 11)
    save("fig2_tree", svg(760, 368, b))


# ---------------------------------------------------------------- 3. how the HER2 drug works
def fig_adc():
    b = ""
    steps = [("1  Bind", "the antibody arm\nlocks onto HER2\nprotein on the surface"),
             ("2  Swallow", "the cell pulls the\nreceptor and drug\ninside"),
             ("3  Release", "the toxic payload is\ncut free inside\nthe cell"),
             ("4  Kill", "the payload damages\nDNA; the cell dies\n(and nearby cells)")]
    for i, (t, s) in enumerate(steps):
        x = 20 + i * 185
        b += box(x, 20, 160, 110, BLUE_L if i == 0 else GREY_L, BLUE if i == 0 else GREY)
        b += text(x + 80, 44, t, 14, 700, fill=BLUE if i == 0 else INK)
        b += text(x + 80, 68, s, 11.5)
        if i < 3:
            b += arrow(x + 162, 75, x + 183, 75)
    b += text(370, 160, "Step 1 needs HER2 protein on the cell surface. RNA alone cannot be grabbed.", 12.5, 700, fill=BLUE)
    b += text(370, 180, "That is why the protein layer (step 8) is the deciding test for the drug claim.", 12)
    save("fig3_adc", svg(760, 195, b))


# ---------------------------------------------------------------- 4. the logical gap
def fig_gap():
    b = ""
    chain = [("What IBCGA saw", "a group of TNBC whose RNA\nlooks a bit like HER2+ cancer\n(144 genes)"),
             ("What they named it", "\"HER2-like\" or\n\"transcriptionally\nHER2-active\""),
             ("What that implies", "these patients might\nrespond to HER2\ndrugs")]
    for i, (t, s) in enumerate(chain):
        x = 20 + i * 250
        b += box(x, 20, 210, 100, GREY_L, GREY)
        b += text(x + 105, 44, t, 13, 700)
        b += text(x + 105, 66, s, 11.5)
        if i < 2:
            b += arrow(x + 212, 70, x + 248, 70)
    qs = [(20, "Q1  Is the HER2 part real?", "Remove the 9 genes that sit\nnext to ERBB2. Does the\nscore still track ERBB2?"),
          (270, "Q2  What else could it be?", "Does the score follow AR\nand the LAR subtype more\nthan it follows HER2?"),
          (520, "Q3  Is there HER2 protein?", "Do high scorers carry more\nHER2 protein, or active\nHER2, than low scorers?")]
    for x, t, s in qs:
        b += arrow(x + 105, 158, x + 105, 124, color=ORANGE, col="ao", dash="4 3")
        b += box(x, 160, 210, 92, "white", ORANGE)
        b += text(x + 105, 182, t, 12.5, 700, fill=ORANGE)
        b += text(x + 105, 202, s, 11.5)
    save("fig4_gap", svg(750, 265, b))


# ---------------------------------------------------------------- 5. singscore by ranks
def fig_singscore():
    import random
    random.seed(4)
    b = ""
    L, R = 60, 700
    def row(y, label, up_pos, dn_pos, score, col):
        s = text(L - 10, y - 18, label, 12.5, 700, anchor="start")
        s += f'<line x1="{L}" y1="{y}" x2="{R}" y2="{y}" stroke="{GREY}" stroke-width="10" stroke-linecap="round" opacity=".18"/>'
        for p in up_pos:
            s += f'<circle cx="{L + p * (R - L)}" cy="{y - 7}" r="4.5" fill="{BLUE}"/>'
        for p in dn_pos:
            s += f'<circle cx="{L + p * (R - L)}" cy="{y + 7}" r="4.5" fill="none" stroke="{ORANGE}" stroke-width="1.8"/>'
        s += text(R + 12, y + 5, f"score {score}", 13, 700, anchor="start", fill=col)
        return s
    hi_up = [random.uniform(.6, 1) for _ in range(12)]
    hi_dn = [random.uniform(0, .45) for _ in range(14)]
    lo_up = [random.uniform(0, 1) for _ in range(12)]
    lo_dn = [random.uniform(0, 1) for _ in range(14)]
    b += text(L, 22, "Inside ONE tumour, all ~18,000 genes are sorted from lowest RNA level (left) to highest (right).", 12.5, 400, anchor="start")
    b += row(80, "Tumour A (HER2-like pattern switched on)", hi_up, hi_dn, "+0.15", BLUE)
    b += row(160, "Tumour B (pattern switched off)", lo_up, lo_dn, "-0.12", GREY)
    b += text(L, 205, "lowest in this tumour", 11, anchor="start", fill=GREY)
    b += text(R, 205, "highest in this tumour", 11, anchor="end", fill=GREY)
    b += arrow(L + 140, 201, R - 140, 201)
    b += f'<circle cx="{L + 6}" cy="230" r="4.5" fill="{BLUE}"/>' + text(L + 18, 234, "a gene the signature says should be UP (33 genes)", 11.5, anchor="start")
    b += f'<circle cx="{L + 386}" cy="230" r="4.5" fill="none" stroke="{ORANGE}" stroke-width="1.8"/>' + text(L + 398, 234, "a gene that should be DOWN (67 genes)", 11.5, anchor="start")
    b += text(L, 262, "Score = how far right the up genes sit + how far left the down genes sit, each scaled to run from -0.5 to +0.5.", 12, anchor="start")
    b += text(L, 280, "Only the order inside one tumour matters, so RNA-seq (TCGA) and microarray (METABRIC) can be compared.", 12, anchor="start")
    save("fig5_singscore", svg(800, 292, b))


# ---------------------------------------------------------------- 6. confounder logic
def fig_confound():
    b = ""
    b += box(270, 15, 200, 58, ORANGE_L, ORANGE)
    b += text(370, 38, "AR / LAR biology", 14, 700, fill=ORANGE)
    b += text(370, 57, "luminal, hormone-driven tumours", 11)
    b += box(40, 165, 200, 58, "#fff", INK)
    b += text(140, 190, "HER2-like score", 14, 700)
    b += text(140, 208, "(100 genes, no amplicon)", 11)
    b += box(500, 165, 200, 58, BLUE_L, BLUE)
    b += text(600, 190, "ERBB2 RNA", 14, 700, fill=BLUE)
    b += text(600, 208, "(the HER2 gene)", 11)
    b += arrow(320, 75, 175, 162, col="ao", color=ORANGE, sw=2.2)
    b += arrow(420, 75, 565, 162, col="ao", color=ORANGE, sw=2.2)
    b += f'<line x1="244" y1="194" x2="496" y2="194" stroke="{GREY}" stroke-width="2" stroke-dasharray="6 4"/>'
    b += text(370, 186, "the link we see", 11.5, 600, fill=GREY)
    b += text(370, 214, "is any of it left once AR is accounted for?", 11.5, italic=True)
    b += text(370, 250, "If AR pushes both up, the score and ERBB2 will move together even with no HER2 biology.", 12)
    b += text(370, 268, "Fixes we used: take LAR and ER-low tumours out, or hold AR fixed with a partial correlation.", 12)
    save("fig6_confound", svg(740, 280, b))


# ---------------------------------------------------------------- 7. HER2 protein ruler
def fig_ruler():
    rows = [("TNBC, bottom score third", -0.43, GREY), ("TNBC, middle third", -0.21, GREY),
            ("TNBC, top score third", -0.13, ORANGE), ("HR+ HER2-negative", 0.14, GREY),
            ("HER2-positive", 0.88, BLUE)]
    L, R, lo, hi = 230, 720, -0.6, 1.0
    X = lambda v: L + (v - lo) / (hi - lo) * (R - L)
    b = text(20, 22, "Median HER2 protein (TCGA RPPA, log scale units). Further right = more HER2 on the cells.", 12.5, anchor="start")
    for i, (lab, v, c) in enumerate(rows):
        y = 55 + i * 34
        b += text(L - 12, y + 4, lab, 12, 600 if c != GREY else 400, anchor="end", fill=c if c != GREY else INK)
        b += f'<line x1="{L}" y1="{y}" x2="{R}" y2="{y}" stroke="{GREY_L}" stroke-width="2"/>'
        b += f'<circle cx="{X(v):.1f}" cy="{y}" r="7" fill="{c}"/>'
        b += text(X(v) + 13, y + 4, f"{v:+.2f}", 11.5, 600, anchor="start", fill=c)
    y = 55 + 5 * 34 + 8
    x1, x2, x3 = X(-0.43), X(-0.13), X(0.88)
    b += f'<line x1="{x1:.1f}" y1="{y}" x2="{x3:.1f}" y2="{y}" stroke="{GREY}" stroke-width="1.5"/>'
    b += f'<line x1="{x1:.1f}" y1="{y}" x2="{x2:.1f}" y2="{y}" stroke="{ORANGE}" stroke-width="5"/>'
    b += text((x1 + x3) / 2, y + 22, "the top third moves only 23% of the way from the bottom third toward HER2-positive", 12, 700, fill=ORANGE)
    save("fig7_ruler", svg(760, y + 34, b))


# ---------------------------------------------------------------- 8. protein forest from real results
def fig_protein():
    p = pd.read_csv("results/protein_test.csv")
    p = p[(p.cohort == "TCGA") & (p.subset == "all TNBC")]
    pick = [("HER2", "HER2 total", BLUE), ("HER2 pY1248", "HER2 active (pY1248)", BLUE),
            ("HER3", "HER3", BLUE), ("EGFR", "EGFR", BLUE),
            ("AKT pS473", "AKT active (downstream)", GREY), ("ERK pT202/Y204", "ERK active (downstream)", GREY),
            ("AR (control)", "AR (luminal control)", ORANGE), ("INPP4B (control)", "INPP4B (luminal control)", ORANGE)]
    L, R, lo, hi = 240, 700, -0.4, 0.65
    X = lambda v: L + (v - lo) / (hi - lo) * (R - L)
    b = text(20, 20, "Correlation of the HER2-like score with each protein, 131 TCGA TNBC. Grey bar = range reached by 1,000 random gene sets.", 12, anchor="start")
    top = 45
    for i, (k, lab, c) in enumerate(pick):
        r = p[p.protein == k].iloc[0]
        y = top + 14 + i * 30
        b += text(L - 12, y + 4, lab, 12, 400, anchor="end")
        b += f'<rect x="{X(r.null_lo):.1f}" y="{y - 7}" width="{X(r.null_hi) - X(r.null_lo):.1f}" height="14" rx="3" fill="{GREY}" opacity=".22"/>'
        b += f'<circle cx="{X(r.rho_score):.1f}" cy="{y}" r="6.5" fill="{c}"/>'
        inside = r.null_lo <= r.rho_score <= r.null_hi
        b += (f'<text x="{X(r.rho_score) + 12:.1f}" y="{y + 4}" font-size="11" font-weight="600" fill="{c}" stroke="white" '
              f'stroke-width="3" paint-order="stroke">{r.rho_score:.2f}' + ("  inside random range" if inside else "  outside") + '</text>')
    yb = top + 14 + len(pick) * 30
    b += f'<line x1="{X(0):.1f}" y1="{top}" x2="{X(0):.1f}" y2="{yb - 12}" stroke="{INK}" stroke-width="1"/>'
    for v in (-0.4, -0.2, 0, 0.2, 0.4, 0.6):
        b += text(X(v), yb + 4, f"{v:g}", 10.5, fill=GREY)
    b += text((L + R) / 2, yb + 22, "Spearman correlation with the score", 11, fill=GREY)
    save("fig8_protein", svg(760, yb + 32, b))


# ---------------------------------------------------------------- 9. the verdict by layer
def fig_verdict():
    rows = [("DNA", "No extra ERBB2 copies, almost no breakpoints (9 of 155).\nControlling for copy number leaves the RNA link unchanged.", "not HER2", GREY),
            ("RNA", "A real but modest link to ERBB2 RNA (0.30 to 0.48 in basal tumours),\nmostly carried by LAR and ER-low tumours (every LAR tumour sits in the top third).", "mostly AR", ORANGE),
            ("Protein", "HER2 0.21 and active HER2 0.18, both inside the random range.\nAR 0.50 and INPP4B 0.55, far outside it.", "not HER2", GREY),
            ("Normal tissue", "Paired with their own normal tissue, TNBC do not raise ERBB2\n(6 of 13 go up); HER2+ tumours do (18 of 24).", "not HER2", GREY)]
    b = ""
    for i, (lay, ev, v, c) in enumerate(rows):
        y = 10 + i * 66
        b += box(10, y, 110, 56, GREY_L, GREY)
        b += text(65, y + 33, lay, 13, 700)
        b += box(130, y, 480, 56, "white", GREY_L)
        b += text(142, y + 23, ev, 11.5, anchor="start")
        b += box(620, y, 120, 56, ORANGE_L if c == ORANGE else "white", c, sw=2)
        b += text(680, y + 33, v, 13, 700, fill=c if c == ORANGE else INK)
    save("fig9_verdict", svg(750, 10 + 4 * 66 + 4, b))


# ---------------------------------------------------------------- 10. full pipeline
def fig_pipeline():
    b = ""
    def node(x, y, w, h, title, sub, f=GREY_L, c=GREY):
        s = box(x, y, w, h, f, c)
        s += text(x + w / 2, y + 18, title, 12, 700, fill=c if c != GREY else INK)
        s += text(x + w / 2, y + 34, sub, 10.5)
        return s
    # column 1: data
    b += text(85, 18, "DATA", 12, 700, fill=GREY)
    data = [("cBioPortal: TCGA", "RNA, clinical, RPPA,\ncopy number, normals"), ("cBioPortal: METABRIC", "microarray RNA,\nclinical, copy number"),
            ("CPTAC 2020", "mass-spec protein"), ("IBCGA preprint", "144 genes, pathway\ntables S21 to S23")]
    for i, (t, s) in enumerate(data):
        b += node(10, 30 + i * 78, 150, 66, t, s)
    b += text(85, 345, "00_download_data.sh", 10.5, 600, fill=GREY)
    # column 2: build
    b += text(300, 18, "BUILD THE SCORE", 12, 700, fill=GREY)
    build = [("01 define TNBC", "ER, PR, HER2 negative:\n159 TCGA, 250 METABRIC"), ("01b signature", "144 genes; 109 measured\nin both; 100 off-amplicon"),
             ("02 singscore", "one score per tumour\n(full, no_amp, amp_only)"), ("02b / 02c annotate", "add PAM50, AR, ERBB2,\nTNBCtype tables")]
    for i, (t, s) in enumerate(build):
        b += node(220, 30 + i * 78, 160, 66, t, s, BLUE_L, BLUE)
        if i < 3:
            b += arrow(300, 96 + i * 78, 300, 106 + i * 78)
    for i in range(4):
        b += f'<line x1="162" y1="{63 + i * 78}" x2="190" y2="{63 + i * 78}" stroke="{GREY}" stroke-width="1.8"/>'
    b += f'<line x1="190" y1="63" x2="190" y2="297" stroke="{GREY}" stroke-width="1.8"/>'
    b += arrow(190, 63, 218, 63) + arrow(190, 141, 218, 141)
    # column 3: tests by question
    b += text(560, 18, "TEST THE CLAIM", 12, 700, fill=GREY)
    q = [("Is the HER2 link real?", "03 basal only   04 HER2 IHC\n05 / 09 pathways   10 vs HER2+"),
         ("What else explains it?", "06 confounders   07 TNBCtype\n12 ER-low sensitivity   15 clusters"),
         ("Is HER2 protein there?", "08 RPPA and CPTAC\n08b robustness"),
         ("DNA and normal tissue?", "11 DNA layer   13 normal vs\ntumour   14 admixture")]
    for i, (t, s) in enumerate(q):
        b += node(450, 30 + i * 78, 220, 66, t, s, "white", ORANGE)
    b += f'<line x1="382" y1="270" x2="420" y2="270" stroke="{GREY}" stroke-width="1.8"/>'
    b += f'<line x1="420" y1="63" x2="420" y2="297" stroke="{GREY}" stroke-width="1.8"/>'
    for i in range(4):
        b += arrow(420, 63 + i * 78, 448, 63 + i * 78)
    # column 4: outputs
    b += text(800, 18, "OUTPUTS", 12, 700, fill=GREY)
    outs = [("results/*.csv", "every number in\nthe slides"), ("results/*.png", "heatmaps, pairs,\ngroup plots"), ("deck + notes", "the 23 slides\nyou present")]
    for i, (t, s) in enumerate(outs):
        b += node(730, 30 + i * 100, 140, 66, t, s)
    for i in range(4):
        b += f'<line x1="672" y1="{63 + i * 78}" x2="700" y2="{63 + i * 78}" stroke="{GREY}" stroke-width="1.8"/>'
    b += f'<line x1="700" y1="63" x2="700" y2="297" stroke="{GREY}" stroke-width="1.8"/>'
    b += arrow(700, 63, 728, 63)
    b += arrow(800, 96, 800, 128) + arrow(800, 196, 800, 228)
    b += text(560, 345, "run_all.sh runs every step in order", 10.5, 600, fill=GREY)
    save("fig10_pipeline", svg(880, 355, b))


for f in [fig_bigpicture, fig_tree, fig_adc, fig_gap, fig_singscore, fig_confound, fig_ruler, fig_protein, fig_verdict, fig_pipeline]:
    f()
print("ok")
