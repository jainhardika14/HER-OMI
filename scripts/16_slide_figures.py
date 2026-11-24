"""HER-OMI step 16: figures for the slides, drawn straight from the data (no hand-built charts).

Colour code used in every figure: pink = HER2 / ERBB2, teal = AR / LAR, grey = random gene sets or neutral.
Outputs results/fig_*.png. Uses the same nulls and seeds as steps 3, 8, 9 and 11, so the numbers on the
figures match the tables.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, "scripts")
from heromi_lib import load_expression, singscore, partial_spearman, random_set_scores

RAW = "data/raw"
PINK, PINK_L, TEAL, TEAL_L = "#C2185B", "#F3B6CE", "#0097A7", "#9FD8DF"
GREY, GREY_L, INK, MUTED, BG = "#8A7D85", "#E4DADF", "#3A1A2A", "#6B5A63", "#FFFBFC"
plt.rcParams.update({"font.family": "Inter", "font.size": 15, "axes.edgecolor": "#CDBFC6", "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "axes.titleweight": "bold",
                     "axes.titlesize": 16.5, "axes.titlecolor": INK, "axes.titlelocation": "left"})
FILES = {"TCGA": f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt",
         "METABRIC": f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt"}

sig = pd.read_csv("data/signature_coverage.csv")
core = sig[sig.TCGA & sig.METABRIC & ~sig.on_her2_amplicon]
n_up, n_dn = int(core.direction.eq("Upregulated").sum()), int(core.direction.eq("Downregulated").sum())
sig_ids = set(sig.loc[sig.entrez_id.notna(), "entrez_id"].astype(int))
ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")

# TNBCtype calls, as in step 7
calls = []
for c in FILES:
    r = pd.read_csv(f"data/tnbctype/{c}_result.csv", index_col=0)
    r.index = r.index.str.replace(".", "-", regex=False)
    calls.append(r[["subtype"]])
ann = ann.join(pd.concat(calls))
ann["subtype"] = ann.subtype.fillna("ER-low")
expr = {c: load_expression(p, set(ann.index[ann.cohort.eq(c)])) for c, p in FILES.items()}
nulls = {c: pd.DataFrame(random_set_scores(e, sig_ids, n_up, n_dn).T, index=e.columns) for c, e in expr.items()}


def save(fig, name):
    fig.savefig(f"results/{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------- 1. replication: real signature vs 1,000 random gene sets ----------------
fig, axes = plt.subplots(1, 4, figsize=(19, 5), sharey=True)
for ax, (c, sub) in zip(axes, [("TCGA", "all TNBC"), ("TCGA", "basal only"), ("METABRIC", "all TNBC"), ("METABRIC", "basal only")]):
    a = ann.loc[expr[c].columns]
    m = np.ones(len(a), bool) if sub == "all TNBC" else a.pam50.eq("Basal").to_numpy()
    erbb2 = expr[c].loc[2064][m]
    obs = stats.spearmanr(a.score_no_amp[m], erbb2)[0]
    null = np.array([stats.spearmanr(nulls[c].iloc[m.nonzero()[0], i], erbb2)[0] for i in range(1000)])
    q95 = np.quantile(null, .95)
    ax.hist(null, bins=40, color=GREY_L, edgecolor=BG, linewidth=.5)
    ax.axvline(q95, color=GREY, lw=1.5, ls="--")
    ax.axvline(obs, color=PINK, lw=3)
    ax.text(obs, ax.get_ylim()[1] * .92, f" signature\n {obs:.2f}", color=PINK, fontweight="bold", va="top", fontsize=13)
    ax.text(-.43, ax.get_ylim()[1] * .97, f"dashed: 95% of\nrandom sets\nbelow {q95:.2f}", color=GREY, ha="left", va="top", fontsize=12.5)
    ax.set_title(f"{c}, {sub.replace(' only', '')} (n={m.sum()})", fontsize=15.5)
    ax.set_xlabel("rho, score vs ERBB2 RNA")
    ax.set_xlim(-.45, .7)
axes[0].set_ylabel("number of random gene sets")
save(fig, "fig_rep_nulls")

# ---------------- 2. score by TNBCtype subtype ----------------
ORDER = ["BL1", "BL2", "IM", "M", "MSL", "UNS", "LAR", "ER-low"]
fig, axes = plt.subplots(1, 2, figsize=(16, 5.4), sharey=True)
rng = np.random.default_rng(607)
for ax, c in zip(axes, FILES):
    a = ann[ann.cohort.eq(c)]
    cut = a.score_no_amp.quantile(2 / 3)
    ax.axhspan(cut, a.score_no_amp.max() + .05, color=PINK_L, alpha=.25, lw=0)
    ax.axhline(cut, color=GREY, ls="--", lw=1.2)
    ax.text(-.45, cut + .006, "top third of scores", ha="left", va="bottom", color=MUTED, fontsize=12.5)
    for i, s in enumerate(ORDER):
        v = a.score_no_amp[a.subtype.eq(s)].to_numpy()
        if not len(v):
            continue
        col = TEAL if s == "LAR" else (TEAL_L if s == "ER-low" else GREY)
        ax.scatter(i + rng.uniform(-.22, .22, len(v)), v, s=26, color=col, alpha=.85, edgecolor="none")
        ax.plot([i - .3, i + .3], [np.median(v)] * 2, color=INK, lw=2.2)
        ax.text(i, ax.get_ylim()[0] if False else -0.33, f"n={len(v)}", ha="center", color=MUTED, fontsize=10.5)
    lar = a.subtype.eq("LAR")
    ax.set_title(f"{c}: all {lar.sum()} LAR tumours sit in the top third")
    ax.set_xticks(range(len(ORDER)), ORDER)
    for t in ax.get_xticklabels():
        if t.get_text() in ("LAR", "ER-low"):
            t.set_color(TEAL); t.set_fontweight("bold")
    ax.set_ylim(-.36, .3)
axes[0].set_ylabel("HER2-like score (100 genes)")
fig.text(.5, -.03, "TNBCtype subtype (Lehmann). BL = basal-like, IM = immunomodulatory, M = mesenchymal, MSL = mesenchymal stem-like, "
         "UNS = unstable, LAR = luminal androgen receptor, ER-low = flagged as possibly ER-positive", ha="center", color=MUTED, fontsize=11)
save(fig, "fig_lar_subtypes")

# ---------------- 3. protein: scatter HER2 and AR, plus all proteins against the random range ----------------
rppa = pd.read_csv(f"{RAW}/brca_tcga/data_rppa.txt", sep="\t").set_index("Composite.Element.REF").T
e = expr["TCGA"]
ids = [s for s in e.columns if s in rppa.index]
d = ann.loc[ids].join(rppa.loc[ids].astype(float))
nul = nulls["TCGA"].loc[ids]
pt = pd.read_csv("results/protein_test.csv")
pt = pt[(pt.cohort == "TCGA") & (pt.subset == "all TNBC")].set_index("protein")
fig = plt.figure(figsize=(17, 5.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.35], wspace=.32)
for k, (col, lab, colr, key) in enumerate([("ERBB2|HER2", "HER2 protein", PINK, "HER2"), ("AR|AR", "AR protein", TEAL, "AR (control)")]):
    ax = fig.add_subplot(gs[k])
    ax.scatter(d.score_no_amp, d[col], s=30, color=colr, alpha=.75, edgecolor="none")
    r = pt.loc[key]
    ax.set_title(f"{lab}: rho = {r.rho_score:.2f}")
    ax.text(.03, .97, f"random gene sets reach {r.null_hi:.2f}", transform=ax.transAxes, va="top", color=MUTED, fontsize=11.5)
    ax.set_xlabel("HER2-like score")
    ax.set_ylabel(f"{lab} (RPPA)")
ax = fig.add_subplot(gs[2])
show = [("HER2", "HER2", PINK), ("HER2 pY1248", "HER2 active (pY1248)", PINK), ("HER3", "HER3", PINK), ("EGFR", "EGFR", PINK),
        ("AKT pS473", "AKT active", GREY), ("ERK pT202/Y204", "ERK active", GREY),
        ("AR (control)", "AR", TEAL), ("INPP4B (control)", "INPP4B", TEAL)]
for i, (k, lab, colr) in enumerate(show[::-1]):
    r = pt.loc[k]
    ax.plot([r.null_lo, r.null_hi], [i, i], color=GREY_L, lw=11, solid_capstyle="round")
    ax.scatter(r.rho_score, i, s=90, color=colr, zorder=3)
    ax.text(r.rho_score + .03, i, f"{r.rho_score:.2f}", va="center", color=colr, fontsize=11.5, fontweight="bold")
ax.set_yticks(range(len(show)), [s[1] for s in show[::-1]])
ax.axvline(0, color="#CDBFC6", lw=1)
ax.set_xlim(-.4, .7)
ax.set_xlabel("correlation with the score (grey bar = 95% of random gene sets)")
ax.set_title("All proteins, 131 TCGA TNBC")
save(fig, "fig_protein")

# ---------------- 4. HER2 protein scale and the HER2-low tumours ----------------
clin = pd.read_csv(f"{RAW}/brca_tcga/data_clinical_patient.txt", sep="\t", comment="#").set_index("PATIENT_ID")
ref = rppa[["ERBB2|HER2"]].astype(float).copy()
ref["pid"] = ref.index.str[:12]
ref = ref.join(clin[["ER_STATUS_BY_IHC", "IHC_HER2", "HER2_FISH_STATUS", "HER2_IHC_SCORE"]], on="pid")
her2pos = ref.IHC_HER2.eq("Positive") | ref.HER2_FISH_STATUS.eq("Positive")
hrneg = ref.ER_STATUS_BY_IHC.eq("Positive") & ~her2pos & (ref.IHC_HER2.eq("Negative") | ref.HER2_FISH_STATUS.eq("Negative"))
d["third"] = pd.qcut(d.score_no_amp, 3, labels=["low", "mid", "high"])
groups = [("TNBC\nlow third", d.loc[d.third.eq("low"), "ERBB2|HER2"], GREY), ("TNBC\nmid third", d.loc[d.third.eq("mid"), "ERBB2|HER2"], GREY),
          ("TNBC\ntop third", d.loc[d.third.eq("high"), "ERBB2|HER2"], TEAL), ("HR+\nHER2-neg", ref.loc[hrneg, "ERBB2|HER2"], GREY),
          ("HER2+", ref.loc[her2pos, "ERBB2|HER2"], PINK)]
fig = plt.figure(figsize=(20, 6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1, 1], wspace=.55)
ax = fig.add_subplot(gs[0])
for i, (lab, v, colr) in enumerate(groups):
    ax.scatter(i + rng.uniform(-.25, .25, len(v)), v, s=12 if len(v) > 100 else 22, color=colr, alpha=.55, edgecolor="none")
    ax.plot([i - .32, i + .32], [v.median()] * 2, color=INK, lw=2.5)
    ax.text(i, 4.6, f"{v.median():+.2f}", ha="center", fontsize=14, color=colr if colr != GREY else INK, fontweight="bold")
ax.set_xticks(range(len(groups)), [g[0] for g in groups])
ax.set_ylabel("HER2 protein (RPPA)")
ax.set_ylim(-2, 5.1)
ax.set_title("HER2 protein by group (median above)")
low = d[pd.to_numeric(d.index.str[:12].map(clin.HER2_IHC_SCORE), errors="coerce").isin([1, 2])]
for k, (col, lab, colr) in enumerate([("ERBB2|HER2", "HER2 protein", PINK), ("AR|AR", "AR protein", TEAL)]):
    ax = fig.add_subplot(gs[k + 1])
    rho = stats.spearmanr(low.score_no_amp, low[col])[0]
    ax.scatter(low.score_no_amp, low[col], s=34, color=colr, alpha=.8, edgecolor="none")
    ax.set_title(f"HER2-low TNBC: {lab}\nrho = {rho:.2f}")
    ax.set_xlabel("HER2-like score")
    ax.set_ylabel(f"{lab} (RPPA)")
    ax.text(.03, .97, f"n = {len(low)} (IHC 1+ or 2+)", transform=ax.transAxes, va="top", color=MUTED, fontsize=11.5)
save(fig, "fig_her2_scale")

# ---------------- 5. the circular pathway evidence ----------------
gmt = {f[0]: f[2:] for f in (l.rstrip("\n").split("\t") for l in open("data/ref/c2.cp.reactome.gmt"))
       if f[0] == "REACTOME_GRB7_EVENTS_IN_ERBB2_SIGNALING"}
pw = gmt["REACTOME_GRB7_EVENTS_IN_ERBB2_SIGNALING"]
SIG = set(sig.gene)
fig, axes = plt.subplots(1, 4, figsize=(20, 4.6), gridspec_kw={"wspace": .5})
for j, c in enumerate(FILES):
    sym = pd.read_csv(FILES[c], sep="\t", usecols=[0, 1]).dropna().drop_duplicates("Hugo_Symbol")
    s2e = dict(zip(sym.Hugo_Symbol, sym.Entrez_Gene_Id.astype(int)))
    idf = lambda genes: [s2e[x] for x in genes if x in s2e and s2e[x] in expr[c].index]
    a = ann.loc[expr[c].columns]
    for k, (genes, lab, colr) in enumerate([(pw, "pathway, all genes", PINK), ([x for x in pw if x not in SIG], "without ERBB2 and GRB7", GREY)]):
        ps = singscore(expr[c], idf(genes), []).total
        rho = stats.spearmanr(ps, a.score_full)[0]
        ax = axes[2 * j + k]
        ax.scatter(a.score_full, ps, s=16, color=colr, alpha=.65, edgecolor="none")
        ax.set_title(f"{c}: {lab}\nrho = {rho:.2f}", fontsize=13.5)
        ax.set_xlabel("signature score")
        if k == 0:
            ax.set_ylabel("GRB7 events in ERBB2\nsignalling (pathway score)")
save(fig, "fig_circular")

# ---------------- 6. DNA: copy number explains ERBB2 RNA, not the score ----------------
ERBB2_POS = 37_856_254
seg = pd.read_csv(f"{RAW}/brca_tcga/data_cna_hg19.seg", sep="\t")
seg.columns = ["ID", "chrom", "start", "end", "nmark", "mean"]
seg = seg[seg.chrom.astype(str).eq("17") & (seg.start <= ERBB2_POS) & (seg.end >= ERBB2_POS)]
cn = seg.assign(sid=seg.ID.str[:15]).groupby("sid")["mean"].mean()
t = ann.loc[expr["TCGA"].columns].assign(sid=lambda x: x.index.str[:15])
t["cn"] = t.sid.map(cn)
t = t.dropna(subset=["cn"])
t["erbb2"] = np.log2(t.ERBB2_expr + 1)
fig, axes = plt.subplots(1, 3, figsize=(19, 5.2), gridspec_kw={"wspace": .55})
for ax, (y, lab, colr) in zip(axes[:2], [("erbb2", "ERBB2 RNA (log2)", PINK), ("score_no_amp", "HER2-like score", GREY)]):
    rho = stats.spearmanr(t.cn, t[y])[0]
    ax.scatter(t.cn, t[y], s=22, color=colr, alpha=.7, edgecolor="none")
    ax.set_xlabel("ERBB2 DNA copy number (log2 ratio)")
    ax.set_ylabel(lab)
    ax.set_title(f"copy number vs {lab.split(' (')[0]}\nrho = {rho:.2f}")
dl = pd.read_csv("results/dna_layer_test.csv")
dl = dl[dl.subset.eq("all TNBC")].set_index("cohort")
ax = axes[2]
for i, c in enumerate(["TCGA", "METABRIC"]):
    ax.bar(i - .18, dl.loc[c, "score_vs_mRNA"], .34, color=PINK_L)
    ax.bar(i + .18, dl.loc[c, "score_vs_mRNA_given_cn"], .34, color=PINK)
    ax.text(i - .18, dl.loc[c, "score_vs_mRNA"] + .01, f"{dl.loc[c, 'score_vs_mRNA']:.2f}", ha="center", fontsize=11.5)
    ax.text(i + .18, dl.loc[c, "score_vs_mRNA_given_cn"] + .01, f"{dl.loc[c, 'score_vs_mRNA_given_cn']:.2f}", ha="center", fontsize=11.5, fontweight="bold")
ax.set_xticks([0, 1], ["TCGA", "METABRIC"])
ax.set_ylim(0, .7)
ax.set_ylabel("score vs ERBB2 RNA (rho)")
ax.set_title("before (light) and after (dark)\naccounting for copy number", fontsize=13.5)
save(fig, "fig_dna")

# ---------------- 7. protein: what we would expect if the score measured HER2, vs what we see ----------------
pt_all = pd.read_csv("results/protein_test.csv")
pt_all = pt_all[(pt_all.cohort == "TCGA") & (pt_all.subset == "all TNBC")].set_index("protein")
fig = plt.figure(figsize=(20, 6.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.45], wspace=.62)
x = d.score_no_amp.to_numpy()
ylim = (d["ERBB2|HER2"].min() - .3, d["ERBB2|HER2"].max() + .3)
# A: illustration only. Same tumours and same HER2 values, re-paired so that HER2 rises with the score
#    at the strength a real link would need (Spearman about 0.6, like AR's 0.50 or stronger).
gen = np.random.default_rng(607)
latent = stats.rankdata(x) / len(x) + gen.normal(0, .42, len(x))
yexp = np.sort(d["ERBB2|HER2"].to_numpy())[stats.rankdata(latent, method="ordinal").astype(int) - 1]
ax = fig.add_subplot(gs[0])
ax.scatter(x, yexp, s=30, facecolor="none", edgecolor=PINK, lw=1.3, alpha=.85)
ax.set_title(f"EXPECTED if the score measured HER2\nrho about {stats.spearmanr(x, yexp)[0]:.1f}, beyond chance", color=PINK)
ax.text(.03, .97, "illustration, not data:\nthe real HER2 values, re-paired\nto show a genuine upward trend", transform=ax.transAxes,
        va="top", color=MUTED, fontsize=12, style="italic")
ax.set_xlabel("HER2-like score"); ax.set_ylabel("HER2 protein (RPPA)"); ax.set_ylim(*ylim)
# B: the real data
ax = fig.add_subplot(gs[1])
ax.scatter(x, d["ERBB2|HER2"], s=30, color=PINK, alpha=.8, edgecolor="none")
r = pt_all.loc["HER2"]
ax.set_title(f"WHAT WE SEE (real data)\nrho = {r.rho_score:.2f}; chance reaches {r.null_hi:.2f}", color=INK)
ax.set_xlabel("HER2-like score"); ax.set_ylabel("HER2 protein (RPPA)"); ax.set_ylim(*ylim)
# C: every protein, with the zone a real link would have to reach
ax = fig.add_subplot(gs[2])
rows = [("HER2", "HER2", PINK), ("HER2 pY1248", "HER2 active (pY1248)", PINK), ("HER3", "HER3", PINK), ("EGFR", "EGFR", PINK),
        ("AR (control)", "AR", TEAL), ("INPP4B (control)", "INPP4B", TEAL)][::-1]
for i, (k, lab, colr) in enumerate(rows):
    r = pt_all.loc[k]
    if colr == PINK:
        ax.add_patch(plt.Rectangle((r.null_hi, i - .32), .75 - r.null_hi, .64, color=PINK_L, alpha=.45, lw=0))
    ax.plot([r.null_lo, r.null_hi], [i, i], color=GREY_L, lw=12, solid_capstyle="butt")
    ax.scatter(r.rho_score, i, s=110, color=colr, zorder=3, edgecolor="white", lw=1.5)
    ax.text(r.rho_score, i + .36, f"{r.rho_score:.2f}", ha="center", color=colr, fontsize=12.5, fontweight="bold")
ax.set_yticks(range(len(rows)), [x_[1] for x_ in rows])
ax.axvline(0, color="#CDBFC6", lw=1)
ax.set_xlim(-.4, .75); ax.set_ylim(-.6, len(rows) - .3)
ax.set_xlabel("correlation with the score")
ax.set_title("Should be in the pink zone.\nHER2 family is not; AR is", loc="left")
ax.text(.66, 4.05, "expected\nzone", color=PINK, fontsize=12, ha="center", fontweight="bold")
fig.text(.985, -.04, "grey bar = what random gene sets reach (chance)   ·   pink zone = where a real HER2 link would have to fall",
         fontsize=12, color=MUTED, ha="right")
save(fig, "fig_protein_expected")
