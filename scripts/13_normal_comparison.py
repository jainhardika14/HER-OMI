"""HER-OMI step 13: normal breast vs TNBC vs HER2-positive, side by side (TCGA, RNA).

Asked for by Prof. Choudhary: put TNBC next to normal breast tissue and next to HER2-positive tumours.
Data: TCGA PanCancer Atlas RNA-seq, tumours and adjacent normal breast processed by the same pipeline
(the legacy matrix used in steps 2 to 12 has no normal samples). All values log2(RSEM + 1).

Two views
  A. Groups. Normal breast (all 114 adjacent normals), TNBC basal-like (bottom score third),
     TNBC HER2-like (top score third), HER2-positive HR-negative (same definition as step 10).
     Thirds are the step 10 thirds, so the groups match the heatmap. Mann-Whitney tests between pairs.
  B. Pairs. Patients whose own normal tissue was also sequenced: 13 of our TNBC and 24 HER2-positive
     patients (any HR status; only 4 are HR-negative). Tumour minus own normal per gene, Wilcoxon signed rank.
The HER2-like score (amplicon-free, 100 genes) is recomputed here on this matrix for every sample,
normals included; singscore is rank based within a sample, so this is legitimate.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, "scripts")
from heromi_lib import load_expression, singscore

P = "data/raw/brca_tcga_pan_can_atlas_2018"
GENES = {"ERBB2": 2064, "GRB7": 2886, "ESR1": 2099, "PGR": 5241, "AR": 367, "FOXA1": 3169, "GATA3": 2625,
         "KRT5": 3852, "EGFR": 1956, "MKI67": 4288}

ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")
tn = ann[ann.cohort.eq("TCGA")].copy()
tn["third"] = pd.qcut(tn.score_no_amp, 3, labels=["low", "mid", "high"])   # identical to step 10
clin = pd.read_csv("data/raw/brca_tcga/data_clinical_patient.txt", sep="\t", comment="#").set_index("PATIENT_ID")
her2 = clin.IHC_HER2.eq("Positive") | clin.HER2_FISH_STATUS.eq("Positive")
hrneg = clin.ER_STATUS_BY_IHC.eq("Negative") & clin.PR_STATUS_BY_IHC.eq("Negative")

tum_cols = pd.read_csv(f"{P}/data_mrna_seq_v2_rsem.txt", sep="\t", nrows=0).columns[2:]
nor_cols = pd.read_csv(f"{P}/normals/data_mrna_seq_v2_rsem_normal_samples.txt", sep="\t", nrows=0).columns[2:]
tum_by_pid = {c[:12]: c for c in tum_cols}
groups = {
    "Normal breast": list(nor_cols),
    "TNBC basal-like": [s for s in tn.index[tn.third.eq("low")] if s in tum_cols],
    "TNBC HER2-like": [s for s in tn.index[tn.third.eq("high")] if s in tum_cols],
    "HER2+ (HR-negative)": [tum_by_pid[p] for p in clin.index[her2 & hrneg] if p in tum_by_pid],
}
pairs = {"TNBC": [p for p in tn.pid.unique() if p in tum_by_pid and f"{p}-11" in nor_cols],
         "HER2+": [p for p in clin.index[her2] if p in tum_by_pid and f"{p}-11" in nor_cols]}
need_t = set(sum([g for k, g in groups.items() if k != "Normal breast"], [])) | \
         {tum_by_pid[p] for v in pairs.values() for p in v}
et = load_expression(f"{P}/data_mrna_seq_v2_rsem.txt", need_t)
en = load_expression(f"{P}/normals/data_mrna_seq_v2_rsem_normal_samples.txt", set(nor_cols))
expr = et.join(en, how="inner")                      # genes measured in every tumour and normal
print(f"matrix: {expr.shape[0]} genes x {expr.shape[1]} samples ({et.shape[1]} tumours, {en.shape[1]} normals)")

sig = pd.read_csv("data/signature_coverage.csv")
core = sig[sig.TCGA & sig.METABRIC & ~sig.on_her2_amplicon]
sym = pd.read_csv(f"{P}/data_mrna_seq_v2_rsem.txt", sep="\t", usecols=[0, 1]).dropna().drop_duplicates("Hugo_Symbol")
s2e = dict(zip(sym.Hugo_Symbol, sym.Entrez_Gene_Id.astype(int)))
ids = lambda genes: [s2e[g] for g in genes if g in s2e and s2e[g] in expr.index]
up, dn = ids(core.gene[core.direction.eq("Upregulated")]), ids(core.gene[core.direction.eq("Downregulated")])
print(f"signature genes found: {len(up)} up, {len(dn)} down")

F = np.log2(expr.loc[list(GENES.values())] + 1)
F.index = list(GENES)
F = F.T
F["HER2-like score"] = singscore(expr, up, dn).total
check = F.loc[[s for s in tn.index if s in F.index], "HER2-like score"]
print(f"PanCan score vs legacy score in our TNBC: rho = {stats.spearmanr(check, tn.loc[check.index, 'score_no_amp'])[0]:.3f}")
FEAT = ["HER2-like score"] + list(GENES)

# ---------------- A. groups ----------------
long = pd.concat([F.loc[g].assign(group=k) for k, g in groups.items()])
long.to_csv("results/normal_comparison_values.csv")
med = long.groupby("group", sort=False)[FEAT].median().T
med.loc["n"] = long.group.value_counts().reindex(med.columns)
med.to_csv("results/normal_comparison_medians.csv")
rows = []
for a, b in [("TNBC HER2-like", "Normal breast"), ("TNBC basal-like", "Normal breast"),
             ("TNBC HER2-like", "HER2+ (HR-negative)"), ("TNBC HER2-like", "TNBC basal-like")]:
    for f in FEAT:
        x, y = F.loc[groups[a], f], F.loc[groups[b], f]
        rows.append(dict(comparison=f"{a} vs {b}", feature=f, median_a=x.median(), median_b=y.median(),
                         diff=x.median() - y.median(), p=stats.mannwhitneyu(x, y).pvalue))
grp = pd.DataFrame(rows)
def bh(p):
    """Benjamini-Hochberg q values, same as step 8b."""
    p = np.asarray(p); o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out


grp["q_BH"] = grp.groupby("comparison").p.transform(bh)
grp.to_csv("results/normal_comparison_groups.csv", index=False)

# ---------------- B. matched pairs ----------------
prow = []
for kind, pids in pairs.items():
    T = F.loc[[tum_by_pid[p] for p in pids]].reset_index(drop=True)
    N = F.loc[[f"{p}-11" for p in pids]].reset_index(drop=True)
    pd.concat([T.assign(pid=pids, tissue="tumour"), N.assign(pid=pids, tissue="normal")]).assign(kind=kind) \
      .to_csv(f"results/normal_pairs_{kind.replace('+', 'pos')}.csv", index=False)
    for f in FEAT:
        d = T[f] - N[f]
        prow.append(dict(pairs=kind, n=len(pids), feature=f, median_tumour=T[f].median(), median_normal=N[f].median(),
                         median_change=d.median(), n_up=int((d > 0).sum()), p_wilcoxon=stats.wilcoxon(T[f], N[f]).pvalue))
pr = pd.DataFrame(prow)
pr.to_csv("results/normal_comparison_pairs.csv", index=False)

pd.set_option("display.width", 220)
print("\nA. Group medians (log2 RSEM + 1; score in singscore units)")
print(med.round(2).to_string())
print("\nA. Tests")
print(grp.round(4).to_string(index=False))
print("\nB. Matched pairs, tumour minus own normal")
print(pr.round(4).to_string(index=False))

# ---------------- figures ----------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BG, INK, MUTED, GRID = "#FFFBFC", "#3A1A2A", "#6B5A63", "#F1E6EB"
COL = {"Normal breast": "#0097A7", "TNBC basal-like": "#8A7D85", "TNBC HER2-like": "#E06C9F",
       "HER2+ (HR-negative)": "#9C1450"}            # validated categorical set (dataviz palette, adjacent order)
LAB = {"Normal breast": "Normal\nbreast", "TNBC basal-like": "TNBC\nbasal-like", "TNBC HER2-like": "TNBC\nHER2-like",
       "HER2+ (HR-negative)": "HER2+\nHR-neg"}
plt.rcParams.update({"font.size": 15, "axes.edgecolor": MUTED, "axes.labelcolor": MUTED, "xtick.color": INK,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(607)


def tidy(ax):
    ax.set_facecolor(BG)
    ax.yaxis.grid(True, color=GRID, lw=1)
    ax.set_axisbelow(True)


# A. groups
panels = [("HER2-like score", "HER2-like score"), ("ERBB2", "ERBB2 mRNA"), ("AR", "AR mRNA"),
          ("FOXA1", "FOXA1 mRNA"), ("ESR1", "ESR1 mRNA"), ("MKI67", "MKI67 mRNA (proliferation)")]
fig, axes = plt.subplots(2, 3, figsize=(18, 11), facecolor=BG)
for ax, (f, title) in zip(axes.flat, panels):
    tidy(ax)
    for i, (g, members) in enumerate(groups.items()):
        v = F.loc[members, f].to_numpy()
        ax.scatter(i + rng.uniform(-0.22, 0.22, len(v)), v, s=22, color=COL[g], alpha=0.55, lw=0)
        q1, m, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([i - 0.3, i + 0.3], [m, m], color=INK, lw=3, solid_capstyle="round")
        ax.plot([i, i], [q1, q3], color=INK, lw=1.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels([LAB[g] for g in groups], fontsize=13.5)
    ax.set_title(title, fontsize=17, color=INK, loc="left", pad=12, fontweight="bold")
    ax.set_xlim(-0.6, 3.6)
axes[0, 0].set_ylabel("singscore")
for ax in list(axes[0, 1:]) + list(axes[1]):
    ax.set_ylabel("log2(RSEM + 1)")
fig.text(0.01, 0.015, f"TCGA PanCancer RNA-seq. n: normal {len(groups['Normal breast'])}, TNBC basal-like "
         f"{len(groups['TNBC basal-like'])}, TNBC HER2-like {len(groups['TNBC HER2-like'])}, HER2+ HR-negative "
         f"{len(groups['HER2+ (HR-negative)'])}. Dots = tumours or normal samples; black bar = median, line = middle half.",
         fontsize=12.5, color=MUTED)
fig.tight_layout(rect=(0, 0.035, 1, 1), w_pad=2.5, h_pad=3)
fig.savefig("results/normal_comparison_groups.png", dpi=150, facecolor=BG)

# B. matched pairs
fig, axes = plt.subplots(2, 3, figsize=(18, 10.5), facecolor=BG)
for r, (kind, colr) in enumerate([("TNBC", COL["TNBC HER2-like"]), ("HER2+", COL["HER2+ (HR-negative)"])]):
    d = pd.read_csv(f"results/normal_pairs_{kind.replace('+', 'pos')}.csv")
    for c, (f, title) in enumerate([("ERBB2", "ERBB2 mRNA"), ("HER2-like score", "HER2-like score"), ("AR", "AR mRNA")]):
        ax = axes[r, c]
        tidy(ax)
        w = d.pivot(index="pid", columns="tissue", values=f)
        for _, row in w.iterrows():
            ax.plot([0, 1], [row.normal, row.tumour], color="#D6C8CF", lw=1.4, marker="o", ms=5, mfc="#D6C8CF", mec="none")
        ax.plot([0, 1], [w.normal.median(), w.tumour.median()], color=colr, lw=4, marker="o", ms=11, mec=BG, mew=2)
        s = pr[(pr.pairs == kind) & (pr.feature == f)].iloc[0]
        ax.set_title(f"{kind}: {title}", fontsize=17, color=INK, loc="left", fontweight="bold", pad=34)
        ax.text(0, 1.03, f"{s.n_up} of {s.n} patients go up · Wilcoxon p = {s.p_wilcoxon:.2g}", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=13, color=MUTED)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["own normal", "tumour"])
        ax.set_xlim(-0.3, 1.3)
fig.text(0.01, 0.01, f"TCGA patients whose adjacent normal breast was also sequenced: {len(pairs['TNBC'])} TNBC, "
         f"{len(pairs['HER2+'])} HER2-positive (any HR status). Grey = one patient; coloured = median.",
         fontsize=12.5, color=MUTED)
fig.tight_layout(rect=(0, 0.04, 1, 1), h_pad=2.5, w_pad=2.5)
fig.savefig("results/normal_comparison_pairs.png", dpi=150, facecolor=BG)
print("saved results/normal_comparison_groups.png and results/normal_comparison_pairs.png")
