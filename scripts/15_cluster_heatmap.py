"""HER-OMI step 15: a clustered heatmap of the signature genes across all TCGA TNBC (a presentation figure).

Rows: the 100 amplicon-free signature genes (33 up, 67 down), z-scored across tumours and clustered
within each block (average linkage, correlation distance). Columns: 159 TCGA TNBC, clustered by
Ward linkage on the same genes. Tracks above the heatmap show the HER2-like score, PAM50, TNBCtype
(LAR and the ER-low flag), ERBB2 mRNA and AR mRNA, so the reader can see what the tumour clusters are.
The script also cuts the tumour tree into two clusters and prints what each contains.
Input: data/tnbctype/TCGA_result.csv (TNBCtype output, see step 7) and the files from steps 2 and 2c.
"""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster, dendrogram
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
sys.path.insert(0, "scripts")
from heromi_lib import load_expression

ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")
a = ann[ann.cohort.eq("TCGA")].copy()
sig = pd.read_csv("data/signature_coverage.csv")
core = sig[sig.TCGA & sig.METABRIC & ~sig.on_her2_amplicon].copy()
expr = np.log2(load_expression("data/raw/brca_tcga/data_mrna_seq_v2_rsem.txt", set(a.index)) + 1)
sym = pd.read_csv("data/raw/brca_tcga/data_mrna_seq_v2_rsem.txt", sep="\t", usecols=[0, 1]).dropna()
s2e = dict(zip(sym.Hugo_Symbol, sym.Entrez_Gene_Id.astype(int)))
core["entrez"] = core.entrez_id.fillna(core.gene.map(s2e))   # Entrez first, as in step 2 (catches renamed genes)
core = core.dropna(subset=["entrez"])
core = core[core.entrez.astype(int).isin(expr.index)]
a = a.loc[expr.columns]

tt = pd.read_csv("data/tnbctype/TCGA_result.csv", index_col=0)
tt.index = tt.index.str.replace(".", "-", regex=False)
a["tnbctype"] = [tt.subtype.get(s, "ER-low flag") for s in a.index]
a["basal"] = a.pam50.eq("Basal")

blocks = []
for d in ["Upregulated", "Downregulated"]:
    g = core[core.direction.eq(d)]
    m = expr.loc[g.entrez.astype(int)]
    m.index = g.gene.values
    z = m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1), axis=0).clip(-2.5, 2.5)
    order = leaves_list(linkage(z.values, "average", metric="correlation"))
    blocks.append((d, z.iloc[order]))
Z = pd.concat([b for _, b in blocks])

col_link = linkage(Z.T.values, "ward")
col_order = leaves_list(col_link)
a["cluster"] = fcluster(col_link, 2, "maxclust")
summary = a.groupby("cluster").agg(n=("score_no_amp", "size"), median_score=("score_no_amp", "median"),
                                    LAR=("tnbctype", lambda s: int((s == "LAR").sum())),
                                    ER_low=("tnbctype", lambda s: int((s == "ER-low flag").sum())),
                                    basal=("basal", "sum"), median_ERBB2=("ERBB2_expr", lambda s: np.log2(s + 1).median()),
                                    median_AR=("AR", "median"))
print(summary.round(3).to_string())
print("score differs between clusters: Mann-Whitney p =",
      f"{stats.mannwhitneyu(*[g.score_no_amp for _, g in a.groupby('cluster')]).pvalue:.2g}")
summary.to_csv("results/cluster_heatmap_summary.csv")

# ---------------- figure ----------------
BG, INK, MUTED = "#FFFBFC", "#3A1A2A", "#6B5A63"
div = LinearSegmentedColormap.from_list("bc", ["#0097A7", "#FFFBFC", "#C2185B"])
seq = LinearSegmentedColormap.from_list("seq", ["#FBEFF4", "#8E1244"])
plt.rcParams.update({"font.size": 15, "text.color": INK})
A = a.iloc[col_order]
cl = A.cluster.to_numpy()
b = int(np.flatnonzero(cl != cl[0])[0])                     # first column of the second cluster
assert (cl[b:] == cl[b]).all(), "clusters are not contiguous in the dendrogram order"
nL, nR = b, len(cl) - b
tracks = [("Cluster", (cl != cl[0]).astype(int), ListedColormap(["#0097A7", "#CDBFC6"]), None),
          ("HER2-like score", A.score_no_amp.rank().to_numpy(), seq, None),
          ("ERBB2 mRNA", A.ERBB2_expr.rank().to_numpy(), seq, None),
          ("AR mRNA", A.AR.rank().to_numpy(), seq, None),
          ("PAM50 basal", A.basal.astype(int).to_numpy(), ListedColormap(["#EDE3E8", "#7A6A73"]), None),
          ("TNBCtype LAR", np.select([A.tnbctype.eq("LAR"), A.tnbctype.eq("ER-low flag")], [2, 1], 0),
           ListedColormap(["#F6EEF2", "#9FD8DF", "#0097A7"]), None)]
fig = plt.figure(figsize=(18, 10.5), facecolor=BG)
gs = fig.add_gridspec(3 + len(tracks), 1, height_ratios=[1.6, 0.55] + [0.28] * (len(tracks) - 1) + [0.15, 6.5], hspace=0.08,
                      left=0.17, right=0.985, top=0.95, bottom=0.07)
ax = fig.add_subplot(gs[0])
dendrogram(col_link, ax=ax, color_threshold=0, above_threshold_color="#A08F98", no_labels=True)
ax.axis("off")
ax.axvline(10 * b, color=INK, lw=2.5, ymax=0.75)

ax.set_title("159 TCGA triple-negative tumours, clustered on the 100 signature genes", loc="left", fontsize=19,
             color=INK, fontweight="bold", x=-0.155)
for i, (name, vals, cmap, _) in enumerate(tracks):
    t = fig.add_subplot(gs[1 + i])
    t.imshow(vals[None, :], aspect="auto", cmap=cmap, interpolation="nearest")
    t.set_xticks([]); t.set_yticks([])
    for s in t.spines.values():
        s.set_visible(False)
    t.text(-0.008, 0.5, name, transform=t.transAxes, ha="right", va="center", fontsize=14, color=INK,
           fontweight="bold" if i == 0 else "normal")
    t.axvline(b - 0.5, color=INK, lw=2.5)
    if i == 0:
        t.text(nL / 2 - 0.5, 0, f"LEFT: {nL}", ha="center", va="center", fontsize=14, color="white", fontweight="bold")
        t.text(b + nR / 2 - 0.5, 0, f"RIGHT CLUSTER: {nR} tumours, mostly basal, low score, little AR",
               ha="center", va="center", fontsize=14, color=INK, fontweight="bold")
h = fig.add_subplot(gs[-1])
h.imshow(Z.iloc[:, col_order].values, aspect="auto", cmap=div, vmin=-2.5, vmax=2.5, interpolation="nearest")
h.set_xticks([]); h.set_yticks([])
n_up = len(blocks[0][1])
h.axhline(n_up - 0.5, color=INK, lw=1.2)
h.axvline(b - 0.5, color=INK, lw=2.5)
h.add_patch(plt.Rectangle((-0.5, -0.5), nL, len(Z), fill=False, edgecolor="#0097A7", lw=3.5))
h.text(-0.01, 1 - (n_up / 2) / len(Z), f"{n_up} genes up\nin HER2-like", transform=h.transAxes, ha="right", va="center",
       fontsize=15, color=INK)
h.text(-0.01, (len(Z) - n_up) / 2 / len(Z), f"{len(Z) - n_up} genes down\nin HER2-like", transform=h.transAxes,
       ha="right", va="center", fontsize=15, color=INK)
for s in h.spines.values():
    s.set_visible(False)
h.set_xlabel("Each column is one tumour. Teal = lower than average, pink = higher than average (z-score).",
             fontsize=13.5, color=MUTED, labelpad=10)
fig.text(0.985, 0.012, "Tracks: darker = higher rank. TNBCtype: dark teal = LAR, light teal = ER-low flag. Step 15.",
         ha="right", fontsize=12.5, color=MUTED)
fig.savefig("results/cluster_heatmap_tcga.png", dpi=150, facecolor=BG)
print("saved results/cluster_heatmap_tcga.png")
