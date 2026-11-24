"""HER-OMI step 10: HER2-like TNBC next to true HER2-positive tumours, layer by layer (TCGA).

Groups (TCGA legacy study, primary tumours with RNA-seq):
  HER2-positive, HR-negative : HER2 IHC positive or FISH positive, ER and PR negative by IHC
  TNBC, HER2-like            : top third of the amplicon-free score among our 159 TNBC
  TNBC, basal-like           : bottom third
Features, grouped by omics layer: ERBB2 copy number (DNA); ERBB2, GRB7, the amplicon-free HER2-like score,
AR, FOXA1, MKI67 (RNA); HER2, HER2 pY1248, AR, INPP4B, AKT pS473 (RPPA protein).
Each row is z-scored across the shown tumours and clipped at +/-2. Output: a heatmap for the slides and a
table of group medians.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, "scripts")
from heromi_lib import load_expression, singscore

RAW = "data/raw"
ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")
tn = ann[ann.cohort.eq("TCGA")].copy()
tn["third"] = pd.qcut(tn.score_no_amp, 3, labels=["low", "mid", "high"])

clin = pd.read_csv(f"{RAW}/brca_tcga/data_clinical_patient.txt", sep="\t", comment="#").set_index("PATIENT_ID")
cols = pd.read_csv(f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", sep="\t", nrows=0).columns[2:]
prim = pd.Series([c for c in cols if c[13:15] == "01"])
pc = clin.reindex(prim.str[:12])
pc.index = prim.values
her2pos = (pc.IHC_HER2.eq("Positive") | pc.HER2_FISH_STATUS.eq("Positive")) & \
          pc.ER_STATUS_BY_IHC.eq("Negative") & pc.PR_STATUS_BY_IHC.eq("Negative")
groups = {"HER2-positive\n(HR-negative)": list(pc.index[her2pos]),
          "TNBC\nHER2-like (top third)": list(tn.index[tn.third.eq("high")]),
          "TNBC\nbasal-like (bottom third)": list(tn.index[tn.third.eq("low")])}
allids = sum(groups.values(), [])

expr = load_expression(f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", set(allids))
sig = pd.read_csv("data/signature_coverage.csv")
core = sig[sig.TCGA & sig.METABRIC & ~sig.on_her2_amplicon]
sym = pd.read_csv(f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", sep="\t", usecols=[0, 1]).dropna().drop_duplicates("Hugo_Symbol")
s2e = dict(zip(sym.Hugo_Symbol, sym.Entrez_Gene_Id.astype(int)))
ids = lambda genes: [s2e[g] for g in genes if g in s2e and s2e[g] in expr.index]
score = singscore(expr, ids(core.gene[core.direction.eq("Upregulated")]),
                  ids(core.gene[core.direction.eq("Downregulated")])).total
lg = np.log2(expr + 1)

seg = pd.read_csv(f"{RAW}/brca_tcga/data_cna_hg19.seg", sep="\t")
seg.columns = ["ID", "chrom", "start", "end", "nmark", "mean"]
seg = seg[seg.chrom.astype(str).eq("17") & (seg.start <= 37_865_000) & (seg.end >= 37_865_000)]
cn = seg.assign(sid=seg.ID.str[:15]).groupby("sid")["mean"].first()

rppa = pd.read_csv(f"{RAW}/brca_tcga/data_rppa.txt", sep="\t").set_index("Composite.Element.REF").T.astype(float)

F = pd.DataFrame(index=expr.columns)
F["ERBB2 copy number"] = [cn.get(s[:15], np.nan) for s in F.index]
for g in ["ERBB2", "GRB7"]:
    F[f"{g} mRNA"] = lg.loc[s2e[g]]
F["HER2-like score (no amplicon genes)"] = score
for g in ["AR", "FOXA1", "MKI67"]:
    F[f"{g} mRNA"] = lg.loc[s2e[g]]
for ab, lab in [("ERBB2|HER2", "HER2 protein"), ("ERBB2|HER2_pY1248", "HER2 pY1248 (active)"),
                ("AR|AR", "AR protein"), ("INPP4B|INPP4B", "INPP4B protein"), ("AKT1 AKT2 AKT3|Akt_pS473", "AKT pS473")]:
    F[lab] = [rppa.at[s, ab] if s in rppa.index else np.nan for s in F.index]
LAYERS = [("DNA", ["ERBB2 copy number"]),
          ("RNA", ["ERBB2 mRNA", "GRB7 mRNA", "HER2-like score (no amplicon genes)", "AR mRNA", "FOXA1 mRNA", "MKI67 mRNA"]),
          ("Protein", ["HER2 protein", "HER2 pY1248 (active)", "AR protein", "INPP4B protein", "AKT pS473"])]
rows = [r for _, rs in LAYERS for r in rs]

# order columns: by group, then by score within group
order = []
for g, members in groups.items():
    m = [s for s in members if s in F.index]
    order += list(F.loc[m].sort_values("HER2-like score (no amplicon genes)", ascending=False).index)
X = F.loc[order, rows]
Z = ((X - X.mean()) / X.std()).clip(-2, 2).T

med = pd.DataFrame({g.replace("\n", " "): X.loc[[s for s in groups[g] if s in X.index]].median() for g in groups})
nn = {g.replace("\n", " "): sum(s in X.index for s in groups[g]) for g in groups}
med.loc["n tumours"] = pd.Series(nn)
med.to_csv("results/her2pos_comparison_medians.csv")
pd.set_option("display.width", 200)
print(med.round(2).to_string())
print("protein available:", {g.replace(chr(10), ' '): int(X.loc[[s for s in groups[g] if s in X.index], 'HER2 protein'].notna().sum()) for g in groups})

# ---- figure ----
cmap = LinearSegmentedColormap.from_list("div", ["#00505A", "#0097A7", "#9FD8DF", "#FBF3F6", "#F3B6CE", "#C2185B", "#7A0E3A"])
cmap.set_bad("#D9CFD4")
fig = plt.figure(figsize=(13, 6.6), facecolor="#FFFBFC")
AX_B, AX_H = 0.05, 0.70
ax = fig.add_axes([0.25, AX_B, 0.66, AX_H])
ax.set_facecolor("#FFFBFC")
im = ax.imshow(np.ma.masked_invalid(Z.values), aspect="auto", cmap=cmap, vmin=-2, vmax=2, interpolation="nearest")
ax.set_yticks(range(len(rows)))
ax.set_yticklabels(rows, fontsize=10.5, color="#3A1A2A")
ax.set_xticks([])
for s in ax.spines.values():
    s.set_visible(False)
x0 = 0
for g, members in groups.items():
    n = sum(s in X.index for s in members)
    if x0:
        ax.axvline(x0 - 0.5, color="#FFFBFC", lw=4)
    ax.text(x0 + n / 2 - 0.5, -1.0, f"{g}\nn = {n}", ha="center", va="bottom", fontsize=11, color="#3A1A2A",
            fontweight="bold")
    x0 += n
y0 = 0
for layer, rs in LAYERS:
    if y0:
        ax.axhline(y0 - 0.5, color="#FFFBFC", lw=4)
    fig.text(0.015, AX_B + AX_H - AX_H * (y0 + len(rs) / 2 - 0.5 + 0.5) / len(rows), layer, fontsize=12, fontweight="bold",
             color="#6B5A63", va="center")
    y0 += len(rs)
cax = fig.add_axes([0.93, 0.22, 0.012, 0.36])
cb = fig.colorbar(im, cax=cax)
cb.set_ticks([-2, 0, 2])
cb.set_ticklabels(["low", "average", "high"])
cb.ax.tick_params(labelsize=9, colors="#6B5A63")
cb.outline.set_visible(False)
fig.text(0.25, 0.962, "HER2-like TNBC share HER2-positive tumours' AR/luminal programme, not their HER2",
         fontsize=13.5, fontweight="bold", color="#3A1A2A")
fig.text(0.25, 0.93, "TCGA primary tumours. Each row scaled across all tumours shown; grey = no RPPA protein data. "
         "Columns sorted by HER2-like score within each group.", fontsize=9.5, color="#6B5A63")
fig.savefig("results/her2pos_comparison_heatmap.png", dpi=200, facecolor=fig.get_facecolor())
print("saved results/her2pos_comparison_heatmap.png")
