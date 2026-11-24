"""HER-OMI step 10b: the HER2-positive comparison repeated in METABRIC (DNA and RNA only; no protein data).

Groups, defined with the same tests as our METABRIC TNBC rule (Decision 1, option 3):
  HER2-positive, ER-negative : ERBB2 GAIN on SNP6, ER negative by IHC, PR negative
  TNBC, HER2-like            : top third of the amplicon-free score among our 250 TNBC
  TNBC, basal-like           : bottom third
Features: ERBB2 copy number (GISTIC call, -2 to 2), ERBB2, GRB7 mRNA, amplicon-free score, AR, FOXA1, MKI67 mRNA.
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
tn = ann[ann.cohort.eq("METABRIC")].copy()
tn["third"] = pd.qcut(tn.score_no_amp, 3, labels=["low", "mid", "high"])
mp = pd.read_csv(f"{RAW}/brca_metabric/data_clinical_patient.txt", sep="\t", comment="#")
ms = pd.read_csv(f"{RAW}/brca_metabric/data_clinical_sample.txt", sep="\t", comment="#")
m = mp.merge(ms, on="PATIENT_ID").set_index("PATIENT_ID")
her2pos = m.HER2_SNP6.eq("GAIN") & m.ER_IHC.eq("Negative") & m.PR_STATUS.eq("Negative")
groups = {"HER2-positive\n(ER-negative)": list(m.index[her2pos]),
          "TNBC\nHER2-like (top third)": list(tn.index[tn.third.eq("high")]),
          "TNBC\nbasal-like (bottom third)": list(tn.index[tn.third.eq("low")])}
allids = sum(groups.values(), [])
path = f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt"
expr = load_expression(path, set(allids))
sig = pd.read_csv("data/signature_coverage.csv")
core = sig[sig.TCGA & sig.METABRIC & ~sig.on_her2_amplicon]
sym = pd.read_csv(path, sep="\t", usecols=[0, 1]).dropna().drop_duplicates("Hugo_Symbol")
s2e = dict(zip(sym.Hugo_Symbol, sym.Entrez_Gene_Id.astype(int)))
ids = lambda genes: [s2e[g] for g in genes if g in s2e and s2e[g] in expr.index]
score = singscore(expr, ids(core.gene[core.direction.eq("Upregulated")]),
                  ids(core.gene[core.direction.eq("Downregulated")])).total
cna = pd.read_csv(f"{RAW}/brca_metabric/data_cna.txt", sep="\t")
cna = cna[cna.Hugo_Symbol.eq("ERBB2")].drop(columns=["Hugo_Symbol", "Entrez_Gene_Id"]).iloc[0]

F = pd.DataFrame(index=expr.columns)
F["ERBB2 copy number (GISTIC)"] = cna.reindex(F.index).astype(float)
for g in ["ERBB2", "GRB7"]:
    F[f"{g} mRNA"] = expr.loc[s2e[g]]
F["HER2-like score (no amplicon genes)"] = score
for g in ["AR", "FOXA1", "MKI67"]:
    F[f"{g} mRNA"] = expr.loc[s2e[g]]
LAYERS = [("DNA", ["ERBB2 copy number (GISTIC)"]),
          ("RNA", ["ERBB2 mRNA", "GRB7 mRNA", "HER2-like score (no amplicon genes)", "AR mRNA", "FOXA1 mRNA", "MKI67 mRNA"])]
rows = [r for _, rs in LAYERS for r in rs]
order = []
for g, members in groups.items():
    mm = [s for s in members if s in F.index]
    order += list(F.loc[mm].sort_values("HER2-like score (no amplicon genes)", ascending=False).index)
X = F.loc[order, rows]
Z = ((X - X.mean()) / X.std()).clip(-2, 2).T
med = pd.DataFrame({g.replace("\n", " "): X.loc[[s for s in groups[g] if s in X.index]].median() for g in groups})
med.loc["n tumours"] = pd.Series({g.replace("\n", " "): sum(s in X.index for s in groups[g]) for g in groups})
med.to_csv("results/her2pos_comparison_medians_metabric.csv")
pd.set_option("display.width", 200)
print(med.round(2).to_string())

cmap = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#b7d3f6", "#f0efec", "#f4bdbb", "#e34948", "#9c2121"])
fig = plt.figure(figsize=(13, 4.9), facecolor="#fcfcfb")
AX_B, AX_H = 0.06, 0.62
ax = fig.add_axes([0.25, AX_B, 0.66, AX_H])
im = ax.imshow(Z.values, aspect="auto", cmap=cmap, vmin=-2, vmax=2, interpolation="nearest")
ax.set_yticks(range(len(rows)))
ax.set_yticklabels(rows, fontsize=10.5, color="#0b0b0b")
ax.set_xticks([])
for s in ax.spines.values():
    s.set_visible(False)
x0 = 0
for g, members in groups.items():
    n = sum(s in X.index for s in members)
    if x0:
        ax.axvline(x0 - 0.5, color="#fcfcfb", lw=4)
    ax.text(x0 + n / 2 - 0.5, -0.9, f"{g}\nn = {n}", ha="center", va="bottom", fontsize=11, color="#0b0b0b", fontweight="bold")
    x0 += n
y0 = 0
for layer, rs in LAYERS:
    if y0:
        ax.axhline(y0 - 0.5, color="#fcfcfb", lw=4)
    fig.text(0.015, AX_B + AX_H - AX_H * (y0 + len(rs) / 2) / len(rows), layer, fontsize=12, fontweight="bold",
             color="#52514e", va="center")
    y0 += len(rs)
cax = fig.add_axes([0.93, 0.15, 0.012, 0.4])
cb = fig.colorbar(im, cax=cax)
cb.set_ticks([-2, 0, 2])
cb.set_ticklabels(["low", "average", "high"])
cb.ax.tick_params(labelsize=9, colors="#52514e")
cb.outline.set_visible(False)
fig.text(0.25, 0.95, "METABRIC repeats it: HER2-like TNBC share the AR/luminal programme, not ERBB2",
         fontsize=13.5, fontweight="bold", color="#0b0b0b")
fig.text(0.25, 0.905, "METABRIC primary tumours (no protein data). Each row scaled across all tumours shown. "
         "Columns sorted by HER2-like score within each group.", fontsize=9.5, color="#52514e")
fig.savefig("results/her2pos_comparison_heatmap_metabric.png", dpi=200, facecolor=fig.get_facecolor())
print("saved results/her2pos_comparison_heatmap_metabric.png")
