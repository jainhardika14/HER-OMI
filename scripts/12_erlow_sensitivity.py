"""HER-OMI step 12: ER-low sensitivity reruns of step 4 (HER2 IHC trend) and step 5 (pathway agreement with S22).

Decision 3 kept the TNBCtype ER-flagged (ER-low) tumours in the main analysis, with this rerun as the check.
Flagged = in the original TNBCtype upload but not the cleaned one (12 TCGA, 11 METABRIC).
"""
import numpy as np
import pandas as pd
from scipy import stats

ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")
up = "results/tnbctype_upload"
FLAG = set()
for c in ("TCGA", "METABRIC"):
    FLAG |= set(pd.read_csv(f"{up}/{c}_TNBC_expression.csv", nrows=0).columns) - \
            set(pd.read_csv(f"{up}/{c}_TNBC_expression_noERflag.csv", nrows=0).columns)
ann["er_flag"] = ann.index.isin(FLAG)

# ---------------- step 4: HER2 IHC trend (TCGA) ----------------
clin = pd.read_csv("data/raw/brca_tcga/data_clinical_patient.txt", sep="\t", comment="#").set_index("PATIENT_ID")
t = ann[ann.cohort.eq("TCGA")].copy()
t["ihc"] = pd.to_numeric(t.pid.map(clin.HER2_IHC_SCORE), errors="coerce")
t = t[t.ihc.isin([0, 1, 2])]
rows4 = []
for version, g0 in [("with ER-flag", t), ("ER-flag removed", t[~t.er_flag])]:
    for subset, g in [("all TNBC", g0), ("basal", g0[g0.pam50.eq("Basal")])]:
        rho, p = stats.spearmanr(g.ihc, g.score_no_amp)
        low, zero = g[g.ihc.isin([1, 2])].score_no_amp, g[g.ihc.eq(0)].score_no_amp
        rows4.append(dict(version=version, subset=subset, n=len(g), n_her2_0=len(zero), trend_rho=rho, trend_p=p,
                          her2low_vs_0_p=stats.mannwhitneyu(low, zero).pvalue))
r4 = pd.DataFrame(rows4)

# ---------------- step 5: pathway profile agreement with IBCGA S22 ----------------
s22 = pd.read_csv("data/ibcga_S22.csv").set_index("set")
ref, padj = -s22["diff"], s22.padj
SIG = set(pd.read_csv("data/IBCGA_HER2like_signature_144genes.csv").gene)
sets = {}
for f in ["data/ref/h.all.gmt", "data/ref/c2.cp.reactome.gmt"]:
    for line in open(f):
        name, _, *genes = line.rstrip("\n").split("\t")
        sets[name] = [x for x in genes if x not in SIG]


def pathway_scores(expr):
    R = expr.rank(axis=0).to_numpy(); N = R.shape[0]
    row = {g: i for i, g in enumerate(expr.index)}
    out = {}
    for name, genes in sets.items():
        idx = [row[g] for g in genes if g in row]; n = len(idx)
        if 10 <= n <= 500:
            lo, hi = (n + 1) / 2, (2 * N - n + 1) / 2
            out[name] = (R[idx].mean(0) - lo) / (hi - lo) - 0.5
    return pd.DataFrame(out, index=expr.columns)


rows5 = []
for c in ("TCGA", "METABRIC"):
    expr = pd.read_csv(f"{up}/{c}_TNBC_expression.csv", index_col=0)
    ps = pathway_scores(expr)
    meta = ann.loc[ps.index]
    for version, keep in [("with ER-flag", np.ones(len(meta), bool)), ("ER-flag removed", ~meta.er_flag.to_numpy())]:
        for subset, m in [("all TNBC", keep), ("basal", keep & meta.pam50.eq("Basal").to_numpy())]:
            prof = ps[m].corrwith(meta.loc[m, "score_no_amp"], method="spearman")
            j = pd.DataFrame({"ours": prof, "s22": ref, "padj": padj}).dropna()
            sg = j[j.padj < 0.05]
            rows5.append(dict(cohort=c, version=version, subset=subset, n=int(m.sum()),
                              rho_all=j.ours.corr(j.s22, method="spearman"),
                              rho_fdr05=sg.ours.corr(sg.s22, method="spearman"),
                              sign_agreement=(np.sign(sg.ours) == np.sign(sg.s22)).mean()))
r5 = pd.DataFrame(rows5)
r4.to_csv("results/erlow_sensitivity_step4.csv", index=False)
r5.to_csv("results/erlow_sensitivity_step5.csv", index=False)
pd.set_option("display.width", 200)
print("Step 4, HER2 IHC trend (TCGA):"); print(r4.round(3).to_string(index=False))
print("\nStep 5, agreement with IBCGA S22:"); print(r5.round(3).to_string(index=False))
