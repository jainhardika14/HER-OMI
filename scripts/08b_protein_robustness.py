"""HER-OMI step 8b: robustness of the protein result (TCGA RPPA).

  1. Remove the 12 TCGA tumours TNBCtype flagged as possibly ER positive (ER-low): does the pattern hold?
  2. Restrict to HER2-low TNBC (IHC 1+ or 2+, the trastuzumab deruxtecan-eligible group): does the score
     pick out tumours with more HER2 protein there?
  3. Benjamini-Hochberg correction across proteins within each subset.
Same method as step 8: Spearman vs 1000 matched random gene sets (two-sided), partial given ERBB2 mRNA.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, "scripts")
from heromi_lib import load_expression, partial_spearman, random_set_scores

RAW = "data/raw"
sig = pd.read_csv("data/signature_coverage.csv")
core = sig[sig.TCGA & sig.METABRIC & ~sig.on_her2_amplicon]
n_up, n_dn = int(core.direction.eq("Upregulated").sum()), int(core.direction.eq("Downregulated").sum())
ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")
up = "results/tnbctype_upload"
FLAG = set(pd.read_csv(f"{up}/TCGA_TNBC_expression.csv", nrows=0).columns) - \
       set(pd.read_csv(f"{up}/TCGA_TNBC_expression_noERflag.csv", nrows=0).columns)
PROTEINS = {"ERBB2|HER2": "HER2", "ERBB2|HER2_pY1248": "HER2 pY1248", "AR|AR": "AR",
            "INPP4B|INPP4B": "INPP4B", "AKT1 AKT2 AKT3|Akt_pS473": "AKT pS473", "EGFR|EGFR_pY1068": "EGFR pY1068"}

rppa = pd.read_csv(f"{RAW}/brca_tcga/data_rppa.txt", sep="\t").set_index("Composite.Element.REF").T
tc = ann[ann.cohort.eq("TCGA")]
ids = [s for s in tc.index if s in rppa.index]
expr = load_expression(f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", set(ids))
ids = list(expr.columns)
sig_ids = set(sig.loc[sig.entrez_id.notna(), "entrez_id"].astype(int))
nulls = pd.DataFrame(random_set_scores(expr, sig_ids, n_up, n_dn).T, index=ids)
clin = pd.read_csv(f"{RAW}/brca_tcga/data_clinical_patient.txt", sep="\t", comment="#").set_index("PATIENT_ID")
d = tc.loc[ids].join(rppa.loc[ids, list(PROTEINS)].astype(float))
d["ERBB2_log2"] = np.log2(d.ERBB2_expr + 1)
d["ihc"] = pd.to_numeric(d.pid.map(clin.HER2_IHC_SCORE), errors="coerce")
d["er_flag"] = d.index.isin(FLAG)

subsets = {
    "all TNBC (reference)": np.ones(len(d), bool),
    "all TNBC, ER-flag removed": ~d.er_flag.to_numpy(),
    "basal, ER-flag removed": (d.pam50.eq("Basal") & ~d.er_flag).to_numpy(),
    "HER2-low (IHC 1+ or 2+)": d.ihc.isin([1, 2]).to_numpy(),
    "HER2-low, basal": (d.ihc.isin([1, 2]) & d.pam50.eq("Basal")).to_numpy(),
}
rows = []
for name, m in subsets.items():
    g, nul = d[m], nulls[m]
    for ab, lab in PROTEINS.items():
        r = stats.spearmanr(g.score_no_amp, g[ab])[0]
        null = np.array([stats.spearmanr(nul[i], g[ab])[0] for i in nul.columns])
        rows.append(dict(subset=name, n=len(g), protein=lab, rho_score=r,
                         p_emp=(1 + np.sum(np.abs(null) >= abs(r))) / (1 + len(null)),
                         rho_ERBB2_mRNA=stats.spearmanr(g.ERBB2_log2, g[ab])[0],
                         partial_given_mRNA=partial_spearman(g.score_no_amp, g[ab], [g.ERBB2_log2])[0]))
res = pd.DataFrame(rows)


def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out


res["q_BH"] = res.groupby("subset").p_emp.transform(bh)
res.to_csv("results/protein_robustness.csv", index=False)
pd.set_option("display.width", 200)
print(f"ER-flag tumours with RPPA: {d.er_flag.sum()} of {len(FLAG)}; HER2-low with RPPA: {d.ihc.isin([1, 2]).sum()}")
print(res.round(3).to_string(index=False))

# HER2 protein in the top score third, with and without ER-flag tumours
for name, g in [("with ER-flag", d), ("ER-flag removed", d[~d.er_flag])]:
    t = pd.qcut(g.score_no_amp, 3, labels=["low", "mid", "high"])
    print(f"HER2 protein median by score third ({name}): "
          + ", ".join(f"{k} {g.loc[t.eq(k), 'ERBB2|HER2'].median():.2f}" for k in ["low", "mid", "high"]))
