"""HER-OMI step 7: is the HER2-like signal just the LAR subtype?

Input: TNBCtype (Lehmann 2011 centroids; Chen et al. 2012 web tool) results for the cleaned uploads,
TCGA 147 and METABRIC 239 tumours, saved in data/tnbctype/ as <COHORT>_result.csv (call per tumour)
and <COHORT>_coef.csv (correlation with each of the six centroids). The 12 TCGA and 11 METABRIC
tumours TNBCtype flagged as possibly ER positive have no call; they are labelled "ER-flag" here.

Questions
  1. Do our subtype frequencies resemble those the preprint reports for the same cohorts (S15)?
  2. Where does each subtype sit on the amplicon-free score? How well does the score pick out LAR?
  3. Does the score's link to ERBB2 mRNA survive (a) removing LAR tumours, (b) controlling for each
     tumour's continuous LAR centroid correlation, in all TNBC and within PAM50 basal?
     Each against 1000 random gene sets put through the identical test.
"""
import numpy as np
import pandas as pd
from scipy import stats

import sys; sys.path.insert(0, "scripts")
from heromi_lib import load_expression, partial_spearman, random_set_scores

RAW = "data/raw"
ORDER = ["BL1", "BL2", "IM", "M", "MSL", "LAR", "UNS", "ER-flag"]
FILES = {"TCGA": f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt",
         "METABRIC": f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt"}
S15 = {"TCGA": {"BL1": 30, "BL2": 21, "IM": 27, "LAR": 15, "M": 35, "MSL": 13, "UNS": 16},
       "METABRIC": {"BL1": 71, "BL2": 21, "IM": 66, "LAR": 39, "M": 48, "MSL": 38, "UNS": 37}}

ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")
sig = pd.read_csv("data/signature_coverage.csv")
both = sig.TCGA & sig.METABRIC & ~sig.on_her2_amplicon
n_up, n_dn = int((both & sig.direction.eq("Upregulated")).sum()), int((both & sig.direction.eq("Downregulated")).sum())


def fix_id(s):
    return s.replace(".", "-")


calls = []
for c in FILES:
    r = pd.read_csv(f"data/tnbctype/{c}_result.csv", index_col=0)
    k = pd.read_csv(f"data/tnbctype/{c}_coef.csv", index_col=0).add_prefix("coef_")
    d = r.join(k)
    d.index = d.index.map(fix_id)
    calls.append(d)
calls = pd.concat(calls)
t = ann.join(calls[["subtype", "correlation", "coef_LAR"]])
t["subtype"] = t.subtype.fillna("ER-flag")
t.to_csv("results/scores_tnbc_tnbctype.csv")

pd.set_option("display.width", 200)
summary = []
for c in FILES:
    a = t[t.cohort.eq(c)].copy()
    called = a[~a.subtype.eq("ER-flag")]
    print(f"\n==================== {c}: {len(a)} TNBC, {len(called)} with a TNBCtype call")

    # 1. frequencies vs preprint S15
    ours = called.subtype.value_counts().reindex(ORDER[:-1]).fillna(0).astype(int)
    theirs = pd.Series(S15[c]).reindex(ORDER[:-1])
    print("1) Subtype frequencies, ours vs preprint S15 (different TNBC definitions, so descriptive only)")
    print(pd.DataFrame({"ours_n": ours, "ours_%": (100 * ours / ours.sum()).round(1),
                        "S15_n": theirs, "S15_%": (100 * theirs / theirs.sum()).round(1)}).T.to_string())
    chi = stats.chi2_contingency(np.vstack([ours.values, theirs.values]))
    print(f"   chi-square p = {chi[1]:.3f}")

    # 2. score by subtype; tertiles over all TNBC in the cohort
    a["tertile"] = pd.qcut(a.score_no_amp, 3, labels=["low", "mid", "high"])
    print("\n2) Amplicon-free score by subtype (median) and tertile counts")
    tab = pd.crosstab(a.subtype, a.tertile).reindex(ORDER).fillna(0).astype(int)
    tab["median_score"] = a.groupby("subtype").score_no_amp.median().reindex(ORDER).round(3)
    tab["basal_PAM50"] = a[a.pam50.eq("Basal")].subtype.value_counts().reindex(ORDER).fillna(0).astype(int)
    print(tab.to_string())
    hi = a[a.tertile.eq("high")]
    lar = called.subtype.eq("LAR")
    auc = stats.mannwhitneyu(called.score_no_amp[lar], called.score_no_amp[~lar]).statistic / (lar.sum() * (~lar).sum())
    print(f"   top tertile: {len(hi)} tumours; LAR {hi.subtype.eq('LAR').sum()}, ER-flag {hi.subtype.eq('ER-flag').sum()}, "
          f"other {(~hi.subtype.isin(['LAR', 'ER-flag'])).sum()}")
    print(f"   LAR tumours in top tertile: {a[a.subtype.eq('LAR')].tertile.eq('high').sum()} of {a.subtype.eq('LAR').sum()}")
    print(f"   AUC of score for LAR vs other called subtypes: {auc:.2f}")
    print(f"   score vs LAR centroid correlation (all called): rho = {stats.spearmanr(called.score_no_amp, called.coef_LAR)[0]:.2f}")
    print("   PAM50 x TNBCtype:")
    print(pd.crosstab(a.pam50.fillna("NA"), a.subtype).reindex(columns=ORDER).fillna(0).astype(int).to_string())

    # 3. ERBB2 link beyond LAR, against matched random gene sets
    expr = load_expression(FILES[c], set(a.index))
    sig_ids = set(sig.loc[sig.entrez_id.notna(), "entrez_id"].astype(int))
    nulls = pd.DataFrame(random_set_scores(expr, sig_ids, n_up, n_dn).T, index=expr.columns)

    tests = {
        "all TNBC":                          (a.index, []),
        "all TNBC, no LAR, no ER-flag":      (a.index[~a.subtype.isin(["LAR", "ER-flag"])], []),
        "all called, given LAR coef":        (called.index, ["coef_LAR"]),
        "basal":                             (a.index[a.pam50.eq("Basal")], []),
        "basal, no LAR, no ER-flag":         (a.index[a.pam50.eq("Basal") & ~a.subtype.isin(["LAR", "ER-flag"])], []),
        "basal called, given LAR coef":      (a.index[a.pam50.eq("Basal") & ~a.subtype.eq("ER-flag")], ["coef_LAR"]),
        "TNBCtype BL1+BL2":                  (a.index[a.subtype.isin(["BL1", "BL2"])], []),
    }
    print("\n3) Amplicon-free score vs ERBB2 mRNA")
    for name, (idx, cov) in tests.items():
        g = a.loc[idx]
        C = [g[x] for x in cov]
        r, p = partial_spearman(g.score_no_amp, g.ERBB2_expr, C)
        nr = np.array([partial_spearman(nulls.loc[idx, i], g.ERBB2_expr, C)[0] for i in range(nulls.shape[1])])
        p_emp = (1 + (nr >= r).sum()) / (1 + len(nr))
        print(f"   {name:34s} n={len(g):3d}  rho={r:+.2f}  p={p:.3g}  null95={np.quantile(nr, .95):.2f}  p_emp={p_emp:.3f}")
        summary.append(dict(cohort=c, test=name, n=len(g), rho=r, p=p, null_95=np.quantile(nr, .95), p_emp=p_emp))

pd.DataFrame(summary).to_csv("results/tnbctype_tests.csv", index=False)
