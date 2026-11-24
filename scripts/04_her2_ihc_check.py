"""HER-OMI step 4: does the amplicon-free score rise with HER2 IHC score (0, 1+, 2+) in TCGA TNBC?
METABRIC has no HER2 IHC, so this is TCGA only. 2+ tumours here are FISH negative by our TNBC definition."""
import pandas as pd
from scipy import stats

s = pd.read_csv("results/scores_tnbc_annot.csv")
clin = pd.read_csv("data/raw/brca_tcga/data_clinical_patient.txt", sep="\t", comment="#")[
    ["PATIENT_ID", "HER2_IHC_SCORE", "IHC_HER2", "HER2_FISH_STATUS"]]
t = s[s.cohort.eq("TCGA")].merge(clin, left_on="pid", right_on="PATIENT_ID")
t["ihc"] = pd.to_numeric(t.HER2_IHC_SCORE, errors="coerce")
print("IHC 3+ tumours inside TNBC (contradictory records, excluded):",
      t.loc[t.ihc.eq(3), ["pid", "IHC_HER2", "HER2_FISH_STATUS"]].to_dict("records"))
t = t[t.ihc.isin([0, 1, 2])]

rows = []
for subset, g in [("all TNBC", t), ("Basal only", t[t.pam50.eq("Basal")])]:
    for k, gg in g.groupby("ihc"):
        rows.append({"subset": subset, "IHC": f"{int(k)}+", "n": len(gg),
                     "median_score": gg.score_no_amp.median(), "median_ERBB2_log2": (gg.ERBB2_expr + 1).pipe(lambda x: x.apply(lambda v: __import__('math').log2(v))).median()})
    rho, p = stats.spearmanr(g.ihc, g.score_no_amp)
    zero, low = g[g.ihc.eq(0)].score_no_amp, g[g.ihc.isin([1, 2])].score_no_amp
    mw = stats.mannwhitneyu(low, zero, alternative="two-sided")
    print(f"{subset}: n={len(g)} | trend Spearman rho={rho:.2f} p={p:.3f} | "
          f"HER2-low vs HER2-0 Mann-Whitney p={mw.pvalue:.3f} (n {len(low)} vs {len(zero)})")
print(pd.DataFrame(rows).round(3).to_string(index=False))
