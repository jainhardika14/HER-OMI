"""HER-OMI step 8: is the IBCGA "HER2-like" state HER2 at the protein level?

IBCGA (Research Square 2026, Fig. 7) call HER2-like TNBC "transcriptionally HER2-active" and "a candidate
stratum for HER2-directed ADC evaluation", but never measure HER2 protein. ADCs bind HER2 protein, and
HER2 signalling needs phosphorylated HER2. So here, within TNBC:

  A. TCGA RPPA (131 of our TNBC): amplicon-free score vs total HER2, HER2 pY1248, other ERBB family and
     downstream phosphoproteins, plus luminal/AR proteins as positive controls. Each against 1000 random
     gene sets matched for size and direction (two-sided empirical p). Then: does the score add anything
     about HER2 protein beyond ERBB2 mRNA (partial correlation)? And where do high-scoring TNBC sit on the
     HER2 protein scale relative to HER2-positive and HR+ HER2-negative tumours?
  B. CPTAC 2020 (28 TNBC, mass-spec proteome + phosphoproteome): the same question, independent platform.

Uses the 100 amplicon-free signature genes measured in TCGA and METABRIC (33 up, 67 down).
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "scripts")
from heromi_lib import load_expression, singscore, partial_spearman, random_set_scores

RAW = "data/raw"
rng = np.random.default_rng(607)

sig = pd.read_csv("data/signature_coverage.csv")
core = sig[sig.TCGA & sig.METABRIC & ~sig.on_her2_amplicon]
n_up, n_dn = int(core.direction.eq("Upregulated").sum()), int(core.direction.eq("Downregulated").sum())
ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")

PROTEINS = {  # RPPA name -> short label
    "ERBB2|HER2": "HER2", "ERBB2|HER2_pY1248": "HER2 pY1248",
    "EGFR|EGFR": "EGFR", "EGFR|EGFR_pY1068": "EGFR pY1068", "ERBB3|HER3": "HER3", "ERBB3|HER3_pY1289": "HER3 pY1289",
    "MAP2K1|MEK1_pS217_S221": "MEK1 pS217/221", "MAPK1 MAPK3|MAPK_pT202_Y204": "ERK pT202/Y204",
    "AKT1 AKT2 AKT3|Akt_pS473": "AKT pS473", "RPS6KB1|p70S6K_pT389": "p70S6K pT389", "SRC|Src_pY416": "SRC pY416",
    "AR|AR": "AR (control)", "GATA3|GATA3": "GATA3 (control)", "INPP4B|INPP4B": "INPP4B (control)",
    "ESR1|ER-alpha": "ER-alpha (control)",
}


def emp_p(r, null):
    return (1 + np.sum(np.abs(null) >= abs(r))) / (1 + len(null))


# ---------------- A. TCGA RPPA ----------------
rppa = pd.read_csv(f"{RAW}/brca_tcga/data_rppa.txt", sep="\t").set_index("Composite.Element.REF").T
tc = ann[ann.cohort.eq("TCGA")]
ids = [s for s in tc.index if s in rppa.index]
expr = load_expression(f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", set(ids))
ids = list(expr.columns)
sig_ids = set(sig.loc[sig.entrez_id.notna(), "entrez_id"].astype(int))
null_scores = pd.DataFrame(random_set_scores(expr, sig_ids, n_up, n_dn).T, index=ids)
d = tc.loc[ids].join(rppa.loc[ids, list(PROTEINS)].astype(float))
d["ERBB2_log2"] = np.log2(d.ERBB2_expr + 1)
print(f"TCGA TNBC with RPPA: {len(d)} (PAM50 basal {d.pam50.eq('Basal').sum()})")

rows = []
for subset, m in [("all TNBC", np.ones(len(d), bool)), ("basal", d.pam50.eq("Basal").to_numpy())]:
    g, nul = d[m], null_scores[m]
    for ab, lab in PROTEINS.items():
        r = stats.spearmanr(g.score_no_amp, g[ab])[0]
        null = np.array([stats.spearmanr(nul[i], g[ab])[0] for i in nul.columns])
        r_mrna = stats.spearmanr(g.ERBB2_log2, g[ab])[0]
        pr, _ = partial_spearman(g.score_no_amp, g[ab], [g.ERBB2_log2])
        rows.append(dict(cohort="TCGA", subset=subset, n=len(g), protein=lab, rho_score=r,
                         null_lo=np.quantile(null, .025), null_hi=np.quantile(null, .975), p_emp=emp_p(r, null),
                         rho_ERBB2_mRNA=r_mrna, partial_score_given_mRNA=pr))
res = pd.DataFrame(rows)

# Where do high scorers sit on the HER2 protein scale? Reference groups from all TCGA tumours with RPPA.
clin = pd.read_csv(f"{RAW}/brca_tcga/data_clinical_patient.txt", sep="\t", comment="#").set_index("PATIENT_ID")
ref = rppa[["ERBB2|HER2"]].astype(float).copy()
ref["pid"] = ref.index.str[:12]
ref = ref.join(clin[["ER_STATUS_BY_IHC", "IHC_HER2", "HER2_FISH_STATUS"]], on="pid")
her2pos = ref.IHC_HER2.eq("Positive") | ref.HER2_FISH_STATUS.eq("Positive")
hrpos_her2neg = ref.ER_STATUS_BY_IHC.eq("Positive") & ~her2pos & (ref.IHC_HER2.eq("Negative") | ref.HER2_FISH_STATUS.eq("Negative"))
d["tertile"] = pd.qcut(d.score_no_amp, 3, labels=["low", "mid", "high"])
scale = {
    "HER2-positive (IHC 3+ or FISH+)": ref.loc[her2pos, "ERBB2|HER2"],
    "HR+ HER2-negative": ref.loc[hrpos_her2neg, "ERBB2|HER2"],
    "TNBC, top score third": d.loc[d.tertile.eq("high"), "ERBB2|HER2"],
    "TNBC, middle third": d.loc[d.tertile.eq("mid"), "ERBB2|HER2"],
    "TNBC, bottom third": d.loc[d.tertile.eq("low"), "ERBB2|HER2"],
}
scale_tbl = pd.DataFrame({k: dict(n=len(v), median=v.median(), q25=v.quantile(.25), q75=v.quantile(.75))
                          for k, v in scale.items()}).T
hi, lo, pos = scale["TNBC, top score third"], scale["TNBC, bottom third"], scale["HER2-positive (IHC 3+ or FISH+)"]
gap_frac = (hi.median() - lo.median()) / (pos.median() - lo.median())
mw = stats.mannwhitneyu(hi, lo)

# ---------------- B. CPTAC 2020 ----------------
cp = pd.read_csv(f"{RAW}/brca_cptac_2020/data_clinical_patient.txt", sep="\t", comment="#")
cs = pd.read_csv(f"{RAW}/brca_cptac_2020/data_clinical_sample.txt", sep="\t", comment="#")
cc = cp.merge(cs, on="PATIENT_ID").set_index("SAMPLE_ID")
tn = cc.index[cc.TNBC_UPDATED_CLINICAL_STATUS.str.lower().eq("positive")]
fp = pd.read_csv(f"{RAW}/brca_cptac_2020/data_mrna_seq_fpkm.txt", sep="\t").dropna(subset=["Hugo_Symbol"])
fp = fp.groupby("Hugo_Symbol")[list(tn)].mean().dropna()            # genes with any missing value dropped
up = [g for g in core.loc[core.direction.eq("Upregulated"), "gene"] if g in fp.index]
dn = [g for g in core.loc[core.direction.eq("Downregulated"), "gene"] if g in fp.index]
cscore = singscore(fp, up, dn).total
pool = [g for g in fp.index if g not in set(sig.gene)]
R = fp.rank(axis=0).values
N = R.shape[0]
row = {g: i for i, g in enumerate(fp.index)}
norm = lambda mr, k: (mr - (k + 1) / 2) / ((2 * N - k + 1) / 2 - (k + 1) / 2) - 0.5
cnull = []
for _ in range(1000):
    pick = rng.choice(pool, len(up) + len(dn), replace=False)
    iu, idn = [row[x] for x in pick[:len(up)]], [row[x] for x in pick[len(up):]]
    cnull.append(norm(R[iu].mean(0), len(up)) + norm((N + 1 - R[idn]).mean(0), len(dn)))
cnull = np.array(cnull)

prot = pd.read_csv(f"{RAW}/brca_cptac_2020/data_protein_quantification.txt", sep="\t").set_index("Composite.Element.REF")
phos = pd.read_csv(f"{RAW}/brca_cptac_2020/data_phosphoprotein_quantification.txt", sep="\t").set_index("ENTITY_STABLE_ID")
pick_prot = {"ERBB2": "HER2", "EGFR": "EGFR", "ERBB3": "HER3", "AR": "AR (control)", "GATA3": "GATA3 (control)",
             "INPP4B": "INPP4B (control)"}
cvals = {}
for gene, lab in pick_prot.items():
    hits = [i for i in prot.index if i.split("|")[0] == gene]
    if hits:
        cvals[lab] = prot.loc[hits[0], tn].astype(float)
for site, lab in [("ERBB2_Y1248y", "HER2 pY1248"), ("ERBB2_Y1196y", "HER2 pY1196"), ("ERBB2_S1054s", "HER2 pS1054")]:
    if site in phos.index:
        cvals[lab] = pd.to_numeric(phos.loc[site, tn], errors="coerce")
erbb2_rna = fp.loc["ERBB2", tn] if "ERBB2" in fp.index else None
crows = []
for lab, v in cvals.items():
    ok = v.notna().to_numpy()
    r = stats.spearmanr(cscore[ok], v[ok])[0]
    null = np.array([stats.spearmanr(cnull[i][ok], v[ok])[0] for i in range(len(cnull))])
    crows.append(dict(cohort="CPTAC", subset="TNBC", n=int(ok.sum()), protein=lab, rho_score=r,
                      null_lo=np.quantile(null, .025), null_hi=np.quantile(null, .975), p_emp=emp_p(r, null),
                      rho_ERBB2_mRNA=stats.spearmanr(erbb2_rna[ok], v[ok])[0],
                      partial_score_given_mRNA=partial_spearman(cscore[ok], v[ok], [erbb2_rna[ok]])[0]))
res = pd.concat([res, pd.DataFrame(crows)], ignore_index=True)

# ---------------- output ----------------
res.to_csv("results/protein_test.csv", index=False)
scale_tbl.to_csv("results/her2_protein_scale.csv")
pd.set_option("display.width", 220)
print("\nScore (amplicon-free) vs protein. null = 95% range of 1000 random gene sets. "
      "partial = score vs protein given ERBB2 mRNA.")
print(res.round(3).to_string(index=False))
print(f"\nCPTAC TNBC: {len(tn)} tumours, signature genes found {len(up)} up / {len(dn)} down; "
      f"score vs ERBB2 mRNA rho = {stats.spearmanr(cscore, erbb2_rna)[0]:.2f}")
print("\nTCGA HER2 protein (RPPA) scale:")
print(scale_tbl.round(3).to_string())
print(f"Top vs bottom score third, HER2 protein: Mann-Whitney p = {mw.pvalue:.3f}; the top third closes "
      f"{100 * gap_frac:.1f}% of the gap from the bottom third to HER2-positive tumours.")
