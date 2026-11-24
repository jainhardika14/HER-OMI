"""HER-OMI step 5: does the HER2-like gradient carry the same pathway biology IBCGA reported?

IBCGA tables (supplementary workbook):
  S22  HER2-like vs basal-like, 1,741 MSigDB Hallmark + Reactome sets. Sign: positive = higher in BASAL
       (inferred: E2F targets and G2M positive, ERBB2 signalling and oestrogen response negative).
  S23  LAR-like vs basal-like TNBC, 50 Hallmark sets. Sign: positive = higher in LAR
       (fatty acid metabolism and androgen response positive, cell cycle negative).
  Both are re-signed below so that positive always means "higher in the non-basal group".

Our side: score every pathway in every TNBC tumour (singscore, up-only), then correlate each pathway with
the amplicon-free HER2-like score. The 144 signature genes are removed from every pathway first, so a
pathway cannot agree with the score just because it shares genes with it.
"""
import numpy as np
import pandas as pd

MIN_SIZE, MAX_SIZE = 10, 500

s22 = pd.read_csv("data/ibcga_S22.csv").set_index("set")
s23 = pd.read_csv("data/ibcga_S23.csv").set_index("set")
ibcga = pd.DataFrame({"HER2like_vs_basal": -s22["diff"], "LAR_vs_basal": s23["diff"],
                      "padj_S22": s22.padj, "padj_S23": s23.padj})

sig_genes = set(pd.read_csv("data/IBCGA_HER2like_signature_144genes.csv").gene)
sets = {}
for f in ["data/ref/h.all.gmt", "data/ref/c2.cp.reactome.gmt"]:
    for line in open(f):
        name, _, *genes = line.rstrip("\n").split("\t")
        sets[name] = [g for g in genes if g not in sig_genes]

scores = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")


def pathway_scores(expr):
    """Up-only singscore for every gene set: mean within-tumour rank, rescaled to -0.5..0.5."""
    R = expr.rank(axis=0).to_numpy()
    N = R.shape[0]
    row = {g: i for i, g in enumerate(expr.index)}
    out = {}
    for name, genes in sets.items():
        idx = [row[g] for g in genes if g in row]
        n = len(idx)
        if MIN_SIZE <= n <= MAX_SIZE:
            lo, hi = (n + 1) / 2, (2 * N - n + 1) / 2
            out[name] = (R[idx].mean(0) - lo) / (hi - lo) - 0.5
    return pd.DataFrame(out, index=expr.columns)


profiles = {}
for cohort in ["TCGA", "METABRIC"]:
    expr = pd.read_csv(f"results/tnbctype_upload/{cohort}_TNBC_expression.csv", index_col=0)
    ps = pathway_scores(expr)
    meta = scores.loc[ps.index]
    for subset, mask in [("all", np.ones(len(meta), bool)), ("basal", meta.pam50.eq("Basal").to_numpy())]:
        sub = ps[mask]
        profiles[f"{cohort}_{subset}"] = sub.corrwith(meta.loc[mask, "score_no_amp"], method="spearman")
    print(f"{cohort}: {ps.shape[1]} pathways scored in {ps.shape[0]} tumours")

prof = pd.DataFrame(profiles)
prof.to_csv("results/pathway_gradient_profiles.csv")
both = prof.join(ibcga, how="inner")
print(f"\nPathways in both our data and S22: {both.HER2like_vs_basal.notna().sum()}")


def rho(a, b):
    ok = a.notna() & b.notna()
    return a[ok].corr(b[ok], method="spearman"), int(ok.sum())


def partial(x, y, z):
    d = pd.DataFrame({"x": x, "y": y, "z": z}).dropna().rank()
    X = np.c_[np.ones(len(d)), d.z]
    res = lambda v: v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
    return np.corrcoef(res(d.x.to_numpy()), res(d.y.to_numpy()))[0, 1]


print("\n1) Replication between our cohorts (TCGA vs METABRIC profile):")
for subset in ["all", "basal"]:
    r, n = rho(prof[f"TCGA_{subset}"], prof[f"METABRIC_{subset}"])
    print(f"   {subset:5s}  rho = {r:.2f}  over {n} pathways")

print("\n2) Agreement with IBCGA HER2-like vs basal (S22), all pathways / only S22 FDR < 0.05:")
sigS22 = both.padj_S22 < 0.05
for col in prof.columns:
    r_all, n_all = rho(both[col], both.HER2like_vs_basal)
    r_sig, n_sig = rho(both.loc[sigS22, col], both.loc[sigS22, "HER2like_vs_basal"])
    agree = (np.sign(both.loc[sigS22, col]) == np.sign(both.loc[sigS22, "HER2like_vs_basal"])).mean()
    print(f"   {col:15s} rho = {r_all:.2f} (n={n_all}) | FDR<0.05: rho = {r_sig:.2f}, sign agreement {agree:.0%} (n={n_sig})")

h = both.dropna(subset=["LAR_vs_basal", "HER2like_vs_basal"])
r_ref, _ = rho(h.HER2like_vs_basal, h.LAR_vs_basal)
print(f"\n3) The 50 Hallmark sets. Warning: IBCGA's own HER2-like and LAR profiles correlate {r_ref:.2f} here,")
print("   so this comparison can only pick up what is left over after their shared biology.")
for col in prof.columns:
    r_h, _ = rho(h[col], h.HER2like_vs_basal)
    r_l, _ = rho(h[col], h.LAR_vs_basal)
    print(f"   {col:15s} vs HER2-like {r_h:.2f} | vs LAR {r_l:.2f} | "
          f"HER2-like given LAR {partial(h[col], h.HER2like_vs_basal, h.LAR_vs_basal):.2f} | "
          f"LAR given HER2-like {partial(h[col], h.LAR_vs_basal, h.HER2like_vs_basal):.2f}")

key = ["HALLMARK_E2F_TARGETS", "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_ESTROGEN_RESPONSE_EARLY",
       "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "REACTOME_SIGNALING_BY_ERBB2"]
print("\n4) Landmark pathways (correlation with the HER2-like score; IBCGA values re-signed):")
print(both.loc[[k for k in key if k in both.index], list(prof.columns) + ["HER2like_vs_basal", "LAR_vs_basal"]]
      .round(3).to_string())
