"""HER-OMI step 6: within PAM50 basal TNBC, does the amplicon-free score's link to ERBB2 mRNA
survive controlling for proliferation, oestrogen response, ESR1 and LAR?

Covariates (every one of the 144 signature genes removed from every covariate gene set):
  prolif11 : PAM50 11-gene proliferation index (Nielsen 2010), singscore up-only
  prolifH  : Hallmark E2F targets + G2M checkpoint, singscore up-only
  estH     : Hallmark oestrogen response early + late, singscore up-only
  ESR1     : ESR1 mRNA
  LAR      : mean z of AR, FOXA1, SPDEF, XBP1, GATA3 within the subset (same proxy as step 3)
Null: 1000 random gene sets matched for size and direction (33 up, 67 down) from non-signature
genes, each run through the identical model.
ER-low sensitivity: drop the tumours TNBCtype flagged as possibly ER positive (12 TCGA, 11 METABRIC).
"""
import numpy as np
import pandas as pd
from scipy import stats

import sys; sys.path.insert(0, "scripts")
from heromi_lib import load_expression, singscore, partial_spearman

RAW = "data/raw"
rng = np.random.default_rng(607)

ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")
sig = pd.read_csv("data/signature_coverage.csv")
SIG = set(sig.gene)
up_n = int((sig.TCGA & sig.METABRIC & sig.direction.eq("Upregulated") & ~sig.on_her2_amplicon).sum())
dn_n = int((sig.TCGA & sig.METABRIC & sig.direction.eq("Downregulated") & ~sig.on_her2_amplicon).sum())


def gmt(path, names):
    return {f[0]: f[2:] for f in (l.rstrip("\n").split("\t") for l in open(path)) if f[0] in names}


H = gmt("data/ref/h.all.gmt", {"HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
                               "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE"})
SETS = {
    "prolif11": ["BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS", "UBE2C"],
    "prolifH": sorted(set(H["HALLMARK_E2F_TARGETS"]) | set(H["HALLMARK_G2M_CHECKPOINT"])),
    "estH": sorted(set(H["HALLMARK_ESTROGEN_RESPONSE_EARLY"]) | set(H["HALLMARK_ESTROGEN_RESPONSE_LATE"])),
}
SETS = {k: [g for g in v if g not in SIG] for k, v in SETS.items()}
LAR = ["AR", "FOXA1", "SPDEF", "XBP1", "GATA3"]

up_tbl = "results/tnbctype_upload"
FLAG = {c: sorted(set(pd.read_csv(f"{up_tbl}/{c}_TNBC_expression.csv", nrows=0).columns) -
                  set(pd.read_csv(f"{up_tbl}/{c}_TNBC_expression_noERflag.csv", nrows=0).columns))
        for c in ("TCGA", "METABRIC")}

FILES = {"TCGA": f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt",
         "METABRIC": f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt"}


def spearman(a, b):
    return stats.spearmanr(a, b)[0]


partial = partial_spearman


MODELS = {"none": [], "prolif11": ["prolif11"], "prolifH": ["prolifH"], "ESR1": ["ESR1"],
          "estH": ["estH"], "LAR": ["LAR"], "prolifH+estH": ["prolifH", "estH"],
          "prolifH+ESR1": ["prolifH", "ESR1"], "prolifH+estH+ESR1+LAR": ["prolifH", "estH", "ESR1", "LAR"]}

rows, covs_out = [], []
for cohort, path in FILES.items():
    a = ann[ann.cohort.eq(cohort)]
    expr = load_expression(path, set(a.index))
    sym = pd.read_csv(path, sep="\t", usecols=[0, 1]).dropna().drop_duplicates("Hugo_Symbol")
    s2e = dict(zip(sym.Hugo_Symbol, sym.Entrez_Gene_Id.astype(int)))
    ids = lambda genes: [s2e[g] for g in genes if g in s2e and s2e[g] in expr.index]

    cov = pd.DataFrame(index=expr.columns)
    for k, genes in SETS.items():
        cov[k] = singscore(expr, ids(genes), []).total
    cov["ESR1"] = expr.loc[s2e["ESR1"]]
    cov = cov.join(a[["score_no_amp", "ERBB2_expr", "pam50"] + LAR])
    cov["er_flag"] = cov.index.isin(FLAG[cohort])
    cov["cohort"] = cohort
    covs_out.append(cov)

    b = cov[cov.pam50.eq("Basal")].copy()
    print(f"\n=== {cohort} basal n={len(b)}; TNBCtype ER-flagged in basal: {b.er_flag.sum()} of {len(FLAG[cohort])}")
    for c in ["prolif11", "prolifH", "ESR1", "estH"]:
        print(f"  {c:9s} vs score {spearman(b.score_no_amp, b[c]):+.2f}   vs ERBB2 {spearman(b.ERBB2_expr, b[c]):+.2f}")

    sig_rows = set(ids(SIG))
    pool = np.array([i for i, g in enumerate(expr.index) if g not in sig_rows])
    N = expr.shape[0]

    for subset, bb in (("basal", b), ("basal_no_erflag", b[~b.er_flag])):
        bb = bb.copy()
        bb["LAR"] = bb[LAR].apply(lambda x: (x - x.mean()) / x.std()).mean(axis=1)   # z within subset
        if subset == "basal":
            print(f"  LAR       vs score {spearman(bb.score_no_amp, bb.LAR):+.2f}   vs ERBB2 {spearman(bb.ERBB2_expr, bb.LAR):+.2f}")
        R = expr[bb.index].rank(axis=0).values
        norm = lambda mr, k: (mr - (k + 1) / 2) / ((2 * N - k + 1) / 2 - (k + 1) / 2) - 0.5
        nulls = []
        for _ in range(1000):
            pick = rng.choice(pool, up_n + dn_n, replace=False)
            nulls.append(norm(R[pick[:up_n]].mean(0), up_n) + norm((N + 1 - R[pick[up_n:]]).mean(0), dn_n))
        for m, cs in MODELS.items():
            C = [bb[c] for c in cs]
            r, p = partial(bb.score_no_amp, bb.ERBB2_expr, C)
            nr = np.array([partial(ns, bb.ERBB2_expr, C)[0] for ns in nulls])
            rows.append(dict(cohort=cohort, subset=subset, model=m, n=len(bb), k=len(cs), rho=r, p=p,
                             null_95=np.quantile(nr, 0.95), p_emp=(1 + (nr >= r).sum()) / 1001))

out = pd.DataFrame(rows)
se = 1 / np.sqrt(out.n - 3 - out.k)
out["ci_lo"], out["ci_hi"] = np.tanh(np.arctanh(out.rho) - 1.96 * se), np.tanh(np.arctanh(out.rho) + 1.96 * se)
base = out[out.model.eq("none")].set_index(["cohort", "subset"]).rho
out["retained"] = out.rho / [base[(c, s)] for c, s in zip(out.cohort, out.subset)]
out.to_csv("results/confounder_test.csv", index=False)
pd.concat(covs_out).to_csv("results/covariates_tnbc.csv")

t = out[out.subset.eq("basal")].set_index(["model", "cohort"])
print("\nCohort difference (Fisher z), basal:")
for m in MODELS:
    x, y = t.loc[(m, "TCGA")], t.loc[(m, "METABRIC")]
    z = (np.arctanh(y.rho) - np.arctanh(x.rho)) / np.hypot(1 / np.sqrt(x.n - 3 - x.k), 1 / np.sqrt(y.n - 3 - y.k))
    print(f"  {m:24s} p = {2 * stats.norm.sf(abs(z)):.3f}")
pd.set_option("display.width", 200)
print(out.round(3).to_string(index=False))
