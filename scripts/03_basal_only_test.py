"""HER-OMI step 3: does the HER2-like signature mean anything once LAR is (mostly) out of the picture?

Idea: LAR tumours are almost never PAM50 Basal. So restrict to PAM50 Basal TNBC and ask whether the
amplicon-free score still tracks ERBB2 mRNA. Two extra guards:
  * null model: score 1000 random gene sets of the same size and direction mix, to see what correlation
    with ERBB2 an arbitrary 100-gene score reaches by chance
  * chromosome 17 check: broad 17q gains (not focal amplification) also raise ERBB2, so we rescore
    without any chromosome 17 gene
"""
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, "scripts")
from heromi_lib import load_expression, singscore

RAW = "data/raw"
rng = np.random.default_rng(607)
N_RANDOM = 1000

s = pd.read_csv("results/scores_tnbc_annot.csv")
sig = pd.read_csv("data/signature_coverage.csv")
hgnc = pd.read_csv("data/ref/hgnc_complete_set.txt", sep="\t", usecols=["symbol", "location"], low_memory=False)
sig = sig[sig.TCGA & sig.METABRIC].merge(hgnc.rename(columns={"symbol": "gene"}), on="gene", how="left")
sym2ent = pd.read_csv(f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", sep="\t", usecols=[0, 1]).dropna()
sym2ent = sym2ent.drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol").Entrez_Gene_Id
sig["entrez_id"] = sig.entrez_id.fillna(sig.gene.map(sym2ent)).astype(int)
sig["chr17"] = sig.location.fillna("").str.startswith("17")
noamp = sig[~sig.on_her2_amplicon]
no17 = sig[~sig.chr17]
print(f"Signature genes on chromosome 17: {sig.chr17.sum()} "
      f"(amplicon {sig.on_her2_amplicon.sum()}, other chr17 {(sig.chr17 & ~sig.on_her2_amplicon).sum()})")
print("  non-amplicon chr17 genes:", noamp[noamp.chr17].gene.tolist())


def ids(rows, d):
    return rows.loc[rows.direction.eq(d), "entrez_id"].tolist()


def spearman(a, b):
    return pd.Series(a).corr(pd.Series(b), method="spearman")


def partial_spearman(x, y, z):
    r = pd.DataFrame({"x": x, "y": y, "z": z}).rank()
    X = np.c_[np.ones(len(r)), r.z]
    res = lambda v: v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
    return np.corrcoef(res(r.x.values), res(r.y.values))[0, 1]


paths = {"TCGA": f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt",
         "METABRIC": f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt"}
out = []
for cohort, path in paths.items():
    for subset in ["all TNBC", "Basal only"]:
        g = s[s.cohort.eq(cohort)].copy()
        if subset == "Basal only":
            g = g[g.pam50.eq("Basal")]
        expr = load_expression(path, set(g["sample"]))
        g = g.set_index("sample").loc[expr.columns]
        erbb2 = expr.loc[2064]
        lar = g[["AR", "FOXA1", "SPDEF", "XBP1", "GATA3"]].apply(lambda x: (x - x.mean()) / x.std()).mean(axis=1)

        obs = singscore(expr, ids(noamp, "Upregulated"), ids(noamp, "Downregulated")).total
        no17_score = singscore(expr, ids(no17, "Upregulated"), ids(no17, "Downregulated")).total

        # null: random gene sets with the same number of up and down genes, never containing ERBB2
        pool = expr.index.drop(2064).to_numpy()
        n_up, n_dn = len(ids(noamp, "Upregulated")), len(ids(noamp, "Downregulated"))
        # same maths as singscore(), but ranking once instead of 1000 times (identical result, much faster)
        R = expr.rank(axis=0).to_numpy(); N = R.shape[0]; row = {gid: i for i, gid in enumerate(expr.index)}
        def fast(up, dn):
            u = R[[row[x] for x in up]].mean(0); d = (N + 1 - R[[row[x] for x in dn]]).mean(0)
            sc = lambda m, n: (m - (n + 1) / 2) / ((2 * N - n + 1) / 2 - (n + 1) / 2) - 0.5
            return sc(u, len(up)) + sc(d, len(dn))
        assert np.allclose(fast(ids(noamp, "Upregulated"), ids(noamp, "Downregulated")), obs.values)
        null = []
        for _ in range(N_RANDOM):
            pick = rng.choice(pool, n_up + n_dn, replace=False)
            null.append(spearman(fast(pick[:n_up], pick[n_up:]), erbb2.values))
        null = np.array(null)
        rho = spearman(obs, erbb2)

        out.append({
            "cohort": cohort, "subset": subset, "n": expr.shape[1],
            "score_SD": obs.std(),
            "rho_ERBB2": rho,
            "null_95th": np.quantile(null, 0.95),
            "empirical_p": (np.sum(np.abs(null) >= abs(rho)) + 1) / (N_RANDOM + 1),
            "rho_ERBB2_no_chr17": spearman(no17_score, erbb2),
            "rho_AR": spearman(obs, g.AR),
            "rho_LARavg": spearman(obs, lar),
            "partial_rho_ERBB2_given_LAR": partial_spearman(obs, erbb2, lar),
        })
        print(f"done {cohort} {subset}")

res = pd.DataFrame(out)
res.to_csv("results/basal_only_test.csv", index=False)
pd.set_option("display.width", 200)
print(res.round(3).to_string(index=False))
