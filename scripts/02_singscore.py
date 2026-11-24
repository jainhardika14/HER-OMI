"""HER-OMI step 2: score every TNBC tumour with the IBCGA HER2-like signature (singscore, bidirectional).

singscore (Foroutan et al. 2018, BMC Bioinformatics) in four steps, done separately for each tumour:
  1. Rank all genes in that tumour by expression, lowest = 1, highest = N.
  2. Up genes: take their mean rank. If the signature is "on", up genes sit near the top.
  3. Down genes: rank in reverse (highest = 1) and take the mean. If "on", down genes sit near the bottom.
  4. Rescale each mean rank so its theoretical min to max maps onto -0.5 to +0.5, then add up and down.
Total score runs from -1 to +1; higher = more HER2-like. Because it only uses ranks inside one tumour,
it does not care whether the numbers came from RNA-seq or a microarray, which is why we use it here.
"""
import numpy as np
import pandas as pd

import sys; sys.path.insert(0, "scripts")
exec(open("scripts/01_define_tnbc.py").read())          # gives t (TCGA) and m (METABRIC) with TNBC flags
RAW = "data/raw"


from heromi_lib import load_expression, singscore


# ---- signature: the 109 genes measured in both cohorts ----
sig = pd.read_csv("data/signature_coverage.csv")
sig = sig[sig.TCGA & sig.METABRIC].copy()
# a gene HGNC no longer lists under its old symbol (FAM136A) still has an Entrez ID in the TCGA matrix
sym2ent = pd.read_csv(f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", sep="\t", usecols=[0, 1]).dropna()
sym2ent = sym2ent.drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol").Entrez_Gene_Id
sig["entrez_id"] = sig.entrez_id.fillna(sig.gene.map(sym2ent)).astype(int)

def gene_ids(rows, direction):
    return rows.loc[rows.direction.eq(direction), "entrez_id"].tolist()

versions = {
    "full":      sig,                                   # 109 genes
    "no_amp":    sig[~sig.on_her2_amplicon],            # 100 genes, the key test
    "amp_only":  sig[sig.on_her2_amplicon],             # 9 genes, all upregulated
}

# ---- expression for TNBC tumours only ----
tcga_ids = {c for c in pd.read_csv(f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", sep="\t", nrows=0).columns
            if c[13:15] == "01" and c[:12] in set(t.loc[t.TNBC, "PATIENT_ID"])}
cohorts = {
    "TCGA": load_expression(f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", tcga_ids),
    "METABRIC": load_expression(f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt",
                                set(m.loc[m.TNBC_cn, "PATIENT_ID"])),
}

rows = []
for cohort, expr in cohorts.items():
    print(f"{cohort}: {expr.shape[1]} TNBC tumours, {expr.shape[0]} genes ranked")
    res = pd.DataFrame(index=expr.columns)
    for name, s in versions.items():
        res[f"score_{name}"] = singscore(expr, gene_ids(s, "Upregulated"), gene_ids(s, "Downregulated")).total
    res["ERBB2_expr"] = expr.loc[2064]                    # Entrez 2064 = ERBB2, for sanity checks
    res["cohort"] = cohort
    rows.append(res)

scores = pd.concat(rows).rename_axis("sample").reset_index()
scores.to_csv("results/scores_tnbc.csv", index=False)
print(scores.groupby("cohort")[[c for c in scores if c.startswith("score")]].describe().T.round(3))
