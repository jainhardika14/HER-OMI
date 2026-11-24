"""HER-OMI step 14: is the HER2-like score just normal tissue mixed into the tumour sample?

Step 13 found that normal breast scores as high as HER2-positive tumours on the HER2-like score. So a tumour
sample with more normal tissue in it could score higher for that reason alone. TCGA's clinical tables here
have no purity estimate, so we use an expression proxy: adjacent normal breast is mostly fat, and adipocyte
genes are almost absent from tumour cells. Proxy = mean z score of 8 adipocyte genes (ADIPOQ, PLIN1, LEP,
CIDEA, FABP4, LPL, GPD1, CD36), computed within each cohort's TNBC.
Questions: does the score track the proxy? Do the score's links to ERBB2 mRNA (step 3) and to AR and INPP4B
protein (step 8) survive controlling for it?
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, "scripts")
from heromi_lib import load_expression, partial_spearman

RAW = "data/raw"
ADIPO = {"ADIPOQ": 9370, "PLIN1": 5346, "LEP": 3952, "CIDEA": 1149, "FABP4": 2167, "LPL": 4023, "GPD1": 2819, "CD36": 948}
ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")
FILES = {"TCGA": (f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", True),
         "METABRIC": (f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt", False)}
rows = []
for c, (path, log) in FILES.items():
    a = ann[ann.cohort.eq(c)]
    e = load_expression(path, set(a.index))
    e = np.log2(e + 1) if log else e
    g = e.loc[[i for i in ADIPO.values() if i in e.index]]
    proxy = g.sub(g.mean(axis=1), axis=0).div(g.std(axis=1), axis=0).mean(axis=0)
    d = a.loc[e.columns].assign(adipose=proxy)
    for subset, m in [("all TNBC", np.ones(len(d), bool)), ("basal", d.pam50.eq("Basal").to_numpy())]:
        x = d[m]
        rows.append(dict(cohort=c, subset=subset, n=len(x), test="score vs adipose proxy",
                         rho=stats.spearmanr(x.score_no_amp, x.adipose)[0], p=stats.spearmanr(x.score_no_amp, x.adipose)[1]))
        r0 = stats.spearmanr(x.score_no_amp, x.ERBB2_expr)[0]
        r1, p1 = partial_spearman(x.score_no_amp, x.ERBB2_expr, [x.adipose])
        rows.append(dict(cohort=c, subset=subset, n=len(x), test="score vs ERBB2 mRNA", rho=r0, p=np.nan))
        rows.append(dict(cohort=c, subset=subset, n=len(x), test="score vs ERBB2 mRNA, given adipose proxy", rho=r1, p=p1))
    if c == "TCGA":
        rppa = pd.read_csv(f"{RAW}/brca_tcga/data_rppa.txt", sep="\t").set_index("Composite.Element.REF").T
        ids = [s for s in d.index if s in rppa.index]
        x = d.loc[ids].join(rppa.loc[ids, ["AR|AR", "INPP4B|INPP4B", "ERBB2|HER2"]].astype(float))
        for ab, lab in [("AR|AR", "AR protein"), ("INPP4B|INPP4B", "INPP4B protein"), ("ERBB2|HER2", "HER2 protein")]:
            rows.append(dict(cohort=c, subset="all TNBC", n=len(x), test=f"score vs {lab}", rho=stats.spearmanr(x.score_no_amp, x[ab])[0], p=np.nan))
            r1, p1 = partial_spearman(x.score_no_amp, x[ab], [x.adipose])
            rows.append(dict(cohort=c, subset="all TNBC", n=len(x), test=f"score vs {lab}, given adipose proxy", rho=r1, p=p1))
out = pd.DataFrame(rows)
out.to_csv("results/normal_admixture_test.csv", index=False)
pd.set_option("display.width", 200)
print(out.round(3).to_string(index=False))
