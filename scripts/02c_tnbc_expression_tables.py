"""HER-OMI step 2c: write one expression table per cohort, TNBC tumours only, for TNBCtype and step 5.

Format TNBCtype accepts: genes as rows (gene symbols), tumours as columns, log2 scale.
  TCGA     : log2(RSEM + 1)
  METABRIC : array values as given (already log2)
Duplicated symbols keep the row with the highest mean expression.

Two versions per cohort:
  <COHORT>_TNBC_expression.csv           all our TNBC with RNA (TCGA 159, METABRIC 250)
  <COHORT>_TNBC_expression_noERflag.csv  without the tumours TNBCtype's quality check flagged as possibly
                                         ER positive (data/tnbctype/er_flagged_samples.txt); this is the
                                         version TNBCtype accepted (TCGA 147, METABRIC 239)
Later steps find the ER-flagged tumours as the difference between the two files.
"""
import os
import numpy as np
import pandas as pd

RAW = "data/raw"
OUT = "results/tnbctype_upload"
os.makedirs(OUT, exist_ok=True)
s = pd.read_csv("results/scores_tnbc.csv")
flag = set(open("data/tnbctype/er_flagged_samples.txt").read().split())

for c, f, log in [("TCGA", f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", True),
                  ("METABRIC", f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt", False)]:
    e = pd.read_csv(f, sep="\t").dropna(subset=["Hugo_Symbol"])
    ids = [x for x in s.loc[s.cohort == c, "sample"] if x in e.columns]
    e = e[["Hugo_Symbol"] + ids].dropna()
    e = e.loc[e[ids].mean(axis=1).groupby(e.Hugo_Symbol).idxmax()].set_index("Hugo_Symbol")
    if log:
        e = np.log2(e + 1)
    e = e.round(4)
    e.to_csv(f"{OUT}/{c}_TNBC_expression.csv")
    keep = [x for x in e.columns if x not in flag]
    e[keep].to_csv(f"{OUT}/{c}_TNBC_expression_noERflag.csv")
    print(f"{c}: {e.shape[0]} genes; {e.shape[1]} TNBC; {e.shape[1] - len(keep)} ER-flagged removed -> {len(keep)}")
