"""HER-OMI step 1: define TNBC in TCGA-BRCA and METABRIC from cBioPortal clinical tables.

Data: cBioPortal datahub (github.com/cBioPortal/datahub), studies brca_tcga (legacy, has IHC fields)
and brca_metabric. Run from the folder that holds data/raw.
"""
import pandas as pd

RAW = "data/raw"

def read(path):
    return pd.read_csv(path, sep="\t", comment="#")

# TCGA: receptor status by IHC, HER2 resolved by FISH when IHC is equivocal or missing (ASCO/CAP logic)
t = read(f"{RAW}/brca_tcga/data_clinical_patient.txt")
er = t.ER_STATUS_BY_IHC.eq("Negative")
pr = t.PR_STATUS_BY_IHC.eq("Negative")
ihc, fish = t.IHC_HER2, t.HER2_FISH_STATUS
her2_neg = (ihc.eq("Negative")
            | (ihc.isin(["Equivocal", "Indeterminate"]) & fish.eq("Negative"))
            | (ihc.isin(["[Not Available]"]) & fish.eq("Negative")))
t["TNBC"] = er & pr & her2_neg
print("TCGA TNBC (IHC, FISH for equivocal or missing IHC):", t.TNBC.sum())   # 160

# METABRIC: no PR or HER2 IHC exists. ER_IHC is IHC; PR_STATUS and HER2_STATUS are not.
m = read(f"{RAW}/brca_metabric/data_clinical_patient.txt").merge(
    read(f"{RAW}/brca_metabric/data_clinical_sample.txt"), on="PATIENT_ID")
er_ihc = m.ER_IHC.eq("Negative")                 # note: positives are spelled "Positve" in the source
pr_neg = m.PR_STATUS.eq("Negative")
# Chosen definition (option 3): ER by IHC, PR by expression (only option), HER2 not gained by SNP6 copy number.
# Rejected alternatives, kept for the write-up: preprint's all-status-column rule gives 320 (51 are ER+ by IHC);
# ER_IHC plus HER2_STATUS gives 276.
m["TNBC_cn"] = er_ihc & pr_neg & m.HER2_SNP6.isin(["NEUTRAL", "LOSS"])
print("METABRIC TNBC (ER IHC, HER2 copy number):", m.TNBC_cn.sum())   # 250
