"""HER-OMI step 2b: add PAM50 calls and LAR / luminal marker genes to the scores.

Writes results/scores_tnbc_annot.csv, the input for steps 3 to 6. (In Session 1 this file was made
interactively; scripted in Session 2 so the pipeline rebuilds end to end.)

PAM50: TCGA from the PanCancer Atlas SUBTYPE column (the legacy brca_tcga study has none);
METABRIC from CLAUDIN_SUBTYPE. Markers: TCGA as log2(RSEM + 1); METABRIC array values as given.
Output matches the Session 1 file exactly, except DHCR24 is now filled in for TCGA (it was blank
there; DHCR24 is not used in any analysis).
"""
import numpy as np
import pandas as pd

import sys; sys.path.insert(0, "scripts")
from heromi_lib import load_expression

RAW = "data/raw"
MARKERS = {"APOD": 347, "AR": 367, "FOXA1": 3169, "GATA3": 2625, "KRT18": 3875,
           "PIP": 5304, "SPDEF": 25803, "TFF1": 7031, "XBP1": 7494, "DHCR24": 1718}

s = pd.read_csv("results/scores_tnbc.csv")
s["pid"] = np.where(s.cohort.eq("TCGA"), s["sample"].str[:12], s["sample"])

pc = pd.read_csv(f"{RAW}/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", sep="\t", comment="#")
mb = pd.read_csv(f"{RAW}/brca_metabric/data_clinical_patient.txt", sep="\t", comment="#")
s["SUBTYPE"] = s.pid.map(dict(zip(pc.PATIENT_ID, pc.SUBTYPE))).where(s.cohort.eq("TCGA"))
s["CLAUDIN_SUBTYPE"] = s.pid.map(dict(zip(mb.PATIENT_ID, mb.CLAUDIN_SUBTYPE))).where(s.cohort.eq("METABRIC"))
s["pam50"] = (s.SUBTYPE.str.replace("BRCA_", "", regex=False).fillna(s.CLAUDIN_SUBTYPE)
              .replace({"Her2": "HER2E"}))

ids = set(s["sample"])
expr = {
    "TCGA": np.log2(load_expression(f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", ids) + 1),
    "METABRIC": load_expression(f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt", ids),
}
for g, eid in MARKERS.items():
    s[g] = np.nan
    for c, e in expr.items():
        if eid in e.index:
            m = s.cohort.eq(c)
            s.loc[m, g] = s.loc[m, "sample"].map(e.loc[eid]).values

s.to_csv("results/scores_tnbc_annot.csv", index=False)
print(s.groupby("cohort").pam50.value_counts(dropna=False).unstack(0))
