"""HER-OMI step 9: are IBCGA's "HER2-active" pathway scores circular?

IBCGA Fig. 7A shows HER2-like tumours scoring high on three Reactome ERBB2 programmes (GRB7-mediated ERBB2
signalling, constitutive signalling by overexpressed ERBB2, ERBB2 ECD mutant signalling) and reads this as
transcriptionally active HER2. These sets are small and contain ERBB2 and GRB7, which are themselves
signature and amplicon genes. Here each set is scored in TCGA and METABRIC TNBC three ways (all genes; ERBB2
and other amplicon genes removed; every signature gene removed) and correlated with the full signature score.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, "scripts")
from heromi_lib import load_expression, singscore

RAW = "data/raw"
SETS = ["REACTOME_GRB7_EVENTS_IN_ERBB2_SIGNALING", "REACTOME_CONSTITUTIVE_SIGNALING_BY_OVEREXPRESSED_ERBB2",
        "REACTOME_SIGNALING_BY_ERBB2_ECD_MUTANTS", "REACTOME_SIGNALING_BY_ERBB2"]
gmt = {f[0]: f[2:] for f in (l.rstrip("\n").split("\t") for l in open("data/ref/c2.cp.reactome.gmt")) if f[0] in SETS}
sig = pd.read_csv("data/signature_coverage.csv")
SIG, AMP = set(sig.gene), set(sig.gene[sig.on_her2_amplicon])
ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")
FILES = {"TCGA": f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt",
         "METABRIC": f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt"}
for s in SETS:
    g = gmt[s]
    print(f"{s}: {len(g)} genes; signature genes inside: {sorted(set(g) & SIG)}")
rows = []
for c, path in FILES.items():
    a = ann[ann.cohort.eq(c)]
    expr = load_expression(path, set(a.index))
    sym = pd.read_csv(path, sep="\t", usecols=[0, 1]).dropna().drop_duplicates("Hugo_Symbol")
    s2e = dict(zip(sym.Hugo_Symbol, sym.Entrez_Gene_Id.astype(int)))
    ids = lambda genes: [s2e[x] for x in genes if x in s2e and s2e[x] in expr.index]
    for s in SETS:
        for version, genes in [("all genes", gmt[s]), ("amplicon genes removed", [x for x in gmt[s] if x not in AMP]),
                               ("all signature genes removed", [x for x in gmt[s] if x not in SIG])]:
            k = ids(genes)
            ps = singscore(expr, k, []).total
            rows.append(dict(cohort=c, gene_set=s.replace("REACTOME_", ""), version=version, n_genes=len(k),
                             rho_full_score=stats.spearmanr(ps, a.loc[expr.columns, "score_full"])[0],
                             rho_ERBB2_mRNA=stats.spearmanr(ps, a.loc[expr.columns, "ERBB2_expr"])[0]))
out = pd.DataFrame(rows)
out.to_csv("results/pathway_overlap_check.csv", index=False)
pd.set_option("display.width", 200)
print(out.round(2).to_string(index=False))
