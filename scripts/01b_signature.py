"""HER-OMI step 1b: extract the IBCGA signature and pathway tables, and check gene coverage.

Inputs : data/ibcga/ibcga_supp_tables.xlsx (medRxiv supplementary tables), data/ref/hgnc_complete_set.txt,
         the TCGA and METABRIC expression matrices.
Outputs: data/IBCGA_HER2like_signature_144genes.csv  the 144 genes of Table S21
         data/signature_coverage.csv                  each gene, its HGNC type and Entrez ID, and whether
                                                      it is measured in TCGA and in METABRIC
         data/ibcga_S22.csv, data/ibcga_S23.csv       pathway tables (HER2-like vs basal; LAR vs basal)
A gene counts as measured if its symbol or its Entrez ID appears in the matrix (Entrez IDs catch genes
renamed since the platform was annotated).
"""
import pandas as pd

RAW = "data/raw"
XLSX = "data/ibcga/ibcga_supp_tables.xlsx"

# ---- Table S21: the signature ----
sig = pd.read_excel(XLSX, "Table S21", header=3).iloc[:, 1:]
sig.columns = ["gene", "overlap_count", "categories", "direction", "on_her2_amplicon"]
sig = sig.dropna(subset=["gene"])
sig["on_her2_amplicon"] = sig.on_her2_amplicon.eq("Yes")
sig.to_csv("data/IBCGA_HER2like_signature_144genes.csv", index=False)
print(f"S21: {len(sig)} genes, {sig.gene.nunique()} unique; {sig.direction.value_counts().to_dict()}; "
      f"amplicon genes: {sig.loc[sig.on_her2_amplicon, 'gene'].tolist()}")

# ---- coverage in TCGA and METABRIC ----
h = pd.read_csv("data/ref/hgnc_complete_set.txt", sep="\t", low_memory=False,
                usecols=["symbol", "entrez_id", "locus_group", "locus_type", "prev_symbol", "alias_symbol"])
cov = sig.merge(h.rename(columns={"symbol": "gene"}), on="gene", how="left")
print("not found in HGNC under this symbol:", cov.loc[cov.entrez_id.isna(), "gene"].tolist())
for name, f in [("TCGA", f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt"),
                ("METABRIC", f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt")]:
    e = pd.read_csv(f, sep="\t", usecols=[0, 1])
    ent, sym = set(e.Entrez_Gene_Id.dropna().astype(int)), set(e.Hugo_Symbol.dropna())
    cov[name] = cov.gene.isin(sym) | cov.entrez_id.fillna(-1).astype(int).isin(ent)
    coding = cov.locus_group.eq("protein-coding gene")
    print(f"{name}: {cov[name].sum()}/144 measured; protein-coding {cov.loc[coding, name].sum()}/{coding.sum()}")
both = cov.TCGA & cov.METABRIC
print(f"measured in both: {both.sum()} ({cov[both].direction.value_counts().to_dict()}; "
      f"amplicon {cov[both].on_her2_amplicon.sum()})")
cov.to_csv("data/signature_coverage.csv", index=False)


# ---- Tables S22 and S23: pathways ----
def pathway_table(sheet):
    d = pd.read_excel(XLSX, sheet, header=None).dropna(how="all").dropna(axis=1, how="all")
    start = d.index[d.iloc[:, 0].eq("Gene_set")][0]
    d.columns = ["set", "diff", "p", "padj"]
    return d.loc[start + 1:].reset_index(drop=True).astype({"diff": float, "p": float, "padj": float})


for sheet, out in [("Table S22", "data/ibcga_S22.csv"), ("Table S23", "data/ibcga_S23.csv")]:
    t = pathway_table(sheet)
    t.to_csv(out, index=False)
    print(f"{sheet}: {len(t)} pathways, {(t.padj < 0.05).sum()} with adjusted p < 0.05")
