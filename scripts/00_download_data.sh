#!/usr/bin/env bash
# HER-OMI step 0: download everything the analysis needs (about 1 GB). Run from the project folder.
set -euo pipefail
B=https://media.githubusercontent.com/media/cBioPortal/datahub/master/public
mkdir -p data/raw/brca_tcga data/raw/brca_metabric data/ibcga data/ref

for f in data_clinical_patient.txt data_clinical_sample.txt data_mrna_seq_v2_rsem.txt; do
  curl -sSL -o data/raw/brca_tcga/$f $B/brca_tcga/$f; done
for f in data_clinical_patient.txt data_clinical_sample.txt data_mrna_illumina_microarray.txt; do
  curl -sSL -o data/raw/brca_metabric/$f $B/brca_metabric/$f; done

# IBCGA supplementary tables (Khanna et al. 2026, medRxiv), Table S21 is the signature
curl -sSL -A "Mozilla/5.0" -o data/ibcga/ibcga_supp_tables.xlsx \
  "https://www.medrxiv.org/content/medrxiv/early/2026/03/30/2026.03.27.26349506/DC2/embed/media-2.xlsx?download=true"

# HGNC gene table, used to match signature symbols to Entrez IDs across old and new gene names
curl -sSL -o data/ref/hgnc_complete_set.txt \
  https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

# MSigDB Hallmark and Reactome gene sets, the collections IBCGA used in Tables S22 and S23
M=https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs
curl -sSL -o data/ref/h.all.gmt $M/h.all.v2024.1.Hs.symbols.gmt
curl -sSL -o data/ref/c2.cp.reactome.gmt $M/c2.cp.reactome.v2024.1.Hs.symbols.gmt

# PanCancer Atlas clinical table, for TCGA PAM50 calls (SUBTYPE column; the legacy study has none).
# Added in Session 2; used by 02b_annotate.py.
mkdir -p data/raw/brca_tcga_pan_can_atlas_2018
curl -sSL -o data/raw/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt \
  $B/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt

# ---- added in Session 3: protein and DNA layers ----
# TCGA RPPA protein (step 8) and hg19 copy-number segments (step 10)
for f in data_rppa.txt data_cna_hg19.seg; do
  curl -sSL -o data/raw/brca_tcga/$f $B/brca_tcga/$f; done
# METABRIC discrete copy-number calls (GISTIC, -2 to 2), for steps 10b and 11
curl -sSL -o data/raw/brca_metabric/data_cna.txt $B/brca_metabric/data_cna.txt
# CPTAC 2020 breast: clinical, RNA (FPKM), mass-spec proteome and phosphoproteome (step 8)
mkdir -p data/raw/brca_cptac_2020
for f in data_clinical_patient.txt data_clinical_sample.txt data_mrna_seq_fpkm.txt \
         data_protein_quantification.txt data_phosphoprotein_quantification.txt; do
  curl -sSL -o data/raw/brca_cptac_2020/$f $B/brca_cptac_2020/$f; done
# IBCGA follow-up preprint (Research Square rs-10980563, Sept 2026): PDF, supplementary text and tables
mkdir -p data/ibcga2
curl -sSL -A "Mozilla/5.0" -o data/ibcga2/main.pdf "https://www.researchsquare.com/article/rs-10980563/v1.pdf?c=1789632336000"
curl -sSL -A "Mozilla/5.0" -o data/ibcga2/supp_text.pdf "https://assets-eu.researchsquare.com/files/rs-10980563/v1/b325e1fadc6a534e9ab46536.pdf"
curl -sSL -A "Mozilla/5.0" -o data/ibcga2/supp_tables.xlsx "https://assets-eu.researchsquare.com/files/rs-10980563/v1/a292d982d7bf710ebc245cfe.xlsx"

# Not downloadable, kept in the project instead:
#   data/signature_coverage.csv, data/ibcga_S22.csv, data/ibcga_S23.csv (made interactively in Session 1)
#   data/tnbctype/<COHORT>_{result,coef,pvalue}.csv (TNBCtype web tool output for the noERflag uploads)

# ---- added for step 13 (normal vs TNBC vs HER2-positive) ----
# TCGA PanCancer RNA-seq for tumours and for the 114 adjacent normal breast samples, processed together
mkdir -p data/raw/brca_tcga_pan_can_atlas_2018/normals
curl -sSL -o data/raw/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt \
  $B/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt
curl -sSL -o data/raw/brca_tcga_pan_can_atlas_2018/normals/data_mrna_seq_v2_rsem_normal_samples.txt \
  $B/brca_tcga_pan_can_atlas_2018/normals/data_mrna_seq_v2_rsem_normal_samples.txt
