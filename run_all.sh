#!/usr/bin/env bash
# Run the whole HER-OMI pipeline in order, from the project folder:   bash run_all.sh
# Each step's printed output is also saved to logs/<step>.log.
# Step 07 needs your TNBCtype result files in data/tnbctype/ (already included in this repository).
set -euo pipefail
mkdir -p logs results
STEPS="01_define_tnbc 01b_signature 02_singscore 02b_annotate 02c_tnbc_expression_tables
       03_basal_only_test 04_her2_ihc_check 05_pathway_comparison 06_confounder_test 07_tnbctype
       08_protein_test 08b_protein_robustness 09_pathway_overlap_check 10_her2pos_comparison
       10b_her2pos_comparison_metabric 11_dna_layer 12_erlow_sensitivity
       13_normal_comparison 14_normal_admixture 15_cluster_heatmap 16_slide_figures"
for s in $STEPS; do
  echo "=== $s"
  python3 scripts/$s.py 2>&1 | tee logs/$s.log
done
