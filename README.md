# HER-OMI

**Is the IBCGA "HER2-like" state in triple-negative breast cancer actually HER2? A DNA, RNA and protein test in Western cohorts.**

DH607 Introduction to Computational Multi-omics, IIT Bombay, Autumn 2026. Hardika Jain (25B2701).

## Question

The Indian Breast Cancer Genome Atlas (IBCGA) reported a 144-gene signature that marks a "HER2-like" transcriptional state in HER2-negative, triple-negative breast cancer (TNBC), and proposed these tumours as candidates for HER2-directed antibody-drug conjugates. The signature had not been tested outside the Indian cohort, and HER2 protein was never measured. This project tests it in TCGA, METABRIC and CPTAC, layer by layer.

## Data

All public: TCGA-BRCA, METABRIC and CPTAC 2020 breast (via the cBioPortal datahub), the IBCGA preprints' supplementary tables (medRxiv 10.64898/2026.03.27.26349506; Research Square rs-10980563), HGNC and MSigDB. `scripts/00_download_data.sh` downloads everything (about 1 GB). TNBCtype subtype calls are in `data/tnbctype/`.

## Main findings

* The signature replicates in TCGA and METABRIC, even inside basal TNBC (score vs ERBB2 RNA 0.30 to 0.55, above 1,000 random gene sets).
* Its top end is the AR-driven LAR subtype: every LAR tumour falls in the top third of the score (AUC 0.97 and 0.96).
* At the protein level the score does not track HER2 or active HER2 (0.21 and 0.18, inside the random range) but does track AR (0.50) and INPP4B (0.55); CPTAC mass spectrometry agrees.
* ERBB2 copy number does not explain the link, and in matched patients TNBC ERBB2 is no higher than in the patient's own normal breast.
* Conclusion: the "HER2-like" state is an AR/luminal programme, so HER2 drug eligibility should still be judged by HER2 protein.

## How to run

See `CHECKLIST.md` for setup on a Mac and a step-by-step walkthrough, including the numbers each step should print. In short:

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
bash scripts/00_download_data.sh
bash run_all.sh
```

## Layout

| Folder | Contents |
|---|---|
| `scripts/` | One script per analysis step, numbered in run order; shared functions in `heromi_lib.py` |
| `data/` | Small inputs kept in the repo (signature, coverage table, pathway tables, TNBCtype calls); large downloads are ignored by git |
| `results/` | Tables and figures produced by the scripts |
| `CHECKLIST.md` | Every step, what it should print, questions and papers |
| `NOTES.md` | My explanations of each step |
| `WALKTHROUGH.md` | A plain-language walkthrough of every script, with the key code |
| `guide/` | The study guide (PDF and Markdown) that explains the whole project from scratch |

## Use of AI tools

Analysis code was developed with the help of Claude (Anthropic). The research question, the choice of signature and cohorts, the TNBC definitions, the decision to test protein and DNA layers, and the interpretation were mine; each step was re-run and checked by me (see `CHECKLIST.md`).
