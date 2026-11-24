# HER-OMI checklist

Every analysis in this project, in the order it runs. For each step: run it yourself, check your numbers match, read the code until you can explain it, read the paper behind it, and question it. Write your answers in `NOTES.md` in your own words.

Each step has four boxes. Tick them in this file as you go (change `[ ]` to `[x]`), so the GitHub history shows your progress.

* **Ran it**: you ran the script and your numbers match "You should see".
* **Understood the code**: you can explain every block of the script to someone else.
* **Read the papers**: you have read the listed papers (at least abstract, methods and the relevant figure).
* **Questions resolved**: you have written answers to the questions in `NOTES.md`, or found a problem and we fixed it.

Items marked **(publication)** are not needed for the talk, but a reviewer would ask for them if this becomes a preprint.

Citations marked *(verify)* were written from memory. Check volume, pages and year against the journal before using them.

---

## Schedule (talk on 24 November)

| When | Session | Steps |
|---|---|---|
| 7 to 9 Nov | 1 | Setup, Steps 0, 1, 1b |
| 10 to 16 Nov (endsems) | Reading only, optional | Singscore paper (Step 2), Lehmann 2011 (Step 3) |
| 17 Nov | 2 | Steps 2, 2b, 2c, and the R cross-check |
| 18 Nov | 3 | Steps 3, 4, 5 |
| 19 Nov | 4 | Steps 6, 7, 12 |
| 20 Nov | 5 | Steps 8, 8b, 9 (the main new result) |
| 21 Nov | 6 | Steps 10, 10b, 11, and GitHub setup |
| 22 to 23 Nov | | Slides and rehearsal (about 20 minutes plus 10 for questions) |

---

## Setup (once, on your Mac)

1. Open **Terminal** (press Cmd + Space, type "Terminal").
2. This Mac has Anaconda (the Terminal line starts with `(base)`), so we use conda to get Python 3.13, the version the scripts were tested with. (Python 3.15 from python.org is also installed but is too new for some packages.)
3. Unzip the project folder somewhere easy, for example your Desktop. Then move into it:
   `cd ~/Desktop/heromi`
4. Make a separate environment for this project and install the packages:
   ```
   conda create -n heromi python=3.13 -y
   conda activate heromi
   pip install -r requirements.txt
   ```
5. **Every time you come back:** open Terminal, then `cd ~/Desktop/heromi` and `conda activate heromi`. You will see `(heromi)` at the start of the line.

- [ ] Setup done: `python -c "import pandas, scipy, matplotlib; print('ok')"` prints `ok`

**Useful:** to read a script, open it in any text editor (VS Code is free and colours the code). Every script starts with a block of plain English between `"""` marks that says what it does. Read that first.

---

## Step 0: Download the data

**What it does:** downloads about 1 GB of public data: TCGA, METABRIC and CPTAC tables from cBioPortal, the two IBCGA preprints' supplementary files, the HGNC gene list, and MSigDB gene sets.

**Run:** `bash scripts/00_download_data.sh` (takes a few minutes; nothing prints if it works)

**You should see:** a `data/raw/` folder containing `brca_tcga`, `brca_metabric`, `brca_tcga_pan_can_atlas_2018` and `brca_cptac_2020`.

**Questions to ask:**
1. What is cBioPortal, and who made these tables? Why use the "legacy" `brca_tcga` study rather than the newer PanCancer Atlas one? (Hint: open both clinical files and look for ER, PR and HER2 columns.)
2. What is the difference between RNA-seq (TCGA) and a microarray (METABRIC)? Why does it matter later?

**Read:** Cerami et al. 2012, *Cancer Discovery* 2(5):401 (cBioPortal) *(verify)*.

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 1: Define TNBC in each cohort

**What it does:** picks the triple-negative tumours. TCGA: ER and PR negative by IHC; HER2 negative by IHC, or by FISH when IHC is equivocal or missing. METABRIC (your option 3): ER negative by IHC, PR negative (mRNA based, the only option), and no ERBB2 copy gain on SNP6 arrays.

**Run:** `python3 scripts/01_define_tnbc.py`

**You should see:** `TCGA TNBC ...: 160` and `METABRIC TNBC ...: 250`

**Questions to ask:**
1. What do ASCO/CAP guidelines say counts as ER negative and HER2 negative? Where does each cohort fall short of that? (See the section "How TNBC gets defined in each cohort" in the Papers, Ideas and Methods doc.)
2. Why is defining HER2 from ERBB2 mRNA circular for *this* project?
3. Change the TCGA rule to "HER2 IHC Negative only" (ignore FISH). How many tumours do you get? (Answer: 116.) Which rule is more defensible?
4. Why did the published METABRIC TNBC set have 320 tumours, and why are 51 of them ER positive by IHC?

**Read:** Allison et al. 2020, *J Clin Oncol* (ASCO/CAP ER and PR guideline); Wolff et al. 2018, *J Clin Oncol* 36(20):2105 (HER2 guideline, reaffirmed 2023) *(verify)*; TCGA Network 2012, *Nature* 490:61; Curtis et al. 2012, *Nature* 486:346 (METABRIC).

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved
- [ ] **(publication)** Find METABRIC documentation confirming that PR_STATUS and HER2_STATUS are mRNA based (we inferred it from the data: each splits perfectly by its own gene's mRNA).

---

## Step 1b: Extract the signature and check which genes exist in Western data

**What it does:** reads Table S21 (the 144 genes), checks which are measured in TCGA and METABRIC (matching by gene symbol or Entrez ID, using HGNC), and extracts the pathway tables S22 and S23.

**Run:** `python3 scripts/01b_signature.py`

**You should see:** `TCGA: 114/144 measured; protein-coding 106/107`, `METABRIC: 112/144 measured; protein-coding 106/107`, `measured in both: 109 ({'Downregulated': 67, 'Upregulated': 42}; amplicon 9)`, and `Table S22: 1741 pathways`.

**Questions to ask:**
1. How did IBCGA build the signature? (Research Square preprint, Methods, "HER2-associated transcriptional program".) Is a gene list from IHC 3+ vs IHC 0 tumours a sensible way to find a "HER2-like" state?
2. Most missing genes are non-coding RNAs and pseudogenes. Why would a microarray lack them? What would we miss if the signal lives in those genes?
3. Why match by Entrez ID as well as symbol? (Hint: FAM136A.)
4. Why use only the 109 genes present in both cohorts, instead of everything each cohort has?

**Read:** Khanna et al. 2026, medRxiv (Methods); Chowdhury, Khanna et al. 2026, Research Square, Methods section on the HER2-associated signature.

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 2: Score every tumour (singscore)

**What it does:** within each tumour, ranks all ~20,000 genes, then asks whether the signature's "up" genes sit near the top and "down" genes near the bottom. Gives a score from −1 to +1. Three versions: all 109 genes ("full"), without the 9 amplicon genes ("no_amp", the main one), and amplicon genes only.

**Run:** `python3 scripts/02_singscore.py` (the singscore maths is in `scripts/heromi_lib.py`)

**You should see:** `TCGA: 159 TNBC tumours, 20505 genes ranked` and `METABRIC: 250 TNBC tumours, 20192 genes ranked`. Mean no_amp score: TCGA −0.078, METABRIC −0.149.

**Questions to ask:**
1. Why a rank-based score instead of averaging expression values? (Think RNA-seq counts vs microarray intensities.)
2. Work through `singscore()` in `heromi_lib.py` with a pen: why are the "down" genes ranked in reverse? Where do the −0.5 and +0.5 limits come from?
3. Why does removing the 9 amplicon genes matter for the project's question?
4. Why is "full vs no_amp correlate 0.99" *not* evidence of anything?

**Read:** Foroutan et al. 2018, *BMC Bioinformatics* 19:404 (singscore) *(verify)*.

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved
- [ ] **Cross-check with R** (you have R): install the official package (`BiocManager::install("singscore")`), score the TCGA table from Step 2c with the same up and down genes, and check it matches our Python scores. We will do this together on 17 Nov. (Talk and publication: it shows our home-made code is correct.)

---

## Step 2b: Add PAM50 subtypes and marker genes

**What it does:** attaches each tumour's PAM50 subtype (Basal, HER2-enriched, luminal, normal-like) and the expression of luminal/AR marker genes (AR, FOXA1, GATA3 and others).

**Run:** `python3 scripts/02b_annotate.py`

**You should see:** TCGA Basal 116, HER2E 12, LumA 4, LumB 1, Normal 9 (17 have no call); METABRIC Basal 141, claudin-low 77, HER2E 17, Normal 10, LumA 5.

**Questions to ask:**
1. What is PAM50? Why can a TNBC tumour be "HER2-enriched" by PAM50 when it is HER2 negative clinically?
2. What is claudin-low, and why does only METABRIC have it?

**Read:** Parker et al. 2009, *J Clin Oncol* 27(8):1160 (PAM50).

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 2c: Make the TNBCtype upload tables

**What it does:** writes one expression table per cohort (genes by tumours, log2) in the format TNBCtype accepts, with and without the 23 tumours TNBCtype flagged as possibly ER positive (listed in `data/tnbctype/er_flagged_samples.txt`).

**Run:** `python3 scripts/02c_tnbc_expression_tables.py`

**You should see:** `TCGA: 20513 genes; 159 TNBC; 12 ER-flagged removed -> 147` and `METABRIC: 20384 genes; 250 TNBC; 11 ER-flagged removed -> 239`.

**Questions to ask:**
1. TNBCtype flagged tumours that are ER negative by IHC but express ESR1 mRNA. What is "ER-low", and why might TCGA's records hide it?

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 3: Inside basal tumours only, does the score still track ERBB2?

**What it does:** restricts to PAM50 basal TNBC (where LAR is rare) and correlates the score with ERBB2 mRNA. Compares against 1,000 random gene sets of the same size, removes all chromosome 17 genes, and controls for a 5-gene LAR marker average.

**Run:** `python3 scripts/03_basal_only_test.py`

**You should see:** basal only, score vs ERBB2: TCGA 0.30 (random-set 95th percentile 0.22, p 0.019); METABRIC 0.48 (95th percentile 0.20, p 0.001).

**Questions to ask:**
1. Why do random gene sets correlate about 0.2 with ERBB2? What would we have wrongly concluded without this comparison?
2. What is a partial correlation? What does "controlling for LAR" do here?
3. Why remove chromosome 17 genes as well as the amplicon genes?

**Read:** Lehmann et al. 2011, *J Clin Invest* 121(7):2750 (TNBC subtypes, LAR); Venet, Dumont and Detours 2011, *PLoS Comput Biol* 7(10):e1002240, "Most random gene expression signatures are significantly associated with breast cancer outcome" *(verify)*.

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 4: Does the score rise with HER2 IHC score (0, 1+, 2+)?

**What it does:** TCGA only (METABRIC has no HER2 IHC). Tests for a trend across IHC categories and compares HER2-low (1+, 2+) with HER2-0.

**Run:** `python3 scripts/04_her2_ihc_check.py`

**You should see:** all TNBC n = 88, trend rho 0.21 (p 0.047), HER2-low vs HER2-0 p 0.025; basal n = 60, rho 0.10 (p 0.43).

**Questions to ask:**
1. Why does HER2 IHC 0 vs 1+ suddenly matter clinically? (Trastuzumab deruxtecan.)
2. With only 17 basal HER2-0 tumours, what size of effect could this test detect? Why is "not significant" not the same as "no effect"?

**Read:** Modi et al. 2022, *N Engl J Med* 387(1):9 (DESTINY-Breast04); Schettini et al. 2021, *npj Breast Cancer* 7:1 (HER2-low biology).

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 5: Does the signature carry IBCGA's pathway biology?

**What it does:** scores ~1,350 Hallmark and Reactome pathways in every tumour (with all 144 signature genes removed from each pathway), correlates each with the score, and compares that pattern with IBCGA's Table S22 (HER2-like vs basal) and S23 (LAR vs basal).

**Run:** `python3 scripts/05_pathway_comparison.py`

**You should see:** TCGA vs METABRIC profiles agree (rho 0.81 all, 0.63 basal); agreement with S22: TCGA all 0.70, basal 0.54; METABRIC all 0.68, basal 0.61.

**Questions to ask:**
1. Why remove the signature genes from every pathway first?
2. S22 never says which group was subtracted from which. How did we infer the direction? How sure are we?
3. Agreement with S22 shows replication. Why does it not show the axis is *specific* to HER2?

**Read:** Liberzon et al. 2015, *Cell Systems* 1(6):417 (Hallmark gene sets) *(verify)*.

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved
- [ ] **(publication)** Confirm the S22/S23 sign convention with the IBCGA authors or their code.

---

## Step 6: Controls for proliferation and oestrogen signalling

**What it does:** within basal tumours, checks whether the score–ERBB2 link survives controlling for proliferation, the Hallmark oestrogen response, ESR1 mRNA and the LAR proxy, each against random gene sets.

**Run:** `python3 scripts/06_confounder_test.py` (slow: a few minutes)

**You should see:** no controls: TCGA 0.30, METABRIC 0.48; all four controls: TCGA 0.18 (borderline), METABRIC 0.40 (p 0.001).

**Questions to ask:**
1. What is over-adjustment? Why might controlling for the oestrogen response remove real signal?
2. ESR1 itself changes nothing, but the oestrogen *response* gene set does. What does that suggest the gene set is really measuring in ER negative tumours?
3. Why is TCGA weaker than METABRIC? List three possible reasons.

**Read:** Nielsen et al. 2010, *Clin Cancer Res* 16(21):5222 (PAM50 proliferation genes) *(verify)*.

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 7: Is the HER2-like group just LAR? (TNBCtype)

**What it does:** uses the Lehmann subtype calls you got from the TNBCtype website. Checks how well the score picks out LAR, and whether the score–ERBB2 link survives removing LAR.

**Before running:** put your TNBCtype output files in `data/tnbctype/`, named `TCGA_result.csv`, `TCGA_coef.csv`, `METABRIC_result.csv`, `METABRIC_coef.csv`.

**Run:** `python3 scripts/07_tnbctype.py`

**You should see:** LAR tumours in the top score third: TCGA 13 of 13, METABRIC 26 of 26; AUC for LAR 0.97 and 0.96; basal with LAR removed: TCGA 0.30, METABRIC 0.47.

**Questions to ask:**
1. If the score finds every LAR tumour, is the "HER2-like group" new at all?
2. Why does controlling for LAR make the TCGA basal result *stronger* (0.30 to 0.40)?

**Read:** Chen et al. 2012, *Cancer Informatics* 11:147 (TNBCtype tool) *(verify)*; Lehmann et al. 2011 again, LAR section.

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 8: The protein test (the main new result)

**What it does:** compares the score with HER2 protein, phosphorylated (active) HER2, and other proteins measured by RPPA in TCGA (131 TNBC) and by mass spectrometry in CPTAC (28 TNBC). Asks whether the score says anything about HER2 protein beyond ERBB2 mRNA.

**Run:** `python3 scripts/08_protein_test.py`

**You should see:** TCGA all TNBC: HER2 0.21 (p 0.09), HER2 pY1248 0.18 (p 0.14), AR 0.50 (p 0.001), INPP4B 0.55 (p 0.001); score vs HER2 given ERBB2 mRNA −0.06. CPTAC: AR 0.75 (p 0.001). The top score third closes only 23% of the gap to HER2-positive tumours.

**Questions to ask:**
1. What is RPPA, and why is it relative rather than an absolute amount like IHC?
2. What does a partial correlation of zero "given ERBB2 mRNA" mean in plain words?
3. Trastuzumab deruxtecan works even at very low HER2. So what exactly does our result say, and not say, about IBCGA's ADC proposal?
4. Why are AR and INPP4B "positive controls"?

**Read:** the IBCGA Research Square preprint, Fig. 7 and the paragraph proposing ADC evaluation; Krug et al. 2020, *Cell* (CPTAC breast proteogenomics) *(verify)*; Li et al. 2013, *Nature Methods*, TCPA (TCGA RPPA) *(verify)*.

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 8b: Is the protein result robust?

**What it does:** repeats Step 8 without the ER-low tumours, and within HER2-low tumours only (the patients eligible for trastuzumab deruxtecan).

**Run:** `python3 scripts/08b_protein_robustness.py`

**You should see:** without ER-low: HER2 0.27 (p 0.027) but **0.00 given ERBB2 mRNA**; HER2-low only (n 61): HER2 −0.02 (p 0.92), AR 0.68 (p 0.001).

**Questions to ask:**
1. Removing ER-low tumours makes the HER2 protein link significant. Why does this not rescue IBCGA's claim?
2. Why is the HER2-low subset the fairest test of the ADC idea?

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 9: Are IBCGA's "HER2-active" pathway scores circular?

**What it does:** scores the three Reactome ERBB2 sets IBCGA used, with and without the genes they share with the signature (ERBB2, GRB7).

**Run:** `python3 scripts/09_pathway_overlap_check.py`

**You should see:** GRB7 events set vs full score: TCGA 0.25 with all genes, 0.10 without ERBB2 and GRB7; METABRIC 0.21 and 0.05. "Constitutive signalling by overexpressed ERBB2" goes negative (−0.35, −0.43).

**Questions to ask:**
1. What is circularity in a gene set test, and why does a 5-gene set make it worse?
2. IBCGA's comparison included HER2-positive tumours; ours is TNBC only. Is our critique fair?

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 10 and 10b: Comparison with true HER2-positive tumours

**What it does:** puts HER2-like TNBC between true HER2-positive (HR negative) tumours and basal TNBC across DNA, RNA and protein, as heatmaps. 10 is TCGA, 10b is METABRIC (no protein).

**Run:** `python3 scripts/10_her2pos_comparison.py` then `python3 scripts/10b_her2pos_comparison_metabric.py`

**You should see:** two PNG heatmaps in `results/`. TCGA medians: ERBB2 mRNA 16.7 (HER2-positive), 12.3 (HER2-like), 11.7 (basal); AR protein 0.49, 0.51, −0.54.

**Questions to ask:**
1. HR negative HER2-positive tumours are AR-rich. What is a "molecular apocrine" tumour, and how does it explain why the signature calls AR-high TNBC "HER2-like"?
2. Some "HER2-positive" tumours on the right of the TCGA heatmap have no amplification. What could they be?

**Read:** Farmer et al. 2005, *Oncogene* 24:4660, "Identification of molecular apocrine breast tumours by microarray analysis" *(verify)*.

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved

---

## Step 11: Does DNA at the ERBB2 locus explain it?

**What it does:** checks ERBB2 copy number (low-level gains and losses) and DNA breakpoints near IKZF3, the rearrangement IBCGA reported in HER2-positive tumours.

**Run:** `python3 scripts/11_dna_layer.py`

**You should see:** score vs ERBB2 mRNA given copy number: TCGA 0.43, METABRIC 0.55 (barely changed). Breakpoint near IKZF3 in 9 of 155 TCGA TNBC, with higher ERBB2 mRNA (p 0.001). The ERBB2 window ranks 2nd of 95 chromosome 17 windows.

**Questions to ask:**
1. Why is a copy-number breakpoint only a *proxy* for a rearrangement?
2. Why compare the ERBB2 window against every other window on chromosome 17?

**Read:** Tarantino et al. 2023, *Nat Commun* 14 (HER2-low genomics); IBCGA Research Square preprint, structural variant section and Fig. 3.

- [ ] Ran it  - [ ] Understood the code  - [ ] Read the papers  - [ ] Questions resolved
- [ ] **(publication)** Replace the breakpoint proxy with real structural variant calls (for example from PCAWG whole genomes).

---

## Step 12: Do the ER-low tumours change Steps 4 and 5?

**Run:** `python3 scripts/12_erlow_sensitivity.py`

**You should see:** Step 4 all TNBC rho 0.19 (p 0.08) without ER-low tumours; Step 5 agreement with S22 stays at 0.54 to 0.68.

- [ ] Ran it  - [ ] Understood the code  - [ ] Questions resolved

---

## Run everything at once

`bash run_all.sh` runs Steps 1 to 12 in order and saves each step's output in `logs/`. Use it only once you have done the steps one by one.

---

## Before the talk

- [ ] GitHub account created and repository pushed (21 Nov)
- [ ] Every citation checked against the journal (all *(verify)* items)
- [ ] LLM use statement written (see the analysis logs in the Claude project for what was done with Claude)
- [ ] Slides built and rehearsed to about 20 minutes

## Towards a preprint (after 24 November)

- [ ] Ask Prof. Choudhary whether he would mentor a short preprint
- [ ] All **(publication)** items above
- [ ] A third cohort for the protein result, if one with HER2 protein and TNBC can be found
