# HER-OMI walkthrough: every step, explained against the code

This file explains the whole pipeline in plain language, one script at a time. Read it next to the code.

**How to set up your screen in VS Code**

1. Open this file, then press **Cmd+Shift+V** to see it formatted (the "preview").
2. In the Explorer on the left, right click the script you are reading (for example `scripts/01_define_tnbc.py`) and choose **Open to the Side**. Now the explanation is on one half of the screen and the code on the other.
3. Every code snippet below is copied exactly from your scripts. To jump to it, click in the script, press **Cmd+F**, and paste a few words from the snippet.

**How each step is laid out**

* **Why this step exists**: the scientific question, in one or two sentences.
* **In and out**: which files it reads and which it writes.
* **The code, piece by piece**: snippet, then what it does.
* **What came out**: the numbers and what they mean.
* **If an examiner asks**: the likely questions, with answers you can say out loud.

You do not need to understand every symbol. Aim to understand *what each block is for* and *why it is there*. That is what your examiners will test.

---

## Part 0. The map

```
 00  download raw data from cBioPortal, IBCGA, HGNC, MSigDB
  |
 01  who is TNBC?  (TCGA 160, METABRIC 250)
 01b which signature genes can we measure?  (109 in both cohorts)
  |
 02  score every TNBC tumour with the signature  (singscore)
 02b add PAM50 subtype and marker genes
 02c write expression tables for TNBCtype and step 5
  |
  |---- RNA layer: is the signal real, and what is it?
  |     03  does the score track ERBB2 inside basal tumours, beyond chance?
  |     04  does it rise with HER2 IHC (0, 1+, 2+)?
  |     05  does it carry the same pathway biology IBCGA reported?
  |     06  does it survive proliferation and oestrogen signalling?
  |     07  is it just the LAR subtype?
  |
  |---- PROTEIN layer: is it actually HER2?
  |     08  score vs HER2 protein and phospho HER2 (TCGA RPPA, CPTAC)
  |     08b robustness: ER-low removed, HER2-low only, multiple testing
  |     09  are IBCGA's "ERBB2 pathway" results circular?
  |
  |---- COMPARISON and DNA layer
  |     10  HER2-like TNBC next to real HER2-positive tumours (heatmap, TCGA)
  |     10b same in METABRIC
  |     11  does DNA at the ERBB2 locus explain it?
  |     12  ER-low sensitivity for steps 4 and 5
```

The story the steps tell, in one paragraph: the signature replicates in Western TNBC, and it tracks ERBB2 mRNA even inside basal tumours, beyond what random gene sets do. But at the protein level it follows AR and INPP4B (luminal androgen markers), not HER2 protein or active HER2. Next to real HER2-positive tumours, "HER2-like" TNBC look like basal TNBC on every HER2 measure and like HER2-positive tumours only on AR and FOXA1. So the state is a luminal/androgen programme, not HER2.

---

## Part 1. A five minute Python primer

These are the only ideas you need to read the scripts.

**`import`** loads a toolbox. `pandas` (called `pd`) handles tables. `numpy` (called `np`) handles numbers and arrays. `scipy.stats` (called `stats`) has the statistical tests.

**A DataFrame** is a table, like an Excel sheet. Rows have an *index* (row names, for us usually tumour IDs or gene IDs), columns have names.

| Code | Meaning |
|---|---|
| `pd.read_csv(path, sep="\t")` | read a file into a table. `sep="\t"` means columns are separated by tabs |
| `comment="#"` | skip lines that start with `#` (cBioPortal puts 4 description lines at the top of clinical files) |
| `t.ER_STATUS_BY_IHC` | one column of table `t` |
| `.eq("Negative")` | for every row: is the value equal to "Negative"? Gives a column of True/False |
| `.isin(["A", "B"])` | is the value one of these? True/False per row |
| `a & b` | AND: True only if both are True |
| `a \| b` | OR: True if either is True |
| `~a` | NOT: flips True and False |
| `t[mask]` | keep only the rows where `mask` is True |
| `.sum()` on True/False | counts the Trues (True counts as 1) |
| `.loc[x]` | pick rows by their name. `.loc[rows, cols]` picks rows and columns |
| `.iloc[:, 1:]` | pick by position. `:` means all rows, `1:` means every column from the second one onwards (Python counts from 0) |
| `merge(other, on="PATIENT_ID")` | join two tables row by row where PATIENT_ID matches |
| `.join(other)` | the same, matching on the index |
| `groupby("x")` | split rows into groups by the value of x, then do something per group |
| `.rank()` | replace each value with its rank (smallest = 1) |
| `.T` | transpose: rows become columns |
| `np.log2(x + 1)` | log base 2. The `+ 1` avoids log(0), which is undefined |
| `pd.qcut(x, 3, labels=["low","mid","high"])` | cut into three equal sized groups (tertiles) |
| `f"text {value}"` | an "f-string": puts the value of a variable inside text, for printing |
| `for a, b in [...]:` | loop: do the indented block once for each item |
| `def name(...):` | define a small reusable function |
| `lambda x: ...` | a one-line function with no name |
| `s[13:15]` | characters 14 and 15 of a string (from position 13 up to, not including, 15) |

**Three statistics you will see everywhere**

* **Spearman correlation (rho)**: does y go up when x goes up? It uses ranks, not raw values, so one extreme tumour cannot drive it and it does not assume a straight line. +1 is perfect agreement, 0 is none, −1 is perfectly opposite.
* **Partial Spearman correlation**: the correlation of x and y *after removing* what a third variable z explains in both. If x and y are only related because both follow z, the partial correlation drops to about 0.
* **Empirical p value from random gene sets**: we score 1,000 random gene sets of the same size and up/down mix as the signature, and ask how often a random set does as well as the real one. If only 1 in 1,000 random sets matches it, p is about 0.001. This protects against the trap that *almost any* large gene set correlates with *something* in tumours.

---

## Part 2. The shared toolbox: `scripts/heromi_lib.py`

Every later script imports four functions from here. Understand these once and the rest becomes much easier.

### 2.1 `load_expression`: read an expression matrix

```python
e = pd.read_csv(path, sep="\t")
e = e.dropna(subset=["Entrez_Gene_Id"])
e["Entrez_Gene_Id"] = e["Entrez_Gene_Id"].astype(int)
```

A cBioPortal expression file has one row per gene. The first two columns are `Hugo_Symbol` (the gene name, like ERBB2) and `Entrez_Gene_Id` (a permanent number for the gene, ERBB2 is 2064). Every other column is one tumour. We drop rows with no Entrez ID and make the IDs whole numbers.

**Why Entrez IDs and not names?** Gene names change over time (HGNC renames genes). The number never changes, so matching by number is safer.

```python
cols = [c for c in e.columns[2:] if c in keep_samples]
e = e[["Entrez_Gene_Id"] + cols].dropna()
```

Keep only the tumours we asked for (for example, only TNBC). `.dropna()` then removes any gene that is missing in even one tumour, because singscore needs every gene ranked in every tumour.

```python
e = e.loc[e[cols].mean(axis=1).groupby(e.Entrez_Gene_Id).idxmax()]
```

Some genes appear twice (especially on microarrays, where one gene can have several probes). For each gene we keep the row with the highest average expression, the usual convention because that probe is most likely to be detecting the gene.

### 2.2 `singscore`: the score itself

This is the method from Foroutan et al. 2018. It answers: *in this one tumour, are the signature's up genes near the top of the expression ranking and its down genes near the bottom?*

```python
N = expr.shape[0]
ranks = expr.rank(axis=0, method="average")
```

`N` is the number of genes. `rank(axis=0)` ranks genes *within each tumour separately* (each column on its own). The lowest expressed gene in a tumour gets rank 1, the highest gets rank N. Ties share an averaged rank.

**This is the key idea of the whole project.** Because the score uses only ranks inside one tumour, it does not care whether the numbers came from RNA-seq (TCGA) or a microarray (METABRIC). Those two platforms produce numbers on completely different scales, but "ERBB2 is in the top 5% of genes in this tumour" means the same thing on both.

```python
if reverse:
    r = N + 1 - r
```

For the down genes, flip the ranks, so a gene with very *low* expression gets a *high* rank. Now "signature on" means high ranks for both halves.

```python
lo, hi = (n + 1) / 2, (2 * N - n + 1) / 2
return (r.mean(axis=0) - lo) / (hi - lo) - 0.5
```

Take the mean rank of the n genes in the set. The smallest possible mean rank is `lo` (all n genes at the very bottom) and the largest is `hi` (all at the very top). This line rescales the mean so that `lo` becomes −0.5 and `hi` becomes +0.5.

*A toy example.* 10 genes in a tumour (N = 10). Two up genes, ranked 9 and 10. Mean rank 9.5. lo = (2+1)/2 = 1.5, hi = (20−2+1)/2 = 9.5. Score = (9.5 − 1.5)/(9.5 − 1.5) − 0.5 = **+0.5**, the maximum. If they had been ranked 1 and 2, the score would be −0.5.

```python
out["total"] = out.up + out.down
```

Up half (−0.5 to +0.5) plus down half (−0.5 to +0.5) gives a total from −1 to +1. Higher means more HER2-like. Real tumours sit in a narrow band (our scores have a spread of about 0.08), because no tumour has all 100 genes at the extremes.

### 2.3 `partial_spearman`: correlation after removing a covariate

```python
R = pd.concat([...for v in [x, y] + list(covs)], axis=1).rank().values
Z = np.column_stack([np.ones(len(R)), R[:, 2:]])
res = [R[:, k] - Z @ np.linalg.lstsq(Z, R[:, k], rcond=None)[0] for k in (0, 1)]
r = np.corrcoef(*res)[0, 1]
```

In words:
1. Put x, y and the covariates side by side and rank each column.
2. Build `Z`: a column of 1s (the intercept) plus the covariate ranks.
3. `lstsq` fits a straight line predicting x's ranks from the covariates; subtracting the prediction leaves the *residual*, the part of x the covariates cannot explain. Same for y.
4. Correlate the two residuals. That is the partial Spearman correlation.

The last two lines turn r into a p value with the usual t test, with degrees of freedom reduced by one per covariate. With no covariates the function gives ordinary Spearman.

### 2.4 `random_set_scores`: the null model

```python
rng = np.random.default_rng(seed)
pool = np.array([i for i, g in enumerate(expr.index) if g not in set(exclude_ids)])
...
pick = rng.choice(pool, n_up + n_down, replace=False)
out[i] = norm(R[pick[:n_up]].mean(0), n_up) + norm((N + 1 - R[pick[n_up:]]).mean(0), n_down)
```

* `pool` is every gene *except* the signature genes.
* Each of the 1,000 rounds picks `n_up + n_down` genes at random (no gene twice), calls the first `n_up` of them "up" and the rest "down", and scores them with exactly the singscore formula.
* `seed=607` makes the "random" choices the same every run, so you and anyone else get identical numbers. (607 is just the course number.)

Why it matters: in tumours, thousands of genes move together with proliferation and differentiation. So a random 100-gene score will often correlate with ERBB2 at, say, 0.2 just by riding those big axes. The null tells us what "by chance" looks like *for this exact data*, and we only believe the signature if it beats that.

---

## Step 0. `scripts/00_download_data.sh`: get the data

**Why:** everything must be rebuildable from public sources. This is a shell script (runs in Terminal, not Python).

```bash
set -euo pipefail
B=https://media.githubusercontent.com/media/cBioPortal/datahub/master/public
```

`set -euo pipefail` means "stop immediately if anything fails", so you never continue with a half-downloaded file. `B` is the address of cBioPortal's public data store on GitHub. cBioPortal (Cerami et al. 2012) is a portal that collects cancer genomics studies in one consistent file format.

```bash
curl -sSL -o data/raw/brca_tcga/$f $B/brca_tcga/$f; done
```

`curl` downloads a file. `-o` says where to save it. `-sSL`: quiet, but show errors, and follow redirects.

What it downloads, and what each is for:

| What | Used for |
|---|---|
| `brca_tcga` clinical patient and sample tables | TNBC definition (ER, PR, HER2 IHC, FISH) |
| `brca_tcga` RNA-seq (`data_mrna_seq_v2_rsem.txt`) | scoring |
| `brca_metabric` clinical tables and microarray | TNBC definition and scoring |
| IBCGA supplementary workbook (medRxiv) | Table S21 is the signature; S22, S23 are their pathway results |
| HGNC gene table | matching gene names to Entrez IDs |
| MSigDB Hallmark and Reactome gene sets | pathway analysis (same collections IBCGA used) |
| PanCancer Atlas clinical table | TCGA PAM50 subtype calls |
| TCGA RPPA, TCGA copy-number segments, METABRIC copy number, CPTAC proteome | protein and DNA layers |
| IBCGA follow-up preprint (Research Square) | to check what they did and did not test |

**Why the legacy `brca_tcga` study and not the newer PanCancer one?** Only the legacy study has the IHC receptor fields we need to define TNBC. The PanCancer one has PAM50 calls, so we take PAM50 from there (step 2b).

**If an examiner asks**
* *Why not GDC directly?* cBioPortal gives TCGA and METABRIC in the same format, with clinical fields already cleaned, so one set of code reads both.
* *Is anything not downloadable?* The TNBCtype results (a web tool, step 7) and three small tables made in Session 1; they are in the repo.

---

## Step 1. `scripts/01_define_tnbc.py`: who counts as TNBC?

**Why:** every result depends on which tumours you call triple negative. TNBC means ER negative, PR negative and HER2 negative, but each cohort measured these differently, so the rule must be written out explicitly per cohort.

**In:** clinical tables. **Out:** nothing saved; it creates two tables in memory, `t` (TCGA) and `m` (METABRIC), with a True/False TNBC column. Step 2 runs this file directly.

### TCGA

```python
def read(path):
    return pd.read_csv(path, sep="\t", comment="#")
```

A tiny helper so we do not repeat the same options. `comment="#"` skips cBioPortal's four description lines at the top.

```python
er = t.ER_STATUS_BY_IHC.eq("Negative")
pr = t.PR_STATUS_BY_IHC.eq("Negative")
```

For each patient: is ER negative by immunohistochemistry (IHC)? Is PR? Each gives a True/False column. IHC means a pathologist stained a tissue slice with an antibody and judged whether the receptor protein is there (under ASCO/CAP rules, below 1% of cells stained = negative).

```python
her2_neg = (ihc.eq("Negative")
            | (ihc.isin(["Equivocal", "Indeterminate"]) & fish.eq("Negative"))
            | (ihc.isin(["[Not Available]"]) & fish.eq("Negative")))
```

HER2 is negative if **any one** of three things is true (the `|` signs are OR):
1. IHC says negative (score 0 or 1+), or
2. IHC was equivocal (2+) or indeterminate, *and* FISH says not amplified, or
3. IHC is missing, *and* FISH says not amplified.

This copies the clinical ASCO/CAP logic: an equivocal IHC is settled by FISH (a test that counts copies of the ERBB2 gene in cells).

```python
t["TNBC"] = er & pr & her2_neg
print("TCGA TNBC (IHC, FISH for equivocal or missing IHC):", t.TNBC.sum())   # 160
```

TNBC = all three negative (`&` is AND). `.sum()` counts the Trues: **160 patients**. One has no primary tumour RNA, so 159 are scored.

*The bug story, worth knowing:* an early version forgot branch 3 and got 144. Missing IHC with negative FISH is still HER2 negative, so 160 is right. IBCGA reported 157 for TCGA.

### METABRIC

```python
m = read(f"{RAW}/brca_metabric/data_clinical_patient.txt").merge(
    read(f"{RAW}/brca_metabric/data_clinical_sample.txt"), on="PATIENT_ID")
```

METABRIC splits its clinical data into a patient table and a sample table. The columns we need live in both, so we join them on patient ID.

```python
er_ihc = m.ER_IHC.eq("Negative")                 # note: positives are spelled "Positve" in the source
pr_neg = m.PR_STATUS.eq("Negative")
m["TNBC_cn"] = er_ihc & pr_neg & m.HER2_SNP6.isin(["NEUTRAL", "LOSS"])
```

METABRIC has **no PR or HER2 IHC**. Its `ER_STATUS`, `PR_STATUS` and `HER2_STATUS` columns turn out to be derived from gene expression (we showed this: each separates perfectly on its own gene's mRNA, AUC about 1.0, while ER_IHC only reaches 0.96, the signature of a real lab test). So our rule (option 3, your choice) is:
* ER: the only real IHC test available.
* PR: expression-based status, the only option.
* HER2: not gained on SNP6 copy number, a DNA measurement. This is the closest thing to FISH, which also measures DNA copies.

Result: **250 tumours**, all with array data.

The comment in the code records the rejected rules: using all three status columns gives 320 (and 51 of those are ER positive by IHC); ER IHC plus HER2_STATUS gives 276.

**If an examiner asks**
* *Why not just use the status columns?* Because they are mRNA calls, so you would be defining TNBC with expression and then studying expression in the same tumours: circular. And 51 of the 320 tumours they give are ER positive by the actual IHC test.
* *Your counts differ from IBCGA's (157, 320). Why?* IBCGA took those counts from earlier publications and never gave the rule. Our rules are explicit and match the ASCO/CAP logic as closely as each cohort's data allow.
* *What about ER-low (1 to 10%) tumours?* They are TNBC by IHC, so they are in. We later found TNBCtype flags 12 TCGA and 11 METABRIC as possibly ER positive, and we rerun key results without them (steps 6, 8b, 12).

---

## Step 1b. `scripts/01b_signature.py`: which signature genes can we measure?

**Why:** IBCGA's signature has 144 genes, but not every gene is measured on every platform. A signature with missing genes is a different signature, so we must say exactly which genes we used.

**In:** IBCGA workbook, HGNC table, both expression files. **Out:** `data/IBCGA_HER2like_signature_144genes.csv`, `data/signature_coverage.csv`, `data/ibcga_S22.csv`, `data/ibcga_S23.csv`.

```python
sig = pd.read_excel(XLSX, "Table S21", header=3).iloc[:, 1:]
sig.columns = ["gene", "overlap_count", "categories", "direction", "on_her2_amplicon"]
```

Read the sheet "Table S21". `header=3` because the real column titles are on the 4th row (rows above are a title). `.iloc[:, 1:]` drops an empty first column. Then we give the columns short names. Every gene has a direction (Upregulated or Downregulated in HER2-like tumours) and a flag for whether it sits on the HER2 amplicon (the stretch of chromosome 17 that gets amplified in HER2-positive cancer).

```python
cov = sig.merge(h.rename(columns={"symbol": "gene"}), on="gene", how="left")
```

Attach HGNC information (Entrez ID, gene type) to each signature gene. `how="left"` keeps all 144 genes even if HGNC does not know one.

```python
cov[name] = cov.gene.isin(sym) | cov.entrez_id.fillna(-1).astype(int).isin(ent)
```

A gene counts as measured in a cohort if its **name or** its **Entrez ID** appears in that cohort's expression file. The ID route rescues genes renamed since the platform was annotated. (`fillna(-1)` puts a dummy number where the ID is missing so the line does not crash.)

What came out:
* TCGA measures 114 of 144, METABRIC 112. Almost all protein-coding genes are covered (106 of 107 in each); most missing genes are noncoding RNAs and pseudogenes that neither platform measures.
* **109 genes are measured in both: 42 up, 67 down, 9 on the amplicon.** These 109 are "the signature" for the rest of the project.

The pathway part (`pathway_table`) finds the row that says `Gene_set` and reads everything below it as IBCGA's pathway results, saved for step 5.

**If an examiner asks**
* *Doesn't losing 35 genes change the signature?* Yes, and we say so as a limitation. But the lost genes are mostly noncoding, and the scoring method is rank based, so it degrades gracefully. Using the same 109 in both cohorts keeps the two cohorts comparable.

---

## Step 2. `scripts/02_singscore.py`: score every TNBC tumour

**Why:** turn each tumour's whole expression profile into one number, "how HER2-like is it".

**In:** step 1's tables, coverage table, expression files. **Out:** `results/scores_tnbc.csv` (one row per tumour, three scores, ERBB2 expression).

```python
exec(open("scripts/01_define_tnbc.py").read())
```

Runs step 1 inside this script, so the tables `t` and `m` (with the TNBC flags) exist here.

```python
sig = sig[sig.TCGA & sig.METABRIC].copy()
...
sig["entrez_id"] = sig.entrez_id.fillna(sig.gene.map(sym2ent)).astype(int)
```

Keep the 109 genes measured in both. One gene (FAM136A) has no Entrez ID in HGNC under that name, so `fillna` fills it from the TCGA file, where the name is listed with its ID.

```python
versions = {
    "full":      sig,                                   # 109 genes
    "no_amp":    sig[~sig.on_her2_amplicon],            # 100 genes, the key test
    "amp_only":  sig[sig.on_her2_amplicon],             # 9 genes, all upregulated
}
```

Three versions of the score:
* **full** (109 genes): the signature as IBCGA meant it, minus unmeasured genes.
* **no_amp** (100 genes: 33 up, 67 down): the amplicon genes removed. **This is the one used almost everywhere.** Reason: amplicon genes like ERBB2 and GRB7 sit next to each other on chr17 and go up together whenever that region has extra copies. If they are in the score, the score will trivially correlate with ERBB2. Removing them asks the real question: does the *rest* of the HER2-like programme track ERBB2?
* **amp_only** (9 genes): for comparison.

```python
tcga_ids = {c for c in pd.read_csv(..., nrows=0).columns
            if c[13:15] == "01" and c[:12] in set(t.loc[t.TNBC, "PATIENT_ID"])}
```

This reads only the header row (`nrows=0`) of the TCGA file and keeps tumour columns that pass two tests. A TCGA sample ID looks like `TCGA-A2-A04P-01`:
* `c[:12]` is the first 12 characters, `TCGA-A2-A04P`, the patient. It must be one of our TNBC patients.
* `c[13:15]` is characters 14 and 15, the sample type. `01` means primary solid tumour (as opposed to `11`, normal tissue, or `06`, metastasis).

```python
res[f"score_{name}"] = singscore(expr, gene_ids(s, "Upregulated"), gene_ids(s, "Downregulated")).total
res["ERBB2_expr"] = expr.loc[2064]
```

For each version, score every tumour, and also store ERBB2 expression (Entrez ID 2064) next to it, because almost every later test compares the score with ERBB2.

**What came out:** 159 TCGA and 250 METABRIC TNBC scored. Full and amplicon-free scores correlate 0.99, *which is expected arithmetically* (9 of 109 genes removed) and is not evidence of anything. The amplicon-free score correlates with ERBB2 mRNA at **0.47 (TCGA) and 0.55 (METABRIC)**.

**If an examiner asks**
* *Why singscore and not GSVA or ssGSEA?* Singscore scores each tumour on its own, so it is not affected by which other tumours are in the cohort, and it is rank based, so it works across RNA-seq and microarray. It also handles up and down genes in one score, which this signature needs (57 up, 87 down).
* *Why remove the amplicon genes?* To avoid a circular result: amplicon genes rise with ERBB2 by DNA dosage, not because of a HER2 programme.

---

## Step 2b. `scripts/02b_annotate.py`: add subtype and marker genes

**Why:** later steps ask "within basal tumours" and "is it LAR?", so each tumour needs its PAM50 subtype and the expression of a few luminal/androgen marker genes.

**In:** `results/scores_tnbc.csv`, clinical tables, expression. **Out:** `results/scores_tnbc_annot.csv`, the main table for steps 3 onwards.

```python
s["pid"] = np.where(s.cohort.eq("TCGA"), s["sample"].str[:12], s["sample"])
```

Make a patient ID column. `np.where(condition, a, b)` means "a where the condition is True, else b". TCGA patient = first 12 characters of the sample ID; METABRIC samples are already patient IDs.

```python
s["SUBTYPE"] = s.pid.map(dict(zip(pc.PATIENT_ID, pc.SUBTYPE))).where(s.cohort.eq("TCGA"))
s["pam50"] = (s.SUBTYPE.str.replace("BRCA_", "", regex=False).fillna(s.CLAUDIN_SUBTYPE)
              .replace({"Her2": "HER2E"}))
```

`dict(zip(a, b))` makes a lookup from patient to subtype; `.map` applies it. `.where(...)` keeps the value only for TCGA rows. Then TCGA's labels like `BRCA_Basal` are cleaned to `Basal`, METABRIC's `CLAUDIN_SUBTYPE` fills the METABRIC rows, and "Her2" is renamed "HER2E" (HER2-enriched, the PAM50 class).

```python
MARKERS = {"APOD": 347, "AR": 367, "FOXA1": 3169, "GATA3": 2625, ...}
...
"TCGA": np.log2(load_expression(...) + 1),
```

Marker genes by Entrez ID. AR, FOXA1, GATA3, SPDEF and XBP1 are classic luminal androgen receptor (LAR) markers (Lehmann 2011). TCGA RNA-seq values are put on a log2 scale (raw RSEM values are very skewed; logs make them behave). METABRIC array values are already log2.

**If an examiner asks**
* *Where does PAM50 come from?* The published calls: TCGA PanCancer Atlas `SUBTYPE` and METABRIC `CLAUDIN_SUBTYPE`. We did not recompute PAM50.

---

## Step 2c. `scripts/02c_tnbc_expression_tables.py`: tables for TNBCtype

**Why:** TNBCtype (step 7) is a web tool that needs a file with gene names as rows, tumours as columns, on a log2 scale. Step 5 reuses the same files.

```python
e = e.loc[e[ids].mean(axis=1).groupby(e.Hugo_Symbol).idxmax()].set_index("Hugo_Symbol")
if log:
    e = np.log2(e + 1)
```

Same duplicate rule as before (keep the highest mean row), but by gene *name* this time, because the tool wants names. TCGA gets log2(x+1).

```python
keep = [x for x in e.columns if x not in flag]
e[keep].to_csv(f"{OUT}/{c}_TNBC_expression_noERflag.csv")
```

Two files per cohort: all TNBC, and all TNBC minus the tumours TNBCtype refused as possibly ER positive. TNBCtype accepted the second (TCGA 147, METABRIC 239). Later steps find the flagged tumours as *the difference between the two files*, which is why you will see this pattern again:

```python
set(pd.read_csv(f"{up}/TCGA_TNBC_expression.csv", nrows=0).columns) - \
set(pd.read_csv(f"{up}/TCGA_TNBC_expression_noERflag.csv", nrows=0).columns)
```

(`-` between two sets means "in the first but not the second".)

---

# Part 3. The RNA layer: is the signal real, and what is it?

## Step 3. `scripts/03_basal_only_test.py`: does the score track ERBB2 inside basal tumours?

**Why:** in step 2 the score correlates with ERBB2 at about 0.5. But there is an easy explanation that has nothing to do with HER2: TNBC contains a minority of luminal-ish tumours (LAR, ER-low) that have higher ERBB2 *and* higher scores for luminal reasons. If the correlation is only "luminal tumours versus basal tumours", it would vanish inside PAM50 basal tumours, where LAR is almost absent. So we look inside basal only, with two guards: the random gene set null, and a chromosome 17 check.

**In:** `results/scores_tnbc_annot.csv`, coverage, HGNC, expression. **Out:** `results/basal_only_test.csv`.

```python
sig["chr17"] = sig.location.fillna("").str.startswith("17")
noamp = sig[~sig.on_her2_amplicon]
no17 = sig[~sig.chr17]
```

HGNC gives each gene's chromosome location (like `17q12`). We flag genes on chromosome 17. Two gene lists: amplicon-free (100), and no chromosome 17 at all. Why the second? A broad gain of the whole long arm of chr17 (17q) is common in breast cancer and raises *every* gene on it, including ERBB2, without any focal amplification. If the score's link to ERBB2 came from other chr17 signature genes riding the same DNA gain, removing all chr17 genes would kill it.

```python
for cohort, path in paths.items():
    for subset in ["all TNBC", "Basal only"]:
```

The whole analysis runs four times: two cohorts × (all TNBC, basal only).

```python
lar = g[["AR", "FOXA1", "SPDEF", "XBP1", "GATA3"]].apply(lambda x: (x - x.mean()) / x.std()).mean(axis=1)
```

A simple LAR proxy: z-score each of the five LAR markers (subtract the mean, divide by the standard deviation, so they are on the same scale), then average them per tumour.

```python
pool = expr.index.drop(2064).to_numpy()
...
def fast(up, dn):
    ...
assert np.allclose(fast(ids(noamp, "Upregulated"), ids(noamp, "Downregulated")), obs.values)
```

The null model. `fast` is the same singscore formula, but it ranks once instead of 1,000 times (much quicker). The `assert` line is a safety check: it proves `fast` gives *exactly* the same answer as the real `singscore` for the real signature, and stops the script if not. Note this step's pool excludes only ERBB2; later steps use the library function, which excludes all signature genes. Both are fair nulls; the numbers barely differ.

```python
"empirical_p": (np.sum(np.abs(null) >= abs(rho)) + 1) / (N_RANDOM + 1),
```

Empirical p: how many random sets reach a correlation at least as strong as ours (`abs` makes it two-sided), plus 1, divided by 1,001. The `+1` is standard: it counts the real signature as one of the possible sets, so p can never be exactly 0. The smallest possible p with 1,000 sets is 1/1001 ≈ 0.001.

```python
"partial_rho_ERBB2_given_LAR": partial_spearman(obs, erbb2, lar),
```

Correlation with ERBB2 after removing what the LAR proxy explains. (This script has its own small copy of `partial_spearman`, written before it moved to the library; same maths.)

**What came out** (`results/basal_only_test.csv`):

| | TCGA all | TCGA basal | METABRIC all | METABRIC basal |
|---|---|---|---|---|
| n | 159 | 116 | 250 | 141 |
| score vs ERBB2 | 0.47 | **0.30** | 0.55 | **0.48** |
| random sets, 95th percentile | 0.23 | 0.22 | 0.22 | 0.20 |
| empirical p | 0.001 | **0.019** | 0.001 | **0.001** |
| without any chr17 gene | 0.45 | 0.28 | 0.54 | 0.47 |
| given LAR proxy | 0.40 | 0.32 | 0.45 | 0.46 |

Reading it: inside basal tumours, the score still tracks ERBB2 mRNA, more than random gene sets do, in both cohorts. Removing chr17 genes and controlling for LAR markers barely changes it. So there is a real continuous "HER2-like" gradient tied to ERBB2 mRNA, even among basal tumours. (Whether it reaches HER2 *protein* is the question of step 8.)

**If an examiner asks**
* *Why 1,000 random sets?* It gives a resolution of about 0.001 for p, enough here, and runs in seconds.
* *Why is TCGA weaker?* Smaller n and a different platform; it is a consistent theme (steps 6 and 7 show TCGA's basal gradient is harder to separate from differentiation). We report both cohorts rather than pooling.
* *Is a correlation of 0.30 meaningful?* It is modest. The point is that it beats the 95th percentile of random 100-gene scores (0.22), which already absorbs the "everything correlates with everything" effect.

---

## Step 4. `scripts/04_her2_ihc_check.py`: does the score rise with HER2 IHC?

**Why:** clinically, HER2 is graded by IHC as 0, 1+, 2+ or 3+. "HER2-low" (1+, or 2+ without amplification) is now a treatment category, because trastuzumab deruxtecan works there (DESTINY-Breast04). If the HER2-like state is real HER2 biology, the score should rise from 0 to 1+ to 2+. TCGA only; METABRIC has no HER2 IHC.

```python
t = s[s.cohort.eq("TCGA")].merge(clin, left_on="pid", right_on="PATIENT_ID")
t["ihc"] = pd.to_numeric(t.HER2_IHC_SCORE, errors="coerce")
```

Join scores with the clinical table. `to_numeric(..., errors="coerce")` turns text like "2" into the number 2, and anything unreadable (like "[Not Available]") into missing.

```python
t = t[t.ihc.isin([0, 1, 2])]
```

One TNBC tumour has a recorded IHC of 3+ but an IHC call of "Negative": the record contradicts itself, so it is printed and excluded.

```python
rho, p = stats.spearmanr(g.ihc, g.score_no_amp)
mw = stats.mannwhitneyu(low, zero, alternative="two-sided")
```

Two tests: a trend test (Spearman of IHC score vs score), and a two-group test, HER2-low (1+ or 2+) vs HER2-0, with Mann-Whitney (compares two groups without assuming a bell curve). The line that computes `median_ERBB2_log2` is clumsy but just takes the median of log2(ERBB2 + 1) per group.

**What came out:** 88 TNBC have an IHC score. All TNBC: weak upward trend (rho 0.21, p = 0.047); HER2-low vs HER2-0, p = 0.025 (70 vs 18 tumours). Basal only: no clear trend (rho 0.10, p = 0.43; only 17 HER2-0). ERBB2 mRNA itself does rise across IHC 0 to 2+.

**If an examiner asks**
* *So the score does follow HER2-low?* Weakly, in all TNBC, and not detectably inside basal. This test is underpowered (only 18 HER2-0 tumours), and IHC 0 vs 1+ is known to be poorly reproducible between pathologists. That is why the protein test (step 8) matters more: RPPA measures HER2 protein on a continuous scale.

---

## Step 5. `scripts/05_pathway_comparison.py`: same pathway biology as IBCGA?

**Why:** IBCGA compared HER2-like with basal TNBC across about 1,700 pathways (their Table S22). If our gradient is the same biological state, the pathways that go up and down with our score should match theirs. This is a replication check of the *biology*, not just the score.

```python
ibcga = pd.DataFrame({"HER2like_vs_basal": -s22["diff"], "LAR_vs_basal": s23["diff"], ...
```

IBCGA's S22 "mean difference" is positive when a pathway is higher in *basal* (we inferred this: E2F targets and G2M are positive there, which are basal features). So we flip the sign (`-s22["diff"]`) so that positive means "higher in HER2-like". The preprint does not state the sign, so this is an inference, logged as such.

```python
for line in open(f):
    name, _, *genes = line.rstrip("\n").split("\t")
    sets[name] = [g for g in genes if g not in sig_genes]
```

A GMT file has one pathway per line: name, a web link, then the genes, separated by tabs. `name, _, *genes =` takes the first item as the name, throws away the link (`_`), and puts the rest in `genes`. **Then every signature gene is removed from every pathway.** Without this, a pathway that happens to contain signature genes would agree with the score automatically.

```python
if MIN_SIZE <= n <= MAX_SIZE:
    lo, hi = (n + 1) / 2, (2 * N - n + 1) / 2
    out[name] = (R[idx].mean(0) - lo) / (hi - lo) - 0.5
```

Score each pathway in each tumour with up-only singscore (same formula as the up half in Part 2). Pathways with under 10 genes are too noisy; over 500 are too vague. About 1,350 pass.

```python
profiles[f"{cohort}_{subset}"] = sub.corrwith(meta.loc[mask, "score_no_amp"], method="spearman")
```

For every pathway, its Spearman correlation with our score. The list of ~1,350 correlations is a **gradient profile**: a fingerprint of what goes up and down along the HER2-like axis. Made for all TNBC and for basal only, in each cohort.

Then four comparisons are printed: (1) TCGA profile vs METABRIC profile; (2) our profiles vs IBCGA S22; (3) on the 50 Hallmark sets, HER2-like vs LAR; (4) a few landmark pathways.

**What came out:**

* Our two cohorts agree with each other: rho 0.81 (all TNBC), 0.63 (basal).
* Our profiles agree with IBCGA's HER2-like vs basal biology:

| Profile | vs S22, all | vs S22, FDR < 0.05 | same direction |
|---|---|---|---|
| TCGA all | 0.70 | 0.78 | 79% |
| TCGA basal | 0.54 | 0.66 | 68% |
| METABRIC all | 0.68 | 0.76 | 76% |
| METABRIC basal | 0.61 | 0.69 | 72% |

* What the axis actually is: **less proliferation** (E2F, G2M, MYC targets strongly negative), less interferon signalling, **more oestrogen response and lipid metabolism**. And **Reactome "Signaling by ERBB2" does not track the score** (−0.05 to +0.10).

So: the signature reproduces IBCGA's HER2-like biology in Western TNBC, but that biology looks like a low-proliferation, hormone-response, lipid axis, not HER2 signalling. This was the first hint that "HER2-like" might not mean HER2.

**If an examiner asks**
* *Isn't agreement with S22 guaranteed?* Partly: any luminal vs basal axis would agree somewhat. Agreement shows the axis replicates; it does not make it HER2 specific. That is why the protein test exists.
* *Why remove signature genes from pathways?* To stop pathways agreeing with the score just because they share genes with it. Step 9 shows IBCGA's ERBB2 pathway result fails exactly this check.

---

## Step 6. `scripts/06_confounder_test.py`: does the gradient survive proliferation and oestrogen signalling?

**Why:** step 5 showed the score is anti-proliferative and pro-oestrogen-response. Maybe the score and ERBB2 are only correlated because *both* follow proliferation or luminal differentiation. Test: partial correlation of score vs ERBB2 inside basal, controlling for each of those.

```python
SETS = {
    "prolif11": ["BIRC5", "CCNB1", "CDC20", ...],
    "prolifH": sorted(set(H["HALLMARK_E2F_TARGETS"]) | set(H["HALLMARK_G2M_CHECKPOINT"])),
    "estH": sorted(set(H["HALLMARK_ESTROGEN_RESPONSE_EARLY"]) | set(H["HALLMARK_ESTROGEN_RESPONSE_LATE"])),
}
SETS = {k: [g for g in v if g not in SIG] for k, v in SETS.items()}
```

Three covariate scores: the 11-gene PAM50 proliferation index (Nielsen 2010), a larger Hallmark proliferation set, and Hallmark oestrogen response. Again, signature genes are removed from each. Plus ESR1 mRNA itself and the LAR proxy.

```python
MODELS = {"none": [], "prolif11": ["prolif11"], ..., "prolifH+estH+ESR1+LAR": [...]}
```

Each model is a list of covariates to control for, from none up to all four together.

```python
r, p = partial(bb.score_no_amp, bb.ERBB2_expr, C)
nr = np.array([partial(ns, bb.ERBB2_expr, C)[0] for ns in nulls])
```

For each model: the real partial correlation, and the same partial correlation for each of the 1,000 random gene sets. So the null is also adjusted the same way. `p_emp` here is one-sided (`nr >= r`) because we are asking whether a *positive* link survives.

```python
se = 1 / np.sqrt(out.n - 3 - out.k)
out["ci_lo"], out["ci_hi"] = np.tanh(np.arctanh(out.rho) - 1.96 * se), ...
```

95% confidence intervals with the Fisher z transform: `arctanh` stretches correlations onto a scale where the error is roughly normal with standard error 1/√(n − 3 − k) (k = number of covariates), add ±1.96 standard errors, then `tanh` converts back. The "Cohort difference" block uses the same transform to test whether TCGA and METABRIC differ.

```python
out["retained"] = out.rho / [base[(c, s)] for c, s in zip(out.cohort, out.subset)]
```

"Retained" = adjusted correlation divided by the unadjusted one: what fraction of the association survives.

**What came out (basal only):**

| Model | TCGA rho | kept | METABRIC rho | kept |
|---|---|---|---|---|
| none | 0.30 | 100% | 0.48 | 100% |
| proliferation | 0.25 | 84% | 0.40 | 85% |
| ESR1 mRNA | 0.30 | 101% | 0.47 | 98% |
| LAR proxy | 0.32 | 106% | 0.46 | 97% |
| oestrogen response | 0.20 | 67% | 0.46 | 96% |
| all four | 0.18 (CI crosses 0) | 59% | 0.40 | 84% |

Proliferation explains about 15% in both cohorts. ESR1 and LAR explain nothing. The Hallmark oestrogen response set removes a third of the TCGA link but almost nothing in METABRIC. Since ESR1 itself does not matter in these ER negative tumours, "oestrogen response" here is best read as a luminal differentiation programme. Dropping the ER-flagged tumours changes nothing.

**If an examiner asks**
* *Could adjusting for oestrogen response be over-adjustment?* Yes. Over-adjustment means controlling for something that is part of the effect you are studying, which removes real signal. estH itself correlates with ERBB2 in TCGA basal (+0.28), so controlling for it may remove part of what we want to measure. Correlation cannot tell which way the arrow points.

---

## Step 7. `scripts/07_tnbctype.py`: is the HER2-like signal just the LAR subtype?

**Why:** Lehmann 2011 split TNBC into subtypes; one, LAR (luminal androgen receptor), is luminal-like and AR driven. If "HER2-like" is just LAR renamed, the project's interest collapses. TNBCtype is the official tool that assigns each tumour to Lehmann's subtypes by correlating it with six centroids (average profiles): BL1, BL2, IM, M, MSL, LAR. (UNS = unstable, no clear call.)

```python
def fix_id(s):
    return s.replace(".", "-")
```

The TNBCtype tool is written in R, which changes `-` in column names to `.`. This changes TCGA IDs back so they match ours.

```python
t["subtype"] = t.subtype.fillna("ER-flag")
```

The 23 tumours the tool refused have no call; they are labelled "ER-flag" and kept.

```python
chi = stats.chi2_contingency(np.vstack([ours.values, theirs.values]))
```

Check 1: are our subtype frequencies similar to the ones IBCGA reported for the same cohorts (their Table S15)? A chi-square test compares two sets of counts.

```python
auc = stats.mannwhitneyu(called.score_no_amp[lar], called.score_no_amp[~lar]).statistic / (lar.sum() * (~lar).sum())
```

Check 2: how well does the score pick out LAR? This is the AUC (area under the ROC curve), computed from the Mann-Whitney U statistic: U divided by (number of LAR × number of non-LAR) is exactly the probability that a random LAR tumour scores higher than a random non-LAR tumour. 1.0 = perfect separation, 0.5 = coin flip.

```python
tests = {
    "all TNBC":                          (a.index, []),
    "all TNBC, no LAR, no ER-flag":      (..., []),
    "all called, given LAR coef":        (called.index, ["coef_LAR"]),
    ...
```

Check 3: does the link between score and ERBB2 survive (a) removing LAR and ER-flag tumours, (b) controlling for each tumour's continuous similarity to the LAR centroid (`coef_LAR`)? Each against the 1,000 random sets.

**What came out:**

* Frequencies match IBCGA's (chi-square p = 0.86 TCGA, 0.19 METABRIC).
* **The score is an excellent LAR detector**: every LAR tumour is in the top score third (13 of 13 TCGA, 26 of 26 METABRIC), AUC 0.97 and 0.96.
* **But the ERBB2 link is not LAR**: removing LAR and ER-flag tumours, all TNBC: 0.34 (TCGA), 0.45 (METABRIC), both beating the null. Given the LAR centroid, inside basal: 0.40 and 0.47.
* One discordant result: within Lehmann's BL1 + BL2 only, TCGA shows nothing (−0.04, n = 42) and METABRIC shows 0.45 (n = 71). Post hoc and small, but consistent with TCGA's weaker gradient.

So: the high-scoring tail overlaps LAR heavily (not new), but the coupling to ERBB2 mRNA is separate from LAR.

**If an examiner asks**
* *Why were tumours rejected by TNBCtype?* The tool checks ESR1 expression and refuses tumours that look ER positive. Ours are ER negative by IHC (often ER-low, under 10%), so they are TNBC clinically; we keep them in the main analysis and test without them.
* *Didn't IBCGA already show HER2-like overlaps LAR?* They reported subtype frequencies; they did not test whether the signature's ERBB2 link survives removing LAR, or test it against random gene sets.

---

# Part 4. The protein layer: is it actually HER2?

This is the heart of the project and what is new. IBCGA called HER2-like TNBC "transcriptionally HER2-active" and a candidate group for HER2-directed antibody drug conjugates (ADCs, such as trastuzumab deruxtecan). But an ADC binds HER2 **protein** on the cell surface, and HER2 signalling needs HER2 to be **phosphorylated** (switched on). IBCGA measured neither. We can, in two independent protein datasets.

## Step 8. `scripts/08_protein_test.py`: score vs HER2 protein

**The two protein datasets**
* **TCGA RPPA** (reverse phase protein array): antibodies measure about 200 proteins and phosphoproteins in each tumour. 131 of our TNBC have it. It includes total HER2 and HER2 phosphorylated at tyrosine 1248 (pY1248), the classic "HER2 is active" mark.
* **CPTAC 2020** (Krug et al. 2020): mass spectrometry proteomics and phosphoproteomics, a completely different technology. Only 28 TNBC, but independent.

```python
PROTEINS = {  # RPPA name -> short label
    "ERBB2|HER2": "HER2", "ERBB2|HER2_pY1248": "HER2 pY1248",
    "EGFR|EGFR": "EGFR", ...
    "AR|AR": "AR (control)", "GATA3|GATA3": "GATA3 (control)", "INPP4B|INPP4B": "INPP4B (control)",
```

The proteins tested: HER2 and active HER2; the rest of the HER family (EGFR, HER3) and their active forms; the downstream signalling chain (MEK, ERK, AKT, p70S6K, SRC); and **positive controls**: AR, GATA3, INPP4B and ER-alpha, luminal/LAR proteins that we expect the score to track if it is a luminal programme. A positive control proves the method *can* detect a protein signal, so a null result for HER2 is meaningful rather than a failed experiment.

```python
def emp_p(r, null):
    return (1 + np.sum(np.abs(null) >= abs(r))) / (1 + len(null))
```

Two-sided empirical p, as in step 3.

```python
rppa = pd.read_csv(f"{RAW}/brca_tcga/data_rppa.txt", sep="\t").set_index("Composite.Element.REF").T
```

The RPPA file has proteins as rows and tumours as columns. `.T` flips it so tumours are rows, matching our score table.

```python
r = stats.spearmanr(g.score_no_amp, g[ab])[0]
null = np.array([stats.spearmanr(nul[i], g[ab])[0] for i in nul.columns])
r_mrna = stats.spearmanr(g.ERBB2_log2, g[ab])[0]
pr, _ = partial_spearman(g.score_no_amp, g[ab], [g.ERBB2_log2])
```

For each protein, four numbers:
1. `rho_score`: score vs protein.
2. `null_lo`, `null_hi`: the range 95% of random gene sets fall in.
3. `rho_ERBB2_mRNA`: ERBB2 mRNA vs the protein (a sanity check: HER2 mRNA should track HER2 protein).
4. `partial_score_given_mRNA`: does the score tell us anything about the protein **beyond what ERBB2 mRNA already tells us**? If this is about 0, the score is at best a weaker copy of ERBB2 mRNA for that protein.

```python
her2pos = ref.IHC_HER2.eq("Positive") | ref.HER2_FISH_STATUS.eq("Positive")
...
gap_frac = (hi.median() - lo.median()) / (pos.median() - lo.median())
```

A scale check: put TNBC score thirds next to real HER2-positive tumours and HR+ HER2-negative tumours on the same HER2 protein axis. `gap_frac` asks: of the distance from the low-score TNBC to HER2-positive tumours, how much do the high-score TNBC cover?

The CPTAC block (part B) does the same with a different file format: genes are named by symbol, so the score uses symbols, and the random null is built inline with the same formula. TNBC status is CPTAC's own `TNBC_UPDATED_CLINICAL_STATUS` column.

**What came out** (`results/protein_test.csv`):

| TCGA, all TNBC (n = 131) | score rho | random 95% range | p | ERBB2 mRNA rho | score given mRNA |
|---|---|---|---|---|---|
| HER2 protein | 0.21 | −0.23 to 0.24 | 0.09 | **0.53** | **−0.05** |
| HER2 pY1248 (active) | 0.18 | −0.21 to 0.23 | 0.14 | 0.39 | −0.01 |
| **AR** (control) | **0.50** | −0.32 to 0.31 | 0.001 | 0.16 | 0.49 |
| **INPP4B** (control) | **0.55** | −0.30 to 0.31 | 0.001 | 0.31 | 0.48 |
| AKT pS473 | 0.34 | −0.28 to 0.29 | 0.01 | 0.30 | 0.23 |

CPTAC (n = 28): AR 0.75 (p 0.001), INPP4B 0.58; HER2 0.40 (p 0.06), but given ERBB2 mRNA only 0.15; active HER2 (pY1248) 0.12.

HER2 protein scale (TCGA RPPA medians): HER2-positive 0.88; HR+ HER2-negative 0.14; TNBC top score third −0.13; bottom third −0.43. The top third is statistically above the bottom (p = 0.023) but closes only **23% of the gap** to HER2-positive tumours, and sits below ordinary HR+ HER2-negative tumours.

**How to say it:** the score does not track HER2 protein or active HER2 beyond chance, and it adds nothing beyond ERBB2 mRNA (which itself tracks HER2 protein well, 0.53, so the measurement works). What the score does track, strongly and in both technologies, is AR and INPP4B: luminal androgen proteins.

**If an examiner asks**
* *Absence of evidence is not evidence of absence?* True in general, which is why we have positive controls (AR, INPP4B pass easily with the same method and n) and why we show ERBB2 mRNA *does* track HER2 protein here. The method has the power to see a protein signal; it sees one for AR, not for HER2.
* *RPPA is semi-quantitative.* Yes, so we only use ranks (Spearman) and relative comparisons within the same platform, and CPTAC replicates the pattern with mass spectrometry.
* *Only 28 CPTAC TNBC.* Yes: CPTAC is confirmation of direction, not a standalone test.

---

## Step 8b. `scripts/08b_protein_robustness.py`: does the protein result hold up?

**Why:** three objections to step 8. (1) Maybe the ER-low tumours drive it. (2) The drug question is really about HER2-low tumours: *within* HER2-low, does the score pick tumours with more HER2? (3) Many proteins were tested, so some correlations could be flukes.

```python
subsets = {
    "all TNBC (reference)": np.ones(len(d), bool),
    "all TNBC, ER-flag removed": ~d.er_flag.to_numpy(),
    "basal, ER-flag removed": ...,
    "HER2-low (IHC 1+ or 2+)": d.ihc.isin([1, 2]).to_numpy(),
    "HER2-low, basal": ...,
}
```

Five subsets of tumours, each a True/False mask.

```python
def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; ...
```

The Benjamini-Hochberg correction. When you test many things, some will look significant by luck. BH sorts the p values, multiplies each by (number of tests ÷ its rank), and makes sure the results never decrease. The output, q, is the expected fraction of false discoveries if you call that result significant. Done per subset across the six proteins.

**What came out** (`results/protein_robustness.csv`):

* **ER-flag removed (n = 120):** HER2 protein 0.27 (p 0.03), so it looks better at first. But given ERBB2 mRNA it is **0.00**. AR 0.44 and INPP4B 0.50 stay strong.
* **Within HER2-low (n = 61):** HER2 protein **−0.02** (nothing). AR **0.68**, INPP4B 0.58. So among the very tumours that might get an ADC, the score sorts them by AR, not by HER2.
* After BH correction, AR and INPP4B stay significant in every subset except AR inside basal with ER-flag removed (q 0.17). HER2 protein passes correction only in the ER-flag-removed set, and there its partial given mRNA is 0.00. Active HER2 (pY1248) never passes.

**If an examiner asks**
* *The ER-flag-removed HER2 correlation is significant. Doesn't that support IBCGA?* It shows the score carries some ERBB2 mRNA information (known from step 3), and mRNA carries HER2 protein information. The partial of 0.00 says the score adds nothing on top of ERBB2 mRNA itself, so for choosing patients you would just measure HER2 (IHC) directly.

---

## Step 9. `scripts/09_pathway_overlap_check.py`: is IBCGA's ERBB2 pathway evidence circular?

**Why:** IBCGA's main evidence for "HER2-active" is that HER2-like tumours score high on Reactome ERBB2 pathways (their Fig. 7A). But some of those pathways are tiny and contain ERBB2 and GRB7, which are signature genes and amplicon genes. A gene set that contains ERBB2 will correlate with ERBB2 almost by definition.

```python
for version, genes in [("all genes", gmt[s]), ("amplicon genes removed", [x for x in gmt[s] if x not in AMP]),
                       ("all signature genes removed", [x for x in gmt[s] if x not in SIG])]:
```

Score each ERBB2 pathway three ways: with all its genes, without amplicon genes, without any signature gene. Then correlate with the full signature score (the version IBCGA used) and with ERBB2 mRNA.

**What came out:** "GRB7 events in ERBB2 signalling" has only 4 measured genes in TCGA (5 in METABRIC), two of them ERBB2 and GRB7. With them: correlation with the score 0.25 (TCGA), 0.21 (METABRIC). Without them: **0.10 and 0.05**. The bigger "Signaling by ERBB2" set barely correlates at all (0.17 and 0.01 with everything in). "Amplicon removed" and "all signature removed" give identical numbers because the only signature genes in these sets are the amplicon ones.

**If an examiner asks**
* *Is this a criticism of IBCGA?* It is a methodological point: their pathway result is mostly ERBB2 and GRB7 counting twice, once in the signature and once in the pathway. It is consistent with our protein result: no independent evidence of HER2 pathway activity.

---

# Part 5. Comparison and DNA

## Step 10. `scripts/10_her2pos_comparison.py`: next to real HER2-positive tumours (the heatmap)

**Why:** your professor suggested a heatmap comparing HER2-like TNBC with true HER2-positive tumours. It puts every layer (DNA, RNA, protein) side by side in one picture.

```python
her2pos = (pc.IHC_HER2.eq("Positive") | pc.HER2_FISH_STATUS.eq("Positive")) & \
          pc.ER_STATUS_BY_IHC.eq("Negative") & pc.PR_STATUS_BY_IHC.eq("Negative")
groups = {"HER2-positive\n(HR-negative)": ..., "TNBC\nHER2-like (top third)": ..., "TNBC\nbasal-like (bottom third)": ...}
```

Three groups. The HER2-positive group is also ER and PR negative, so it differs from TNBC **only** in HER2; that makes it the fairest comparison. (`\n` in a label is a line break on the figure.)

```python
seg = seg[seg.chrom.astype(str).eq("17") & (seg.start <= 37_865_000) & (seg.end >= 37_865_000)]
```

The copy number file lists DNA segments with their average copy number (log2 ratio, 0 = normal two copies). Keep the segment that covers position 37,865,000 on chromosome 17, where ERBB2 sits (hg19 coordinates). The underscores in the number are just for readability.

```python
Z = ((X - X.mean()) / X.std()).clip(-2, 2).T
```

Each feature is z-scored across all tumours shown (so different units become comparable: red = above average, blue = below), clipped at ±2 so a few extreme values do not wash out the colours, and transposed so features are rows.

The rest is drawing (matplotlib): the heatmap, white gaps between groups and between layers, labels, colour bar, title.

**What came out** (`results/her2pos_comparison_medians.csv`, medians):

| | HER2-positive | TNBC HER2-like | TNBC basal-like |
|---|---|---|---|
| ERBB2 copy number | **2.14** | −0.05 | −0.21 |
| ERBB2 mRNA | **16.7** | 12.3 | 11.7 |
| HER2 protein | **2.38** | −0.09 | −0.43 |
| HER2 pY1248 | **1.68** | 0.05 | −0.04 |
| AR protein | 0.49 | **0.51** | −0.54 |
| INPP4B protein | 0.04 | **0.32** | −0.59 |
| FOXA1 mRNA | 12.1 | **11.2** | 5.6 |

On every HER2 measure, HER2-like TNBC sit with basal TNBC, far from HER2-positive. On AR and FOXA1 they sit with HER2-positive tumours. The figure title says it: *HER2-like TNBC share HER2-positive tumours' AR/luminal programme, not their HER2.* (HR-negative HER2-positive tumours are often "molecular apocrine", AR-high, which is likely why a signature trained against them picks up AR.)

## Step 10b. `scripts/10b_her2pos_comparison_metabric.py`

Same idea in METABRIC, DNA and RNA only (no protein). The HER2-positive group uses the same tests as our METABRIC TNBC rule, flipped: ERBB2 GAIN on SNP6, ER negative by IHC, PR negative. Same pattern.

---

## Step 11. `scripts/11_dna_layer.py`: does DNA at the ERBB2 locus explain it?

**Why:** if not amplification, maybe a small copy gain, or a structural break near ERBB2 (IBCGA reported an IKZF3-associated rearrangement in HER2-positive tumours), raises ERBB2 in HER2-like TNBC.

```python
dna = pd.DataFrame({
    "erbb2_cn": per.apply(lambda s: s.loc[(s.start <= ERBB2_POS) & (s.end >= ERBB2_POS), "mean"].mean()),
    "bp_window": per.apply(lambda s: int(((s.start > WIN[0]) & (s.start < WIN[1])).sum())),
    "chr17_segments": per.size(),
})
```

For each tumour (`per` is the segments grouped by tumour):
* `erbb2_cn`: copy number of the segment covering ERBB2.
* `bp_window`: how many segments *start* inside chr17:37.6 to 38.1 Mb (CDK12 to IKZF3). A new segment starting means the copy number changes there, a **breakpoint**. TCGA has no rearrangement calls here, so this is a proxy, not proof.
* `chr17_segments`: total segments on chr17, a measure of how scrambled the chromosome is overall (genome instability).

```python
r0 = stats.spearmanr(g.score_no_amp, g.erbb2)[0]
r1 = partial_spearman(g.score_no_amp, g.erbb2, [g[cn_col]])[0]
```

Score vs ERBB2 mRNA, before and after controlling for copy number. If DNA dosage explained the link, r1 would collapse.

```python
partial_rho_given_instability=partial_spearman(gg.bp, gg[var], [gg.chr17_segments])[0]
```

Do tumours with a breakpoint in the window differ, beyond simply having more broken chromosomes?

```python
for w0 in starts:
    has = s17[(s17.start > w0) & (s17.start < w0 + 500_000)].groupby("sid").size()
```

A negative control: repeat the breakpoint test in every 0.5 Mb window along chr17. If many random windows also show an effect, the ERBB2 window is not special.

**What came out** (`results/dna_layer_test.csv`, `results/dna_breakpoint_test.csv`):

* Copy number does track ERBB2 mRNA (0.45 TCGA, 0.20 METABRIC), but controlling for it leaves the link between score and ERBB2 mRNA intact: **0.47 → 0.43** (TCGA), **0.55 → 0.55** (METABRIC); basal similar. So the HER2-like link is transcriptional, not DNA dosage.
* Only **9 of 155** TCGA TNBC have a breakpoint in the window. They have higher ERBB2 mRNA (p < 0.001), higher HER2 protein and a higher score, and this survives the instability control (partial 0.32). Of 94 chr17 windows, only one other (around 25 Mb) gives a stronger ERBB2 mRNA difference. Inside basal (4 tumours) nothing is clear. **Suggestive only, small n**, a lead for future work with real rearrangement data.

**If an examiner asks**
* *Is a segment boundary a rearrangement?* Not necessarily; it is any change in copy number. Real structural variant calls (from whole genome sequencing) would be needed. We present it as a hypothesis.

---

## Step 12. `scripts/12_erlow_sensitivity.py`: ER-low sensitivity for steps 4 and 5

**Why:** Decision 3 kept the ER-low (TNBCtype ER-flagged) tumours in, with a promise to check that results do not depend on them. Steps 6 and 8b already did; this covers steps 4 and 5.

```python
for version, g0 in [("with ER-flag", t), ("ER-flag removed", t[~t.er_flag])]:
```

The step 4 IHC test and the step 5 pathway agreement are rerun twice, with and without the flagged tumours. The code inside is a compact copy of steps 4 and 5.

**What came out:** step 4: trend 0.21 → 0.19, HER2-low vs 0 p 0.03 → 0.04. Step 5: agreement with S22 0.70 → 0.62 (TCGA all), others essentially unchanged (METABRIC 0.68 → 0.68). **Robust.**

---

# Part 6. One page to memorise

**The question:** IBCGA found a "HER2-like" TNBC state in Indian tumours from RNA alone and suggested it as a target for HER2 ADCs. Is it actually HER2?

**The design:** test it in Western TNBC (TCGA 159, METABRIC 250, CPTAC 28) across three layers: RNA, protein, DNA. Every correlation is checked against 1,000 random gene sets.

**The answers:**
1. **It replicates** (RNA): the gradient exists in Western TNBC, matches IBCGA's pathway biology (rho 0.70), and tracks ERBB2 mRNA even in basal tumours (0.30 and 0.48, beating random sets).
2. **It is not just LAR**, but it contains all LAR tumours (AUC 0.97).
3. **It is not HER2 at the protein level**: no link to HER2 protein or active HER2 beyond ERBB2 mRNA (partial about 0), none within HER2-low tumours (−0.02), while AR (0.50 to 0.75) and INPP4B track it strongly in two technologies.
4. **IBCGA's ERBB2 pathway evidence is circular** (0.25 → 0.10 when ERBB2 and GRB7 are removed).
5. **Next to real HER2-positive tumours**, HER2-like TNBC look basal on every HER2 measure and HER2-positive only on AR and FOXA1.
6. **DNA does not explain it**; a rare breakpoint near ERBB2 is a lead, not a finding.

**The conclusion:** "HER2-like" TNBC is a luminal androgen (AR-high, molecular apocrine-like) programme, not a HER2 state. Selecting patients for HER2 ADCs should rely on measuring HER2 protein (IHC), not on this signature.

**Limitations to state yourself, before anyone asks:** Western cohorts only (the IBCGA data are controlled access); 109 of 144 genes measurable; METABRIC TNBC defined partly by mRNA (PR) and copy number (HER2); CPTAC is small; RPPA is semi-quantitative; the breakpoint proxy is not a rearrangement call; everything is correlational.

---

## Running everything

From the project folder, in the VS Code terminal (with your conda environment active):

```bash
bash scripts/00_download_data.sh     # once, about 1 GB
bash run_all.sh                      # every step in order
```

Or run one step at a time, for example `python scripts/03_basal_only_test.py`, and compare what it prints with the numbers in this file. Each step writes its table into `results/`; open any `.csv` there in VS Code or Excel.

---

# Part 7. Added on 23 November: normal tissue (steps 13 and 14)

## Step 13. `scripts/13_normal_comparison.py`: normal breast vs TNBC vs HER2-positive

**Why:** Prof. Choudhary asked for TNBC side by side with normal breast tissue and with HER2-positive tumours. The legacy TCGA matrix used in steps 2 to 12 has no normal samples, so this step uses the TCGA PanCancer release, where tumours and 114 adjacent normal breast samples were processed by the same pipeline (that matters: values from two different pipelines cannot be compared gene by gene).

```python
groups = {
    "Normal breast": list(nor_cols),
    "TNBC basal-like": [s for s in tn.index[tn.third.eq("low")] if s in tum_cols],
    "TNBC HER2-like": [s for s in tn.index[tn.third.eq("high")] if s in tum_cols],
    "HER2+ (HR-negative)": [tum_by_pid[p] for p in clin.index[her2 & hrneg] if p in tum_by_pid],
}
```

Four groups. The TNBC thirds are the same thirds as the step 10 heatmap, so the two figures agree.

```python
pairs = {"TNBC": [p for p in tn.pid.unique() if p in tum_by_pid and f"{p}-11" in nor_cols],
```

The paired view: patients whose own normal tissue (sample code `-11`) was sequenced too. 13 TNBC and 24 HER2-positive patients. Comparing a tumour with the same person's normal tissue removes differences between people.

```python
expr = et.join(en, how="inner")
```

Tumours and normals joined into one matrix, keeping only genes measured in all of them.

```python
F["HER2-like score"] = singscore(expr, up, dn).total
```

The score is recomputed for every sample, normals included. This is fair because singscore ranks genes inside each sample. Check: in our TNBC it agrees with the step 2 score at rho 0.994.

```python
p_wilcoxon=stats.wilcoxon(T[f], N[f]).pvalue
```

For pairs the test is the Wilcoxon signed-rank test: like Mann-Whitney, but for two measurements on the same patient.

**What came out:**
* Normal breast scores 0.07 on the HER2-like score, as high as HER2-positive tumours (0.09), and higher than HER2-like TNBC (−0.01). The signature reads luminal, differentiated tissue.
* ERBB2 mRNA: HER2-like TNBC 12.24, normal breast 12.24, identical. Basal-like TNBC are 0.5 log2 below normal (p 0.001); HER2-positive are 4.6 log2 above (about 24 times).
* AR in HER2-like TNBC (6.8) is close to HER2-positive (8.0), not significantly different after correction (q 0.19).
* Pairs: in TNBC, ERBB2 does not change from own normal to tumour (6 of 13 up, p 0.68), while the score and AR fall in 12 of 13. In HER2-positive patients ERBB2 rises (18 of 24, p 0.001).

**If an examiner asks**
* *Why is normal tissue “HER2-like”?* Normal breast epithelium is luminal and AR/FOXA1-positive; the signature reads that programme. This supports the conclusion that the state is luminal/AR biology, not HER2.
* *Only 13 TNBC pairs?* Yes, that is all TCGA has. The paired result agrees with the group result from 114 normals, which is why it is believable.

## Step 14. `scripts/14_normal_admixture.py`: is it just normal tissue in the sample?

**Why:** step 13 raises an obvious worry: if normal breast scores high, maybe “HER2-like” tumours are simply samples with more normal tissue mixed in. TCGA's tables here have no purity estimate, so this step uses a proxy. Normal breast next to a tumour is mostly fat, and fat-cell genes are almost absent from tumour cells.

```python
ADIPO = {"ADIPOQ": 9370, "PLIN1": 5346, "LEP": 3952, "CIDEA": 1149, "FABP4": 2167, "LPL": 4023, "GPD1": 2819, "CD36": 948}
...
proxy = g.sub(g.mean(axis=1), axis=0).div(g.std(axis=1), axis=0).mean(axis=0)
```

z-score each of the 8 fat-cell genes across tumours, then average per tumour: a higher proxy means more fat (more normal tissue) in the sample.

```python
r1, p1 = partial_spearman(x.score_no_amp, x.ERBB2_expr, [x.adipose])
```

Then the key links are recomputed with the proxy as a covariate.

**What came out:** the score rises only weakly with the proxy (0.22 TCGA, 0.16 METABRIC). Controlling for it changes nothing: score vs ERBB2 0.47 → 0.47 (TCGA), 0.55 → 0.57 (METABRIC); AR protein 0.50 → 0.47; HER2 protein 0.21 → 0.22. Admixture does not explain the results.

**If an examiner asks**
* *Why not a real purity estimate?* None is in the cBioPortal tables used here. ABSOLUTE or ESTIMATE scores would be the better test and are listed as future work.
