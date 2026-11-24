"""HER-OMI step 11: does DNA at the ERBB2 locus explain the HER2-like state in non-amplified TNBC?

Two DNA routes could raise ERBB2 without amplification:
  (a) low-level copy gain or loss of the locus (known to track HER2-low; Tarantino et al. 2023)
  (b) a structural break near ERBB2, e.g. the IKZF3-associated rearrangement IBCGA reported in HER2-positive
      tumours. TCGA has no SV calls here, so the proxy is a copy-number breakpoint (a segment boundary) inside
      chr17:37.6-38.1 Mb (hg19; CDK12 to IKZF3). A breakpoint is not proof of a rearrangement.
Questions, within TNBC:
  1. How much of ERBB2 mRNA does ERBB2 copy number explain?
  2. Does the score's link to ERBB2 mRNA survive controlling for copy number? If yes, the link is transcriptional,
     not DNA dosage. Tested against 1000 matched random gene sets put through the same partial correlation.
  3. Do tumours with a breakpoint in the ERBB2 to IKZF3 window differ in ERBB2 mRNA, HER2 protein or the score,
     beyond what general genome instability (number of chr17 segments) explains?
TCGA: continuous segment means (hg19). METABRIC: discrete GISTIC calls (-2 to 2), no segments available here.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, "scripts")
from heromi_lib import load_expression, partial_spearman, random_set_scores

RAW = "data/raw"
ERBB2_POS, WIN = 37_865_000, (37_600_000, 38_100_000)
sig = pd.read_csv("data/signature_coverage.csv")
core = sig[sig.TCGA & sig.METABRIC & ~sig.on_her2_amplicon]
n_up, n_dn = int(core.direction.eq("Upregulated").sum()), int(core.direction.eq("Downregulated").sum())
sig_ids = set(sig.loc[sig.entrez_id.notna(), "entrez_id"].astype(int))
ann = pd.read_csv("results/scores_tnbc_annot.csv").set_index("sample")
rows = []


def test(cohort, subset, g, nul, cn_col):
    """Score vs ERBB2 mRNA, raw and given copy number, each against the random-set null."""
    r0 = stats.spearmanr(g.score_no_amp, g.erbb2)[0]
    r1 = partial_spearman(g.score_no_amp, g.erbb2, [g[cn_col]])[0]
    n0 = np.array([stats.spearmanr(nul[i], g.erbb2)[0] for i in nul.columns])
    n1 = np.array([partial_spearman(nul[i], g.erbb2, [g[cn_col]])[0] for i in nul.columns])
    rows.append(dict(cohort=cohort, subset=subset, n=len(g),
                     cn_vs_ERBB2_mRNA=stats.spearmanr(g[cn_col], g.erbb2)[0],
                     cn_vs_score=stats.spearmanr(g[cn_col], g.score_no_amp)[0],
                     score_vs_mRNA=r0, p_emp=(1 + np.sum(n0 >= r0)) / 1001,
                     score_vs_mRNA_given_cn=r1, p_emp_given_cn=(1 + np.sum(n1 >= r1)) / 1001))


# ---------------- TCGA ----------------
tc = ann[ann.cohort.eq("TCGA")]
seg = pd.read_csv(f"{RAW}/brca_tcga/data_cna_hg19.seg", sep="\t")
seg.columns = ["ID", "chrom", "start", "end", "nmark", "mean"]
seg["sid"] = seg.ID.str[:15]
seg = seg[seg.chrom.astype(str).eq("17")]
per = seg.groupby("sid")
dna = pd.DataFrame({
    "erbb2_cn": per.apply(lambda s: s.loc[(s.start <= ERBB2_POS) & (s.end >= ERBB2_POS), "mean"].mean()),
    "bp_window": per.apply(lambda s: int(((s.start > WIN[0]) & (s.start < WIN[1])).sum())),
    "chr17_segments": per.size(),
})
tc = tc.assign(sid=tc.index.str[:15]).join(dna, on="sid").dropna(subset=["erbb2_cn"])
expr = load_expression(f"{RAW}/brca_tcga/data_mrna_seq_v2_rsem.txt", set(tc.index))
tc = tc.loc[expr.columns]
tc["erbb2"] = np.log2(tc.ERBB2_expr + 1)
nul = pd.DataFrame(random_set_scores(expr, sig_ids, n_up, n_dn).T, index=expr.columns)
test("TCGA", "all TNBC", tc, nul, "erbb2_cn")
b = tc.pam50.eq("Basal")
test("TCGA", "basal", tc[b], nul[b], "erbb2_cn")

# breakpoint proxy
rppa = pd.read_csv(f"{RAW}/brca_tcga/data_rppa.txt", sep="\t").set_index("Composite.Element.REF").T
tc["HER2_protein"] = tc.index.map(lambda s: float(rppa.at[s, "ERBB2|HER2"]) if s in rppa.index else np.nan)
tc["bp"] = tc.bp_window.gt(0).astype(int)
bp_rows = []
for subset, g in [("all TNBC", tc), ("basal", tc[tc.pam50.eq("Basal")])]:
    for var in ["erbb2", "HER2_protein", "score_no_amp", "erbb2_cn"]:
        gg = g.dropna(subset=[var])
        yes, no = gg.loc[gg.bp.eq(1), var], gg.loc[gg.bp.eq(0), var]
        bp_rows.append(dict(subset=subset, variable=var, n_bp=len(yes), n_no=len(no),
                            median_bp=yes.median(), median_no=no.median(),
                            p_mannwhitney=stats.mannwhitneyu(yes, no).pvalue if len(yes) > 1 else np.nan,
                            partial_rho_given_instability=partial_spearman(gg.bp, gg[var], [gg.chr17_segments])[0]))
bp = pd.DataFrame(bp_rows)
# is the ERBB2 window special? the same breakpoint test in every other 0.5 Mb window on chr17
starts = np.arange(0, 81_000_000, 500_000)
sid_ok = set(tc.sid)
s17 = seg[seg.sid.isin(sid_ok)]
win_p = []
for w0 in starts:
    has = s17[(s17.start > w0) & (s17.start < w0 + 500_000)].groupby("sid").size()
    flag = tc.sid.map(has).fillna(0).gt(0)
    if 5 <= flag.sum() <= len(tc) - 5:
        win_p.append((w0, flag.sum(), stats.mannwhitneyu(tc.erbb2[flag], tc.erbb2[~flag]).pvalue))
win_p = pd.DataFrame(win_p, columns=["window_start", "n_bp", "p_erbb2_mRNA"])
erbb2_win_p = bp.loc[(bp.subset == "all TNBC") & (bp.variable == "erbb2"), "p_mannwhitney"].iloc[0]
rank = (win_p.p_erbb2_mRNA <= erbb2_win_p).sum()

# ---------------- METABRIC ----------------
mb = ann[ann.cohort.eq("METABRIC")]
cna = pd.read_csv(f"{RAW}/brca_metabric/data_cna.txt", sep="\t")
cna = cna[cna.Hugo_Symbol.eq("ERBB2")].drop(columns=["Hugo_Symbol", "Entrez_Gene_Id"]).iloc[0]
mb = mb.assign(erbb2_cn=cna.reindex(mb.index)).dropna(subset=["erbb2_cn"])
mexpr = load_expression(f"{RAW}/brca_metabric/data_mrna_illumina_microarray.txt", set(mb.index))
mb = mb.loc[mexpr.columns]
mb["erbb2"] = mb.ERBB2_expr
mnul = pd.DataFrame(random_set_scores(mexpr, sig_ids, n_up, n_dn).T, index=mexpr.columns)
test("METABRIC", "all TNBC", mb, mnul, "erbb2_cn")
b = mb.pam50.eq("Basal")
test("METABRIC", "basal", mb[b], mnul[b], "erbb2_cn")

res = pd.DataFrame(rows)
res.to_csv("results/dna_layer_test.csv", index=False)
bp.to_csv("results/dna_breakpoint_test.csv", index=False)
pd.set_option("display.width", 200)
print("METABRIC ERBB2 GISTIC calls in our TNBC:", mb.erbb2_cn.value_counts().sort_index().to_dict())
print(f"TCGA ERBB2 log2 copy number in TNBC: median {tc.erbb2_cn.median():.2f}; "
      f"gain > 0.2: {(tc.erbb2_cn > 0.2).sum()}, loss < -0.2: {(tc.erbb2_cn < -0.2).sum()} of {len(tc)}")
print("\n1-2) Copy number vs ERBB2 mRNA and the score; score-ERBB2 link with and without copy number")
print(res.round(3).to_string(index=False))
print("\n3) Breakpoint in chr17:37.6-38.1 Mb (CDK12 to IKZF3)")
print(bp.round(3).to_string(index=False))
print(f"   ERBB2 window vs {len(win_p)} other 0.5 Mb chr17 windows (>=5 tumours with a breakpoint): "
      f"{rank} windows have an ERBB2 mRNA p value at or below the ERBB2 window's ({erbb2_win_p:.4f})")
print(win_p.nsmallest(5, "p_erbb2_mRNA").assign(window_Mb=lambda x: x.window_start / 1e6).round(4).to_string(index=False))
