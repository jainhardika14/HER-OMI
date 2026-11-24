"""Shared helpers for HER-OMI: expression loading and singscore."""
import numpy as np
import pandas as pd


def load_expression(path, keep_samples):
    """Read a cBioPortal matrix, keep chosen samples, index genes by Entrez ID (one row per gene)."""
    e = pd.read_csv(path, sep="\t")
    e = e.dropna(subset=["Entrez_Gene_Id"])
    e["Entrez_Gene_Id"] = e["Entrez_Gene_Id"].astype(int)
    cols = [c for c in e.columns[2:] if c in keep_samples]
    e = e[["Entrez_Gene_Id"] + cols].dropna()             # genes with any missing value are dropped
    # duplicates: keep the row with the highest mean expression (the most informative probe)
    e = e.loc[e[cols].mean(axis=1).groupby(e.Entrez_Gene_Id).idxmax()]
    return e.set_index("Entrez_Gene_Id")[cols]


def singscore(expr, up, down):
    """expr: genes x samples. up, down: lists of gene IDs. Returns a DataFrame of up, down and total scores."""
    N = expr.shape[0]
    ranks = expr.rank(axis=0, method="average")           # step 1, ties share an averaged rank

    def scaled(gene_ids, reverse):
        r = ranks.loc[[g for g in gene_ids if g in ranks.index]]
        if reverse:
            r = N + 1 - r                                  # step 3, flip so low expression = high rank
        n = r.shape[0]
        lo, hi = (n + 1) / 2, (2 * N - n + 1) / 2          # smallest and largest possible mean rank
        return (r.mean(axis=0) - lo) / (hi - lo) - 0.5     # step 4

    out = pd.DataFrame({"up": scaled(up, False)})
    if down:
        out["down"] = scaled(down, True)
        out["total"] = out.up + out.down
    else:
        out["total"] = out.up
    return out


# ---- added in Session 2 ----
def partial_spearman(x, y, covs):
    """Spearman partial correlation of x and y given covariates (rank, regress out, correlate).
    Returns (rho, parametric p). With no covariates this is ordinary Spearman."""
    from scipy import stats
    R = pd.concat([pd.Series(np.asarray(v, dtype=float)) for v in [x, y] + list(covs)], axis=1).rank().values
    Z = np.column_stack([np.ones(len(R)), R[:, 2:]])
    res = [R[:, k] - Z @ np.linalg.lstsq(Z, R[:, k], rcond=None)[0] for k in (0, 1)]
    r = np.corrcoef(*res)[0, 1]
    dof = len(R) - 2 - len(covs)
    return r, 2 * stats.t.sf(abs(r * np.sqrt(dof / (1 - r ** 2))), dof)


def random_set_scores(expr, exclude_ids, n_up, n_down, n_sets=1000, seed=607):
    """Singscores for random gene sets of the given size and direction, drawn from genes not in
    exclude_ids. Returns an array n_sets x samples (columns in expr order)."""
    rng = np.random.default_rng(seed)
    R = expr.rank(axis=0).values
    N = R.shape[0]
    pool = np.array([i for i, g in enumerate(expr.index) if g not in set(exclude_ids)])
    norm = lambda mr, k: (mr - (k + 1) / 2) / ((2 * N - k + 1) / 2 - (k + 1) / 2) - 0.5
    out = np.empty((n_sets, R.shape[1]))
    for i in range(n_sets):
        pick = rng.choice(pool, n_up + n_down, replace=False)
        out[i] = norm(R[pick[:n_up]].mean(0), n_up) + norm((N + 1 - R[pick[n_up:]]).mean(0), n_down)
    return out
