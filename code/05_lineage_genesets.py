"""Step A of the registered plan (doi:10.17605/OSF.IO/HQ4A8): lineage gene sets from the
Gut Cell Atlas. Uses the atlas and GENCODE v19 only; reads no GWAS file.

- Cells: Diagnosis in {Healthy adult, Pediatric healthy}, Region != lymph node, not a predicted doublet.
- Pseudobulk: raw counts summed per donor x lineage (atlas field `category`); pairs with < 20 cells dropped; CPM.
- Lineages with < 200 cells or < 5 donors (donors with a kept pseudobulk) after filtering are dropped.
- Genes: autosomal protein-coding GENCODE v19 genes matched by Ensembl ID; expressed = mean CPM >= 1 in >= 1 lineage.
- Specificity(g, L) = mean over donors of CPM in L / sum over lineages of that mean (Skene 2016; Bryois 2020).
- Gene set for L = top 10% of expressed genes by specificity (also 5% and 20% for the registered sensitivity analysis).
- Stability: 100 donor bootstraps (seed 20261008), median Jaccard overlap of the top-10% set with the main set.

Outputs in data/processed/genesets/ and results/simulations/genesets_summary.txt.
Run from the repo root:  .venv/Scripts/python code/05_lineage_genesets.py
"""
import gzip
import hashlib
import pathlib
import re

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ATLAS = pathlib.Path("data/raw/atlas/Full_obj_raw_counts_nosoupx_v2.h5ad")
GTF = pathlib.Path("data/raw/genes/gencode.v19.annotation.gtf.gz")
OUT = pathlib.Path("data/processed/genesets")
RES = pathlib.Path("results/simulations")
KEEP_DX = {"Healthy adult", "Pediatric healthy"}
MIN_CELLS_PAIR, MIN_CELLS_LINEAGE, MIN_DONORS_LINEAGE = 20, 200, 5
TOPS = {"top05": 0.05, "top10": 0.10, "top20": 0.20}
N_BOOT, SEED = 100, 20261008
summary = []


def say(msg):
    print(msg)
    summary.append(msg)


def read_obs():
    with h5py.File(ATLAS, "r") as h:
        obs, cats = h["obs"], h["obs"]["__categories"]
        lab = lambda k: np.array([c.decode() for c in cats[k][:]], dtype=object)[obs[k][:]]
        o = pd.DataFrame({"donor": lab("Sample name"), "dx": lab("Diagnosis"), "region": lab("Region"),
                          "lineage": lab("category"), "doublet": obs["predicted_doublets"][:]})
        genes = pd.DataFrame({"symbol": [g.decode() for g in h["var"]["_index"][:]],
                              "ensembl": [g.decode() for g in h["var"]["gene_ids"][:]]})
    return o, genes


def gencode_genes():
    rows = []
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1).split(".")[0]
            gtype = re.search(r'gene_type "([^"]+)"', f[8]).group(1)
            name = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
            rows.append((gid, name, f[0].removeprefix("chr"), int(f[3]), int(f[4]), gtype))
    g = pd.DataFrame(rows, columns=["ensembl", "gencode_name", "chr", "start", "end", "gene_type"])
    return g[(g["gene_type"] == "protein_coding") & g["chr"].isin([str(i) for i in range(1, 23)])]


def pseudobulk(o, keep):
    """Sum raw counts over cells into (donor, lineage) groups, reading the CSR matrix in row chunks."""
    groups = o.loc[keep, ["donor", "lineage"]].drop_duplicates().sort_values(["lineage", "donor"])
    gidx = {k: i for i, k in enumerate(map(tuple, groups.values))}
    cell_group = np.full(len(o), -1)
    cell_group[keep] = [gidx[(d, l)] for d, l in zip(o.loc[keep, "donor"], o.loc[keep, "lineage"])]
    with h5py.File(ATLAS, "r") as h:
        X = h["X"]
        indptr = X["indptr"][:]
        n_genes = int(X.attrs["shape"][1])
        sums = np.zeros((len(groups), n_genes))
        step = 20000
        for s in range(0, len(o), step):
            e = min(s + step, len(o))
            rows = cell_group[s:e]
            if (rows < 0).all():
                continue
            a, b = indptr[s], indptr[e]
            chunk = sp.csr_matrix((X["data"][a:b], X["indices"][a:b], indptr[s:e + 1] - a), shape=(e - s, n_genes))
            sel = rows >= 0
            G = sp.csr_matrix((np.ones(sel.sum()), (rows[sel], np.flatnonzero(sel))), shape=(len(groups), e - s))
            sums += (G @ chunk).toarray()
    groups = groups.reset_index(drop=True)
    groups["n_cells"] = np.bincount(cell_group[keep], minlength=len(groups))
    return groups, sums


def specificity(groups, cpm, donors=None):
    """Mean over donors of CPM per lineage, then share of each gene's total across lineages.
    `donors` (optional list, may repeat) implements the donor bootstrap."""
    lin = sorted(groups["lineage"].unique())
    means = []
    for L in lin:
        idx = groups.index[groups["lineage"] == L]
        if donors is not None:
            pos = {d: i for d, i in zip(groups.loc[idx, "donor"], idx)}
            idx = [pos[d] for d in donors if d in pos]
        means.append(cpm[idx].mean(axis=0) if len(idx) else np.full(cpm.shape[1], np.nan))
    m = np.vstack(means)
    return lin, m, m / m.sum(axis=0, keepdims=True)


def top_sets(lin, spec, universe, frac):
    n_top = int(round(frac * len(universe)))
    return {L: set(universe[np.argsort(-spec[i, :], kind="stable")[:n_top]]) for i, L in enumerate(lin)}, n_top


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    o, genes = read_obs()
    keep = o["dx"].isin(KEEP_DX).values & (o["region"] != "lymph node").values & ~o["doublet"].values
    say(f"cells kept: {keep.sum():,} of {len(o):,} (healthy postnatal, not lymph node, not doublet); "
        f"donors: {o.loc[keep, 'donor'].nunique()}")
    groups, sums = pseudobulk(o, keep)
    groups = groups[groups["n_cells"] >= MIN_CELLS_PAIR]
    sums = sums[groups.index]
    groups = groups.reset_index(drop=True)
    lin_stats = groups.groupby("lineage").agg(cells=("n_cells", "sum"), donors=("donor", "nunique"))
    ok = lin_stats[(lin_stats["cells"] >= MIN_CELLS_LINEAGE) & (lin_stats["donors"] >= MIN_DONORS_LINEAGE)].index
    for L, r in lin_stats.iterrows():
        say(f"lineage {L}: {r.cells:,} cells in kept pseudobulks, {r.donors} donors -> "
            f"{'kept' if L in ok else 'DROPPED'}")
    sel = groups["lineage"].isin(ok).values
    groups, sums = groups[sel].reset_index(drop=True), sums[sel]
    say(f"K = {len(ok)} lineages retained")
    cpm = sums / sums.sum(axis=1, keepdims=True) * 1e6

    gc = gencode_genes()
    on = genes["ensembl"].isin(gc["ensembl"]).values
    lin, mean_cpm, spec = specificity(groups, cpm)
    expressed = on & (np.nanmax(mean_cpm, axis=0) >= 1)
    universe = genes["ensembl"].values[expressed]
    say(f"autosomal protein-coding GENCODE v19 genes matched in atlas: {on.sum():,}; expressed (mean CPM >= 1 in a "
        f"lineage): {len(universe):,}")

    gc = gc.set_index("ensembl")
    pd.DataFrame({"ensembl": universe}).join(gc, on="ensembl").to_csv(OUT / "control_all_expressed.tsv",
                                                                       sep="\t", index=False)
    main_sets = None
    for tag, frac in TOPS.items():
        sets, n_top = top_sets(lin, spec[:, expressed], universe, frac)
        rows = [(L, g) for L in lin for g in sorted(sets[L])]
        pd.DataFrame(rows, columns=["lineage", "ensembl"]).join(gc, on="ensembl").to_csv(
            OUT / f"lineage_{tag}.tsv", sep="\t", index=False)
        say(f"{tag}: {n_top} genes per lineage")
        if tag == "top10":
            main_sets = sets

    # donor bootstrap stability (registered section 6) and 10 bootstrap gene sets (exploratory 5)
    rng = np.random.default_rng(SEED)
    all_donors = sorted(groups["donor"].unique())
    jac = {L: [] for L in lin}
    boot_rows = []
    for b in range(N_BOOT):
        d = list(rng.choice(all_donors, size=len(all_donors), replace=True))
        _, m_b, s_b = specificity(groups, cpm, donors=d)
        expr_b = on & (np.nanmax(m_b, axis=0) >= 1)
        sets_b, _ = top_sets(lin, np.nan_to_num(s_b[:, expr_b]), genes["ensembl"].values[expr_b], 0.10)
        for L in lin:
            jac[L].append(len(sets_b[L] & main_sets[L]) / len(sets_b[L] | main_sets[L]))
            if b < 10:
                boot_rows += [(b, L, g) for g in sorted(sets_b[L])]
    pd.DataFrame(boot_rows, columns=["boot", "lineage", "ensembl"]).to_csv(OUT / "lineage_top10_boot10.tsv",
                                                                         sep="\t", index=False)
    stab = pd.DataFrame({L: np.percentile(v, [50, 5, 95]) for L, v in jac.items()},
                        index=["median_jaccard", "p05", "p95"]).T
    stab.to_csv(RES / "geneset_bootstrap_stability.tsv", sep="\t", float_format="%.3f")
    say("donor-bootstrap stability of top-10% sets (median Jaccard [5th-95th pct]): " + "; ".join(
        f"{L}: {r.median_jaccard:.2f} [{r.p05:.2f}-{r.p95:.2f}]" for L, r in stab.iterrows()))

    groups.to_csv(RES / "pseudobulk_groups.tsv", sep="\t", index=False)
    digest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.glob("*.tsv"))}
    (OUT / "MANIFEST.sha256").write_text("".join(f"{h}  {n}\n" for n, h in digest.items()))
    (RES / "genesets_summary.txt").write_text("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()
