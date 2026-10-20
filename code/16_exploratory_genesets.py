"""Gene sets for registered exploratory analyses 1 and 4 (atlas only; reads no GWAS file).

1. Fine cell types (atlas field `Integrated_05`): same cell filters, pseudobulk and specificity as Step A
   (code/05), keeping fine types with >= 100 cells and >= 3 donors (registered rule; donor-type pairs need >= 20
   cells, as in Step A). Top 10% of expressed genes per fine type.
4. Region-specific lineage sets: Step A repeated on small-intestine cells only (Region = SmallInt) and on
   large-intestine cells only (Region in {LargeInt, REC}; appendix excluded). Lineage rule as Step A (>= 200 cells,
   >= 5 donors).
Writes data/processed/genesets_exploratory/*.tsv and results/main_analysis/exploratory_genesets_summary.txt.
Run from the repo root:  .venv/Scripts/python code/16_exploratory_genesets.py
"""
import importlib.util
import pathlib

import h5py
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("stepA", ROOT / "code/05_lineage_genesets.py")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)
OUT = ROOT / "data/processed/genesets_exploratory"
RES = ROOT / "results/main_analysis"
summary = []


def say(s):
    print(s, flush=True)
    summary.append(s)


def fine_labels():
    with h5py.File(ROOT / A.ATLAS, "r") as h:
        obs, cats = h["obs"], h["obs"]["__categories"]
        return np.array([c.decode() for c in cats["Integrated_05"][:]], dtype=object)[obs["Integrated_05"][:]]


def build(o, keep, min_cells, min_donors, genes, gc, tag):
    groups, sums = A.pseudobulk(o, keep)
    groups = groups[groups["n_cells"] >= A.MIN_CELLS_PAIR]
    sums = sums[groups.index]
    groups = groups.reset_index(drop=True)
    st = groups.groupby("lineage").agg(cells=("n_cells", "sum"), donors=("donor", "nunique"))
    ok = st[(st["cells"] >= min_cells) & (st["donors"] >= min_donors)].index
    sel = groups["lineage"].isin(ok).values
    groups, sums = groups[sel].reset_index(drop=True), sums[sel]
    cpm = sums / sums.sum(axis=1, keepdims=True) * 1e6
    on = genes["ensembl"].isin(gc.index).values
    lin, mean_cpm, spec_ = A.specificity(groups, cpm)
    expressed = on & (np.nanmax(mean_cpm, axis=0) >= 1)
    universe = genes["ensembl"].values[expressed]
    sets, n_top = A.top_sets(lin, spec_[:, expressed], universe, 0.10)
    rows = [(L, g) for L in lin for g in sorted(sets[L])]
    pd.DataFrame(rows, columns=["set", "ensembl"]).join(gc, on="ensembl").to_csv(OUT / f"{tag}.tsv", sep="\t",
                                                                                index=False)
    say(f"{tag}: {len(ok)} of {len(st)} groups retained (>= {min_cells} cells, >= {min_donors} donors); "
        f"{len(universe):,} expressed genes; {n_top} genes per set; dropped: {sorted(set(st.index) - set(ok))}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    o, genes = A.read_obs()
    base_keep = o["dx"].isin(A.KEEP_DX).values & (o["region"] != "lymph node").values & ~o["doublet"].values
    gc = A.gencode_genes().set_index("ensembl")
    # 1. fine cell types
    fo = o.assign(lineage=fine_labels())
    build(fo, base_keep, 100, 3, genes, gc, "fine_types_top10")
    # 4. region-specific lineages
    build(o, base_keep & (o["region"] == "SmallInt").values, A.MIN_CELLS_LINEAGE, A.MIN_DONORS_LINEAGE, genes, gc,
          "small_intestine_lineages_top10")
    build(o, base_keep & o["region"].isin(["LargeInt", "REC"]).values, A.MIN_CELLS_LINEAGE, A.MIN_DONORS_LINEAGE,
          genes, gc, "large_intestine_lineages_top10")
    (RES / "exploratory_genesets_summary.txt").write_text("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()
