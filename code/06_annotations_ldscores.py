"""Step B of the registered plan: SNP annotations from gene sets and their LD scores, per ancestry.
Uses gene sets (code/05), GENCODE v19 coordinates and 1000 Genomes reference genotypes only; reads no GWAS file.

Annotation columns (one LDSC prefix per column, so each lineage can enter its own model):
  primary      L_<lineage>_top10   genes in the lineage's top-10% set, +/-100 kb
               control_100kb        all expressed genes, +/-100 kb
  sensitivity  L_<lineage>_top05, L_<lineage>_top20 (+/-100 kb); L_<lineage>_top10_50kb, control_50kb (+/-50 kb)
  phase 4      R<i>_random          8 random sets of the same size as the top-10% sets (seed 20261008), used as
                                    shuffled-label controls in simulations
LD scores: ldsc --l2 --ld-wind-cm 1 on the ancestry's own 1000G PLINK files, printed for exactly the SNPs in that
ancestry's baseline v1.2 LD-score files (so all reference LD-score files share identical rows).

Run from the repo root:  .venv/Scripts/python code/06_annotations_ldscores.py [--chroms 22] [--workers 4]
"""
import argparse
import gzip
import subprocess
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from common import ANC, ANNOT, CHROMS, GENESETS, LDSC_DIR, LDSC_PY, LINEAGES, slug

SEED = 20261008


def gene_sets():
    """Return {column name: (window_bp, DataFrame of genes with chr/start/end)}."""
    ctrl = pd.read_csv(GENESETS / "control_all_expressed.tsv", sep="\t", dtype={"chr": str})
    sets = {"control_100kb": (100_000, ctrl), "control_50kb": (50_000, ctrl)}
    for tag in ["top05", "top10", "top20"]:
        t = pd.read_csv(GENESETS / f"lineage_{tag}.tsv", sep="\t", dtype={"chr": str})
        for L in LINEAGES:
            sets[f"L_{slug(L)}_{tag}"] = (100_000, t[t["lineage"] == L])
            if tag == "top10":
                sets[f"L_{slug(L)}_top10_50kb"] = (50_000, t[t["lineage"] == L])
    n_top = len(pd.read_csv(GENESETS / "lineage_top10.tsv", sep="\t").query("lineage == @LINEAGES[0]"))
    rng = np.random.default_rng(SEED)
    for i in range(1, 9):
        sets[f"R{i}_random"] = (100_000, ctrl.iloc[np.sort(rng.choice(len(ctrl), n_top, replace=False))])
    return sets


def in_windows(bp, genes, window):
    """Boolean: SNP position within [start - window, end + window] of any gene (positions on one chromosome)."""
    if genes.empty:
        return np.zeros(len(bp), dtype=bool)
    iv = np.sort(np.c_[genes["start"].values - window, genes["end"].values + window], axis=0)
    starts, ends = iv[:, 0], np.maximum.accumulate(iv[:, 1])  # merged-interval lookup
    k = np.searchsorted(starts, bp, side="right") - 1
    return (k >= 0) & (bp <= ends[np.clip(k, 0, None)])


def run_chrom(args):
    anc, c, sets = args
    out = ANNOT / anc
    bim = pd.read_csv(f"{ANC[anc]['bfile']}".format(chr=c) + ".bim", sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    cols = {}
    for name, (w, g) in sets.items():
        cols[name] = in_windows(bim["BP"].values, g[g["chr"] == str(c)], w).astype(np.int8)
    ann = pd.concat([bim[["CHR", "BP", "SNP", "CM"]], pd.DataFrame(cols)], axis=1)
    work = out / "_combined"
    work.mkdir(parents=True, exist_ok=True)
    ann.to_csv(work / f"all.{c}.annot.gz", sep="\t", index=False, compression="gzip")
    base = pd.read_csv(f"{ANC[anc]['baseline']}".format(chr=c) + ".l2.ldscore.gz", sep="\t", usecols=["SNP"])
    base["SNP"].to_csv(work / f"printsnps.{c}.txt", index=False, header=False)
    log = subprocess.run([str(LDSC_PY), str(LDSC_DIR / "ldsc.py"), "--l2", "--bfile",
                          f"{ANC[anc]['bfile']}".format(chr=c), "--ld-wind-cm", "1",
                          "--annot", str(work / f"all.{c}.annot.gz"),
                          "--print-snps", str(work / f"printsnps.{c}.txt"), "--out", str(work / f"all.{c}")],
                         capture_output=True, text=True)
    if log.returncode != 0:
        raise RuntimeError(f"ldsc --l2 failed for {anc} chr{c}:\n{log.stdout[-2000:]}\n{log.stderr[-2000:]}")
    # split the combined output into one LDSC prefix per annotation column
    ld = pd.read_csv(work / f"all.{c}.l2.ldscore.gz", sep="\t")
    if not ld["SNP"].equals(base["SNP"]):
        raise RuntimeError(f"{anc} chr{c}: printed SNPs differ from baseline rows")
    names = list(sets)
    M = np.loadtxt(work / f"all.{c}.l2.M", ndmin=1)
    M550 = np.loadtxt(work / f"all.{c}.l2.M_5_50", ndmin=1)
    for i, name in enumerate(names):
        d = out / name
        d.mkdir(exist_ok=True)
        ann[["CHR", "BP", "SNP", "CM", name]].to_csv(d / f"{name}.{c}.annot.gz", sep="\t", index=False,
                                                     compression="gzip")
        ld[["CHR", "SNP", "BP", f"{name}L2"]].to_csv(d / f"{name}.{c}.l2.ldscore.gz", sep="\t", index=False,
                                                    compression="gzip")
        (d / f"{name}.{c}.l2.M").write_text(f"{M[i]:g}\n")
        (d / f"{name}.{c}.l2.M_5_50").write_text(f"{M550[i]:g}\n")
    return anc, c, len(bim), len(ld), {n: int(cols[n].sum()) for n in ("control_100kb", f"L_{slug(LINEAGES[0])}_top10")}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--chroms", type=int, nargs="*", default=list(CHROMS))
    p.add_argument("--ancestries", nargs="*", default=["EAS", "EUR"])
    p.add_argument("--workers", type=int, default=4)
    a = p.parse_args()
    sets = gene_sets()
    jobs = [(anc, c, sets) for c in sorted(a.chroms, key=lambda c: c) for anc in a.ancestries]
    with ProcessPoolExecutor(a.workers) as ex:
        for anc, c, n_ref, n_print, counts in ex.map(run_chrom, jobs):
            print(f"{anc} chr{c}: {n_ref:,} reference SNPs, {n_print:,} printed; annotated SNPs {counts}", flush=True)


if __name__ == "__main__":
    main()
