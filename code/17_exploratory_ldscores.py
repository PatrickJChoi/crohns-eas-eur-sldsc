"""LD scores for the registered exploratory gene sets (reads no GWAS file).

Sets: fine cell types and small/large-intestine lineages (code/16), and the 10 donor-bootstrap top-10% lineage sets
(code/05). Same windows (+/-100 kb, GENCODE v19), reference panels and ldsc --l2 settings as Step B (code/06).
Exploratory models are fitted without --overlap-annot (verified in the simulations to leave the coefficients and jackknife
unchanged), so only LD-score, M and M_5_50 files are written per set.
Run from the repo root:  .venv/Scripts/python code/17_exploratory_ldscores.py [--workers 3]
"""
import argparse
import importlib.util
import pathlib
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
import sys  # noqa: E402

sys.path.insert(0, str(ROOT / "code"))
from common import ANC, ANNOT, CHROMS, GENESETS, LDSC_DIR, LDSC_PY  # noqa: E402

spec = importlib.util.spec_from_file_location("stepB", ROOT / "code/06_annotations_ldscores.py")
B = importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)
EXP = ROOT / "data/processed/genesets_exploratory"
OUT_ROOT = ROOT / "data/processed/annot_exploratory"


def clean(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")


def exploratory_sets():
    sets = {}
    for tag, short in [("fine_types_top10", "F"), ("small_intestine_lineages_top10", "SI"),
                       ("large_intestine_lineages_top10", "LI")]:
        t = pd.read_csv(EXP / f"{tag}.tsv", sep="\t", dtype={"chr": str})
        for s, g in t.groupby("set"):
            sets[f"{short}_{clean(s)}"] = (100_000, g)
    boot = pd.read_csv(GENESETS / "lineage_top10_boot10.tsv", sep="\t")
    coords = pd.read_csv(GENESETS / "control_all_expressed.tsv", sep="\t", dtype={"chr": str})
    boot = boot.merge(coords, on="ensembl")
    for (b, L), g in boot.groupby(["boot", "lineage"]):
        sets[f"BOOT{b}_{clean(L)}"] = (100_000, g)
    return sets


def run_chrom(args):
    anc, c, sets = args
    out = OUT_ROOT / anc
    work = out / "_combined"
    work.mkdir(parents=True, exist_ok=True)
    bim = pd.read_csv(f"{ANC[anc]['bfile']}".format(chr=c) + ".bim", sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    cols = {n: B.in_windows(bim["BP"].values, g[g["chr"] == str(c)], w).astype(np.int8) for n, (w, g) in sets.items()}
    ann = pd.concat([bim[["CHR", "BP", "SNP", "CM"]], pd.DataFrame(cols)], axis=1)
    ann.to_csv(work / f"all.{c}.annot.gz", sep="\t", index=False, compression={"method": "gzip", "compresslevel": 1})
    base = pd.read_csv(f"{ANC[anc]['baseline']}".format(chr=c) + ".l2.ldscore.gz", sep="\t", usecols=["SNP"])
    base["SNP"].to_csv(work / f"printsnps.{c}.txt", index=False, header=False)
    r = subprocess.run([str(LDSC_PY), str(LDSC_DIR / "ldsc.py"), "--l2", "--bfile",
                        f"{ANC[anc]['bfile']}".format(chr=c), "--ld-wind-cm", "1", "--annot",
                        str(work / f"all.{c}.annot.gz"), "--print-snps", str(work / f"printsnps.{c}.txt"),
                        "--out", str(work / f"all.{c}")], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"ldsc --l2 failed {anc} chr{c}: {r.stdout[-1500:]} {r.stderr[-1500:]}")
    ld = pd.read_csv(work / f"all.{c}.l2.ldscore.gz", sep="\t")
    if not ld["SNP"].equals(base["SNP"]):
        raise RuntimeError(f"{anc} chr{c}: printed SNPs differ from baseline rows")
    M = np.loadtxt(work / f"all.{c}.l2.M", ndmin=1)
    M550 = np.loadtxt(work / f"all.{c}.l2.M_5_50", ndmin=1)
    for i, name in enumerate(sets):
        d = out / name
        d.mkdir(exist_ok=True)
        ld[["CHR", "SNP", "BP", f"{name}L2"]].to_csv(d / f"{name}.{c}.l2.ldscore.gz", sep="\t", index=False,
                                                    compression={"method": "gzip", "compresslevel": 1})
        (d / f"{name}.{c}.l2.M").write_text(f"{M[i]:g}\n")
        (d / f"{name}.{c}.l2.M_5_50").write_text(f"{M550[i]:g}\n")
    (work / f"all.{c}.annot.gz").unlink()  # large; the per-set LD scores are what the models use
    return anc, c


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--workers", type=int, default=3)
    a = p.parse_args()
    sets = exploratory_sets()
    print(f"{len(sets)} exploratory gene sets", flush=True)
    jobs = [(anc, c, sets) for c in CHROMS for anc in ["EAS", "EUR"]]
    with ProcessPoolExecutor(a.workers) as ex:
        for anc, c in ex.map(run_chrom, jobs):
            print(f"{anc} chr{c} done", flush=True)


if __name__ == "__main__":
    main()
