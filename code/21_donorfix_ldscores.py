"""Donor-corrected Step B (post hoc correction, reported alongside the registered run).
Same code as the registered Step B (code/06.run_chrom), applied to the 8 donor-corrected top-10% lineage sets and the
donor-corrected all-expressed-genes control, written to data/processed/annot_donorfix/ (registered annotations
in data/processed/annot/ are not touched: each worker process sets the output root itself, because Windows worker
processes re-import modules and would not see an override made in the parent).
Run from the repo root:  .venv/Scripts/python code/21_donorfix_ldscores.py [--workers 3]
"""
import argparse
import importlib.util
import pathlib
from concurrent.futures import ProcessPoolExecutor

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("stepB", ROOT / "code/06_annotations_ldscores.py")
B = importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)
FIX_ANNOT = ROOT / "data/processed/annot_donorfix"
FIX_GENESETS = ROOT / "data/processed/genesets_donorfix"


def primary_sets():
    B.GENESETS = FIX_GENESETS
    sets = B.gene_sets()
    keep = [f"L_{B.slug(L)}_top10" for L in B.LINEAGES] + ["control_100kb"]
    return {k: sets[k] for k in keep}


def job(args):
    B.ANNOT = FIX_ANNOT  # set inside the worker process
    if B.ANNOT != FIX_ANNOT:
        raise RuntimeError("output root not set")
    return B.run_chrom(args)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--workers", type=int, default=3)
    a = p.parse_args()
    sets = primary_sets()
    print(f"{len(sets)} donor-corrected sets: {list(sets)}", flush=True)
    jobs = [(anc, c, sets) for c in B.CHROMS for anc in ["EAS", "EUR"]]
    with ProcessPoolExecutor(a.workers) as ex:
        for anc, c, n_ref, n_print, counts in ex.map(job, jobs):
            print(f"{anc} chr{c}: {n_ref:,} reference SNPs, {n_print:,} printed", flush=True)


if __name__ == "__main__":
    main()
