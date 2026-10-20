"""Donor-corrected Steps D-E (post hoc correction, 2026-09-09; reported ALONGSIDE the
registered result, never replacing it).

Same models and tests as the registered primary run (code/13): baseline v1.2 + control + one lineage per model,
--overlap-annot, primary Wald test EAS vs EUR (de Lange) and validations (i) de Lange vs FinnGen, (ii) EAS vs FinnGen,
on the same prepared summary statistics (data/processed/real/; Step C does not depend on gene sets). Only the gene
sets, control and LD scores differ: donor-corrected (code/20, code/21; data/processed/annot_donorfix/).
Writes results/main_analysis/donorfix/primary_results.json and primary_table.tsv.
Run from the repo root:  .venv/Scripts/python code/24_donorfix_primary.py
"""
import importlib.util
import json
import pathlib
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from common import LINEAGES, slug  # noqa: E402


def load(name):
    spec = importlib.util.spec_from_file_location(name[3:-3], ROOT / "code" / name)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


fit, wtest = load("08_fit_models.py"), load("09_primary_test.py")
FIX = ROOT / "data/processed/annot_donorfix"
REAL = ROOT / "data/processed/real"
OUT = ROOT / "results/main_analysis/donorfix"
DATASETS = {"EAS": "EAS", "EUR_deLange": "EUR", "EUR_FinnGen": "EUR"}
NAMES = [f"L_{slug(L)}_top10" for L in LINEAGES]


def fit_fixed(ds, anc, workers=3):
    fit.ANNOT = FIX  # control annotation also comes from the donor-corrected folder (same process; threads)
    outdir = OUT / f"fits_{ds}"
    outdir.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(workers) as ex:
        fits = list(ex.map(lambda n: fit.run_one(REAL / f"{ds}.sumstats.gz", anc, n, "control_100kb", outdir,
                                                 overlap=True, annot_root=FIX), NAMES))
    for f in fits:  # guard: every model used donor-corrected annotation files
        log = (outdir / f"{anc}.{f['name']}.log").read_text()
        if "annot_donorfix" not in log or "annot\\EAS" in log.replace("annot_donorfix", "") \
                or "annot\\EUR" in log.replace("annot_donorfix", ""):
            raise SystemExit(f"{ds} {f['name']}: model did not use only donor-corrected annotations")
    return {"names": NAMES, "tau_star": [f["tau_star"] for f in fits],
            "delete": np.array([f["delete"] for f in fits]).T.tolist(), "fits": fits}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fits = {ds: fit_fixed(ds, anc) for ds, anc in DATASETS.items()}
    res = {"note": "post hoc donor-corrected rerun; the registered result is results/main_analysis/primary_results.json",
           "primary_EAS_vs_EUR": wtest.wald(fits["EAS"], fits["EUR_deLange"]),
           "validation_i_deLange_vs_FinnGen": wtest.wald(fits["EUR_deLange"], fits["EUR_FinnGen"]),
           "validation_ii_EAS_vs_FinnGen": wtest.wald(fits["EAS"], fits["EUR_FinnGen"]),
           "tau_star": {k: v["tau_star"] for k, v in fits.items()},
           "fits": {k: [{kk: vv for kk, vv in f.items() if kk != "delete"} for f in v["fits"]] for k, v in fits.items()}}
    (OUT / "primary_results.json").write_text(json.dumps(res, indent=1))
    p = res["primary_EAS_vs_EUR"]
    pd.DataFrame({"lineage": LINEAGES, "tau_star_EAS": fits["EAS"]["tau_star"], "se_EAS": p["se_a"],
                  "tau_star_EUR": fits["EUR_deLange"]["tau_star"], "se_EUR": p["se_b"],
                  "tau_star_FinnGen": fits["EUR_FinnGen"]["tau_star"], "D_EAS_minus_EUR": p["D"], "se_D": p["se_D"],
                  "z": p["z"], "p": p["p_lineage"], "p_holm": p["p_holm"]}).to_csv(
        OUT / "primary_table.tsv", sep="\t", index=False, float_format="%.4g")
    print(json.dumps({k: {kk: res[k][kk] for kk in ("W", "K", "p")} for k in res if k.startswith(("primary", "valid"))},
                     indent=1))


if __name__ == "__main__":
    main()
