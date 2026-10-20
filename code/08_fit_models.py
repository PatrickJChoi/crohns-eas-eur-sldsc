"""Step D of the registered plan: one LDSC-SEG model per ancestry and gene set, then tau* and its
block-jackknife delete values.

Model: ldsc --h2 with ref-ld = baseline (v1.2 by default) + control annotation + one gene-set annotation, ancestry-
matched weights and frequency files, --overlap-annot --print-coefficients --print-delete-vals (200 blocks).
tau*_L = tau_L * sd_L * M / h2g (Gazal 2017), with M = number of reference SNPs with MAF >= 5% (M_5_50 of the
baseline 'base' column), sd_L = sqrt(p(1-p)) for the binary annotation with p = M_5_50(L) / M, and h2g = sum_c tau_c
M_5_50(c), all from the same model. Per-block tau* uses the per-block tau_L (.part_delete) and per-block h2g (.delete).

Library use: fit_set(sumstats, anc, names, outdir) -> dict with tau_star (K,), delete (B, K).
CLI: .venv/Scripts/python code/08_fit_models.py SUMSTATS.sumstats.gz {EAS,EUR} OUTDIR [--sets primary|random|...]
"""
import argparse
import json
import pathlib
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from common import ANC, ANNOT, CHROMS, LDSC_DIR, LDSC_PY, LINEAGES, slug

SET_GROUPS = {
    "primary": ([f"L_{slug(L)}_top10" for L in LINEAGES], "control_100kb"),
    "random": ([f"R{i}_random" for i in range(1, 9)], "control_100kb"),
    "top05": ([f"L_{slug(L)}_top05" for L in LINEAGES], "control_100kb"),
    "top20": ([f"L_{slug(L)}_top20" for L in LINEAGES], "control_100kb"),
    "w50kb": ([f"L_{slug(L)}_top10_50kb" for L in LINEAGES], "control_50kb"),
}


def m550(prefix):
    """Sum over chromosomes of the M_5_50 vector for an LDSC prefix (one value per annotation column)."""
    return sum(np.loadtxt(f"{prefix}{c}.l2.M_5_50", ndmin=1) for c in CHROMS)


def run_one(sumstats, anc, name, control, outdir, baseline="baseline", overlap=True, annot_root=ANNOT):
    out = pathlib.Path(outdir) / f"{anc}.{name}"
    base_prefix = str(ANC[anc][baseline]).replace("{chr}", "")
    ctrl_prefix = str(ANNOT / anc / control / f"{control}.")
    set_prefix = str(annot_root / anc / name / f"{name}.")
    cmd = [str(LDSC_PY), str(LDSC_DIR / "ldsc.py"), "--h2", str(sumstats),
           "--ref-ld-chr", ",".join([base_prefix, ctrl_prefix, set_prefix]),
           "--w-ld-chr", str(ANC[anc]["weights"]).replace("{chr}", ""),
           "--print-coefficients", "--print-delete-vals", "--out", str(out)]
    if overlap:
        cmd += ["--overlap-annot", "--frqfile-chr", str(ANC[anc]["frq"]).replace("{chr}", "")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not pathlib.Path(f"{out}.part_delete").exists():
        raise RuntimeError(f"ldsc --h2 failed ({anc} {name}):\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    # full-data coefficients: .results (overlap) or the 'Coefficients:' log line
    if overlap:
        tau = pd.read_csv(f"{out}.results", sep="\t")["Coefficient"].values
    else:
        line = next(l for l in open(f"{out}.log") if l.startswith("Coefficients:"))
        tau = np.array(line.split(":", 1)[1].split(), dtype=float)
    part = np.loadtxt(f"{out}.part_delete")          # (B, n_annot) per-block tau
    tot = np.loadtxt(f"{out}.delete").reshape(-1)    # (B,) per-block total h2
    M = np.concatenate([m550(base_prefix), m550(ctrl_prefix), m550(set_prefix)])
    ts = tau_star_from(tau, M, part, tot)
    log_h2 = float(re.search(r"Total Observed scale h2: (\S+)", open(f"{out}.log").read()).group(1))
    if not np.isclose(ts["h2"], log_h2, rtol=2e-3, atol=1e-6):
        raise RuntimeError(f"h2 check failed for {anc} {name}: {ts['h2']} vs log {log_h2}")
    return {"name": name, "tau": float(tau[-1]), **ts,
            "tau_se_ldsc": float(pd.read_csv(f"{out}.results", sep="\t")["Coefficient_std_error"].values[-1])
            if overlap else None}


def tau_star_from(tau, M, part, tot):
    """tau* of the last annotation and its per-block values (Gazal 2017).
    tau: full-data coefficients (n_annot,); M: M_5_50 per annotation, first = baseline 'base' (all MAF >= 5% SNPs);
    part: per-block coefficients (B, n_annot); tot: per-block total h2 (B,)."""
    tau, M = np.asarray(tau, float), np.asarray(M, float)
    M_all = M[0]
    p = M[-1] / M_all
    sd = np.sqrt(p * (1 - p))                          # SD of a binary annotation over MAF >= 5% SNPs
    h2 = float(tau @ M)
    return {"h2": h2, "sd": float(sd), "M": float(M_all), "tau_star": float(tau[-1] * sd * M_all / h2),
            "delete": (np.asarray(part)[:, -1] * sd * M_all / np.asarray(tot).reshape(-1)).tolist()}


def fit_set(sumstats, anc, group, outdir, baseline="baseline", overlap=True, workers=4):
    names, control = SET_GROUPS[group]
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(workers) as ex:
        fits = list(ex.map(lambda n: run_one(sumstats, anc, n, control, outdir, baseline, overlap), names))
    res = {"anc": anc, "group": group, "baseline": baseline, "names": names,
           "tau_star": [f["tau_star"] for f in fits], "delete": np.array([f["delete"] for f in fits]).T.tolist(),
           "fits": [{k: v for k, v in f.items() if k != "delete"} for f in fits]}
    (pathlib.Path(outdir) / f"{anc}.{group}.{baseline}.json").write_text(json.dumps(res))
    return res


def main():
    p = argparse.ArgumentParser()
    p.add_argument("sumstats")
    p.add_argument("anc", choices=["EAS", "EUR"])
    p.add_argument("outdir")
    p.add_argument("--sets", default="primary", choices=list(SET_GROUPS))
    p.add_argument("--baseline", default="baseline", choices=["baseline", "baselineLD"])
    p.add_argument("--no-overlap", action="store_true")
    a = p.parse_args()
    r = fit_set(a.sumstats, a.anc, a.sets, a.outdir, a.baseline, not a.no_overlap)
    print(json.dumps({"anc": r["anc"], "tau_star": r["tau_star"]}))


if __name__ == "__main__":
    main()
