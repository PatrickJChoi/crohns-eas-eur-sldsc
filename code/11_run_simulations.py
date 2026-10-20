"""Simulation driver: run simulated GWAS through the registered pipeline (Steps C-E, unchanged code).

Two stages, so the simulator's reference genotypes and the LDSC processes are never in memory at the same time:
  gen  simulate EAS and EUR summary statistics (code/10), write them in the real files' formats, prepare and munge
       them (code/07). Keeps <scenario>/rep<r>/{EAS,EUR}.sumstats.gz in data/processed/sim/ (not tracked).
  fit  fit the 8 lineage models per ancestry (code/08; plus 8 random gene sets if the scenario asks), run the joint
       Wald test (code/09), append one JSON line per replicate to results/simulations/sims/<scenario>.jsonl.
Both stages skip replicates already done, so they can be resumed.

Scenario file (results/simulations/scenarios/<name>.json):
  {"EAS": {"h2": ..., "e_ctrl": ..., "e": {...}}, "EUR": {...}, "random_sets": false}
Usage: .venv/Scripts/python code/11_run_simulations.py {gen,fit} SCENARIO [SCENARIO ...] --reps 20 --seed 1
       [--no-overlap] [--workers 4]
"""
import argparse
import importlib.util
import json
import pathlib
import re
import subprocess
import sys
import time
import zlib

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import LDSC_DIR, LDSC_PY, ROOT  # noqa: E402


def load(name):
    spec = importlib.util.spec_from_file_location(name[3:-3], ROOT / "code" / name)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


prep = load("07_prepare_sumstats.py")
fit = load("08_fit_models.py")
wtest = load("09_primary_test.py")
SIM = ROOT / "data/processed/sim"
RES = ROOT / "results/simulations"


def scenario(name):
    return json.loads((RES / "scenarios" / f"{name}.json").read_text())


def munge(prepared, out):
    r = subprocess.run([str(LDSC_PY), str(LDSC_DIR / "munge_sumstats.py"), "--sumstats", str(prepared),
                        "--N-col", "N", "--out", str(out)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stdout[-2000:] + r.stderr[-2000:])


def gen(a):
    sim = load("10_simulate.py")
    hm3 = set(pd.read_csv(prep.HM3, header=None)[0])
    rows = {"EAS": sim.eas_rows(hm3), "EUR": sim.eur_rows(hm3)}
    matched, bim = prep.load_eas_matched(), prep.reference_bim("EUR")
    n_map = rows["EAS"].set_index("rsid")["N"]
    N_of = {"EAS": lambda s: n_map.reindex(s).fillna(n_map.median()).values,
            "EUR": lambda s: np.full(len(s), sim.EUR_NEFF)}
    t0 = time.time()
    refs = {anc: sim.Reference(anc, set(rows[anc]["rsid"])) for anc in ["EAS", "EUR"]}
    print(f"reference genotypes loaded in {time.time() - t0:.0f}s", flush=True)
    for name in a.scenarios:
        scen = scenario(name)
        rng = np.random.default_rng([a.seed, zlib.crc32(name.encode())])  # stable across runs
        for rep in range(a.reps):
            work = SIM / name / f"rep{rep}"
            # draw this replicate's random numbers even when skipping, so resumed runs reproduce the same data
            if all((work / f"{anc}.sumstats.gz").exists() for anc in ["EAS", "EUR"]):
                for anc in ["EAS", "EUR"]:
                    refs[anc].simulate(scen[anc], N_of[anc], rng)
                continue
            t1 = time.time()
            work.mkdir(parents=True, exist_ok=True)
            for anc in ["EAS", "EUR"]:
                z = refs[anc].simulate(scen[anc], N_of[anc], rng).set_index("SNP")["Z"]
                x = rows[anc].assign(z=rows[anc]["rsid"].map(z)).dropna(subset=["z"])
                if anc == "EAS":
                    se = 1 / np.sqrt(x["N"].values)
                    raw = pd.DataFrame({"chromosome": x["chromosome"], "base_pair_location": x["base_pair_location"],
                                        "effect_allele": x["g1"], "other_allele": x["g2"], "beta": x["z"] * se,
                                        "standard_error": se, "Direction": x["Direction"]})
                    raw_path = work / "sim_EAS.tsv"
                    raw.to_csv(raw_path, sep="\t", index=False)
                    prepared = prep.prepare_eas(raw_path, hm3, matched)
                else:
                    se = 1 / np.sqrt(sim.EUR_NEFF)
                    raw = pd.DataFrame({"MarkerName": x["MarkerName"], "Allele1": x["Allele1"],
                                        "Allele2": x["Allele2"], "Effect": x["z"] * se, "StdErr": se,
                                        "Direction": x["Direction"]})
                    raw_path = work / "sim_EUR.tsv.gz"
                    raw.to_csv(raw_path, sep="\t", index=False)
                    prepared = prep.prepare_eur(raw_path, hm3, bim)
                prepared.to_csv(work / f"{anc}.prepared.tsv.gz", sep="\t", index=False)
                munge(work / f"{anc}.prepared.tsv.gz", work / anc)
                raw_path.unlink()
                (work / f"{anc}.prepared.tsv.gz").unlink()
            print(f"gen {name} rep {rep} ({time.time() - t1:.0f}s)", flush=True)


def h2_z(outdir, anc, name):
    txt = open(pathlib.Path(outdir) / f"{anc}.{name}.log").read()
    m = re.search(r"Total Observed scale h2: (\S+) \((\S+)\)", txt)
    return float(m.group(1)) / float(m.group(2)), float(re.search(r"Mean Chi\^2: (\S+)", txt).group(1))


def fit_stage(a):
    (RES / "sims").mkdir(parents=True, exist_ok=True)
    for name in a.scenarios:
        scen = scenario(name)
        out_file = RES / "sims" / f"{name}.jsonl"
        done = {json.loads(l)["rep"] for l in out_file.read_text().splitlines()} if out_file.exists() else set()
        for rep in range(a.reps):
            work = SIM / name / f"rep{rep}"
            if rep in done or not (work / "EUR.sumstats.gz").exists():
                continue
            t1 = time.time()
            groups = ["primary"] + (["random"] if scen.get("random_sets") else [])
            fits, diag = {}, {}
            for anc in ["EAS", "EUR"]:
                for g in groups:
                    fits[(anc, g)] = fit.fit_set(work / f"{anc}.sumstats.gz", anc, g, work / "fits",
                                                 overlap=not a.no_overlap, workers=a.workers)
                diag[anc] = h2_z(work / "fits", anc, fit.SET_GROUPS["primary"][0][0])
            out = {"scenario": name, "rep": rep, "seed": a.seed,
                   "h2_z": {k: v[0] for k, v in diag.items()}, "mean_chi2": {k: v[1] for k, v in diag.items()},
                   "primary": wtest.wald(fits[("EAS", "primary")], fits[("EUR", "primary")]),
                   "tau_star": {anc: fits[(anc, "primary")]["tau_star"] for anc in ["EAS", "EUR"]}}
            if "random" in groups:
                out["random"] = wtest.wald(fits[("EAS", "random")], fits[("EUR", "random")])
                out["tau_star_random"] = {anc: fits[(anc, "random")]["tau_star"] for anc in ["EAS", "EUR"]}
            with open(out_file, "a") as f:
                f.write(json.dumps(out) + "\n")
            print(f"fit {name} rep {rep}: W={out['primary']['W']:.2f} p={out['primary']['p']:.3f} "
                  f"h2 z EAS {out['h2_z']['EAS']:.1f} EUR {out['h2_z']['EUR']:.1f} ({time.time() - t1:.0f}s)",
                  flush=True)
            for f in (work / "fits").iterdir():
                f.unlink()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("stage", choices=["gen", "fit"])
    p.add_argument("scenarios", nargs="+")
    p.add_argument("--reps", type=int, default=1)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--no-overlap", action="store_true")
    p.add_argument("--workers", type=int, default=4)
    a = p.parse_args()
    gen(a) if a.stage == "gen" else fit_stage(a)


if __name__ == "__main__":
    main()
