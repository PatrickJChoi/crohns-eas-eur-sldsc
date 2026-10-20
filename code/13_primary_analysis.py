"""The single run of the registered primary analysis (doi:10.17605/OSF.IO/HQ4A8).

1. Step C on the real files (code/07 with --unblind): EAS (GCST90446793), EUR (de Lange 2017, GCST004132),
   FinnGen R13 K11_CD_STRICT2 (EUR replication).
2. Step D: 8 lineage models per dataset, baseline v1.2 + control + lineage, --overlap-annot, as registered.
3. Step E: primary Wald test EAS vs EUR; validation (i) EUR (de Lange) vs FinnGen, (ii) EAS vs FinnGen.
Writes results/main_analysis/primary_results.json and primary_table.tsv; derived summary statistics go to
data/processed/real/ (not tracked).
Run from the repo root:  .venv/Scripts/python code/13_primary_analysis.py
"""
import importlib.util
import json
import pathlib
import re
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from common import LINEAGES  # noqa: E402

REAL = ROOT / "data/processed/real"
OUT = ROOT / "results/main_analysis"
SOURCES = {  # name: (Step C source, raw file, ancestry for LD references)
    "EAS": ("EAS", ROOT / "data/raw/gwas/GCST90446793.tsv", "EAS"),
    "EUR_deLange": ("EUR", ROOT / "data/raw/gwas/cd_build37_40266_20161107.txt.gz", "EUR"),
    "EUR_FinnGen": ("FIN", ROOT / "data/raw/finngen_r13/finngen_R13_K11_CD_STRICT2.gz", "EUR"),
}


def load(name):
    spec = importlib.util.spec_from_file_location(name[3:-3], ROOT / "code" / name)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


fit = load("08_fit_models.py")
wtest = load("09_primary_test.py")


def diagnostics(log_path):
    txt = open(log_path).read()
    g = lambda pat: re.search(pat, txt)
    h2 = g(r"Total Observed scale h2: (\S+) \((\S+)\)")
    icpt = g(r"Intercept: (\S+) \((\S+)\)")
    return {"h2_obs": float(h2.group(1)), "h2_se": float(h2.group(2)),
            "lambda_gc": float(g(r"Lambda GC: (\S+)").group(1)), "mean_chi2": float(g(r"Mean Chi\^2: (\S+)").group(1)),
            "intercept": float(icpt.group(1)), "intercept_se": float(icpt.group(2))}


def main():
    REAL.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    fits, diag, nsnp = {}, {}, {}
    for name, (src, raw, anc) in SOURCES.items():
        prefix = REAL / name
        if not pathlib.Path(f"{prefix}.sumstats.gz").exists():
            r = subprocess.run([sys.executable, str(ROOT / "code/07_prepare_sumstats.py"), src, str(raw), str(prefix),
                                "--unblind"], capture_output=True, text=True, cwd=ROOT)
            if r.returncode != 0:
                sys.exit(f"Step C failed for {name}:\n{r.stdout}\n{r.stderr}")
            print(r.stdout.strip(), flush=True)
        nsnp[name] = int(re.search(r"Writing summary statistics for (\d+) SNPs",
                                   open(f"{prefix}.log").read()).group(1))
        fits[name] = fit.fit_set(f"{prefix}.sumstats.gz", anc, "primary", OUT / f"fits_{name}",
                                 baseline="baseline", overlap=True, workers=3)
        diag[name] = diagnostics(OUT / f"fits_{name}" / f"{anc}.{fit.SET_GROUPS['primary'][0][0]}.log")
        print(f"{name}: fitted 8 lineage models", flush=True)
    res = {"registration": "doi:10.17605/OSF.IO/HQ4A8", "regression_snps": nsnp, "diagnostics": diag,
           "primary_EAS_vs_EUR": wtest.wald(fits["EAS"], fits["EUR_deLange"]),
           "validation_i_deLange_vs_FinnGen": wtest.wald(fits["EUR_deLange"], fits["EUR_FinnGen"]),
           "validation_ii_EAS_vs_FinnGen": wtest.wald(fits["EAS"], fits["EUR_FinnGen"]),
           "tau_star": {k: v["tau_star"] for k, v in fits.items()},
           "fits": {k: v["fits"] for k, v in fits.items()}}
    (OUT / "primary_results.json").write_text(json.dumps(res, indent=1))
    p = res["primary_EAS_vs_EUR"]
    se = {k: [f["tau_se_ldsc"] * f["sd"] * f["M"] / f["h2"] for f in v["fits"]] for k, v in fits.items()}
    tab = pd.DataFrame({"lineage": LINEAGES,
                        "tau_star_EAS": fits["EAS"]["tau_star"], "se_EAS": p["se_a"],
                        "tau_star_EUR": fits["EUR_deLange"]["tau_star"], "se_EUR": p["se_b"],
                        "tau_star_FinnGen": fits["EUR_FinnGen"]["tau_star"],
                        "se_FinnGen": res["validation_i_deLange_vs_FinnGen"]["se_b"],
                        "D_EAS_minus_EUR": p["D"], "se_D": p["se_D"], "z": p["z"], "p": p["p_lineage"],
                        "p_holm": p["p_holm"]})
    tab.to_csv(OUT / "primary_table.tsv", sep="\t", index=False, float_format="%.4g")
    print(json.dumps({k: {kk: res[k][kk] for kk in ("W", "K", "p")} for k in res if k.startswith(("primary", "valid"))},
                     indent=1))


if __name__ == "__main__":
    main()
