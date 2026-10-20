"""Post hoc check (2026-09-09) of the registered assumption Cov(D) = Cov(A) + Cov(B).

The registered test treats the two datasets' jackknife errors as independent. LDSC's block jackknife also reflects where
causal variants sit in the genome; if datasets share causal variants, errors correlate and the summed covariance
overstates Var(D). Here all three CD datasets are restricted to one common SNP set, so LDSC's 200 jackknife blocks cover
the same genomic stretches in every dataset (same SNPs, same order). Then:
  (a) the registered Wald test on the common SNPs (summed covariances);
  (b) a direct jackknife of D: per-block D_b = tau*_A,b - tau*_B,b, pseudovalues of D, Cov(D) from them;
  (c) the per-lineage correlation between the two datasets' tau* pseudovalues.
Common set = SNPs in all three munged files, in both ancestries' baseline LD-score rows and weight files, with
chi2 <= 80 in every dataset (LDSC's default cut, applied jointly so no dataset drops different SNPs).
Writes results/main_analysis/aligned_block_check.json. Run:  .venv/Scripts/python code/23_aligned_block_check.py
"""
import importlib.util
import json
import pathlib
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from common import ANC, CHROMS  # noqa: E402


def load(name):
    spec = importlib.util.spec_from_file_location(name[3:-3], ROOT / "code" / name)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


fit, wtest = load("08_fit_models.py"), load("09_primary_test.py")
REAL = ROOT / "data/processed/real"
WORK = ROOT / "data/processed/aligned"
OUT = ROOT / "results/main_analysis/aligned"
DATASETS = {"EAS": "EAS", "EUR_deLange": "EUR", "EUR_FinnGen": "EUR"}


def snps_of(prefix_pattern, suffix):
    return set(pd.concat(pd.read_csv(str(prefix_pattern).format(chr=c) + suffix, sep="\t", usecols=["SNP"])["SNP"]
                         for c in CHROMS))


def direct(res_a, res_b):
    ta, tb = np.asarray(res_a["tau_star"]), np.asarray(res_b["tau_star"])
    da, db = np.asarray(res_a["delete"]), np.asarray(res_b["delete"])
    B = da.shape[0]
    D, Db = ta - tb, da - db
    pseudo = B * D[None, :] - (B - 1) * Db
    C = np.atleast_2d(np.cov(pseudo.T, ddof=1)) / B
    W = float(D @ np.linalg.solve(C, D))
    pa, pb = B * ta[None, :] - (B - 1) * da, B * tb[None, :] - (B - 1) * db
    corr = [float(np.corrcoef(pa[:, k], pb[:, k])[0, 1]) for k in range(len(D))]
    return {"W": W, "K": len(D), "p": float(stats.chi2.sf(W, len(D))), "se_D": np.sqrt(np.diag(C)).tolist(),
            "D": D.tolist(), "pseudovalue_corr": corr}


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    ss = {k: pd.read_csv(REAL / f"{k}.sumstats.gz", sep="\t").dropna(subset=["Z", "N"]) for k in DATASETS}
    common = set.intersection(*(set(d["SNP"]) for d in ss.values()))
    for anc in ["EAS", "EUR"]:
        common &= snps_of(ANC[anc]["baseline"], ".l2.ldscore.gz")
        common &= snps_of(ANC[anc]["weights"], ".l2.ldscore.gz")
    high = set().union(*(set(d.loc[d["Z"] ** 2 > 80, "SNP"]) for d in ss.values()))
    common -= high
    print(f"common SNPs: {len(common):,} (removed {len(high)} with chi2 > 80 in any dataset)", flush=True)
    fits = {}
    for k, anc in DATASETS.items():
        path = WORK / f"{k}.sumstats.gz"
        ss[k][ss[k]["SNP"].isin(common)].to_csv(path, sep="\t", index=False, compression="gzip")
        fits[k] = fit.fit_set(path, anc, "primary", OUT / f"fits_{k}", overlap=True, workers=2)
        print(f"{k}: fitted", flush=True)
    # every fit must use exactly the common SNPs, so the 200 blocks coincide across datasets
    used = {}
    for k, anc in DATASETS.items():
        for log in (OUT / f"fits_{k}").glob("*.log"):
            counts = re.findall(r"\((\d+) SNPs remain\)|(\d+) SNPs remain", log.read_text())
            last = [int(a or b) for a, b in counts][-1] if counts else None  # final count LDSC used
            used[f"{k}/{log.stem}"] = last
    if len(set(used.values())) != 1 or next(iter(used.values())) != len(common):
        raise SystemExit(f"fits did not all use the common SNP set: {sorted(set(used.values()))} vs {len(common)}")
    res = {"common_snps": len(common), "removed_high_chi2": len(high)}
    for a, b in [("EAS", "EUR_deLange"), ("EUR_deLange", "EUR_FinnGen"), ("EAS", "EUR_FinnGen")]:
        res[f"{a}_vs_{b}"] = {"summed_covariance": wtest.wald(fits[a], fits[b]), "direct_aligned": direct(fits[a], fits[b])}
        s, d = res[f"{a}_vs_{b}"]["summed_covariance"], res[f"{a}_vs_{b}"]["direct_aligned"]
        print(f"{a} vs {b}: summed W={s['W']:.2f} p={s['p']:.3f} | direct W={d['W']:.2f} p={d['p']:.3f} | "
              f"pseudovalue corr {np.round(d['pseudovalue_corr'], 2).tolist()}", flush=True)
    (ROOT / "results/main_analysis/aligned_block_check.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
