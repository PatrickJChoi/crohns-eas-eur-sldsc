"""Registered exploratory analyses (plan section 7). All results are EXPLORATORY.

E1 fine cell types: per-ancestry tau* for each fine type (code/16 sets) and their cross-ancestry correlation.
E2 primary model excluding regression SNPs within +/-1 Mb of NOD2 (GENCODE v19), both ancestries.
E3 ulcerative colitis with the same pipeline (EAS GCST90446795, EUR GCST004133), to see whether any CD pattern is
   CD-specific. EAS UC variants are matched to 1000G through the identifier table built from the CD file (same cohorts
   and positions); coverage is reported. UC cohort sizes from Liu 2023 Methods.
E4 small-intestine-only and large-intestine-only lineage references (code/16 sets).
E5 tau* and the primary Wald test recomputed on 10 donor-bootstrap gene sets (code/05).
Models are fitted without --overlap-annot (verified in the simulations: coefficients and jackknife unchanged). Each analysis writes
results/main_analysis/exploratory/E*.json and is skipped if already done.
Run from the repo root:  .venv/Scripts/python code/18_exploratory_analyses.py [E1 E2 E3 E4 E5]
"""
import importlib.util
import json
import pathlib
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from common import ANC, LDSC_DIR, LDSC_PY, LINEAGES, slug  # noqa: E402


def load(name):
    spec = importlib.util.spec_from_file_location(name[3:-3], ROOT / "code" / name)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


prep, fit, wtest = load("07_prepare_sumstats.py"), load("08_fit_models.py"), load("09_primary_test.py")
B20 = load("17_exploratory_ldscores.py")
REAL = ROOT / "data/processed/real"
OUT = ROOT / "results/main_analysis/exploratory"
EXP_ANNOT = ROOT / "data/processed/annot_exploratory"
UC_COHORTS = {"SHA1": (2400, 6279), "ICH1": (1124, 3724), "KOR1": (1569, 4419), "JPN1": (1769, 1034)}
UC_EUR_NEFF = 4 / (1 / 12366 + 1 / 33609)
ANCS = {"EAS": REAL / "EAS.sumstats.gz", "EUR": REAL / "EUR_deLange.sumstats.gz"}


def fit_many(sumstats, anc, names, outdir, annot_root=fit.ANNOT, workers=3):
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(workers) as ex:
        fits = list(ex.map(lambda n: fit.run_one(sumstats, anc, n, "control_100kb", outdir, overlap=False,
                                                 annot_root=annot_root), names))
    return {"names": names, "tau_star": [f["tau_star"] for f in fits],
            "delete": np.array([f["delete"] for f in fits]).T.tolist(),
            "se": [float(np.sqrt(np.diag(wtest.jackknife_cov([f["tau_star"]], np.array([f["delete"]]).T)))[0])
                   for f in fits]}


def done(tag):
    return (OUT / f"{tag}.json").exists()


def save(tag, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{tag}.json").write_text(json.dumps(obj, indent=1))
    print(f"{tag} written", flush=True)


def primary_names():
    return [f"L_{slug(L)}_top10" for L in LINEAGES]


def e1():
    names = sorted(n.name for n in (EXP_ANNOT / "EAS").iterdir() if n.name.startswith("F_"))
    r = {a: fit_many(ss, a, names, OUT / f"E1_{a}", EXP_ANNOT) for a, ss in ANCS.items()}
    x, y = np.array(r["EAS"]["tau_star"]), np.array(r["EUR"]["tau_star"])
    save("E1_fine_types", {"names": names, "tau_star_EAS": x.tolist(), "se_EAS": r["EAS"]["se"],
                           "tau_star_EUR": y.tolist(), "se_EUR": r["EUR"]["se"],
                           "pearson": stats.pearsonr(x, y)[:2], "spearman": stats.spearmanr(x, y)[:2]})


def e2():
    A = load("05_lineage_genesets.py")
    g = A.gencode_genes()
    nod2 = g[g["gencode_name"] == "NOD2"].iloc[0]
    lo, hi = int(nod2["start"]) - 1_000_000, int(nod2["end"]) + 1_000_000
    region = set()
    for anc in ANCS:
        bim = pd.read_csv(str(ANC[anc]["bfile"]).format(chr=16) + ".bim", sep=r"\s+", header=None,
                          names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
        region |= set(bim.loc[(bim["BP"] >= lo) & (bim["BP"] <= hi), "SNP"])
    r, counts = {}, {}
    for anc, ss in ANCS.items():
        d = pd.read_csv(ss, sep="\t")
        kept = d[~d["SNP"].isin(region)]
        counts[anc] = [len(d), len(kept)]
        p = OUT / f"E2_{anc}.noNOD2.sumstats.gz"
        OUT.mkdir(parents=True, exist_ok=True)
        kept.to_csv(p, sep="\t", index=False, compression="gzip")
        r[anc] = fit_many(p, anc, primary_names(), OUT / f"E2_{anc}")
    save("E2_exclude_NOD2", {"region_hg19": f"chr16:{lo}-{hi}", "snps_before_after": counts,
                             "wald": wtest.wald(r["EAS"], r["EUR"]),
                             "tau_star": {a: v["tau_star"] for a, v in r.items()}})


def munge(prepared, out):
    rr = subprocess.run([str(LDSC_PY), str(LDSC_DIR / "munge_sumstats.py"), "--sumstats", str(prepared),
                         "--N-col", "N", "--out", str(out)], capture_output=True, text=True)
    if rr.returncode != 0:
        raise RuntimeError(rr.stdout[-2000:] + rr.stderr[-2000:])


def e3():
    hm3 = set(pd.read_csv(prep.HM3, header=None)[0])
    eas_raw = ROOT / "data/raw/gwas/GCST90446795.tsv"
    eur_raw = ROOT / "data/raw/gwas/uc_build37_45975_20161107.txt.gz"
    raw_n = sum(1 for _ in open(eas_raw)) - 1
    eas = prep.prepare_eas(eas_raw, hm3, cohorts=UC_COHORTS)
    eur_dir = pd.read_csv(eur_raw, sep="\t", usecols=["Direction"])["Direction"]
    eur = prep.prepare_eur(eur_raw, hm3, n_eff=UC_EUR_NEFF)
    info = {"EAS_raw_variants": raw_n, "EAS_regression_snps": len(eas), "EUR_regression_snps": len(eur),
            "EUR_variants_missing_a_cohort": int(eur_dir.str.contains(r"\?").sum())}
    r = {}
    for anc, d in [("EAS", eas), ("EUR", eur)]:
        p = OUT / f"E3_UC_{anc}"
        OUT.mkdir(parents=True, exist_ok=True)
        d.to_csv(f"{p}.prepared.tsv.gz", sep="\t", index=False)
        munge(f"{p}.prepared.tsv.gz", p)
        r[anc] = fit_many(f"{p}.sumstats.gz", anc, primary_names(), OUT / f"E3_{anc}")
    save("E3_ulcerative_colitis", {"info": info, "wald": wtest.wald(r["EAS"], r["EUR"]),
                                   "tau_star": {a: v["tau_star"] for a, v in r.items()},
                                   "se": {a: v["se"] for a, v in r.items()}})


def e4():
    res = {}
    for region in ["SI", "LI"]:
        names = sorted(n.name for n in (EXP_ANNOT / "EAS").iterdir() if n.name.startswith(f"{region}_"))
        r = {a: fit_many(ss, a, names, OUT / f"E4_{region}_{a}", EXP_ANNOT) for a, ss in ANCS.items()}
        res[region] = {"names": names, "wald": wtest.wald({**r["EAS"], "names": names}, {**r["EUR"], "names": names}),
                       "tau_star": {a: v["tau_star"] for a, v in r.items()}}
    save("E4_region_specific", res)


def e5():
    res = []
    for b in range(10):
        names = [f"BOOT{b}_{B20.clean(L)}" for L in LINEAGES]
        r = {a: fit_many(ss, a, names, OUT / f"E5_b{b}_{a}", EXP_ANNOT) for a, ss in ANCS.items()}
        w = wtest.wald(r["EAS"], r["EUR"])
        res.append({"boot": b, "W": w["W"], "p": w["p"], "D": w["D"], "se_D": w["se_D"],
                    "tau_star": {a: v["tau_star"] for a, v in r.items()}})
    save("E5_bootstrap_genesets", res)


def main():
    todo = sys.argv[1:] or ["E1", "E2", "E3", "E4", "E5"]
    tags = {"E1": "E1_fine_types", "E2": "E2_exclude_NOD2", "E3": "E3_ulcerative_colitis",
            "E4": "E4_region_specific", "E5": "E5_bootstrap_genesets"}
    for t in todo:
        if done(tags[t]):
            print(f"{tags[t]} already done", flush=True)
            continue
        {"E1": e1, "E2": e2, "E3": e3, "E4": e4, "E5": e5}[t]()


if __name__ == "__main__":
    main()
