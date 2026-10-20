"""Registered sensitivity analyses (plan section 6; they do not change the primary conclusion).

(a) drop every SNP with an ImmunoChip (ICH1) contribution from BOTH ancestries' regression SNPs;
(b) baselineLD v2.2 instead of baseline v1.2;
(c) top 5% and top 20% gene sets;
(d) +/-50 kb windows (with the +/-50 kb control).
Each repeats Steps D-E (EAS vs EUR de Lange) on the summary statistics prepared in code/13.
Writes results/main_analysis/sensitivity.json and sensitivity_table.tsv.
Run from the repo root:  .venv/Scripts/python code/14_sensitivity_analyses.py
"""
import gzip
import importlib.util
import json
import pathlib
import sys

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from common import ANC, CHROMS, LINEAGES  # noqa: E402

REAL = ROOT / "data/processed/real"
OUT = ROOT / "results/main_analysis/sensitivity"
EAS_MATCHED = ROOT / "data/processed/eas_variants_matched_hg19.tsv.gz"


def load(name):
    spec = importlib.util.spec_from_file_location(name[3:-3], ROOT / "code" / name)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


fit = load("08_fit_models.py")
wtest = load("09_primary_test.py")


def same_rows(anc):
    """baselineLD v2.2 must list exactly the SNP rows of baseline v1.2 (our annotations were printed for those)."""
    for c in CHROMS:
        a = pd.read_csv(str(ANC[anc]["baseline"]).format(chr=c) + ".l2.ldscore.gz", sep="\t", usecols=["SNP"])
        b = pd.read_csv(str(ANC[anc]["baselineLD"]).format(chr=c) + ".l2.ldscore.gz", sep="\t", usecols=["SNP"])
        if not a["SNP"].equals(b["SNP"]):
            raise SystemExit(f"{anc} chr{c}: baselineLD rows differ from baseline v1.2")


def drop_ich1():
    m = pd.read_csv(EAS_MATCHED, sep="\t", usecols=["rsid", "ich1"])
    ich1 = set(m.loc[m["ich1"], "rsid"])
    out = {}
    for name in ["EAS", "EUR_deLange"]:
        ss = pd.read_csv(REAL / f"{name}.sumstats.gz", sep="\t")
        kept = ss[~ss["SNP"].isin(ich1)]
        path = OUT / f"{name}.noICH1.sumstats.gz"
        kept.to_csv(path, sep="\t", index=False, compression="gzip")
        out[name] = (path, len(ss), len(kept))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res, rows = {}, []

    def record(tag, a, b, note=""):
        w = wtest.wald(a, b)
        res[tag] = {"wald": w, "tau_star_EAS": a["tau_star"], "tau_star_EUR": b["tau_star"], "note": note}
        rows.append({"analysis": tag, "W": w["W"], "K": w["K"], "p": w["p"], "note": note,
                     **{f"D_{L}": d for L, d in zip(LINEAGES, w["D"])}})
        print(f"{tag}: W={w['W']:.2f} p={w['p']:.3g} {note}", flush=True)

    # (a) symmetric ImmunoChip exclusion
    paths = drop_ich1()
    note = "; ".join(f"{k}: {n0:,} -> {n1:,} SNPs" for k, (_, n0, n1) in paths.items())
    a = fit.fit_set(paths["EAS"][0], "EAS", "primary", OUT / "a_noICH1_EAS", workers=3)
    b = fit.fit_set(paths["EUR_deLange"][0], "EUR", "primary", OUT / "a_noICH1_EUR", workers=3)
    record("a_drop_ImmunoChip_SNPs", a, b, note)
    # (b) baselineLD v2.2
    same_rows("EAS"), same_rows("EUR")
    a = fit.fit_set(REAL / "EAS.sumstats.gz", "EAS", "primary", OUT / "b_baselineLD_EAS", baseline="baselineLD", workers=3)
    b = fit.fit_set(REAL / "EUR_deLange.sumstats.gz", "EUR", "primary", OUT / "b_baselineLD_EUR",
                    baseline="baselineLD", workers=3)
    record("b_baselineLD_v2.2", a, b)
    # (c) top 5% / top 20% and (d) +/-50 kb
    for group in ["top05", "top20", "w50kb"]:
        a = fit.fit_set(REAL / "EAS.sumstats.gz", "EAS", group, OUT / f"{group}_EAS", workers=3)
        b = fit.fit_set(REAL / "EUR_deLange.sumstats.gz", "EUR", group, OUT / f"{group}_EUR", workers=3)
        record({"top05": "c_top5pct", "top20": "c_top20pct", "w50kb": "d_50kb_windows"}[group], a, b)
    (OUT.parent / "sensitivity.json").write_text(json.dumps(res, indent=1))
    pd.DataFrame(rows).to_csv(OUT.parent / "sensitivity_table.tsv", sep="\t", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
