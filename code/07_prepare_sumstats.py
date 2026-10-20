"""Step C of the registered plan: turn raw GWAS files into LDSC-ready summary statistics, then munge.

Blinding guard: effect/SE columns of the REAL GWAS files (anything under data/raw/) are read only when --unblind is
given, which is reserved for the single primary run. In the simulations the same code runs on simulated files written in the
real files' exact formats (code/10_simulate.py).

  EAS (GWAS-SSF, GRCh38): positions lifted to hg19 and matched to 1000G EAS (cached by code/03 from identifier columns
       only); effective N per variant from the Direction string (cohort order SHA1, ICH1, KOR1, JPN1).
  EUR (de Lange 2017, GRCh37): chr:pos matched to 1000G EUR; effective N = 34,005 for every variant.
  FIN (FinnGen R13, GRCh38 with rsIDs): rsID matched to 1000G EUR with allele check; effective N = 10,378.
All: Z = effect / SE; strand-ambiguous SNPs dropped; duplicated rsIDs dropped; regression SNPs = HapMap3 no-MHC;
then munge_sumstats.py with defaults.

Output columns: SNP A1 A2 Z N P (P = two-sided P from Z; munge_sumstats requires it).

Usage: .venv/Scripts/python code/07_prepare_sumstats.py {EAS,EUR,FIN} RAW_FILE OUT_PREFIX [--unblind]
"""
import argparse
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats

from common import ANC, CHROMS, LDSC_DIR, LDSC_PY, ROOT

HM3 = ROOT / "data/raw/ldsc_ref/hm3_no_MHC.list.txt"
EAS_MATCHED = ROOT / "data/processed/eas_variants_matched_hg19.tsv.gz"
EAS_COHORTS = {"SHA1": (2552, 6279), "ICH1": (1611, 3724), "KOR1": (1619, 4419), "JPN1": (1590, 1034)}
EUR_NEFF = 4 / (1 / 12194 + 1 / 28072)
FIN_NEFF = 4 / (1 / 2608 + 1 / 497350)
COMP = str.maketrans("ACGT", "TGCA")


def neff(cases, controls):
    return 4 / (1 / cases + 1 / controls)


def eas_neff_from_direction(direction, cohorts=None):
    """Per-variant effective N: sum of the effective N of cohorts whose Direction character is not '?'."""
    cohorts = cohorts or EAS_COHORTS
    w = np.array([neff(*cohorts[c]) for c in ["SHA1", "ICH1", "KOR1", "JPN1"]])
    present = np.array([[ch != "?" for ch in s] for s in direction], dtype=float)
    return present @ w


def is_ambiguous(a1, a2):
    pair = a1.str.upper() + a2.str.upper()
    return pair.isin(["AT", "TA", "CG", "GC"])


def allele_match(g1, g2, r1, r2):
    """1 = same alleles (allowing strand flip), -1 = swapped, 0 = inconsistent."""
    g1, g2, r1, r2 = (s.str.upper() for s in (g1, g2, r1, r2))
    f1, f2 = g1.str.translate(COMP), g2.str.translate(COMP)
    same = ((g1 == r1) & (g2 == r2)) | ((f1 == r1) & (f2 == r2))
    swap = ((g1 == r2) & (g2 == r1)) | ((f1 == r2) & (f2 == r1))
    return np.where(same, 1, np.where(swap, -1, 0))


def load_eas_matched():
    return pd.read_csv(EAS_MATCHED, sep="	", dtype={"chromosome": str},
                       usecols=["chromosome", "base_pair_location", "g1", "g2", "rsid", "a1", "a2"])


def reference_bim(anc):
    cols = ["chr", "rsid", "cm", "pos", "a1", "a2"]
    return pd.concat(pd.read_csv(f"{ANC[anc]['bfile']}".format(chr=c) + ".bim", sep=r"\s+", header=None,
                                 names=cols, dtype={"chr": str}) for c in CHROMS)


def finish(df, hm3):
    """Common filters; output columns SNP A1 A2 Z N (A1 = effect allele on the reference strand)."""
    df = df[(df["am"] != 0) & ~is_ambiguous(df["a1"], df["a2"])]
    df = df[df["rsid"].isin(hm3)]
    df = df[~df["rsid"].duplicated(keep=False)]
    a1 = np.where(df["am"] == 1, df["a1"], df["a2"])
    a2 = np.where(df["am"] == 1, df["a2"], df["a1"])
    out = pd.DataFrame({"SNP": df["rsid"].values, "A1": a1, "A2": a2, "Z": df["Z"].values, "N": df["N"].values})
    out = out[np.isfinite(out["Z"]) & (out["N"] > 0)]
    # munge_sumstats (this LDSC version) requires a P column; two-sided P from Z, accurate in the tails
    return out.assign(P=2 * stats.norm.sf(np.abs(out["Z"].values)))


def prepare_eas(raw, hm3, matched=None, cohorts=None):
    d = pd.read_csv(raw, sep="\t", dtype={"chromosome": str},
                    usecols=["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta",
                             "standard_error", "Direction"])
    m = matched if matched is not None else load_eas_matched()
    d = d.rename(columns={"effect_allele": "g1", "other_allele": "g2"})
    x = d.merge(m, on=["chromosome", "base_pair_location", "g1", "g2"])
    x["am"] = allele_match(x["g1"], x["g2"], x["a1"], x["a2"])
    x["Z"] = x["beta"] / x["standard_error"]
    x["N"] = eas_neff_from_direction(x["Direction"].values, cohorts)
    return finish(x, hm3)


def prepare_eur(raw, hm3, bim=None, n_eff=EUR_NEFF):
    d = pd.read_csv(raw, sep="\t", usecols=["MarkerName", "Allele1", "Allele2", "Effect", "StdErr"])
    cp = d["MarkerName"].str.extract(r"^(\w+):(\d+)_")
    d["chr"], d["pos"] = cp[0], cp[1].astype(int)
    x = d.merge(bim if bim is not None else reference_bim("EUR"), on=["chr", "pos"])
    x["am"] = allele_match(x["Allele1"], x["Allele2"], x["a1"], x["a2"])
    x["Z"] = x["Effect"] / x["StdErr"]
    x["N"] = n_eff
    return finish(x, hm3)


def prepare_fin(raw, hm3):
    d = pd.read_csv(raw, sep="\t", usecols=["rsids", "ref", "alt", "beta", "sebeta"])
    d = d.assign(rsid=d["rsids"].str.split(",")).explode("rsid")
    x = d.merge(reference_bim("EUR")[["rsid", "a1", "a2"]], on="rsid")
    x["am"] = allele_match(x["alt"], x["ref"], x["a1"], x["a2"])  # FinnGen beta is for the alt allele
    x["Z"] = x["beta"] / x["sebeta"]
    x["N"] = FIN_NEFF
    return finish(x, hm3)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("source", choices=["EAS", "EUR", "FIN"])
    p.add_argument("raw")
    p.add_argument("out")
    p.add_argument("--unblind", action="store_true")
    a = p.parse_args()
    raw = pathlib.Path(a.raw).resolve()
    if (ROOT / "data/raw") in raw.parents and not a.unblind:
        sys.exit("Blinded: real GWAS effect sizes are read only in the primary run with --unblind.")
    hm3 = set(pd.read_csv(HM3, header=None)[0])
    out = {"EAS": prepare_eas, "EUR": prepare_eur, "FIN": prepare_fin}[a.source](raw, hm3)
    out.to_csv(a.out + ".prepared.tsv.gz", sep="\t", index=False)
    print(f"{a.source}: {len(out):,} regression SNPs written", flush=True)
    r = subprocess.run([str(LDSC_PY), str(LDSC_DIR / "munge_sumstats.py"), "--sumstats", a.out + ".prepared.tsv.gz",
                        "--N-col", "N", "--out", a.out], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"munge failed:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")


if __name__ == "__main__":
    main()
