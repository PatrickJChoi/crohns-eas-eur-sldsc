"""Donor-corrected Step A (post hoc correction, 2026-10-09; reported ALONGSIDE the
registered run, never replacing it).

The registered Step A (code/05) used the atlas field `Sample name` as the donor. Two children appear under several
sample names for separately sorted cell fractions (T036, T036NEG, T036POS; T110NEG, T110POS: same age, sex and library
names). Here the donor is the sample name with a trailing NEG/POS removed, giving 15 donors. Everything else is the
registered Step A code, run unchanged into separate output folders:
  data/processed/genesets_donorfix/   results/main_analysis/donorfix/
Run from the repo root:  .venv/Scripts/python code/20_donorfix_genesets.py
"""
import importlib.util
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("stepA", ROOT / "code/05_lineage_genesets.py")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)

_registered_read_obs = A.read_obs


def donor_of(sample_name):
    return re.sub(r"(NEG|POS)$", "", sample_name)


def read_obs_true_donors():
    o, genes = _registered_read_obs()
    o["donor"] = o["donor"].map(donor_of)
    return o, genes


def main():
    A.read_obs = read_obs_true_donors
    A.OUT = ROOT / "data/processed/genesets_donorfix"
    A.RES = ROOT / "results/main_analysis/donorfix"
    o, _ = read_obs_true_donors()
    keep = o["dx"].isin(A.KEEP_DX) & (o["region"] != "lymph node") & ~o["doublet"]
    n = o.loc[keep, "donor"].nunique()
    print(f"true donors among kept cells: {n}", flush=True)
    if n != 15:
        raise SystemExit(f"expected 15 donors after merging NEG/POS samples, got {n}")
    A.main()


if __name__ == "__main__":
    main()
