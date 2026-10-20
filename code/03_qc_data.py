"""Data feasibility checks. No association results are summarised: the
checks use only variant identifiers, alleles, allele frequencies, the meta-analysis
"Direction" strings (which cohorts contributed a variant) and atlas metadata.

Steps
1. Unpack the 1000 Genomes reference archives into data/processed/ref/.
2. East Asian GWAS (GRCh38): which cohorts contributed each variant (Direction order
   SHA1, ICH1, KOR1, JPN1 per the consortium readme); how many variants came from the
   ImmunoChip cohort ICH1; lift positions to hg19 and match to 1000G EAS rsIDs/alleles;
   coverage of the HapMap3 regression SNPs; allele-frequency agreement with 1000G EAS.
3. European GWAS (GRCh37): match chr:pos and alleles to 1000G EUR rsIDs; HapMap3
   coverage; cohort-contribution patterns (IBDseq, IIBDGC, GWAS3).
4. Gut Cell Atlas: donors, ages, regions, diagnoses and cell-type labels (metadata only).
Writes results/data_qc/qc_summary.txt and supporting tables.
Run from the repo root with the analysis environment:  .venv/Scripts/python code/03_qc_data.py
"""
import gzip
import pathlib
import sys
import tarfile

import h5py
import numpy as np
import pandas as pd
from pyliftover import LiftOver

RAW = pathlib.Path("data/raw")
REF = pathlib.Path("data/processed/ref")
PROC = pathlib.Path("data/processed")
OUT = pathlib.Path("results/data_qc")
EAS_COHORTS = ["SHA1", "ICH1", "KOR1", "JPN1"]  # readme.sst.txt (IIBDGC Drive, version 12202022)
COMP = str.maketrans("ACGT", "TGCA")
summary = []


def say(msg):
    print(msg)
    summary.append(msg)


def unpack():
    REF.mkdir(parents=True, exist_ok=True)
    for tgz in sorted((RAW / "ldsc_ref").glob("*.tgz")):
        flag = REF / (tgz.name + ".done")
        if flag.exists():
            continue
        # each archive gets its own folder: the EAS and EUR baseline archives hold identically named files
        with tarfile.open(tgz) as t:
            t.extractall(REF / tgz.name.removesuffix(".tgz"), filter="data")
        flag.touch()


def find_bims(pattern):
    bims = sorted(REF.rglob(pattern), key=lambda p: int(p.stem.split(".")[-1]))
    if len(bims) != 22:
        raise SystemExit(f"expected 22 bim files for {pattern}, found {len(bims)}")
    return bims


def read_bims(pattern):
    cols = ["chr", "rsid", "cm", "pos", "a1", "a2"]
    return pd.concat([pd.read_csv(b, sep=r"\s+", header=None, names=cols, dtype={"chr": str}) for b in find_bims(pattern)])


def is_ambiguous(df):
    """A/T and C/G SNPs: strand cannot be resolved from the alleles alone."""
    pair = df["g1"].str.upper() + df["g2"].str.upper()
    return pair.isin(["AT", "TA", "CG", "GC"])


def allele_match(df):
    """1 = same alleles, -1 = swapped, 0 = mismatch (strand flips counted as matches)."""
    g1, g2, r1, r2 = (df[c].str.upper() for c in ("g1", "g2", "a1", "a2"))
    same = ((g1 == r1) & (g2 == r2)) | ((g1.str.translate(COMP) == r1) & (g2.str.translate(COMP) == r2))
    swap = ((g1 == r2) & (g2 == r1)) | ((g1.str.translate(COMP) == r2) & (g2.str.translate(COMP) == r1))
    return np.where(same, 1, np.where(swap, -1, 0))


def eas_checks(hm3):
    cols = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "Direction"]
    d = pd.read_csv(RAW / "gwas/GCST90446793.tsv", sep="\t", usecols=cols, dtype={"chromosome": str})
    say(f"EAS variants in file: {len(d):,}")
    present = np.array([[c != "?" for c in s] for s in d["Direction"].values])
    pat = pd.Series(["".join("x" if p else "-" for p in row) for row in present]).value_counts()
    say("EAS cohort-contribution patterns (SHA1 ICH1 KOR1 JPN1; x = contributed): "
        + "; ".join(f"{k}: {v:,}" for k, v in pat.head(8).items()))
    d["ich1"] = present[:, 1]
    say(f"EAS variants with an ImmunoChip (ICH1) contribution: {d['ich1'].sum():,} ({d['ich1'].mean():.1%})")
    # lift GRCh38 -> hg19
    lo = LiftOver(str(RAW / "liftover/hg38ToHg19.over.chain.gz"))
    hg19 = []
    for c, p in zip(d["chromosome"].values, d["base_pair_location"].values):
        r = lo.convert_coordinate(f"chr{c}", int(p) - 1)
        hg19.append(r[0][1] + 1 if r and r[0][0] == f"chr{c}" else -1)
    d["pos19"] = hg19
    say(f"EAS variants lifted to hg19 on the same chromosome: {(d['pos19'] > 0).mean():.1%}")
    ref = read_bims("1000G.EAS.QC.*.bim")
    m = d[d["pos19"] > 0].merge(ref, left_on=["chromosome", "pos19"], right_on=["chr", "pos"])
    m = m.rename(columns={"effect_allele": "g1", "other_allele": "g2"})
    m["am"] = allele_match(m)
    m = m[m["am"] != 0]
    m["ambiguous"] = is_ambiguous(m)
    say(f"EAS variants matched to a 1000G EAS rsID with consistent alleles: {m['rsid'].nunique():,} "
        f"(strand-ambiguous A/T or C/G: {m['ambiguous'].mean():.1%})")
    keep = ["chromosome", "base_pair_location", "pos19", "rsid", "g1", "g2", "a1", "a2", "am", "ambiguous",
            "ich1", "Direction"]
    PROC.mkdir(parents=True, exist_ok=True)
    m[keep].to_csv(PROC / "eas_variants_matched_hg19.tsv.gz", sep="\t", index=False)
    in_hm3 = m[m["rsid"].isin(hm3)]
    say(f"HapMap3 regression SNPs covered by EAS GWAS: {in_hm3['rsid'].nunique():,} of {len(hm3):,} "
        f"({in_hm3['rsid'].nunique() / len(hm3):.1%}); of these with an ICH1 contribution: {in_hm3['ich1'].mean():.1%}")
    return m


def eur_checks(hm3):
    d = pd.read_csv(RAW / "gwas/cd_build37_40266_20161107.txt.gz", sep="\t",
                    usecols=["MarkerName", "Allele1", "Allele2", "Direction"])
    say(f"EUR variants in file: {len(d):,}")
    pat = d["Direction"].map(lambda s: "".join("-" if c == "?" else "x" for c in s)).value_counts()
    say("EUR cohort-contribution patterns (IBDseq IIBDGC GWAS3, per column order of per-cohort P values; "
        "x = contributed): " + "; ".join(f"{k}: {v:,}" for k, v in pat.head(8).items()))
    chrpos = d["MarkerName"].str.extract(r"^(\w+):(\d+)_")
    d["chr"], d["pos"] = chrpos[0], chrpos[1].astype(int)
    ref = read_bims("1000G.EUR.QC.*.bim")
    m = d.merge(ref, on=["chr", "pos"]).rename(columns={"Allele1": "g1", "Allele2": "g2"})
    m["am"] = allele_match(m)
    m = m[m["am"] != 0]
    say(f"EUR variants matched to a 1000G EUR rsID with consistent alleles: {m['rsid'].nunique():,}")
    cov = m["rsid"].isin(hm3)
    say(f"HapMap3 regression SNPs covered by EUR GWAS: {m.loc[cov, 'rsid'].nunique():,} of {len(hm3):,} "
        f"({m.loc[cov, 'rsid'].nunique() / len(hm3):.1%})")
    return m


def eas_freq_check(m):
    """Compare GWAS effect-allele frequency with 1000G EAS and EUR frequencies (ancestry sanity check)."""
    def frq(pattern):
        f = pd.concat(pd.read_csv(p, sep=r"\s+") for p in sorted(REF.rglob(pattern)))
        return f[["SNP", "A1", "MAF"]].rename(columns={"SNP": "rsid", "A1": "frq_a1", "MAF": "frq_maf"})
    # effect allele = g1; convert the GWAS frequency to the reference A1 allele using the allele match
    m = m[~m["ambiguous"]]  # orientation of A/T and C/G SNPs cannot be checked
    gwas_f_a1 = np.where(m["am"] == 1, m["effect_allele_frequency"], 1 - m["effect_allele_frequency"])
    g = pd.DataFrame({"rsid": m["rsid"].values, "a1": m["a1"].str.upper().values, "gf": gwas_f_a1,
                      "hm3": m["rsid"].isin(HM3).values})
    for anc, pattern in [("EAS", "1000G.EAS.QC.*.frq"), ("EUR", "1000G.EUR.QC.*.frq")]:
        x = g.merge(frq(pattern), on="rsid")
        ref_f = np.where(x["frq_a1"].str.upper() == x["a1"], x["frq_maf"], 1 - x["frq_maf"])
        h = x["hm3"].values
        say(f"EAS GWAS allele frequency vs 1000G {anc} (non-ambiguous SNPs): all r = "
            f"{np.corrcoef(x['gf'], ref_f)[0, 1]:.3f} (n = {len(x):,}); HapMap3 regression SNPs r = "
            f"{np.corrcoef(x['gf'][h], ref_f[h])[0, 1]:.3f} (n = {h.sum():,})")


def atlas_checks():
    """Read cell metadata only (no expression) from the AnnData file (old 0.1.0 layout:
    category labels in obs/__categories, integer codes in obs/<field>)."""
    path = RAW / "atlas/Full_obj_raw_counts_nosoupx_v2.h5ad"
    with h5py.File(path, "r") as h:
        shape = tuple(h["X"].attrs["shape"])
        obs = h["obs"]
        cats = obs["__categories"]
        o = {}
        for k in ["Sample name", "Diagnosis", "Age", "Age_group", "Gender", "Region", "Region code",
                  "Fraction", "10X", "category", "Integrated_05", "predicted_doublets"]:
            v = obs[k][:]
            if k in cats:
                labels = np.array([c.decode() if isinstance(c, bytes) else c for c in cats[k][:]], dtype=object)
                v = np.where(v >= 0, labels[np.clip(v, 0, None)], None)
            o[k] = v
    o = pd.DataFrame(o)
    say(f"atlas: {shape[0]:,} cells x {shape[1]:,} genes; donors ('Sample name'): {o['Sample name'].nunique()}")
    don = o.groupby("Sample name").agg(diagnosis=("Diagnosis", "first"), n_diag=("Diagnosis", "nunique"),
                                       age_group=("Age_group", "first"), age=("Age", "first"),
                                       sex=("Gender", "first"), cells=("Diagnosis", "size"))
    say(f"donors with >1 diagnosis label: {(don['n_diag'] > 1).sum()}")
    say("donors by diagnosis: " + "; ".join(f"{k}: {v}" for k, v in don["diagnosis"].value_counts().items()))
    say("cells by diagnosis: " + "; ".join(f"{k}: {v:,}" for k, v in o["Diagnosis"].value_counts().items()))
    post = o[o["Diagnosis"].isin(["Healthy adult", "Pediatric healthy"])]
    say(f"postnatal healthy cells: {len(post):,} from {post['Sample name'].nunique()} donors; regions: "
        + "; ".join(f"{k}: {v:,}" for k, v in post["Region"].value_counts().items()))
    lin = post.groupby("category").agg(cells=("Sample name", "size"), donors=("Sample name", "nunique"))
    say("postnatal healthy cells by lineage (cells / donors): "
        + "; ".join(f"{k}: {r.cells:,}/{r.donors}" for k, r in lin.iterrows()))
    fine = post.groupby("Integrated_05").agg(cells=("Sample name", "size"), donors=("Sample name", "nunique"))
    say(f"fine cell types among postnatal healthy cells: {len(fine)}; with >=100 cells and >=3 donors: "
        f"{((fine['cells'] >= 100) & (fine['donors'] >= 3)).sum()}")
    don.to_csv(OUT / "atlas_donors.tsv", sep="\t")
    fine.sort_values("cells", ascending=False).to_csv(OUT / "atlas_celltypes_postnatal_healthy.tsv", sep="\t")
    return o


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if sys.argv[1:] == ["--atlas-only"]:
        atlas_checks()
        (OUT / "qc_summary_atlas.txt").write_text("\n".join(summary) + "\n")
        return
    unpack()
    global HM3
    hm3 = HM3 = set(pd.read_csv(RAW / "ldsc_ref/hm3_no_MHC.list.txt", header=None)[0])
    say(f"HapMap3 regression SNP list (no MHC): {len(hm3):,}")
    m_eas = eas_checks(hm3)
    eas_freq_check(m_eas)
    eur_checks(hm3)
    atlas_checks()
    (OUT / "qc_summary.txt").write_text("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()
