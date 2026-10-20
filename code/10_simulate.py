"""Simulator: GWAS summary statistics with realistic LD and planted (or no) enrichment,
written in the real files' formats so that Steps C-E run unchanged.

Generative model (per ancestry, using that ancestry's 1000G Phase 3 reference genotypes):
- Causal SNPs: a random 2% of reference SNPs with MAF >= 5% (per replicate). Effect of standardized genotype
  beta_j ~ N(0, s2_j), s2_j proportional to 1 + e_ctrl * control_j + sum_L e_L * lineage_L_j, scaled so sum s2_j = h2.
- LD is approximated as block-diagonal with 2 cM blocks. Within block b, with standardized reference genotypes X
  (n individuals) and u_b ~ N(0, I_n):   z_j = sqrt(N_j) * x_j' (X_c beta_c) / n  +  x_j' u_b / sqrt(n)
  for each HapMap3 regression SNP j, so z has the reference sample's LD and E[z_j^2] = 1 + N_j * sum_k r_jk^2 s2_k.
- N_j: EAS per-variant effective N from the real Direction strings (coverage metadata, not outcomes); EUR 34,005.
Only identifier and Direction columns of the real GWAS files are read. Effect/SE columns are never read.

Usage: .venv/Scripts/python code/10_simulate.py {EAS,EUR} SCENARIO.json OUTDIR --reps 1 --seed 1
SCENARIO.json: {"h2": 0.1, "e_ctrl": 1.0, "e": {"Myeloid": 2.0, ...}}  (lineages absent from "e" get 0)
"""
import argparse
import json
import pathlib

import numpy as np
import pandas as pd

from common import ANC, ANNOT, CHROMS, LINEAGES, ROOT, slug

HM3 = ROOT / "data/raw/ldsc_ref/hm3_no_MHC.list.txt"
EAS_MATCHED = ROOT / "data/processed/eas_variants_matched_hg19.tsv.gz"
EUR_RAW = ROOT / "data/raw/gwas/cd_build37_40266_20161107.txt.gz"
EAS_COHORTS = [(2552, 6279), (1611, 3724), (1619, 4419), (1590, 1034)]
EUR_NEFF = 4 / (1 / 12194 + 1 / 28072)
BLOCK_CM, CAUSAL_FRAC = 2.0, 0.02


def read_bed(prefix, rows, n):
    """Genotype dosage (int8, 0/1/2 of the bim A1 allele; missing -> -1) for the given bim row indices."""
    per = (n + 3) // 4
    raw = np.memmap(prefix + ".bed", dtype=np.uint8, mode="r", offset=3).reshape(-1, per)[rows]
    codes = np.stack([(raw >> s) & 3 for s in (0, 2, 4, 6)], axis=2).reshape(len(rows), -1)[:, :n]
    lut = np.array([2, -1, 1, 0], dtype=np.int8)  # 00 hom A1, 01 missing, 10 het, 11 hom A2
    return lut[codes]


def standardize(g):
    g = g.astype(np.float32)
    miss = g < 0
    if miss.any():
        g[miss] = np.nan
        mu = np.nanmean(g, axis=1, keepdims=True)
        g = np.where(miss, mu, g)
    mu = g.mean(axis=1, keepdims=True)
    sd = g.std(axis=1, keepdims=True)
    sd[sd == 0] = 1
    return (g - mu) / sd


class Reference:
    """Per-chromosome blocks of reference genotypes for regression SNPs and candidate causal SNPs."""

    def __init__(self, anc, regression_snps):
        self.anc, self.blocks = anc, []
        lin_cols = [f"L_{slug(L)}_top10" for L in LINEAGES]
        for c in CHROMS:
            prefix = str(ANC[anc]["bfile"]).format(chr=c)
            bim = pd.read_csv(prefix + ".bim", sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
            n = len(pd.read_csv(prefix + ".fam", sep=r"\s+", header=None))
            frq = pd.read_csv(str(ANC[anc]["frq"]).format(chr=c) + ".frq", sep=r"\s+")
            maf = bim[["SNP"]].merge(frq[["SNP", "MAF"]], on="SNP", how="left")["MAF"].values
            ann = {"control": pd.read_csv(ANNOT / anc / "control_100kb" / f"control_100kb.{c}.annot.gz",
                                          sep="\t", usecols=["control_100kb"]).values[:, 0]}
            for L, col in zip(LINEAGES, lin_cols):
                ann[L] = pd.read_csv(ANNOT / anc / col / f"{col}.{c}.annot.gz", sep="\t", usecols=[col]).values[:, 0]
            ann = pd.DataFrame(ann)
            common = np.flatnonzero(maf >= 0.05)
            reg = np.flatnonzero(bim["SNP"].isin(regression_snps).values)
            block_id = np.floor((bim["CM"].values - bim["CM"].min()) / BLOCK_CM).astype(int)
            for b in np.unique(block_id[reg]):
                r_idx = reg[block_id[reg] == b]
                c_idx = common[block_id[common] == b]
                self.blocks.append({"chr": c, "reg_snps": bim["SNP"].values[r_idx],
                                    "Xr": standardize(read_bed(prefix, r_idx, n)).astype(np.float16),  # memory
                                    "common_idx": c_idx, "common_ann": ann.iloc[c_idx].reset_index(drop=True),
                                    "prefix": prefix, "n": n})
        self.n = self.blocks[0]["n"]

    def simulate(self, scen, N_of, rng):
        """Return DataFrame SNP, Z for all regression SNPs."""
        e = scen.get("e", {})
        # per-SNP prior variance weights for candidate causal SNPs, then scale to h2 over the chosen causal set
        draws = []
        for blk in self.blocks:
            a = blk["common_ann"]
            w = 1 + scen.get("e_ctrl", 0) * a["control"].values + sum(e.get(L, 0) * a[L].values for L in LINEAGES)
            pick = rng.random(len(w)) < CAUSAL_FRAC
            draws.append((pick, w))
        total_w = sum(w[p].sum() for p, w in draws)
        scale = scen["h2"] / total_w
        out_snp, out_z = [], []
        for blk, (pick, w) in zip(self.blocks, draws):
            Xr, n = blk["Xr"].astype(np.float32), self.n
            u = rng.standard_normal(n).astype(np.float32)
            noise = Xr @ u / np.sqrt(n)
            if pick.any():
                Xc = standardize(read_bed(blk["prefix"], blk["common_idx"][pick], n))
                beta = rng.standard_normal(pick.sum()) * np.sqrt(w[pick] * scale)
                g = Xc.T @ beta.astype(np.float32)          # (n,) genetic values in the reference sample
                signal = Xr @ g / n                          # x_j' X_c beta / n = sum_k r_jk beta_k
            else:
                signal = np.zeros(len(Xr), np.float32)
            N = N_of(blk["reg_snps"])
            out_snp.append(blk["reg_snps"])
            out_z.append(np.sqrt(N) * signal + noise)
        return pd.DataFrame({"SNP": np.concatenate(out_snp), "Z": np.concatenate(out_z)})


def eas_rows(hm3):
    m = pd.read_csv(EAS_MATCHED, sep="\t", dtype={"chromosome": str})
    m = m[m["rsid"].isin(hm3) & ~m["ambiguous"] & ~m["rsid"].duplicated(keep=False)]
    w = np.array([4 / (1 / a + 1 / b) for a, b in EAS_COHORTS])
    m["N"] = np.array([[ch != "?" for ch in s] for s in m["Direction"].values], dtype=float) @ w
    return m


def eur_rows(hm3):
    d = pd.read_csv(EUR_RAW, sep="\t", usecols=["MarkerName", "Allele1", "Allele2", "Direction"])
    cp = d["MarkerName"].str.extract(r"^(\w+):(\d+)_")
    d["chr"], d["pos"] = cp[0], cp[1].astype(int)
    bim = pd.concat(pd.read_csv(str(ANC["EUR"]["bfile"]).format(chr=c) + ".bim", sep=r"\s+", header=None,
                                names=["chr", "rsid", "cm", "pos", "a1", "a2"], dtype={"chr": str}) for c in CHROMS)
    x = d.merge(bim, on=["chr", "pos"])
    x = x[x["rsid"].isin(hm3) & ~x["rsid"].duplicated(keep=False)]
    return x


def main():
    p = argparse.ArgumentParser()
    p.add_argument("anc", choices=["EAS", "EUR"])
    p.add_argument("scenario")
    p.add_argument("outdir")
    p.add_argument("--reps", type=int, default=1)
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()
    scen = json.loads(pathlib.Path(a.scenario).read_text())
    out = pathlib.Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    hm3 = set(pd.read_csv(HM3, header=None)[0])
    rows = eas_rows(hm3) if a.anc == "EAS" else eur_rows(hm3)
    if a.anc == "EAS":
        n_map = rows.set_index("rsid")["N"]
        N_of = lambda snps: n_map.reindex(snps).fillna(n_map.median()).values
    else:
        N_of = lambda snps: np.full(len(snps), EUR_NEFF)
    ref = Reference(a.anc, set(rows["rsid"]))
    rng = np.random.default_rng(a.seed)
    for r in range(a.reps):
        sim = ref.simulate(scen, N_of, rng).set_index("SNP")["Z"]
        if a.anc == "EAS":
            x = rows.assign(z=rows["rsid"].map(sim)).dropna(subset=["z"])
            se = 1 / np.sqrt(x["N"].values)
            df = pd.DataFrame({"chromosome": x["chromosome"], "base_pair_location": x["base_pair_location"],
                               "effect_allele": x["g1"], "other_allele": x["g2"], "beta": x["z"] * se,
                               "standard_error": se, "Direction": x["Direction"]})
            df.to_csv(out / f"sim_EAS_rep{r}.tsv", sep="\t", index=False)
        else:
            x = rows.assign(z=rows["rsid"].map(sim)).dropna(subset=["z"])
            se = 1 / np.sqrt(EUR_NEFF)
            df = pd.DataFrame({"MarkerName": x["MarkerName"], "Allele1": x["Allele1"], "Allele2": x["Allele2"],
                               "Effect": x["z"] * se, "StdErr": se, "Direction": x["Direction"]})
            df.to_csv(out / f"sim_EUR_rep{r}.tsv.gz", sep="\t", index=False)
        print(f"{a.anc} rep {r}: {len(df):,} SNPs, mean chi2 {np.mean(x['z'] ** 2):.3f}", flush=True)


if __name__ == "__main__":
    main()
