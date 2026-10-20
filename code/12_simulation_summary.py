"""Summarize the pipeline-validation simulations.

Reads results/simulations/sims/<scenario>.jsonl and writes results/simulations/simulation_summary.txt plus tables:
- false-positive rate of the primary Wald test where H0 holds (null, shared), with exact binomial 95% CI, and a
  Kolmogorov-Smirnov test of p-value uniformity;
- the same for the 8 random gene sets (shuffled-label control);
- mean tau* per lineage and ancestry per scenario (does the pipeline recover planted enrichment, and only there?);
- calibration of the jackknife: empirical SD across replicates of D = tau*_EAS - tau*_EUR vs the mean jackknife SE;
- power by planted difference (diff scenarios) and the minimum detectable myeloid difference at 80% power
  (linear interpolation of power against the mean observed myeloid difference).
Run from the repo root:  .venv/Scripts/python code/12_simulation_summary.py
"""
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

RES = pathlib.Path("results/simulations")
SHORT = ["B", "Endo", "Epi", "Mes", "Myel", "Neur", "Plas", "T"]
out = []


def say(s=""):
    print(s)
    out.append(s)


def load(name):
    p = RES / "sims" / f"{name}.jsonl"
    return [json.loads(l) for l in p.read_text().splitlines()] if p.exists() else []


def rate(ps, alpha=0.05):
    k, n = int((np.asarray(ps) < alpha).sum()), len(ps)
    ci = stats.binomtest(k, n).proportion_ci(method="exact") if n else (np.nan, np.nan)
    return k, n, ci


def main():
    say("PIPELINE-VALIDATION SIMULATIONS (registered doi:10.17605/OSF.IO/HQ4A8)")
    say("Simulated data only; no real GWAS association result is used.\n")
    rows = []
    for scen in ["null", "shared", "diff4", "diff6", "diff10"]:
        R = load(scen)
        if not R:
            continue
        ps = [r["primary"]["p"] for r in R]
        k, n, ci = rate(ps)
        ks = stats.kstest(ps, "uniform").pvalue
        D = np.array([r["primary"]["D"] for r in R])
        seD = np.array([r["primary"]["se_D"] for r in R])
        h2z = {a: np.mean([r["h2_z"][a] for r in R]) for a in ["EAS", "EUR"]}
        holm = np.array([r["primary"]["p_holm"] for r in R]) < 0.05
        myel_power, other_any = holm[:, 4].mean(), np.delete(holm, 4, axis=1).any(axis=1).mean()
        rows.append({"scenario": scen, "reps": n, "reject_at_0.05": k / n, "ci_low": ci[0], "ci_high": ci[1],
                     "holm_myeloid": myel_power, "holm_any_other": other_any,
                     "ks_uniform_p": ks, "mean_D_myeloid": D[:, 4].mean(), "mean_h2z_EAS": h2z["EAS"],
                     "mean_h2z_EUR": h2z["EUR"]})
        say(f"[{scen}] {n} replicates; mean LDSC h2 z: EAS {h2z['EAS']:.1f}, EUR {h2z['EUR']:.1f}")
        say(f"  primary Wald test rejects at 0.05: {k}/{n} = {k / n:.3f} (95% CI {ci[0]:.3f}-{ci[1]:.3f}); "
            f"KS test of p-value uniformity p = {ks:.2f}")
        say(f"  secondary per-lineage tests (Holm, 0.05): myeloid flagged in {myel_power:.2f} of reps; "
            f"any other lineage flagged in {other_any:.2f}")
        for anc in ["EAS", "EUR"]:
            T = np.array([r["tau_star"][anc] for r in R])
            say(f"  mean tau* {anc}: " + "  ".join(f"{s} {m:+.2f}" for s, m in zip(SHORT, T.mean(0))))
        say("  jackknife calibration, SD of D across reps / mean jackknife SE of D: " + "  ".join(
            f"{s} {sd / se:.2f}" for s, sd, se in zip(SHORT, D.std(0, ddof=1), seD.mean(0))))
        if "random" in R[0]:
            pr = [r["random"]["p"] for r in R]
            k2, n2, ci2 = rate(pr)
            say(f"  random gene sets (shuffled labels): Wald rejects {k2}/{n2} = {k2 / n2:.3f} "
                f"(95% CI {ci2[0]:.3f}-{ci2[1]:.3f}); KS p = {stats.kstest(pr, 'uniform').pvalue:.2f}")
            for anc in ["EAS", "EUR"]:
                T = np.array([r["tau_star_random"][anc] for r in R])
                z = T.mean(0) / (T.std(0, ddof=1) / np.sqrt(len(T)))
                say(f"  random sets mean tau* {anc}: " + "  ".join(f"R{i + 1} {m:+.2f}" for i, m in enumerate(T.mean(0)))
                    + f"   (largest |mean/SE| = {np.abs(z).max():.1f})")
        say()
    h0 = load("null") + load("shared")
    if h0:
        k, n, ci = rate([r["primary"]["p"] for r in h0])
        say(f"Pooled scenarios where H0 is true (null + shared): primary test rejects {k}/{n} = {k / n:.3f} "
            f"(95% CI {ci[0]:.3f}-{ci[1]:.3f})")
    df = pd.DataFrame(rows)
    df.to_csv(RES / "simulation_summary_table.tsv", sep="\t", index=False, float_format="%.4f")
    pw = df[df["scenario"].str.startswith("diff") | (df["scenario"] == "shared")].sort_values("mean_D_myeloid")
    if len(pw) >= 2:
        pw = pw.assign(absD=pw["mean_D_myeloid"].abs()).sort_values("absD")
        x, y = pw["absD"].values, pw["reject_at_0.05"].values
        say("Power of the primary test by mean myeloid difference |D|: " + "; ".join(
            f"{s}: |D| {a:.2f} -> power {b:.2f}" for s, a, b in zip(pw["scenario"], x, y)))
        if (y >= 0.8).any() and (y < 0.8).any():
            i = int(np.argmax(y >= 0.8))  # first planted size reaching 80%; interpolate from the one below
            mde = np.interp(0.8, y[i - 1:i + 1], x[i - 1:i + 1])
            say(f"Minimum detectable myeloid tau* difference at 80% power (interpolated): {mde:.2f}")
        elif (y >= 0.8).all():
            say("All planted differences reached 80% power; the minimum detectable difference is below "
                f"{x.min():.2f}.")
        else:
            say(f"No planted difference reached 80% power; the minimum detectable difference exceeds {x.max():.2f}.")
    (RES / "simulation_summary.txt").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
