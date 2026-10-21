"""Build the Markdown tables of the results report directly from the result files, so every number in
docs/results.md traces to results/main_analysis/. Writes results/main_analysis/report_tables.md.
Run from the repo root:  .venv/Scripts/python code/19_report_tables.py
"""
import json
import pathlib
import re

import numpy as np
import pandas as pd
from scipy import stats

RES = pathlib.Path("results/main_analysis")
E = RES / "exploratory"
L8 = ["B cells", "Endothelial", "Epithelial", "Mesenchymal", "Myeloid", "Neuronal", "Plasma cells", "T cells"]
md = []


def table(df):
    md.append("| " + " | ".join(df.columns) + " |")
    md.append("|" + "---|" * len(df.columns))
    for _, r in df.iterrows():
        md.append("| " + " | ".join(str(v) for v in r.values) + " |")
    md.append("")


def f(x, d=2):
    return f"{x:+.{d}f}"


def main():
    p = json.load(open(RES / "primary_results.json"))
    t = pd.read_csv(RES / "primary_table.tsv", sep="\t")
    md.append("### T1. Run record and LDSC diagnostics (results/main_analysis/primary_results.json; LDSC logs)")
    rows = []
    for k, dg in p["diagnostics"].items():
        anc = "EAS" if k == "EAS" else "EUR"
        log = open(RES / f"fits_{k}" / f"{anc}.L_B_cells_top10.log").read()
        m = re.search(r"Removed (\d+) SNPs with chi\^2 > \S+ \((\d+) SNPs remain\)", log)
        used = f"{int(m.group(2)):,} ({int(m.group(1))} removed, chi2 > 80)" if m else "see log"
        rows.append([k, f"{p['regression_snps'][k]:,}", used, f"{dg['h2_obs']:.3f} ({dg['h2_se']:.3f})",
                     f"{dg['h2_obs'] / dg['h2_se']:.1f}", f"{dg['mean_chi2']:.3f}", f"{dg['lambda_gc']:.3f}",
                     f"{dg['intercept']:.3f} ({dg['intercept_se']:.3f})"])
    table(pd.DataFrame(rows, columns=["dataset", "SNPs after munge", "SNPs used by LDSC", "h2 observed, effective-N scale (SE)", "h2 z",
                                      "mean chi2", "lambda GC", "LDSC intercept (SE)"]))

    md.append("### T2. Primary analysis: tau* per lineage and the East Asian minus European difference "
              "(results/main_analysis/primary_table.tsv)")
    rows = []
    for _, r in t.iterrows():
        lo, hi = r.D_EAS_minus_EUR - 1.96 * r.se_D, r.D_EAS_minus_EUR + 1.96 * r.se_D
        zs = stats.norm.ppf(1 - 0.05 / (2 * len(t)))  # Bonferroni-simultaneous 95% over the 8 lineages
        slo, shi = r.D_EAS_minus_EUR - zs * r.se_D, r.D_EAS_minus_EUR + zs * r.se_D
        rows.append([r.lineage, f"{f(r.tau_star_EAS)} ({r.se_EAS:.2f})", f"{f(r.tau_star_EUR)} ({r.se_EUR:.2f})",
                     f"{f(r.tau_star_FinnGen)} ({r.se_FinnGen:.2f})", f(r.D_EAS_minus_EUR),
                     f"{lo:+.2f} to {hi:+.2f}", f"{slo:+.2f} to {shi:+.2f}", f"{r.z:+.2f}", f"{r.p:.2f}",
                     f"{r.p_holm:.2f}"])
    table(pd.DataFrame(rows, columns=["lineage", "tau* EAS (SE)", "tau* EUR de Lange (SE)", "tau* EUR FinnGen (SE)",
                                      "D = EAS - EUR", "95% CI of D", "simultaneous 95% CI (Bonferroni, 8)", "z", "p",
                                      "p Holm"]))

    md.append("### T3. Joint Wald tests, primary and validation (results/main_analysis/primary_results.json)")
    rows = [[n, f"{p[k]['W']:.2f}", p[k]["K"], f"{p[k]['p']:.3f}"] for n, k in
            [("Primary: EAS vs EUR (de Lange)", "primary_EAS_vs_EUR"),
             ("Validation (i): EUR de Lange vs EUR FinnGen", "validation_i_deLange_vs_FinnGen"),
             ("Validation (ii): EAS vs EUR FinnGen", "validation_ii_EAS_vs_FinnGen")]]
    table(pd.DataFrame(rows, columns=["test", "W", "df", "p"]))

    s = json.load(open(RES / "sensitivity.json"))
    md.append("### T4. Registered sensitivity analyses (results/main_analysis/sensitivity.json)")
    rows = [[k, f"{v['wald']['W']:.2f}", v["wald"]["K"], f"{v['wald']['p']:.3f}",
             f"{max(abs(z) for z in v['wald']['z']):.2f}", f"{min(v['wald']['p_holm']):.2f}", v.get("note", "")]
            for k, v in s.items()]
    table(pd.DataFrame(rows, columns=["analysis", "W", "df", "p", "largest per-lineage |z|", "smallest p Holm", "note"]))

    md.append("### T4b. Per-lineage tau* in each sensitivity analysis, EAS / EUR (results/main_analysis/sensitivity.json)")
    rows = [["primary (registered)"] + [f"{a:+.2f} / {b:+.2f}" for a, b in
                                        zip(p["tau_star"]["EAS"], p["tau_star"]["EUR_deLange"])]]
    rows += [[k] + [f"{a:+.2f} / {b:+.2f}" for a, b in zip(v["tau_star_EAS"], v["tau_star_EUR"])] for k, v in s.items()]
    table(pd.DataFrame(rows, columns=["analysis"] + L8))

    md.append("### T5. Exploratory analyses (results/main_analysis/exploratory/) — EXPLORATORY")
    e2 = json.load(open(E / "E2_exclude_NOD2.json"))
    e3 = json.load(open(E / "E3_ulcerative_colitis.json"))
    e4 = json.load(open(E / "E4_region_specific.json"))
    e5 = json.load(open(E / "E5_bootstrap_genesets.json"))
    e1 = json.load(open(E / "E1_fine_types.json"))
    rows = [["E2 primary model without +/-1 Mb around NOD2", f"{e2['wald']['W']:.2f}", e2["wald"]["K"],
             f"{e2['wald']['p']:.3f}", f"region {e2['region_hg19']}; SNPs {e2['snps_before_after']}"],
            ["E3 ulcerative colitis, EAS vs EUR", f"{e3['wald']['W']:.2f}", e3["wald"]["K"], f"{e3['wald']['p']:.3f}",
             f"regression SNPs prepared: EAS {e3['info']['EAS_regression_snps']:,}, EUR {e3['info']['EUR_regression_snps']:,}"]]
    for reg, v in e4.items():
        rows.append([f"E4 {'small' if reg == 'SI' else 'large'}-intestine-only reference", f"{v['wald']['W']:.2f}",
                     v["wald"]["K"], f"{v['wald']['p']:.3f}", "neuronal lineage below threshold, dropped"])
    ps = [b["p"] for b in e5]
    rows.append(["E5 10 donor-bootstrap gene sets", "-", 8, f"{min(ps):.3f} to {max(ps):.3f}",
                 "sets of 1,414-1,489 genes after coordinate matching; SD of D across bootstrap sets: " + ", ".join(f"{x:.3f}" for x in
                                                              np.array([b['D'] for b in e5]).std(0, ddof=1))])
    rows.append(["E1 51 fine cell types: correlation of tau* EAS vs EUR", "-", "-", "not reported",
                 f"Pearson r = {e1['pearson'][0]:.2f}, Spearman rho = {e1['spearman'][0]:.2f} (gene sets overlap, so "
                 "fine types are not independent and no p-value is given)"])
    table(pd.DataFrame(rows, columns=["exploratory analysis", "W", "df", "p", "detail"]))

    md.append("### T6. E3 ulcerative colitis tau* per lineage (exploratory)")
    rows = [[L, f"{f(e3['tau_star']['EAS'][i])} ({e3['se']['EAS'][i]:.2f})",
             f"{f(e3['tau_star']['EUR'][i])} ({e3['se']['EUR'][i]:.2f})"] for i, L in enumerate(L8)]
    table(pd.DataFrame(rows, columns=["lineage", "UC tau* EAS (SE)", "UC tau* EUR (SE)"]))

    x, y = np.array(e1["tau_star_EAS"]), np.array(e1["tau_star_EUR"])
    sx, sy = np.array(e1["se_EAS"]), np.array(e1["se_EUR"])
    z = (x - y) / np.sqrt(sx ** 2 + sy ** 2)
    order = np.argsort(-(x + y))[:8]
    md.append("### T7. E1 fine cell types with the highest average tau* (exploratory; 51 types)")
    rows = [[e1["names"][i][2:].replace("_", " "), f"{f(x[i])} ({sx[i]:.2f})", f"{f(y[i])} ({sy[i]:.2f})",
             f"{z[i]:+.2f}"] for i in order]
    table(pd.DataFrame(rows, columns=["fine type", "tau* EAS (SE)", "tau* EUR (SE)", "z of difference"]))
    md.append(f"Largest |z| of an EAS-EUR difference across all 51 fine types: {np.abs(z).max():.2f}.")
    # Post hoc follow-ups after independent review (2026-10-09)
    fx = RES / "donorfix" / "primary_results.json"
    if fx.exists():
        q = json.load(open(fx))
        md.append("### T8. POST HOC: donor-corrected rerun (15 true donors; results/main_analysis/donorfix/) — "
                  "reported alongside, not replacing, the registered result")
        rows = [[n, f"{q[k]['W']:.2f}", q[k]["K"], f"{q[k]['p']:.3f}", f"{p[k]['W']:.2f}", f"{p[k]['p']:.3f}"]
                for n, k in [("Primary: EAS vs EUR (de Lange)", "primary_EAS_vs_EUR"),
                             ("Validation (i): de Lange vs FinnGen", "validation_i_deLange_vs_FinnGen"),
                             ("Validation (ii): EAS vs FinnGen", "validation_ii_EAS_vs_FinnGen")]]
        table(pd.DataFrame(rows, columns=["test", "W corrected", "df", "p corrected", "W registered", "p registered"]))
        qp = q["primary_EAS_vs_EUR"]
        rows = [[L, f"{f(q['tau_star']['EAS'][i])}", f"{f(q['tau_star']['EUR_deLange'][i])}",
                 f"{f(q['tau_star']['EUR_FinnGen'][i])}", f"{f(qp['D'][i])}", f"{qp['se_D'][i]:.2f}",
                 f"{qp['p_holm'][i]:.2f}"] for i, L in enumerate(L8)]
        table(pd.DataFrame(rows, columns=["lineage", "tau* EAS", "tau* EUR de Lange", "tau* EUR FinnGen", "D", "se(D)",
                                          "p Holm"]))
    al = RES / "aligned_block_check.json"
    if al.exists():
        a = json.load(open(al))
        md.append(f"### T9. POST HOC: aligned-block check of Cov(D) = Cov(A) + Cov(B) "
                  f"(results/main_analysis/aligned_block_check.json; {a['common_snps']:,} common SNPs)")
        rows = []
        for k, lab in [("EAS_vs_EUR_deLange", "EAS vs EUR (de Lange)"), ("EUR_deLange_vs_EUR_FinnGen",
                                                                        "de Lange vs FinnGen"),
                       ("EAS_vs_EUR_FinnGen", "EAS vs FinnGen")]:
            s_, d_ = a[k]["summed_covariance"], a[k]["direct_aligned"]
            ratio = np.array(d_["se_D"]) / np.array(s_["se_D"])
            rows.append([lab, f"{s_['W']:.2f}", f"{s_['p']:.3f}", f"{d_['W']:.2f}", f"{d_['p']:.3f}",
                         f"{min(d_['pseudovalue_corr']):.2f} to {max(d_['pseudovalue_corr']):.2f}",
                         f"{ratio.min():.2f} to {ratio.max():.2f}"])
        table(pd.DataFrame(rows, columns=["comparison", "W summed (registered formula)", "p summed", "W direct (aligned)",
                                          "p direct", "error correlation across lineages", "se(D) direct / summed"]))
        d_ = a["EAS_vs_EUR_deLange"]["direct_aligned"]
        zs = stats.norm.ppf(1 - 0.05 / 16)
        rows = [[L, f"{f(d_['D'][i])}", f"{d_['D'][i] - 1.96 * d_['se_D'][i]:+.2f} to {d_['D'][i] + 1.96 * d_['se_D'][i]:+.2f}",
                 f"{d_['D'][i] - zs * d_['se_D'][i]:+.2f} to {d_['D'][i] + zs * d_['se_D'][i]:+.2f}"]
                for i, L in enumerate(L8)]
        md.append("EAS vs EUR on the common SNPs, direct (aligned-block) intervals (post hoc):")
        md.append("")
        table(pd.DataFrame(rows, columns=["lineage", "D (common SNPs)", "95% CI", "simultaneous 95% CI (Bonferroni, 8)"]))
    (RES / "report_tables.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
