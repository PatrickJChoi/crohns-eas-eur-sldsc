"""Claims audit for the manuscript (docs/manuscript.md).

Every quantitative claim is listed with the exact text it appears as in the manuscript, the source file, and code that
recomputes the value from that file. The audit FAILS if a claim's text is missing from the manuscript or the recomputed
value differs from the text. Qualitative claims are listed with their supporting result and a flag for overreach.
Writes docs/claims_audit.md. Run from the repo root:  .venv/Scripts/python code/27_claims_audit.py
"""
import json
import pathlib
import re

import numpy as np
import pandas as pd
from scipy import stats

MS = pathlib.Path("docs/manuscript.md").read_text(encoding="utf-8")
R5 = pathlib.Path("results/main_analysis")
P = json.load(open(R5 / "primary_results.json"))
T = pd.read_csv(R5 / "primary_table.tsv", sep="\t").set_index("lineage")
S = json.load(open(R5 / "sensitivity.json"))
E = {k: json.load(open(R5 / "exploratory" / f"{k}.json")) for k in
     ["E1_fine_types", "E2_exclude_NOD2", "E3_ulcerative_colitis", "E4_region_specific", "E5_bootstrap_genesets"]}
FIX = json.load(open(R5 / "donorfix" / "primary_results.json"))
AL = json.load(open(R5 / "aligned_block_check.json"))
SIM = {s: [json.loads(l) for l in open(f"results/simulations/sims/{s}.jsonl")] for s in
       ["null", "shared", "diff4", "diff6", "diff10"]}
GS = pathlib.Path("results/simulations/genesets_summary.txt").read_text(encoding="utf-8")
MAN = pd.read_csv("data/raw/finngen_r13/finngen_R13_manifest.tsv", sep="\t").set_index("phenocode")
CODE10 = pathlib.Path("code/07_prepare_sumstats.py").read_text(encoding="utf-8")
ZS = stats.norm.ppf(1 - 0.05 / 16)


def f2(x):
    return f"{x:+.2f}".replace("-", "−") if x < 0 else f"+{x:.2f}"


def ldsc_used(ds):
    anc = "EAS" if ds == "EAS" else "EUR"
    log = (R5 / f"fits_{ds}" / f"{anc}.L_B_cells_top10.log").read_text()
    return int(re.search(r"\((\d+) SNPs remain\)", log).group(1))


ci = lambda L: (T.loc[L, "D_EAS_minus_EUR"] - 1.96 * T.loc[L, "se_D"], T.loc[L, "D_EAS_minus_EUR"] + 1.96 * T.loc[L, "se_D"])
sci = lambda L: (T.loc[L, "D_EAS_minus_EUR"] - ZS * T.loc[L, "se_D"], T.loc[L, "D_EAS_minus_EUR"] + ZS * T.loc[L, "se_D"])
rate = lambda R, key="p": np.mean([r["primary"][key] < 0.05 for r in R])
D_al, se_al = np.array(AL["EAS_vs_EUR_deLange"]["direct_aligned"]["D"]), np.array(AL["EAS_vs_EUR_deLange"]["direct_aligned"]["se_D"])
se_sum = np.array(AL["EAS_vs_EUR_deLange"]["summed_covariance"]["se_D"])
h0 = SIM["null"] + SIM["shared"]
k0 = sum(r["primary"]["p"] < 0.05 for r in h0)
ci0 = stats.binomtest(k0, len(h0)).proportion_ci(method="exact")
sens_snps = S["a_drop_ImmunoChip_SNPs"]["note"]
eas_n = [int(x.replace(",", "")) for x in re.findall(r"EAS: ([\d,]+) -> ([\d,]+)", sens_snps)[0]]
eur_n = [int(x.replace(",", "")) for x in re.findall(r"EUR_deLange: ([\d,]+) -> ([\d,]+)", sens_snps)[0]]
fin_ratio = (T["se_FinnGen"] / T["se_EUR"]).mean()

# (claim id, text as it appears in the manuscript, recomputed text, source)
CLAIMS = [
    ("Q1", "7,372 cases, 15,456 controls", "{:,} cases, {:,} controls".format(
        sum(int(a) for a, b in re.findall(r'"\w+": \((\d+), (\d+)\)', CODE10.split("EAS_COHORTS =")[1].split("\n")[0])),
        sum(int(b) for a, b in re.findall(r'"\w+": \((\d+), (\d+)\)', CODE10.split("EAS_COHORTS =")[1].split("\n")[0]))),
     "code/07_prepare_sumstats.py EAS_COHORTS (per-cohort counts from Liu 2023 Methods)"),
    ("Q2", "12,194\ncases, 28,072 controls", "12,194\ncases, 28,072 controls" if "1 / 12194 + 1 / 28072" in CODE10 else "?",
     "code/07_prepare_sumstats.py EUR_NEFF (de Lange 2017)"),
    ("Q3", "2,608 cases, 497,350 controls", f"{MAN.loc['K11_CD_STRICT2', 'num_cases']:,} cases, "
                                           f"{MAN.loc['K11_CD_STRICT2', 'num_controls']:,} controls",
     "data/raw/finngen_r13/finngen_R13_manifest.tsv"),
    ("Q4", "W = 4.04, df = 8, p = 0.853", f"W = {P['primary_EAS_vs_EUR']['W']:.2f}, df = {P['primary_EAS_vs_EUR']['K']}, "
                                         f"p = {P['primary_EAS_vs_EUR']['p']:.3f}", "results/main_analysis/primary_results.json"),
    ("Q5", "all adjusted p = 1.00", "all adjusted p = 1.00" if (T["p_holm"].round(2) == 1).all() else "?",
     "results/main_analysis/primary_table.tsv p_holm"),
    ("Q6", "T cells (D = −0.15, 95% CI −0.40 to +0.10)",
     f"T cells (D = {f2(T.loc['T cells', 'D_EAS_minus_EUR'])}, 95% CI {f2(ci('T cells')[0])} to {f2(ci('T cells')[1])})",
     "primary_table.tsv"),
    ("Q7", "T cells (EAS +0.20, EUR +0.35)", f"T cells (EAS {f2(T.loc['T cells', 'tau_star_EAS'])}, EUR "
                                            f"{f2(T.loc['T cells', 'tau_star_EUR'])})", "primary_table.tsv"),
    ("Q8", "(EAS +0.23, EUR +0.20)", f"(EAS {f2(T.loc['Myeloid', 'tau_star_EAS'])}, EUR {f2(T.loc['Myeloid', 'tau_star_EUR'])})",
     "primary_table.tsv (myeloid)"),
    ("Q9", "spanned −0.40 to +0.27", f"spanned {f2(min(ci(L)[0] for L in T.index))} to {f2(max(ci(L)[1] for L in T.index))}",
     "primary_table.tsv (marginal 95% CIs)"),
    ("Q10", "spanned −0.50 to +0.36", f"spanned {f2(min(sci(L)[0] for L in T.index))} to "
                                     f"{f2(max(sci(L)[1] for L in T.index))}", "primary_table.tsv (Bonferroni 8)"),
    ("Q11", "(−0.50 to +0.21)", f"({f2(sci('T cells')[0])} to {f2(sci('T cells')[1])})", "primary_table.tsv (T cells, simultaneous)"),
    ("Q12", "W = 7.60, p = 0.473", f"W = {P['validation_i_deLange_vs_FinnGen']['W']:.2f}, p = "
                                   f"{P['validation_i_deLange_vs_FinnGen']['p']:.3f}", "primary_results.json validation (i)"),
    ("Q13", "W = 7.91, p = 0.442", f"W = {P['validation_ii_EAS_vs_FinnGen']['W']:.2f}, p = "
                                   f"{P['validation_ii_EAS_vs_FinnGen']['p']:.3f}", "primary_results.json validation (ii)"),
    ("Q14", "roughly twice as imprecise", "roughly twice as imprecise" if 1.6 <= fin_ratio <= 2.4 else f"ratio {fin_ratio:.2f}",
     f"primary_table.tsv mean se_FinnGen / se_EUR = {fin_ratio:.2f}"),
    ("Q15", "(p = 0.521 to 0.995", f"(p = {min(v['wald']['p'] for v in S.values()):.3f} to "
                                  f"{max(v['wald']['p'] for v in S.values()):.3f}", "results/main_analysis/sensitivity.json"),
    ("Q16", "(EAS +0.20 to −0.03; EUR +0.35 to +0.14)",
     f"(EAS {f2(T.loc['T cells', 'tau_star_EAS'])} to {f2(S['a_drop_ImmunoChip_SNPs']['tau_star_EAS'][7])}; EUR "
     f"{f2(T.loc['T cells', 'tau_star_EUR'])} to {f2(S['a_drop_ImmunoChip_SNPs']['tau_star_EUR'][7])})",
     "sensitivity.json a_drop_ImmunoChip_SNPs (T cells)"),
    ("Q17", "about 10% of regression SNPs", "about 10% of regression SNPs" if all(
        0.08 <= 1 - b / a <= 0.12 for a, b in [eas_n, eur_n]) else "?",
     f"sensitivity.json note: EAS {1 - eas_n[1] / eas_n[0]:.1%}, EUR {1 - eur_n[1] / eur_n[0]:.1%} removed"),
    ("Q18", "NOD2 region (p = 0.822)", f"NOD2 region (p = {E['E2_exclude_NOD2']['wald']['p']:.3f})", "exploratory/E2"),
    ("Q19", "UC (p = 0.632)", f"UC (p = {E['E3_ulcerative_colitis']['wald']['p']:.3f})", "exploratory/E3"),
    ("Q20", "small-intestine-only (p = 0.837)", f"small-intestine-only (p = {E['E4_region_specific']['SI']['wald']['p']:.3f})",
     "exploratory/E4 SI"),
    ("Q21", "large-intestine-only (p = 0.995)", f"large-intestine-only (p = {E['E4_region_specific']['LI']['wald']['p']:.3f})",
     "exploratory/E4 LI"),
    ("Q22", "bootstrap gene sets (p = 0.814 to 0.990)", "bootstrap gene sets (p = {:.3f} to {:.3f})".format(
        min(b["p"] for b in E["E5_bootstrap_genesets"]), max(b["p"] for b in E["E5_bootstrap_genesets"])), "exploratory/E5"),
    ("Q23", "(Pearson r = 0.73)", f"(Pearson r = {E['E1_fine_types']['pearson'][0]:.2f})", "exploratory/E1"),
    ("Q24", "|z| = 1.79", "|z| = {:.2f}".format(np.abs((np.array(E["E1_fine_types"]["tau_star_EAS"]) - np.array(
        E["E1_fine_types"]["tau_star_EUR"])) / np.sqrt(np.array(E["E1_fine_types"]["se_EAS"]) ** 2 + np.array(
        E["E1_fine_types"]["se_EUR"]) ** 2)).max()), "exploratory/E1"),
    ("Q25", "51 fine cell types", f"{len(E['E1_fine_types']['names'])} fine cell types", "exploratory/E1"),
    ("Q26", "(W = 2.16, p = 0.976)", f"(W = {FIX['primary_EAS_vs_EUR']['W']:.2f}, p = {FIX['primary_EAS_vs_EUR']['p']:.3f})",
     "results/main_analysis/donorfix/primary_results.json"),
    ("Q27", "On 950,331 SNPs", f"On {AL['common_snps']:,} SNPs", "results/main_analysis/aligned_block_check.json"),
    ("Q28", "(0.24 to 0.40 across lineages)", "({:.2f} to {:.2f} across lineages)".format(
        min(AL["EAS_vs_EUR_deLange"]["direct_aligned"]["pseudovalue_corr"]),
        max(AL["EAS_vs_EUR_deLange"]["direct_aligned"]["pseudovalue_corr"])), "aligned_block_check.json"),
    ("Q29", "gave W = 7.22, p = 0.513", "gave W = {:.2f}, p = {:.3f}".format(
        AL["EAS_vs_EUR_deLange"]["direct_aligned"]["W"], AL["EAS_vs_EUR_deLange"]["direct_aligned"]["p"]),
     "aligned_block_check.json"),
    ("Q30", "13–22% smaller", "{:.0f}–{:.0f}% smaller".format(100 * (1 - (se_al / se_sum).max()),
                                                             100 * (1 - (se_al / se_sum).min())), "aligned_block_check.json"),
    ("Q31", "intervals spanned −0.46 to +0.32", f"intervals spanned {f2((D_al - ZS * se_al).min())} to "
                                               f"{f2((D_al + ZS * se_al).max())}", "aligned_block_check.json"),
    ("Q32", "4 of 80 replicates (5.0%; 95% CI\n1.4–12.3%)", f"{k0} of {len(h0)} replicates ({k0 / len(h0):.1%}; 95% CI\n"
                                                           f"{100 * ci0[0]:.1f}–{100 * ci0[1]:.1f}%)",
     "results/simulations/sims/null.jsonl + shared.jsonl"),
    ("Q33", "no rejections in 60 replicates", "no rejections in {} replicates".format(len(SIM["shared"])) if
     sum(r["random"]["p"] < 0.05 for r in SIM["shared"]) == 0 else "?", "results/simulations/sims/shared.jsonl (random sets)"),
    ("Q34", "0%, 25% and 75% for\nmean differences of 0.25, 0.46 and 0.72",
     "{:.0f}%, {:.0f}% and {:.0f}% for\nmean differences of {:.2f}, {:.2f} and {:.2f}".format(
         *[100 * rate(SIM[s]) for s in ["diff4", "diff6", "diff10"]],
         *[abs(np.mean([r["primary"]["D"][4] for r in SIM[s]])) for s in ["diff4", "diff6", "diff10"]]),
     "results/simulations/sims/diff*.jsonl"),
    ("Q35", "myeloid lineage in 95% of replicates and another lineage in 5%",
     "myeloid lineage in {:.0f}% of replicates and another lineage in {:.0f}%".format(
         100 * np.mean([r["primary"]["p_holm"][4] < 0.05 for r in SIM["diff10"]]),
         100 * np.mean([any(p < 0.05 for i, p in enumerate(r["primary"]["p_holm"]) if i != 4) for r in SIM["diff10"]])),
     "results/simulations/sims/diff10.jsonl"),
    ("Q36", "LDSC used 986,470 SNPs for the EAS GWAS, 1,149,713 for the European GWAS and 1,175,811 for FinnGen",
     f"LDSC used {ldsc_used('EAS'):,} SNPs for the EAS GWAS, {ldsc_used('EUR_deLange'):,} for the European GWAS and "
     f"{ldsc_used('EUR_FinnGen'):,} for FinnGen", "results/main_analysis/fits_*/ *.log"),
    ("Q37", "0.847 (s.e. 0.070), 0.510 (0.036) and 0.477 (0.082)", "{:.3f} (s.e. {:.3f}), {:.3f} ({:.3f}) and {:.3f} ({:.3f})".format(
        *[P["diagnostics"][d][k] for d in ["EAS", "EUR_deLange", "EUR_FinnGen"] for k in ["h2_obs", "h2_se"]]),
     "primary_results.json diagnostics"),
    ("Q38", "intercepts\nwere 0.997, 1.070 and 1.024", "intercepts\nwere {:.3f}, {:.3f} and {:.3f}".format(
        *[P["diagnostics"][d]["intercept"] for d in ["EAS", "EUR_deLange", "EUR_FinnGen"]]), "primary_results.json diagnostics"),
    ("Q39", "(142,797 cells)", "({:,} cells)".format(int(re.search(r"cells kept: ([\d,]+)", GS).group(1).replace(",", ""))),
     "results/simulations/genesets_summary.txt"),
    ("Q40", "(14,927 genes)", "({:,} genes)".format(int(re.search(r"expressed \(mean CPM >= 1 in a lineage\): ([\d,]+)", GS)
                                                       .group(1).replace(",", ""))), "genesets_summary.txt"),
    ("Q41", "(1,493 genes)", "({:,} genes)".format(int(re.search(r"top10: (\d+) genes", GS).group(1))), "genesets_summary.txt"),
    ("Q42", "intercept of 1.070", f"intercept of {P['diagnostics']['EUR_deLange']['intercept']:.3f}", "primary_results.json"),
    ("Q43", "(−0.15 in the registered\nmodel, −0.17 without those SNPs)",
     f"({f2(T.loc['T cells', 'D_EAS_minus_EUR'])} in the registered\nmodel, "
     f"{f2(S['a_drop_ImmunoChip_SNPs']['wald']['D'][7])} without those SNPs)", "primary_table.tsv; sensitivity.json (a)"),
]

QUAL = [
    ("Q-A", "No evidence that CD heritability is distributed differently ... (title, abstract, results)", "Q4; registered "
     "interpretation rule (plan section 5): p >= 0.05 means no evidence, not equality", "supported; worded as no evidence"),
    ("Q-B", "Largest tau* in both ancestries for T cells and myeloid cells (registered model)", "Q7, Q8; Table 1",
     "descriptive estimate, not a registered test; depends on the analysis choice (T4b) — stated in Results"),
    ("Q-C", "Much of the T-cell signal lies in ImmunoChip-covered regions", "Q16, Q17", "interpretation of a sensitivity "
     "analysis; wording kept at 'much of' / 'a large part'"),
    ("Q-D", "The registered test is conservative", "Q28, Q29 (post hoc)", "supported by a post hoc check; labelled post hoc"),
    ("Q-E", "Differences as large as the lineage enrichments cannot be excluded", "Q10, Q11, Q7; registered minimum "
     "detectable difference > 0.72 (Q34)", "supported"),
    ("Q-F", "In line with recent single-cell studies (Zhao 2026; Krzak 2026)", "their abstracts (T cells, neutrophils; "
     "monocytes and macrophages)", "consistency statement only; different methods and data"),
    ("Q-G", "Consistent with Amariuta 2020; does not test or contradict Shi 2021", "abstracts of both papers",
     "framing; does not claim our result confirms either"),
]


def main():
    rows, fails = [], 0
    for cid, text, recomputed, src in CLAIMS:
        in_ms = text in MS
        match = text == recomputed
        ok = in_ms and match
        fails += not ok
        rows.append(f"| {cid} | {text.replace(chr(10), ' ')} | {src} | {'yes' if in_ms else 'NO'} | "
                    f"{'yes' if match else 'NO: ' + recomputed.replace(chr(10), ' ')} | {'PASS' if ok else 'FAIL'} |")
    md = ["# Claims audit — manuscript (generated by code/27_claims_audit.py)", "",
          f"Quantitative claims: {len(CLAIMS)}; failures: {fails}.", "",
          "| id | claim as written | source | text found in manuscript | value recomputed from source | status |",
          "|---|---|---|---|---|---|"] + rows + [
          "", "## Qualitative claims and overreach flags", "",
          "| id | claim | support | assessment |", "|---|---|---|---|"] + [
          f"| {a} | {b} | {c} | {d} |" for a, b, c, d in QUAL]
    pathlib.Path("docs/claims_audit.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"{len(CLAIMS)} quantitative claims, {fails} failures")
    for r in rows:
        if "FAIL" in r:
            print(r)


if __name__ == "__main__":
    main()
