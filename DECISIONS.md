# DECISIONS

Dated, one line per methodological decision: what was chosen and why.

## Software and data
- 2026-09-06 — Python 3.11.5; versions pinned in requirements-analysis.txt and requirements-ldsc.txt (see ENVIRONMENT.md).
- 2026-09-07 — FinnGen data freeze R13 (K11_CD_STRICT2, 2,608 cases / 497,350 controls, from finngen_R13_manifest.tsv) used for European validation; it is the latest release.
- 2026-09-07 — Blinding rule: no heritability, genetic correlation or per-variant association result is read before the plan is registered.
- 2026-09-07 — Data: EAS CD GWAS GCST90446793 (Liu 2023), EUR CD GWAS GCST004132 (de Lange 2017), 1000G Phase 3 EAS/EUR S-LDSC references (Zenodo 10515792; baselineLD v2.2 and baseline v1.2), Gut Cell Atlas full object, UCSC hg38ToHg19 chain.
- 2026-09-07 — GWAS Catalog's listed MD5 for GCST90446793.tsv-meta.yaml is stale (three downloads agree with each other, not the list); the data file's MD5 matches. Yaml recorded without a published checksum.
- 2026-09-07 — TLS for downloads verified against certifi's CA bundle (Windows store lacked an intermediate for UCSC); verification was never disabled.
- 2026-09-07 — LDSC: CBIIT/ldsc (Python 3 rewrite linked from bulik/ldsc README) pinned at commit 1f09cf0, in .venv-ldsc with numpy 1.23.5, pandas 1.5.3, scipy 1.9.3, bitarray 2.9.3 (nearest versions to the authors' pins with Windows wheels for Python 3.11).
- 2026-09-07 — Analysis environment .venv (numpy 2.1.3, pandas 2.2.3, scipy 1.17.1, h5py 3.16.0, pyliftover 0.4.1, certifi 2026.7.22); pins in requirements-analysis.txt.
- 2026-09-07 — Strand-ambiguous (A/T, C/G) SNPs excluded (15.5% of matched EAS variants).
- 2026-09-07 — Feasibility (docs/data_feasibility.md): data usable; the ImmunoChip cohort (ICH1) in the EAS file must be handled in the plan; no independent EAS validation cohort exists (stated as a limitation).

## Analysis plan (docs/analysis_plan.md)
- 2026-09-07 — Primary = one joint Wald test (df = K lineages, expected 8) of EAS−EUR differences in τ*, jackknife covariance; chosen so there is exactly one primary test and no lineage is picked by looking at data.
- 2026-09-07 — Lineage level (atlas `category`), not fine cell types: fine myeloid types are too sparse (e.g., 245 monocytes); fine types are exploratory.
- 2026-09-07 — ImmunoChip cohort handled by per-variant effective N (primary), the standard way LDSC handles varying N; symmetric exclusion of ICH1-contributed SNPs is a sensitivity analysis.
- 2026-09-07 — LDSC-SEG model with baseline v1.2 + all-genes control (reference readme recommends v1.2 for cell-type tests); baselineLD v2.2 is a sensitivity analysis.
- 2026-09-07 — Healthy postnatal gut cells only; lymph node excluded (not intestinal mucosa); pseudobulk per donor so donors, not cells, are the unit.
- 2026-09-07 — Validation: FinnGen R13 K11_CD_STRICT2 for the EUR side; no EAS validation (limitation).
- 2026-09-08 — Plan registered on OSF before any association result was read: doi:10.17605/OSF.IO/HQ4A8 (https://osf.io/hq4a8, 2026-09-08 07:16 UTC); project https://osf.io/e9mny. Registered files are byte-identical (SHA-256) to docs/osf_registration/. Deviations after this date are logged below.

## Pipeline validation (docs/pipeline_validation.md)
- 2026-09-08 — Built: code/04 (GENCODE v19, MD5-verified), 05 (Step A gene sets: K = 8 lineages, 1,493 genes each, from 14,927 expressed genes; bootstrap median Jaccard 0.66–0.87), 06 (Step B annotations/LD scores), 07–09 (Steps C–E, blinding guard on real GWAS), 10 (LD-aware simulator), tests/test_pipeline.py.
- 2026-09-08 — Step A lineage drop rule applied as: cells and donors counted within kept pseudobulks (>= 20 cells per donor-lineage). All 8 lineages pass (neuronal: 376 cells, 7 donors).
- 2026-09-08 — Step C format clarification (not an analysis change): this LDSC version's munge_sumstats requires a P column, so the prepared files add P = 2*norm.sf(|Z|) (scipy, accurate in the tails); munge then rebuilds chi-square from P and sign from Z, as LDSC does.
- 2026-09-08 — Simulations fit models without --overlap-annot: verified on a pilot that per-block delete values are bit-identical and full-data tau agrees to ~5 significant digits (log rounding); about 7x faster. The real-data primary run keeps --overlap-annot as registered.
- 2026-09-08 — Simulation calibration (pilot): simulated observed-scale h2 EAS 0.35, EUR 0.15 so that LDSC h2 z-scores approximate the published 7.9 (EAS) and 9.8 (EUR). Scenarios: null (no lineage enrichment, 20 reps); shared (T cells and myeloid enrichment e = 2 in both, 60 reps, plus 8 random gene sets as shuffled-label controls); diff4/diff6/diff10 (EUR myeloid e = 4/6/10 vs EAS 2, 20 reps each). Seeds: --seed 1 with a CRC32 of the scenario name.
- 2026-09-08 — Simulations complete (140 replicates). H0 rejection rate 4/80 = 5.0% (95% CI 1.4–12.3%); random gene sets 0/60; planted myeloid difference flagged in the right lineage; jackknife calibrated (0.93–1.08 in shared). Primary-test power for a one-lineage difference: 0%, 25%, 75% at mean |delta tau*| 0.25, 0.46, 0.72; minimum detectable difference at 80% power > 0.72.
- 2026-09-08 — Simulation runs: generate and fit stages separated after an out-of-memory error; 6 diff10 fits rerun with 3 workers after a second out-of-memory error (same simulated data, resumable driver; no partial results saved; 140/140 replicates recorded once each).

## Primary analysis (docs/results.md)
- 2026-09-08 — Primary run with code/13_primary_analysis.py (script committed before it was run; first and only unblinded run of the registered primary analysis).
- 2026-09-08 — PRIMARY RESULT (single registered run): joint Wald test EAS vs EUR, W = 4.04, K = 8, p = 0.853 -> no evidence of an ancestry difference (registered interpretation: not proof of equality). Validation: (i) de Lange vs FinnGen W = 7.60, p = 0.473; (ii) EAS vs FinnGen W = 7.91, p = 0.442. No per-lineage difference significant after Holm (all p_holm = 1). Per-lineage 95% CIs of D span about -0.40 to +0.27. Results: results/main_analysis/primary_results.json, primary_table.tsv.
- 2026-09-08 — Precision note (post hoc, approximate): real se(D) are 0.50–0.79 of the simulated ones, so the simulated minimum detectable difference (> 0.72 tau*) is conservative; scaled, about 0.5 for myeloid.
- 2026-09-09 — Registered sensitivity analyses done, all with --overlap-annot. Sensitivity (d) used +/-50 kb windows for both lineage and control annotations (control_50kb).
- 2026-09-09 — Exploratory analyses were fitted without --overlap-annot (coefficients and jackknife unchanged, verified in simulations; full-data tau from LDSC's log, ~5 significant digits).
- 2026-09-09 — Exploratory 1: registered text says "the 65 fine cell types"; applying the registered rule (>= 100 cells, >= 3 donors) after Step A's filters gives 51 (the feasibility count of 65 was made without Step A's filters).
- 2026-09-09 — Exploratory 3 (UC): GCST90446795 (EAS UC, 1,375,998,802 bytes, MD5 matches the Catalog) and GCST004133 (EUR UC, de Lange 2017, no published MD5) via code/15. UC cohort sizes from Liu 2023 Methods: SHA1 2,400, ICH1 1,124, KOR1 1,569, JPN1 1,769 UC cases; controls as for CD (6,279; 3,724; 4,419; 1,034). EUR UC effective N from 12,366 / 33,609.
- 2026-09-09 — Exploratory 4: large intestine = atlas regions LargeInt + REC (appendix excluded; the plan did not define it). Neuronal fell below the lineage threshold in both regions (K = 7).
- 2026-09-09 — Exploratory models used the primary control annotation (all expressed genes, main atlas); ~8% of fine-type genes lie outside that universe. E5 bootstrap gene sets kept 1,414–1,489 genes after coordinate matching.
- 2026-09-09 — Documented filters: munge (registered defaults) removed 72,526 EAS SNPs with effective N < 9,669; LDSC removed chi2 > 80 SNPs (93 EAS, 384 de Lange, 0 FinnGen).

## Independent review and post hoc analyses
- 2026-09-09 — Independent review (plan, code and results only) reproduced W = 4.045, p = 0.853 from raw LDSC files and found no leakage. Reporting and test fixes made without new analysis. Two issues led to post hoc analyses, reported alongside (never replacing) the registered result: (1) Step A used atlas "Sample name" as donor; T036/T036NEG/T036POS and T110NEG/T110POS are single children (same age, sex and library names), so 15 donors, not 18; (2) the registered Cov(D) = Cov(EAS) + Cov(EUR) assumption was checked with block-aligned jackknife values.
- 2026-09-09 — Archived the IIBDGC readme (data/raw/gwas/IIBDGC_IBD_EAS_sumstats_readme.sst.txt, 1,852 bytes, SHA-256 7de4e2c4bbab25bd231276275e4a1604ef4f348452114a76c6b134f19ef2ef3a; code/22). It states 'order for EAS: SHA1, ICH1, KOR1, JPN1', consistent with the ImmunoChip cohort in column 2 of the Direction field. Donor-corrected gene sets (code/20): 15 donors, K = 8, top-10% sets 1,492 genes, Jaccard 0.85–0.95 vs registered.
- 2026-09-10 — POST HOC results. Donor-corrected rerun (code/24): primary W = 2.16, df 8, p = 0.976; validation (i) p = 0.477; (ii) p = 0.233; per-lineage tau* changed by <= 0.044 (EAS, de Lange) and <= 0.10 (FinnGen). Aligned-block check (code/23; 950,331 common SNPs): EAS-EUR jackknife errors correlated 0.24–0.40 across lineages, so the registered summed-covariance test is conservative; direct aligned W = 7.22, p = 0.513 (summed on same SNPs W = 4.39, p = 0.820); se(D) 13–22% smaller; all intervals still include 0. Conclusion unchanged.
