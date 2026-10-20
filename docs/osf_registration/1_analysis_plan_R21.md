# Analysis plan (pre-registration) — Study R2-1

**Do the gut cell types that carry Crohn's disease heritability differ between East Asian and European ancestry?**
Version frozen for OSF registration, 2026-09-08.
No genetic association result (heritability, enrichment, per-variant statistic) has been computed or viewed by us.

## 1. Hypothesis
H0: for every postnatal gut cell lineage, the standardized heritability effect of that lineage's specifically expressed
genes (τ\*) is the same for Crohn's disease (CD) in East Asian (EAS) and European (EUR) ancestry. H1: it differs for at
least one lineage. Two-sided; no direction is predicted. Prior work points both ways: cell-type regulatory annotations
capture consistent heritability across EAS and EUR for many traits (Amariuta 2020), while causal effect sizes near
immune-specific genes are more population-specific (Shi 2021) and CD genetics is more ancestry-dependent than UC
(Liu 2023).

## 2. Data (fixed; provenance in data/raw/r21/PROVENANCE.tsv)
- EAS CD GWAS: GCST90446793 (Liu 2023), 7,372 cases / 15,456 controls; case fraction 32.3%.
- EUR CD GWAS: GCST004132 (de Lange 2017), 12,194 / 28,072; case fraction 30.3%.
- EUR replication: FinnGen R13 K11_CD_STRICT2, 2,608 / 497,350; case fraction 0.52%.
- Cell reference: Gut Cell Atlas full object (Elmentaite 2021). LD: 1000G Phase 3 EAS and EUR baseline v1.2,
  weights, frequencies, HapMap3 no-MHC SNP list (Zenodo 10515792). Genes: GENCODE v19 (hg19) protein-coding.

## 3. Primary analysis (exact)
**Step A: lineage gene sets (atlas only, no GWAS data).** Keep cells with Diagnosis ∈ {Healthy adult, Pediatric
healthy}, Region ≠ lymph node, predicted_doublets = False. Lineages = atlas field `category`. Sum raw counts per donor ×
lineage (pseudobulk, donor as the unit; Squair 2021); keep donor-lineage pairs with ≥20 cells; convert to CPM. Genes:
autosomal protein-coding GENCODE v19 genes matched by Ensembl ID, with mean CPM ≥ 1 in at least one lineage. Specificity
of gene g for lineage L = (mean over donors of CPM in L) / (sum of that over all lineages) (Skene 2016; Bryois 2020).
Gene set for L = top 10% of genes by specificity. Lineages with < 200 cells or < 5 donors after filtering are dropped
before any GWAS data are processed; K = the number retained (expected 8).

**Step B: annotations and LD scores, per ancestry.** Annotation for L: 1000G SNPs within ±100 kb of a gene in L's set
(GENCODE v19 gene start–end). Control annotation: the same for all expressed genes. LD scores computed with `ldsc --l2
--ld-wind-cm 1`, separately with the EAS and EUR 1000G PLINK files, printed for HapMap3 no-MHC SNPs.

**Step C: summary statistics.** EAS: lift GRCh38→hg19 (UCSC chain), match position and alleles to 1000G EAS; per-variant
effective N = sum over contributing cohorts (Direction string; order SHA1, ICH1, KOR1, JPN1) of 4/(1/cases + 1/controls)
using 2,552/6,279, 1,611/3,724, 1,619/4,419, 1,590/1,034. EUR: match hg19 position and alleles to 1000G EUR; effective
N = 34,005 for every variant (all variants present in all cohorts). Both: Z = effect/SE aligned to the reference
allele; drop strand-ambiguous A/T and C/G SNPs; regression SNPs = HapMap3 no-MHC; then `munge_sumstats.py` with
defaults.

**Step D: model.** For each ancestry a and lineage L: `ldsc --h2` with baseline v1.2 + control + annotation L
(LDSC-SEG; Finucane 2018), ancestry-matched weights and frequency files, `--overlap-annot --print-coefficients
--print-delete-vals`. τ\*_{L,a} = τ_L · sd(annotation L) · M / h²_g, with M = reference SNPs with MAF ≥ 5% and h²_g from the
same model (Gazal 2017). The block-jackknife delete values (200 blocks; identical SNP blocks across the K fits within an
ancestry) give each ancestry's K×K covariance of τ\*.

**Step E: primary test.** D = τ\*_EAS − τ\*_EUR (vector over K lineages); Cov(D) = Cov_EAS + Cov_EUR (independent
people). Wald W = Dᵀ Cov(D)⁻¹ D, compared with χ²_K. **Primary outcome: the p-value of W; α = 0.05.**
Secondary (reported, not primary): per-lineage z = D_L / se(D_L), Holm-adjusted across K.

## 4. Validation (independent patients)
EUR side: rerun Steps C–E with FinnGen R13 K11_CD_STRICT2 replacing de Lange 2017 (FinnGen rsIDs matched to 1000G EUR;
effective N 10,378). Two tests: (i) consistency, W for τ\*_deLange − τ\*_FinnGen; (ii) replication, W for
τ\*_EAS − τ\*_FinnGen. **EAS side: no independent cohort exists** (PI accepted at Gate 2); this is stated as a limitation.

## 5. Success criteria and interpretation
| Result | Interpretation |
|---|---|
| Primary p < 0.05 and replication p < 0.05 | CD heritability is distributed differently across gut cell lineages in EAS vs EUR; EUR side robust across two independent cohorts; EAS side rests on one dataset. Holm-significant lineages name where. |
| Primary p < 0.05, replication p ≥ 0.05, same sign for all Holm-significant lineages | Suggestive; replication underpowered or EUR-dataset-specific. |
| Primary p < 0.05, replication disagrees in sign | Not supported; likely dataset-specific. |
| Primary p ≥ 0.05 | No evidence of an ancestry difference. Not proof of equality: we report per-lineage difference CIs and the minimum detectable difference from the Phase 4 simulations. |
| Consistency test (i) p < 0.05 | The two EUR datasets disagree (e.g., registry vs clinical ascertainment); all cross-ancestry results are read with extra caution. |
No accuracy metrics are used, so class-balance baselines do not apply; case fractions are listed in Section 2.

## 6. Stability and robustness (pre-specified; do not change the primary conclusion)
- Gene-set stability: 100 donor bootstraps of Step A; median Jaccard overlap with the main gene set per lineage.
- Sensitivity: (a) drop every SNP with an ImmunoChip (ICH1) contribution from **both** ancestries' regression SNPs;
  (b) baselineLD v2.2 instead of baseline v1.2; (c) top 5% and top 20% gene sets; (d) ±50 kb windows.
- Pipeline checks before unblinding (Phase 4): null and planted-signal simulations; minimum detectable difference.

## 7. Exploratory analyses (labelled exploratory)
1. The 65 fine cell types with ≥100 cells and ≥3 donors: per-ancestry τ\* and their cross-ancestry correlation.
2. Primary model excluding ±1 Mb around NOD2 (absent EAS risk variants).
3. Same pipeline for ulcerative colitis (would need the EAS and EUR UC files), to see whether any difference is CD-specific.
4. Small-intestine-only and large-intestine-only references.
5. τ\* recomputed on 10 donor-bootstrap gene sets.

## 8. Positioning (methods literature; all DOIs verified)
LDSC (Bulik-Sullivan 2015, doi:10.1038/ng.3211); stratified LDSC (Finucane 2015, doi:10.1038/ng.3404); LDSC-SEG
(Finucane 2018, doi:10.1038/s41588-018-0081-4); τ\* (Gazal 2017, doi:10.1038/ng.3954); top-10% specificity gene sets
(Skene 2016, doi:10.3389/fnins.2016.00016; Bryois 2020, doi:10.1038/s41588-020-0610-9); single-cell–GWAS integration
(Jagadeesh 2022, doi:10.1038/s41588-022-01187-9); pseudobulk with donor as unit (Squair 2021,
doi:10.1038/s41467-021-25960-2); trans-ethnic genetic correlation (Brown 2016, doi:10.1016/j.ajhg.2016.05.001);
annotation-stratified trans-ethnic correlation, S-LDXR (Shi 2021, doi:10.1038/s41467-021-21286-1); consistent
annotation heritability across populations (Amariuta 2020, doi:10.1038/s41588-020-00740-8); CD across ancestries
(Liu 2023, doi:10.1038/s41588-023-01384-0). Gap: S-LDXR and Liu 2023 ask whether **effect sizes** are shared within
annotations. This study asks whether the **location of heritability across gut cell lineages** differs, for CD
specifically, with single-cell-defined lineages.
