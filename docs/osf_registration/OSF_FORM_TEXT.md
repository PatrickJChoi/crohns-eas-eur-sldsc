# Text to paste into the OSF Preregistration form (Study R2-1)

Step 1: "No" (no existing project). Step 2: "OSF Preregistration". Then fill the fields below.
Upload the three files in this folder (1_, 2_, 3_) where the form allows files (Study design / Analysis plan sections or "Other").

## Title
Do the gut cell types that carry Crohn's disease heritability differ between East Asian and European ancestry? A pre-registered stratified LD score regression study

## Description (summary)
Crohn's disease (CD) genetics differs between East Asian and European ancestry (e.g., NOD2 risk variants are rare in East Asians). We will test whether CD heritability is concentrated near genes specific to the same gut cell lineages in both ancestries. We use public GWAS summary statistics (East Asian: GCST90446793, Liu et al. 2023; European: GCST004132, de Lange et al. 2017; European replication: FinnGen R13 K11_CD_STRICT2) and lineage gene sets from healthy postnatal cells of the Gut Cell Atlas (Elmentaite et al. 2021), analysed with stratified LD score regression (LDSC-SEG). The primary outcome is one joint Wald test of East Asian minus European differences in standardized heritability effect (τ*) across cell lineages. The data were downloaded and checked for provenance and format, but no association result (heritability, enrichment, or per-variant statistic) has been computed or viewed.

## Hypotheses
H0: for every postnatal gut cell lineage, τ* for that lineage's specifically expressed genes is equal for CD in East Asian and European ancestry. H1: it differs for at least one lineage. Two-sided; no direction predicted.

## Study type
Observational study (secondary analysis of existing public data).

## Blinding
No blinding of participants applies. Analyst blinding: no association results were computed or viewed before registration.

## Study design
Cross-ancestry comparison of cell-type heritability enrichment using GWAS summary statistics and single-cell reference data. Full steps A–E: see 1_analysis_plan_R21.md, section 3.

## Randomization
Not applicable.

## Existing data
Registration prior to analysis of the data.

## Explanation of existing data
All data are public and were downloaded on 2026-10-07 (provenance and checksums in 3_code_and_environment.zip). Before registration we only checked variant identifiers, alleles, allele frequencies, which cohorts contributed each variant, and single-cell metadata (2_data_feasibility_R21.md). Published SNP-heritability values from Liu et al. 2023 were read in the paper. No association statistic was summarised by us.

## Data collection procedures
No new data are collected. Sources, accessions and download scripts: 2_data_feasibility_R21.md and 3_code_and_environment.zip.

## Sample size
East Asian CD GWAS: 7,372 cases / 15,456 controls (4 cohorts). European CD GWAS: 12,194 cases / 28,072 controls. European replication (FinnGen): 2,608 cases / 497,350 controls. Cell reference: healthy postnatal gut cells from 18 donors (153,961 cells before excluding lymph node).

## Sample size rationale
Fixed by the public data. Published liability-scale SNP-heritability for CD is 0.213 (s.e. 0.027) in East Asians and 0.196 (s.e. 0.020) in Europeans (Liu et al. 2023). The minimum detectable difference will be estimated by simulation before the real data are analysed and will be reported.

## Stopping rule
Not applicable (fixed existing data; the primary analysis is run once).

## Manipulated variables
None.

## Measured variables
GWAS summary statistics (effect, standard error, alleles, contributing cohorts); single-cell gene expression (raw counts) with donor, diagnosis, region and lineage labels.

## Indices
Gene specificity per lineage (mean donor-level CPM in the lineage divided by the sum across lineages); top 10% of genes define each lineage gene set; τ* = τ · sd(annotation) · M / h²g (Gazal et al. 2017). See plan section 3.

## Statistical models
For each ancestry and lineage: LDSC --h2 with baseline v1.2 + all-expressed-genes control + lineage annotation, ancestry-matched 1000 Genomes LD scores and weights, block jackknife (200 blocks). Primary test: Wald W = Dᵀ Cov(D)⁻¹ D, D = τ*(East Asian) − τ*(European) across K lineages (expected K = 8), Cov(D) = sum of the two ancestries' jackknife covariances, compared with χ²(K). Secondary: per-lineage z, Holm-adjusted.

## Transformations
East Asian positions lifted GRCh38→hg19 (UCSC chain); per-variant effective N for East Asian variants from contributing cohorts; Z = effect/SE aligned to the reference allele; strand-ambiguous SNPs removed; HapMap3 no-MHC regression SNPs; munge_sumstats defaults.

## Inference criteria
α = 0.05 for the primary Wald test. Replication: the same test with FinnGen in place of de Lange 2017. Interpretation table: plan section 5. A non-significant result is reported as "no evidence of a difference", with confidence intervals and the minimum detectable difference.

## Data exclusion
Cells: only Healthy adult and Pediatric healthy; lymph node excluded; predicted doublets excluded; donor–lineage pairs with < 20 cells excluded; lineages with < 200 cells or < 5 donors dropped before GWAS data are processed. Variants: unmatched to 1000 Genomes, allele-inconsistent, strand-ambiguous, outside HapMap3 no-MHC, or removed by munge_sumstats defaults.

## Missing data
Variants missing from a cohort contribute through per-variant effective N; variants absent from either GWAS are not regression SNPs for that ancestry.

## Exploratory analyses
(1) 65 fine cell types; (2) excluding ±1 Mb around NOD2; (3) ulcerative colitis with the same pipeline; (4) small- and large-intestine-only references; (5) τ* on 10 donor-bootstrap gene sets. Sensitivity analyses (pre-specified, do not change the primary conclusion): symmetric exclusion of ImmunoChip-contributed SNPs; baselineLD v2.2; top 5% / 20% gene sets; ±50 kb windows.

## Other
Limitation accepted before registration: no independent East Asian CD GWAS with public summary statistics exists, so the East Asian side cannot be replicated. The East Asian meta-analysis includes an ImmunoChip cohort (ICH1); handled by per-variant effective N (primary) and symmetric exclusion (sensitivity).
