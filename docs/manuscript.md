# No evidence that Crohn's disease heritability is distributed differently across gut cell lineages in East Asian and European ancestry: a pre-registered stratified LD score regression study

Patrick Choi¹

¹ Department of Integrative Biology, University of South Florida, Tampa, FL, USA

ORCID: 0009-0004-1051-2074

Correspondence: patrickchoi@usf.edu

Pre-registration: OSF, doi:10.17605/OSF.IO/HQ4A8 (registered 8 September 2026)

---

## Abstract

**Background.** The genetic architecture of Crohn's disease (CD) differs between East Asian and European ancestry; for
example, major NOD2 risk variants are rare in East Asians. Whether CD genetic risk nevertheless acts through the same
intestinal cell types in both ancestries has not been tested.

**Methods.** In a pre-registered analysis, we built gene sets specific to eight postnatal gut cell lineages from
healthy cells of the Gut Cell Atlas and estimated CD heritability enrichment near those genes with stratified LD score
regression (LDSC-SEG), separately in an East Asian GWAS (7,372 cases, 15,456 controls) and a European GWAS (12,194
cases, 28,072 controls), each with ancestry-matched linkage-disequilibrium references. The primary outcome was a joint
Wald test of East Asian minus European differences in the standardized effect τ\* across the eight lineages. European
results were replicated in FinnGen (2,608 cases). The pipeline was first validated on simulated data.

**Results.** The primary test gave no evidence of an ancestry difference (W = 4.04, df = 8, p = 0.853), and no single
lineage differed after Holm correction. In the registered model, the largest τ\* values in both ancestries were for T
cells and myeloid cells. No registered sensitivity analysis, exploratory analysis or post hoc check gave evidence of a
difference. Precision was limited: simultaneous 95% confidence intervals for the per-lineage differences spanned −0.50
to +0.36 τ\* units, so differences as large as the lineage enrichments themselves cannot be excluded.

**Conclusions.** We found no evidence that the gut cell lineages carrying CD heritability differ between East Asian and
European ancestry, but the available East Asian data cannot rule out differences of the size of the enrichments. Larger
East Asian GWAS and an independent East Asian cohort are needed.

---

## Introduction

Crohn's disease is a chronic inflammatory disease of the gastrointestinal tract with a substantial heritable component.
Most CD genetic associations were discovered in individuals of European ancestry. The largest East Asian study to date
found that, although inflammatory bowel disease genetic effects are broadly consistent across ancestries, "genetics
underlying CD appears more ancestry dependent than UC", driven by allele frequency at NOD2 and effect size at TNFSF15
(Liu et al. 2023).

A separate line of work maps disease heritability onto cell types by testing whether heritability is enriched near genes
specifically expressed in each cell type (Finucane et al. 2018; Jagadeesh et al. 2022). Recent single-cell studies of
the intestine linked CD heritability mainly to immune cells: adaptive immune cells, in particular T cells, and
neutrophils (Zhao et al. 2026), and inflammatory monocytes and macrophages (Krzak et al. 2026). To our knowledge,
whether the same cell types carry CD heritability in East Asian ancestry has not been tested. Two prior results frame the question. Across many traits, heritability captured by
cell-type-specific regulatory annotations is consistent between East Asian and European populations (Amariuta et al.
2020). On the other hand, causal effect sizes near immune-specific genes are more population-specific than elsewhere in
the genome (Shi et al. 2021). Liu et al. (2023) reported that cross-ancestry genetic correlation appeared similar across
functional annotations when partitioned; that analysis did not use cell-type gene sets.

We therefore asked, in a pre-registered study, whether CD heritability is distributed differently across gut cell
lineages in East Asian and European ancestry. We prespecified a single joint test across lineages, so that no lineage
was chosen after seeing the data, and we validated the pipeline on simulated data before analysing the real data.

## Methods

### Pre-registration and analysis order
The hypothesis, data, gene-set construction, model, primary test, validation, interpretation rules, sensitivity
analyses and exploratory analyses were registered on OSF before any heritability or enrichment estimate was computed
(doi:10.17605/OSF.IO/HQ4A8). The analysis pipeline was then tested on simulated and shuffled data, and the primary
analysis was run once. Deviations are listed below and in the study's decision log.

### Data
**GWAS summary statistics.** East Asian (EAS) CD: fixed-effect meta-analysis of four cohorts (China SHA1, 2,552 cases;
ImmunoChip cohort ICH1, 1,611; Korea KOR1, 1,619; Japan JPN1, 1,590; 15,456 controls in total; GWAS Catalog
GCST90446793; Liu et al. 2023). European (EUR) CD: 12,194 cases and 28,072 controls (GCST004132; de Lange et al. 2017).
European replication: FinnGen release 13, endpoint K11_CD_STRICT2 (2,608 cases, 497,350 controls; Kurki et al. 2023).
Ulcerative colitis (UC) GWAS from the same two consortia (GCST90446795, GCST004133) were used in one exploratory analysis.

**Single-cell reference.** The Gut Cell Atlas full object (Elmentaite et al. 2021). We used cells from healthy adult and
healthy pediatric donors, excluding lymph node cells and predicted doublets (142,797 cells), grouped into the atlas's
eight lineage categories.

**Reference panels and genes.** 1000 Genomes Phase 3 East Asian and European genotypes (1000 Genomes Project Consortium
2015) with the matching baseline v1.2 LD scores, regression weights and allele frequencies distributed for stratified
LD score regression (Finucane et al. 2015; Gazal et al. 2017); GENCODE v19 gene coordinates (Harrow et al. 2012).

### Lineage gene sets
Raw counts were summed per donor and lineage (pseudobulk; donor-lineage pairs with at least 20 cells), converted to
counts per million, and averaged across donors (Squair et al. 2021). For each autosomal protein-coding gene expressed at
a mean of at least 1 count per million in some lineage (14,927 genes), specificity for a lineage was its mean expression
in that lineage divided by the sum across lineages (Skene and Grant 2016; Bryois et al. 2020). Each lineage's gene set
was its top 10% of genes (1,493 genes). All eight lineages met the registered minimum of 200 cells and 5 donors.

### Annotations, summary statistics and model
For each lineage, the annotation contained 1000 Genomes SNPs within 100 kb of a gene in the set; a control annotation
contained SNPs within 100 kb of any expressed gene. LD scores were computed separately for each ancestry with its own
reference panel (1 cM windows). EAS positions were lifted from GRCh38 to hg19 and matched to the reference by position
and alleles; per-variant effective sample size was rebuilt from the cohorts contributing to each variant. Strand-
ambiguous SNPs were removed and HapMap3 SNPs outside the MHC were used for regression.

For each ancestry and lineage we fitted LDSC-SEG (Finucane et al. 2018): stratified LD score regression with the
baseline v1.2 annotations, the control annotation and the lineage annotation (Bulik-Sullivan et al. 2015; Finucane et
al. 2015), using the LDSC Python 3 implementation (Breeze et al. 2025). The lineage effect was expressed as τ\*, the
coefficient scaled by the annotation's standard deviation and the number of common SNPs and divided by SNP-heritability
(Gazal et al. 2017). Uncertainty came from LDSC's 200-block jackknife.

### Primary test, validation and interpretation
For the eight lineages, D = τ\*(EAS) − τ\*(EUR). Following the registered plan, the covariance of D was the sum of the
two ancestries' jackknife covariances, and the joint Wald statistic W = Dᵀ Cov(D)⁻¹ D was compared with a χ² distribution
on 8 degrees of freedom (α = 0.05). Per-lineage tests used Holm's step-down correction (Holm 1979). The European side was
replicated by repeating the analysis with FinnGen. The registered interpretation rules state that a non-significant
primary result means "no evidence of an ancestry difference", not equality, and is reported with confidence intervals
and the minimum detectable difference.

### Sensitivity and exploratory analyses
Registered sensitivity analyses: (a) removing every regression SNP with an ImmunoChip contribution from both
ancestries; (b) baselineLD v2.2 instead of baseline v1.2; (c) top-5% and top-20% gene sets; (d) 50 kb windows.
Registered exploratory analyses: 51 fine cell types; excluding 1 Mb around NOD2; UC; small- and large-intestine-only
references; and ten donor-bootstrap gene sets.

### Pipeline validation on simulated data
Before analysing the real data, we simulated GWAS summary statistics from each ancestry's 1000 Genomes genotypes with
block-diagonal LD, planted known enrichment, and ran them through the unchanged pipeline (140 replicates). We checked
the false-positive rate, behaviour with random gene sets, recovery of planted enrichment and power.

### Post hoc analyses after independent review
An independent review of the plan, code and results identified two issues that led to post hoc analyses, both specified
before being run and reported alongside the registered result. First, the registered gene-set step used the atlas
sample name as the donor identifier; two children contributed separately sorted cell fractions under several sample
names, so the reference has 15 rather than 18 donors. We reran the analysis with true donors. Second, the registered
covariance treats the two ancestries' jackknife errors as independent, although shared causal architecture can correlate
them. We restricted all datasets to a common SNP set, so that jackknife blocks covered the same genomic segments, and
computed Cov(D) directly from block-aligned values.

### Data and code availability
All input data are public (accessions above). Code, logs and result files are in the study repository (https://github.com/PatrickJChoi/crohns-eas-eur-sldsc).
The registered plan and data checks are on OSF (doi:10.17605/OSF.IO/HQ4A8).

## Results

### Pipeline validation
In simulations without a true ancestry difference, the primary test rejected in 4 of 80 replicates (5.0%; 95% CI
1.4–12.3%). Random gene sets in place of lineage gene sets produced no rejections in 60 replicates. At the largest planted
myeloid difference, the per-lineage tests flagged the myeloid lineage in 95% of replicates and another lineage in 5%. Power to detect a one-lineage difference was 0%, 25% and 75% for
mean differences of 0.25, 0.46 and 0.72 τ\* units, so the registered minimum detectable difference at 80% power was
above 0.72.

### Data
LDSC used 986,470 SNPs for the EAS GWAS, 1,149,713 for the European GWAS and 1,175,811 for FinnGen. SNP-heritability on
the effective-sample-size scale was 0.847 (s.e. 0.070), 0.510 (0.036) and 0.477 (0.082), respectively. LDSC intercepts
were 0.997, 1.070 and 1.024.

### Primary analysis
The joint test gave no evidence that CD heritability is distributed differently across the eight lineages in East
Asian and European ancestry (W = 4.04, df = 8, p = 0.853; Table 1, Figure 2). No lineage differed after Holm correction
(all adjusted p = 1.00). The largest absolute difference was for T cells (D = −0.15, 95% CI −0.40 to +0.10). In the
registered model the largest τ\* values in both ancestries were for T cells (EAS +0.20, EUR +0.35) and myeloid cells
(EAS +0.23, EUR +0.20) (Figure 1).

The confidence intervals show what the data can and cannot exclude. Marginal 95% intervals for D spanned −0.40 to +0.27,
and simultaneous (Bonferroni, eight lineages) 95% intervals spanned −0.50 to +0.36. The simultaneous interval for T cells
(−0.50 to +0.21) includes a difference as large as the European T-cell τ\* itself.

### Replication and sensitivity analyses
Repeating the analysis with FinnGen gave no evidence that the two European datasets differ (W = 7.60, p = 0.473); this
is not evidence of agreement, because FinnGen estimates were roughly twice as imprecise. With FinnGen as the European
dataset there was again no evidence of an ancestry difference (W = 7.91, p = 0.442). None of the registered sensitivity
analyses gave evidence of a difference (p = 0.521 to 0.995; Table 2). Per-lineage estimates did change, however:
removing ImmunoChip-covered SNPs (about 10% of regression SNPs, concentrated in known immune loci) lowered T-cell τ\* in
both ancestries (EAS +0.20 to −0.03; EUR +0.35 to +0.14), while myeloid τ\* was stable. Much of the T-cell signal
therefore lies in those regions.

### Exploratory analyses
Excluding the NOD2 region (p = 0.822), analysing UC (p = 0.632), using small-intestine-only (p = 0.837) or
large-intestine-only (p = 0.995) references, and using bootstrap gene sets (p = 0.814 to 0.990) each gave no evidence of
an ancestry difference. Across 51 fine cell types, EAS and EUR τ\* were correlated (Pearson r = 0.73); because the fine
gene sets overlap, no p-value is reported. The highest-ranked fine types were T-cell and NK-cell types in both ancestries,
and no fine type differed by more than |z| = 1.79.

### Post hoc analyses
With true donors (15), the primary test again gave no evidence of a difference (W = 2.16, p = 0.976). On 950,331 SNPs
common to all datasets, EAS and EUR jackknife errors were positively correlated (0.24 to 0.40 across lineages), so the
registered test is conservative. Computing the covariance directly from aligned blocks gave W = 7.22, p = 0.513, and
standard errors of D 13–22% smaller than with the registered formula. All intervals still included zero; simultaneous
intervals spanned −0.46 to +0.32.

## Discussion

In this pre-registered study we found no evidence that the gut cell lineages carrying Crohn's disease heritability
differ between East Asian and European ancestry. The same conclusion was reached with a European replication cohort,
in the registered sensitivity and exploratory analyses, and in two post hoc checks prompted by an independent review. In both ancestries the
largest enrichments in the registered model were for T cells and myeloid cells, in line with recent single-cell
studies that linked CD heritability to T cells (Zhao et al. 2026) and to monocytes and macrophages (Krzak et al. 2026).

The main limitation is precision. Our study could not exclude ancestry differences as large as the lineage
enrichments themselves: the minimum detectable difference was above 0.72 τ\* units in simulation, and the simultaneous
intervals in the real data spanned −0.50 to +0.36. A null result here is therefore weak evidence of similarity. It is
consistent with the general finding that functional annotations capture similar heritability across populations
(Amariuta et al. 2020), but it does not test, and does not contradict, the finding that causal effect sizes near immune
genes are more population-specific (Shi et al. 2021), which concerns effect sizes within annotations rather than where
heritability is located.

The T-cell signal depended heavily on SNPs in regions covered by the ImmunoChip, which are dense in known immune loci.
The drop was similar in both ancestries and did not change the estimated T-cell difference (−0.15 in the registered
model, −0.17 without those SNPs), but it means the T-cell enrichment should not be over-interpreted: a large part of it
comes from known association regions rather than from polygenic signal spread across the genome.

Other limitations follow. No independent East Asian CD GWAS with public summary statistics exists, so the East Asian
side rests on one dataset. That dataset includes an ImmunoChip cohort, which we handled through per-variant effective
sample size and, in sensitivity analysis, by removing affected SNPs. The registered covariance assumption made the
primary test conservative; a direct, block-aligned estimate gave the same conclusion. The cell-type reference came from
15 healthy donors, and the gene-set construction initially miscounted donors; the corrected analysis gave the same
conclusion. The analysis concerns common-variant SNP-heritability; rare variants, including NOD2 coding variants that
are absent in East Asians, are not directly represented. Finally, the de Lange LDSC intercept of 1.070 indicates some
inflation not due to polygenic signal; LDSC estimates the intercept freely.

Larger East Asian GWAS, an independent East Asian cohort, and cell-type references from East Asian donors would allow
smaller differences to be detected.

## Acknowledgements
We want to acknowledge the participants and investigators of the FinnGen study. We thank the International IBD Genetics
Consortium and the authors of the East Asian and European GWAS and of the Gut Cell Atlas for making their data public.

## References
- Amariuta T, et al. Improving the trans-ancestry portability of polygenic risk scores by prioritizing variants in predicted cell-type-specific regulatory elements. Nat Genet. 2020. doi:10.1038/s41588-020-00740-8
- Breeze CE, et al. LDscore: a scalable, Python 3-powered web platform for LD score regression analysis. bioRxiv. 2025. doi:10.64898/2025.12.19.695639
- Bryois J, et al. Genetic identification of cell types underlying brain complex traits yields insights into the etiology of Parkinson's disease. Nat Genet. 2020. doi:10.1038/s41588-020-0610-9
- Bulik-Sullivan BK, et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet. 2015. doi:10.1038/ng.3211
- de Lange KM, et al. Genome-wide association study implicates immune activation of multiple integrin genes in inflammatory bowel disease. Nat Genet. 2017. doi:10.1038/ng.3760
- Elmentaite R, et al. Cells of the human intestinal tract mapped across space and time. Nature. 2021. doi:10.1038/s41586-021-03852-1
- Finucane HK, et al. Partitioning heritability by functional annotation using genome-wide association summary statistics. Nat Genet. 2015. doi:10.1038/ng.3404
- Finucane HK, et al. Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nat Genet. 2018. doi:10.1038/s41588-018-0081-4
- Gazal S, et al. Linkage disequilibrium-dependent architecture of human complex traits shows action of negative selection. Nat Genet. 2017. doi:10.1038/ng.3954
- Harrow J, et al. GENCODE: the reference human genome annotation for The ENCODE Project. Genome Res. 2012. doi:10.1101/gr.135350.111
- Holm S. A simple sequentially rejective multiple test procedure. Scand J Stat. 1979;6(2):65–70. JSTOR: https://www.jstor.org/stable/4615733 (no DOI)
- Jagadeesh KA, et al. Identifying disease-critical cell types and cellular processes by integrating single-cell RNA-sequencing and human genetics. Nat Genet. 2022. doi:10.1038/s41588-022-01187-9
- Krzak M, et al. Single-cell RNA sequencing of terminal ileal biopsies identifies signatures of Crohn's disease pathogenesis. Nat Genet. 2026. doi:10.1038/s41588-026-02634-7
- Kurki MI, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature. 2023. doi:10.1038/s41586-022-05473-8
- Liu Z, et al. Genetic architecture of the inflammatory bowel diseases across East Asian and European ancestries. Nat Genet. 2023. doi:10.1038/s41588-023-01384-0
- Shi H, et al. Population-specific causal disease effect sizes in functionally important regions impacted by selection. Nat Commun. 2021. doi:10.1038/s41467-021-21286-1
- Skene NG, Grant SGN. Identification of vulnerable cell types in major brain disorders using single cell transcriptomes and expression weighted cell type enrichment. Front Neurosci. 2016. doi:10.3389/fnins.2016.00016
- Squair JW, et al. Confronting false discoveries in single-cell differential expression. Nat Commun. 2021. doi:10.1038/s41467-021-25960-2
- The 1000 Genomes Project Consortium. A global reference for human genetic variation. Nature. 2015. doi:10.1038/nature15393
- Zhao Y, et al. Single-cell analysis of chromatin accessibility in the human intestine identifies regulatory programs and clarifies genetic associations in Crohn's disease. Nat Genet. 2026. doi:10.1038/s41588-026-02755-z

## Tables and figures
**Table 1.** Primary analysis: τ\* per lineage (s.e.) in each dataset, East Asian minus European difference D with
marginal and simultaneous 95% CIs, and Holm-adjusted p. Source: results/main_analysis/report_tables.md, Table T2.

**Table 2.** Joint Wald tests for the primary analysis, validation and registered sensitivity analyses, with per-lineage
τ\* for each sensitivity analysis. Source: results/main_analysis/report_tables.md, Tables T3, T4 and T4b.

**Figure 1.** τ\* with 95% CIs by gut cell lineage for the East Asian GWAS, the European GWAS (de Lange) and FinnGen,
registered model. File: results/manuscript/fig1_tau_star.pdf.

**Figure 2.** East Asian minus European difference D by lineage: registered marginal (thick) and simultaneous (thin)
95% CIs, and post hoc block-aligned 95% CIs (open circles; common SNP set). File: results/manuscript/fig2_difference.pdf.
