# Phase 2 — Study R2-1: feasibility and provenance (2026-10-07)

Question: is Crohn's disease (CD) heritability enriched in the same gut cell types in East Asian (EAS) and European
(EUR) ancestry?

**Recommendation at Gate 2: CONDITIONAL GO.** The data are usable and every technical problem found has a workable fix.
One limitation cannot be fixed with public data: **no independent East Asian cohort exists to validate the EAS side**
(red flag 2). Go only if the PI accepts that limitation for the pre-registration; otherwise no-go.

No association results were summarised in this phase. Checks used only variant IDs, alleles, allele frequencies, the
meta-analysis "Direction" strings (which cohorts contributed each variant) and atlas metadata. Heritability values
quoted below come from the published paper, not from our analysis.

## Inventory
Downloads, URLs, UTC times, sizes and MD5s: `data/raw/r21/PROVENANCE.tsv` (+ `PROVENANCE_NOTES.md`). All 11 files that
have a published checksum match it. Checks: `code/06_qc_r21.py` → `results/phase2_R21_qc/qc_summary.txt`.

| Dataset | Accession / source | Who | Labels and where they trace to | Build / format |
|---|---|---|---|---|
| EAS CD GWAS | GWAS Catalog GCST90446793; Liu 2023 Nat Genet doi:10.1038/s41588-023-01384-0 | 7,372 CD / 15,456 controls from 4 cohorts: SHA1 (Chinese; 2,552 CD / 6,279 ctrl), ICH1 (Hong Kong, Korea, Japan, **ImmunoChip**; 1,611 / 3,724), KOR1 (Korean; 1,619 / 4,419), JPN1 (Japanese; 1,590 / 1,034) | CD diagnosed by "clinical characteristics, radiological and endoscopic examination, and histological features" per European consensus / ECCO-ESGAR guidelines (paper Methods, SHA1); the paper states SHA1, KOR1 and JPN1 "followed similar clinical diagnostic criteria"; no separate definition is given for ICH1 (described in an earlier paper) | GRCh38; no rsIDs, no per-variant N; METAL meta-analysis; 13.5 M variants |
| EUR CD GWAS | GCST004132; de Lange 2017 Nat Genet doi:10.1038/ng.3760 | 12,194 CD / 28,072 controls (new UK GWAS, IBD sequencing cohort, earlier IIBDGC GWAS) | "diagnosed using accepted endoscopic, histopathological and radiological criteria" (Online Methods) | GRCh37; no rsIDs, no N, no allele frequency; 9.57 M variants, every one present in all 3 cohorts |
| EUR replication | FinnGen R13 K11_CD_STRICT2 (already downloaded) | 2,608 / 497,350 | Registry: IBD drug reimbursement + >1 CD hospital record | GRCh38 |
| Cell-type reference | Gut Cell Atlas, Elmentaite 2021 Nature doi:10.1038/s41586-021-03852-1 | 428,469 cells, 41 donors; postnatal healthy: 153,961 cells from 18 donors (11 pediatric, 7 adult) | GEO-style obs fields `Diagnosis`, `Age_group`, `Region`, `category` (lineage), `Integrated_05` (fine type) | h5ad (AnnData 0.1.0), 33,538 genes |
| LD references | Zenodo 10515792 (S-LDSC reference files), 1000 Genomes Phase 3 EAS and EUR | — | — | hg19; HapMap3 list 1,217,311 SNPs (no MHC) |

Liu 2023 also included FinnGen R7 in its EUR meta-analysis, but this study uses only its EAS file. FinnGen is not part
of de Lange 2017, so FinnGen is an independent EUR replication.

## Checks
- **Variant matching (needed because neither GWAS has rsIDs).** EAS: 97.3% of variants lifted GRCh38→hg19 (UCSC chain);
  8.08 M matched a 1000G EAS rsID with consistent alleles; **87.6%** of HapMap3 regression SNPs covered. EUR: 8.53 M
  matched; **94.7%** of HapMap3 SNPs covered. 15.5% of matched EAS variants are strand-ambiguous (A/T, C/G) and will be
  dropped, as LDSC's own pre-processing does by default.
- **Ancestry/allele sanity.** EAS GWAS allele frequency vs 1000G EAS: r = 0.982 (non-ambiguous SNPs; 0.976 on HapMap3
  SNPs); vs 1000G EUR: r = 0.581. The file is EAS and alleles are oriented correctly. (EUR GWAS has no frequency column.)
- **Duplicates / overlap.** EAS and EUR samples are different people. Liu 2023 removed cross-cohort overlap among SHA1,
  KOR1 and ICH1 with identity-by-descent (Methods; JPN1 handled separately). FinnGen is not in de Lange 2017. The Korean CD GWAS GCST90103738 (1,621 cases) is
  almost certainly the KOR1 cohort (identical 4,419 controls; 1,621 vs 1,619 cases; same research group), so it is
  **not** independent.
- **Multiple samples per patient.** Not applicable to summary statistics. In the atlas, donors contribute several
  samples (regions, cell fractions); Phase 3 should compute cell-type expression per donor first so that no single donor
  dominates.
- **Per-variant sample size (EAS).** Cohort coverage varies a lot: e.g., 5.86 M variants come from SHA1+KOR1+JPN1,
  2.16 M from SHA1 alone, 1.46 M from JPN1 alone. Per-variant N can be rebuilt from the Direction string and the cohort
  sizes above (cohort order SHA1, ICH1, KOR1, JPN1 per the consortium readme).
- **Power.** Liu 2023 reports liability-scale SNP heritability for CD of 0.213 ± 0.027 (EAS, excluding ImmunoChip) and
  0.196 ± 0.020 (EUR), i.e., about 8 and 10 standard errors from zero. That is enough for single-ancestry enrichment
  estimates. Power for the **difference** between ancestries is not known; Phase 3 must set it by simulation, or frame
  the primary result as an estimate with a confidence interval.
- **Reference resolution.** All 8 postnatal lineages are well represented (440–85,521 cells; 17–18 donors each; myeloid
  2,803 cells). Fine myeloid types are thin (e.g., monocytes 245 cells, macrophages 746). 65 of 89 fine types have
  ≥100 cells and ≥3 donors.
- **Software.** LDSC runs under Python 3: NIH's maintained rewrite (github.com/CBIIT/ldsc, commit 1f09cf0, linked from
  the original authors' README), in its own pinned environment (`requirements-ldsc.txt`). It has not yet been validated
  on known inputs (Phase 4).

## Red flags
1. **ImmunoChip cohort inside the EAS file.** ICH1 contributes to 6.7% of all EAS variants and **10.4% of the HapMap3
   regression SNPs**, concentrated in known immune loci. Liu 2023 excluded ImmunoChip from its own heritability
   analyses "because of its lack of genome-wide coverage". Left unhandled, it gives immune regions extra sample size and
   could inflate immune-cell enrichment in EAS only, the very contrast under test. Fixes for Phase 3: (a) rebuild
   per-variant N from the Direction string; or (b) drop ICH1-contributed SNPs from both ancestries' regression sets
   (symmetric). No ICH1-free EAS file is published (checked the IIBDGC download folder).
2. **No independent EAS validation (cannot be fixed with public data).** No other public full summary statistics
   exist for EAS CD: GWAS Catalog 2026-10-07 lists GCST90103738 (= KOR1) and studies without full summary statistics.
   EUR results can be replicated in FinnGen. A difference between ancestries would rest on a single EAS dataset.
3. **Novelty: moderate.** Liu 2023 already showed cross-ancestry genetic correlation is similar across broad
   functional annotations (not cell types). The unpublished IIBDGC multi-ancestry paper could include a cell-type
   comparison.
4. **Sparse myeloid subtypes** in the only openly downloadable atlas, so a myeloid hypothesis should be at lineage
   level.
5. Both GWAS lack rsIDs and N; handled above, but every mapping step must be scripted and checked in Phase 4.
