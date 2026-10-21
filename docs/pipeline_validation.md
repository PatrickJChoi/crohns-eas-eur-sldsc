# Pipeline validation on simulated data (2026-10-08)

Registered plan: doi:10.17605/OSF.IO/HQ4A8. Everything here uses **simulated** GWAS results; no real-GWAS association
result has been computed. Numbers: `results/simulations/simulation_summary.txt` (from `code/12_simulation_summary.py`), per-
replicate results in `results/simulations/sims/*.jsonl`.

## What was built and checked
| Check | Result |
|---|---|
| Pipeline Steps A–E as registered | code/05 (gene sets), 06 (annotations, LD scores), 07 (summary statistics), 08 (LDSC-SEG fits, τ\*), 09 (joint Wald test) |
| Step A: K lineages | 8 retained; 1,493 genes per lineage from 14,927 expressed protein-coding genes |
| Gene-set stability (registered §6) | donor bootstrap, median Jaccard: epithelial 0.87, mesenchymal 0.80, endothelial 0.79, B 0.78, T 0.76, myeloid 0.75, plasma 0.75, neuronal 0.66 |
| Leakage / blinding tests (`tests/test_pipeline.py`, 9/9 pass) | Step A opens no GWAS file and reproduces the gene sets byte for byte; Step B for EAS passes only EAS reference files to LDSC and touches no GWAS file; Step C refuses the real GWAS files without `--unblind`; arithmetic of effective N, allele matching, jackknife covariance, Wald and Holm |
| Our jackknife vs LDSC's | standard errors from the delete values reproduce LDSC's to ~1e-15 |
| `--overlap-annot` | dropped in simulations only (≈7× faster); verified delete values bit-identical and full τ equal to ~5 significant digits. The real-data analysis keeps it as registered. |

## Simulation design
Summary statistics were generated from 1000 Genomes genotypes of each ancestry (block-diagonal LD, 2 cM blocks; 2% of
common SNPs causal; East Asian per-variant N from the real cohort-coverage strings). They were written in the real files'
formats and run through Steps C–E unchanged. Simulated heritability precision (LDSC h² z) averaged about 9 in both
ancestries (published: 7.9 East Asian, 9.8 European).

| Scenario | Planted | Reps |
|---|---|---|
| null | no lineage enrichment in either ancestry | 20 |
| shared | T cells and myeloid enrichment, the same in both ancestries | 60 (plus 8 random gene sets each) |
| diff4 / diff6 / diff10 | European myeloid enrichment raised (EAS unchanged) | 20 each |

## Results
**1. No false findings when there is no difference.** Primary test rejection rate where H0 is true: null 2/20, shared
2/60; pooled **4/80 = 5.0% (95% CI 1.4–12.3%)**, p-values consistent with uniform (KS p = 0.87 and 0.47).

**2. Shuffled labels find nothing.** Replacing lineage gene sets with 8 random gene sets of the same size: Wald test
rejected **0/60** (95% CI 0–6%); random-set τ\* averaged −0.01 to +0.05, against +0.30 to +0.32 for the planted
lineages. The small positive values are expected: random sets include some planted T-cell/myeloid genes.

**3. Planted signal is recovered, in the right place.** In "shared", T cells and myeloid τ\* averaged +0.30 to +0.32 in
both ancestries; other lineages −0.05 to +0.05. In "null", all lineages −0.04 to +0.05. When the myeloid difference was
planted, the per-lineage (secondary) tests flagged **myeloid** in 95% of diff10 replicates and another lineage in 5%.

**4. The jackknife is honest.** Spread of the East Asian minus European difference across replicates divided by the
average jackknife standard error: 0.93–1.08 in "shared" (60 reps), 0.80–1.15 in "null". With a large planted difference
it was below 1 for myeloid (0.62), i.e., conservative.

**5. Power and minimum detectable difference.**
| Mean myeloid τ\* difference | 0.25 | 0.46 | 0.72 |
|---|---|---|---|
| Primary joint test (8 lineages) | 0% | 25% | 75% |
| Secondary myeloid test (Holm) | 0% | 30% | 95% |
For a difference in one lineage, the **primary test reaches 80% power only for differences above about 0.72 τ\* units**,
more than twice the size of the planted enrichment itself (≈0.3). The single-lineage secondary test is more sensitive
(80% somewhere between 0.46 and 0.72) because the joint test spreads its attention over 8 lineages.

## What this means for the real-data analysis
- The pipeline behaved as intended in these simulations: a false-positive rate consistent with 5% (4/80; 95% CI
  1.4–12.3%, so rates up to about 12% are not excluded), no signal from random gene sets (a shuffled-label analogue:
  random gene sets rather than permuted cell labels), planted signal found in the right lineage, jackknife SEs close
  to the empirical spread, and no leakage in the tests.
- **Power is limited.** A significant primary result would be credible. The primary test has about 75–80% power only
  for one-lineage differences of about 0.72 τ\* or more; power is not the same as ruling out, and smaller differences
  would often be missed. Per the registered interpretation table, a null is reported as "no evidence of a difference",
  with these detectable sizes stated.

*Revised 2026-10-09 after the independent review of the results (wording only; numbers unchanged).*

## Limitations of the simulations
Block-diagonal LD from ~500 reference genomes per ancestry; enrichment planted only through the lineage and control gene
sets (no other functional annotations); the ImmunoChip cohort is represented only through per-variant N, not through
its selection of immune regions (the registered sensitivity analysis that drops ImmunoChip-contributed SNPs addresses
this in the real-data analysis); simulated h² precision matched the published values only approximately (~9 vs 7.9 / 9.8).
Causal variants and effects were drawn **independently** for each ancestry (no shared genetic architecture), so the
simulations could not show whether the registered Cov(D) = Cov(EAS) + Cov(EUR) is conservative when ancestries share
causal variants. Gene sets were built with atlas sample names as donors (two children appear under several sample
names; 15 true donors, not 18), so the donor-bootstrap stability figures above are somewhat optimistic.
Two technical incidents, both resolved: the simulator and LDSC ran out of memory when run together (fixed by running
them in separate stages), and 6 diff10 fits were rerun after an out-of-memory error (resumed run, same simulated data).
