### T1. Run record and LDSC diagnostics (results/main_analysis/primary_results.json; LDSC logs)
| dataset | SNPs after munge | SNPs used by LDSC | h2 observed, effective-N scale (SE) | h2 z | mean chi2 | lambda GC | LDSC intercept (SE) |
|---|---|---|---|---|---|---|---|
| EAS | 994,292 | 986,470 (93 removed, chi2 > 80) | 0.847 (0.070) | 12.1 | 1.227 | 1.143 | 0.997 (0.009) |
| EUR_deLange | 1,152,942 | 1,149,713 (384 removed, chi2 > 80) | 0.510 (0.036) | 14.1 | 1.401 | 1.253 | 1.070 (0.012) |
| EUR_FinnGen | 1,178,735 | 1,175,811 (0 removed, chi2 > 80) | 0.477 (0.082) | 5.8 | 1.114 | 1.093 | 1.024 (0.009) |

### T2. Primary analysis: tau* per lineage and the East Asian minus European difference (results/main_analysis/primary_table.tsv)
| lineage | tau* EAS (SE) | tau* EUR de Lange (SE) | tau* EUR FinnGen (SE) | D = EAS - EUR | 95% CI of D | simultaneous 95% CI (Bonferroni, 8) | z | p | p Holm |
|---|---|---|---|---|---|---|---|---|---|
| B cells | +0.04 (0.07) | +0.08 (0.07) | -0.00 (0.13) | -0.04 | -0.24 to +0.15 | -0.31 to +0.23 | -0.44 | 0.66 | 1.00 |
| Endothelial | -0.10 (0.06) | -0.07 (0.06) | +0.03 (0.11) | -0.03 | -0.19 to +0.14 | -0.26 to +0.20 | -0.34 | 0.74 | 1.00 |
| Epithelial | -0.06 (0.05) | -0.04 (0.06) | -0.14 (0.10) | -0.03 | -0.18 to +0.13 | -0.24 to +0.19 | -0.33 | 0.74 | 1.00 |
| Mesenchymal | -0.01 (0.06) | -0.03 (0.06) | +0.19 (0.10) | +0.02 | -0.16 to +0.19 | -0.23 to +0.26 | +0.19 | 0.85 | 1.00 |
| Myeloid | +0.23 (0.08) | +0.20 (0.09) | +0.33 (0.16) | +0.03 | -0.21 to +0.27 | -0.30 to +0.36 | +0.24 | 0.81 | 1.00 |
| Neuronal | -0.17 (0.05) | -0.11 (0.05) | -0.07 (0.11) | -0.06 | -0.20 to +0.09 | -0.26 to +0.14 | -0.78 | 0.43 | 1.00 |
| Plasma cells | +0.06 (0.07) | +0.18 (0.07) | -0.01 (0.14) | -0.13 | -0.32 to +0.07 | -0.40 to +0.14 | -1.29 | 0.20 | 1.00 |
| T cells | +0.20 (0.09) | +0.35 (0.10) | +0.38 (0.18) | -0.15 | -0.40 to +0.10 | -0.50 to +0.21 | -1.15 | 0.25 | 1.00 |

### T3. Joint Wald tests, primary and validation (results/main_analysis/primary_results.json)
| test | W | df | p |
|---|---|---|---|
| Primary: EAS vs EUR (de Lange) | 4.04 | 8 | 0.853 |
| Validation (i): EUR de Lange vs EUR FinnGen | 7.60 | 8 | 0.473 |
| Validation (ii): EAS vs EUR FinnGen | 7.91 | 8 | 0.442 |

### T4. Registered sensitivity analyses (results/main_analysis/sensitivity.json)
| analysis | W | df | p | largest per-lineage |z| | smallest p Holm | note |
|---|---|---|---|---|---|---|
| a_drop_ImmunoChip_SNPs | 7.15 | 8 | 0.521 | 1.37 | 1.00 | EAS: 994,292 -> 885,845 SNPs; EUR_deLange: 1,152,942 -> 1,046,019 SNPs |
| b_baselineLD_v2.2 | 3.99 | 8 | 0.858 | 1.48 | 1.00 |  |
| c_top5pct | 1.34 | 8 | 0.995 | 0.64 | 1.00 |  |
| c_top20pct | 2.48 | 8 | 0.963 | 0.87 | 1.00 |  |
| d_50kb_windows | 6.62 | 8 | 0.578 | 1.63 | 0.82 |  |

### T4b. Per-lineage tau* in each sensitivity analysis, EAS / EUR (results/main_analysis/sensitivity.json)
| analysis | B cells | Endothelial | Epithelial | Mesenchymal | Myeloid | Neuronal | Plasma cells | T cells |
|---|---|---|---|---|---|---|---|---|
| primary (registered) | +0.04 / +0.08 | -0.10 / -0.07 | -0.06 / -0.04 | -0.01 / -0.03 | +0.23 / +0.20 | -0.17 / -0.11 | +0.06 / +0.18 | +0.20 / +0.35 |
| a_drop_ImmunoChip_SNPs | -0.04 / +0.11 | -0.09 / +0.03 | -0.01 / +0.06 | +0.04 / -0.01 | +0.27 / +0.21 | -0.18 / -0.11 | +0.01 / +0.20 | -0.03 / +0.14 |
| b_baselineLD_v2.2 | +0.04 / +0.09 | -0.09 / -0.09 | -0.06 / -0.07 | +0.00 / -0.02 | +0.16 / +0.17 | -0.13 / -0.11 | +0.02 / +0.17 | +0.15 / +0.34 |
| c_top5pct | -0.02 / -0.01 | -0.01 / -0.02 | +0.07 / +0.01 | +0.00 / -0.03 | +0.08 / +0.12 | -0.12 / -0.08 | +0.15 / +0.12 | +0.34 / +0.42 |
| c_top20pct | +0.17 / +0.16 | -0.03 / -0.11 | -0.08 / -0.07 | -0.00 / +0.00 | +0.21 / +0.20 | -0.26 / -0.20 | +0.17 / +0.19 | +0.24 / +0.34 |
| d_50kb_windows | +0.03 / +0.10 | -0.08 / -0.11 | -0.08 / -0.05 | -0.03 / -0.05 | +0.25 / +0.18 | -0.18 / -0.11 | +0.01 / +0.14 | +0.17 / +0.39 |

### T5. Exploratory analyses (results/main_analysis/exploratory/) — EXPLORATORY
| exploratory analysis | W | df | p | detail |
|---|---|---|---|---|
| E2 primary model without +/-1 Mb around NOD2 | 4.37 | 8 | 0.822 | region chr16:49727514-51766988; SNPs {'EAS': [994292, 993587], 'EUR': [1152942, 1152087]} |
| E3 ulcerative colitis, EAS vs EUR | 6.14 | 8 | 0.632 | regression SNPs prepared: EAS 1,066,794, EUR 1,153,071 |
| E4 small-intestine-only reference | 3.49 | 7 | 0.837 | neuronal lineage below threshold, dropped |
| E4 large-intestine-only reference | 0.96 | 7 | 0.995 | neuronal lineage below threshold, dropped |
| E5 10 donor-bootstrap gene sets | - | 8 | 0.814 to 0.990 | sets of 1,414-1,489 genes after coordinate matching; SD of D across bootstrap sets: 0.028, 0.024, 0.018, 0.035, 0.042, 0.020, 0.028, 0.034 |
| E1 51 fine cell types: correlation of tau* EAS vs EUR | - | - | not reported | Pearson r = 0.73, Spearman rho = 0.76 (gene sets overlap, so fine types are not independent and no p-value is given) |

### T6. E3 ulcerative colitis tau* per lineage (exploratory)
| lineage | UC tau* EAS (SE) | UC tau* EUR (SE) |
|---|---|---|
| B cells | +0.09 (0.08) | +0.10 (0.08) |
| Endothelial | -0.13 (0.08) | -0.02 (0.07) |
| Epithelial | +0.05 (0.10) | +0.21 (0.08) |
| Mesenchymal | -0.04 (0.08) | +0.01 (0.07) |
| Myeloid | +0.22 (0.11) | +0.16 (0.09) |
| Neuronal | -0.15 (0.08) | -0.13 (0.06) |
| Plasma cells | -0.01 (0.09) | +0.06 (0.07) |
| T cells | +0.42 (0.10) | +0.21 (0.08) |

### T7. E1 fine cell types with the highest average tau* (exploratory; 51 types)
| fine type | tau* EAS (SE) | tau* EUR (SE) | z of difference |
|---|---|---|---|
| Activated CD4 T | +0.18 (0.08) | +0.39 (0.09) | -1.79 |
| Treg | +0.25 (0.08) | +0.32 (0.07) | -0.73 |
| Activated CD8 T | +0.20 (0.09) | +0.32 (0.09) | -0.96 |
| NK cell | +0.20 (0.08) | +0.23 (0.08) | -0.30 |
| CD8 Tmem | +0.14 (0.08) | +0.27 (0.09) | -1.06 |
| Th17 | +0.14 (0.08) | +0.25 (0.07) | -1.02 |
| SELL CD8 T | +0.17 (0.08) | +0.19 (0.08) | -0.15 |
| gdT | +0.09 (0.08) | +0.27 (0.08) | -1.48 |

Largest |z| of an EAS-EUR difference across all 51 fine types: 1.79.
### T8. POST HOC: donor-corrected rerun (15 true donors; results/main_analysis/donorfix/) — reported alongside, not replacing, the registered result
| test | W corrected | df | p corrected | W registered | p registered |
|---|---|---|---|---|---|
| Primary: EAS vs EUR (de Lange) | 2.16 | 8 | 0.976 | 4.04 | 0.853 |
| Validation (i): de Lange vs FinnGen | 7.56 | 8 | 0.477 | 7.60 | 0.473 |
| Validation (ii): EAS vs FinnGen | 10.48 | 8 | 0.233 | 7.91 | 0.442 |

| lineage | tau* EAS | tau* EUR de Lange | tau* EUR FinnGen | D | se(D) | p Holm |
|---|---|---|---|---|---|---|
| B cells | +0.08 | +0.11 | +0.10 | -0.03 | 0.11 | 1.00 |
| Endothelial | -0.11 | -0.06 | +0.05 | -0.05 | 0.09 | 1.00 |
| Epithelial | -0.03 | -0.04 | -0.11 | +0.01 | 0.08 | 1.00 |
| Mesenchymal | -0.01 | -0.02 | +0.23 | +0.01 | 0.09 | 1.00 |
| Myeloid | +0.21 | +0.19 | +0.31 | +0.01 | 0.12 | 1.00 |
| Neuronal | -0.15 | -0.11 | -0.06 | -0.04 | 0.08 | 1.00 |
| Plasma cells | +0.10 | +0.14 | +0.05 | -0.04 | 0.11 | 1.00 |
| T cells | +0.21 | +0.34 | +0.39 | -0.13 | 0.13 | 1.00 |

### T9. POST HOC: aligned-block check of Cov(D) = Cov(A) + Cov(B) (results/main_analysis/aligned_block_check.json; 950,331 common SNPs)
| comparison | W summed (registered formula) | p summed | W direct (aligned) | p direct | error correlation across lineages | se(D) direct / summed |
|---|---|---|---|---|---|---|
| EAS vs EUR (de Lange) | 4.39 | 0.820 | 7.22 | 0.513 | 0.24 to 0.40 | 0.78 to 0.87 |
| de Lange vs FinnGen | 7.66 | 0.468 | 7.92 | 0.441 | 0.01 to 0.40 | 0.82 to 1.00 |
| EAS vs FinnGen | 7.25 | 0.510 | 8.46 | 0.389 | 0.09 to 0.33 | 0.87 to 0.96 |

EAS vs EUR on the common SNPs, direct (aligned-block) intervals (post hoc):

| lineage | D (common SNPs) | 95% CI | simultaneous 95% CI (Bonferroni, 8) |
|---|---|---|---|
| B cells | -0.06 | -0.22 to +0.10 | -0.28 to +0.17 |
| Endothelial | -0.00 | -0.14 to +0.13 | -0.19 to +0.19 |
| Epithelial | -0.04 | -0.18 to +0.10 | -0.23 to +0.16 |
| Mesenchymal | +0.01 | -0.13 to +0.14 | -0.18 to +0.20 |
| Myeloid | +0.03 | -0.18 to +0.24 | -0.27 to +0.32 |
| Neuronal | -0.04 | -0.17 to +0.09 | -0.22 to +0.14 |
| Plasma cells | -0.12 | -0.27 to +0.03 | -0.33 to +0.09 |
| T cells | -0.18 | -0.38 to +0.02 | -0.46 to +0.09 |

