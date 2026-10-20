# Crohn's disease heritability by gut cell lineage in East Asian and European ancestry

Code, logs and results for a pre-registered stratified LD score regression (S-LDSC / LDSC-SEG) study. It asks
whether Crohn's disease heritability is distributed differently across gut cell lineages in East Asian and European
ancestry.

- **Archived code (v1.0.0):** doi:[10.5281/zenodo.23491800](https://doi.org/10.5281/zenodo.23491800)
- **Pre-registration:** doi:[10.17605/OSF.IO/HQ4A8](https://doi.org/10.17605/OSF.IO/HQ4A8). The plan was registered
  before any association result was read.
- **Primary result:** joint Wald test of East Asian minus European differences in τ\* across 8 lineages,
  W = 4.04, df = 8, p = 0.853. This is no evidence of an ancestry difference; it is not proof of equality.
- **Manuscript:** [docs/manuscript.md](docs/manuscript.md). Every number in it is checked against the result files by
  [code/27_claims_audit.py](code/27_claims_audit.py), with the output in [docs/claims_audit.md](docs/claims_audit.md).

Author: Patrick Choi, Department of Integrative Biology, University of South Florida
(ORCID [0009-0004-1051-2074](https://orcid.org/0009-0004-1051-2074)).

## Data

All inputs are public. Download URLs, dates and checksums are recorded in `data/raw/PROVENANCE.tsv` and
`data/raw/finngen_r13/PROVENANCE.tsv`.

| Input | Source |
|---|---|
| East Asian Crohn's disease GWAS (Liu 2023) | GWAS Catalog GCST90446793 |
| European Crohn's disease GWAS (de Lange 2017) | GWAS Catalog GCST004132 |
| European validation GWAS | FinnGen R13, endpoint K11_CD_STRICT2 |
| Gut single-cell reference | Gut Cell Atlas (Elmentaite 2021), full object |
| LD reference | 1000 Genomes Phase 3 EAS and EUR S-LDSC files, Zenodo 10515792 |
| Gene coordinates | GENCODE v19 |
| Exploratory ulcerative colitis GWAS | GCST90446795 (East Asian), GCST004133 (European) |

The large downloaded and derived data (~35 GB) are not in this repository; the scripts rebuild them.

## Environment

See [ENVIRONMENT.md](ENVIRONMENT.md). In short:

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements-analysis.txt
python -m venv .venv-ldsc
.venv-ldsc/Scripts/python -m pip install -r requirements-ldsc.txt
git clone https://github.com/CBIIT/ldsc.git tools/ldsc
git -C tools/ldsc checkout 1f09cf0
```

Paths above are for Windows; on Linux or macOS use `.venv/bin/python`. Scripts are run from the repository root.

## Scripts, in run order

| Script | What it does | Writes |
|---|---|---|
| `01_download_finngen.py` | FinnGen R13 summary statistics | `data/raw/finngen_r13/` |
| `02_download_data.py` | GWAS, atlas, LD references, liftover chain | `data/raw/` |
| `03_qc_data.py` | Data feasibility checks (no association results) | `results/data_qc/` |
| `04_download_gencode.py` | GENCODE v19 annotation | `data/raw/genes/` |
| `05_lineage_genesets.py` | Step A: lineage-specific gene sets from donor pseudobulks | `data/processed/genesets/` |
| `06_annotations_ldscores.py` | Step B: SNP annotations and LD scores | `data/processed/annot/` |
| `07_prepare_sumstats.py` | Step C: LDSC-ready summary statistics | `data/processed/real/` |
| `08_fit_models.py` | Step D: LDSC-SEG models and τ\* | (called by later scripts) |
| `09_primary_test.py` | Step E: joint Wald test | (called by later scripts) |
| `10_simulate.py`, `11_run_simulations.py`, `12_simulation_summary.py` | Pipeline validation on simulated data | `results/simulations/` |
| `13_primary_analysis.py` | The single registered primary run | `results/main_analysis/` |
| `14_sensitivity_analyses.py` | Registered sensitivity analyses | `results/main_analysis/sensitivity/` |
| `15_download_uc.py` … `18_exploratory_analyses.py` | Registered exploratory analyses | `results/main_analysis/exploratory/` |
| `19_report_tables.py` | Tables for [docs/results.md](docs/results.md) | `results/main_analysis/report_tables.md` |
| `20_donorfix_genesets.py`, `21_donorfix_ldscores.py`, `24_donorfix_primary.py` | Post hoc: rerun with true donor IDs | `results/main_analysis/donorfix/` |
| `22_archive_iibdgc_readme.py` | Archives the readme that fixes the cohort order | `data/raw/gwas/` |
| `23_aligned_block_check.py` | Post hoc: block-aligned jackknife covariance check | `results/main_analysis/aligned_block_check.json` |
| `25_verify_references.py` | Checks every manuscript reference on Crossref | `results/manuscript/reference_check.tsv` |
| `26_figures.py` | Manuscript figures | `results/manuscript/` |
| `27_claims_audit.py` | Checks every number in the manuscript against its source file | `docs/claims_audit.md` |

Tests: `.venv/Scripts/python -m unittest tests/test_pipeline.py`. They check correctness and that no test data leak
into fitted steps.

## Documents

| File | Content |
|---|---|
| [docs/analysis_plan.md](docs/analysis_plan.md) | The analysis plan |
| [docs/data_feasibility.md](docs/data_feasibility.md) | Data inventory, provenance and checks done before the plan |
| [docs/pipeline_validation.md](docs/pipeline_validation.md) | Simulation results: false-positive rate, planted-signal recovery, power |
| [docs/results.md](docs/results.md) | Full results: primary, validation, sensitivity, exploratory, post hoc |
| [docs/manuscript.md](docs/manuscript.md) | Manuscript |
| [DECISIONS.md](DECISIONS.md) | Dated log of methodological decisions and deviations from the plan |
| [docs/osf_registration/](docs/osf_registration/) | Byte-exact copies of the files registered on OSF |

## Names in the registered files

The files in `docs/osf_registration/` and the run logs in `results/` were written before the repository was tidied,
so they use older names. The code and results are otherwise unchanged.

| Older name | Current name |
|---|---|
| Study R2-1 | this study |
| `code/05_…` to `code/30_…` | renumbered `code/02_…` to `code/27_…` (same order; e.g. `16_phase5_primary.py` → `13_primary_analysis.py`) |
| `code/r21common.py`, `tests/test_r21.py` | `code/common.py`, `tests/test_pipeline.py` |
| `data/raw/r21/`, `data/processed/r21/` | `data/raw/`, `data/processed/` |
| `results/phase2_R21_qc/`, `phase4_R21/`, `phase5_R21/`, `phase6/` | `results/data_qc/`, `simulations/`, `main_analysis/`, `manuscript/` |
| `docs/ANALYSIS_PLAN_R21.md`, `PHASE2_R21_FEASIBILITY.md`, `PHASE4_R21_REPORT.md`, `PHASE5_R21_RESULTS.md` | `docs/analysis_plan.md`, `data_feasibility.md`, `pipeline_validation.md`, `results.md` |
| Phase 2 / 4 / 5 | data feasibility / pipeline validation / main analysis |

## License

Code is released under the [MIT License](LICENSE). Input data remain under the terms of their original sources.
