"""Shared paths and constants (registered plan doi:10.17605/OSF.IO/HQ4A8)."""
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
REF = ROOT / "data/processed/ref"
GENESETS = ROOT / "data/processed/genesets"
ANNOT = ROOT / "data/processed/annot"
LDSC_PY = ROOT / ".venv-ldsc/Scripts/python.exe"
LDSC_DIR = ROOT / "tools/ldsc"
CHROMS = range(1, 23)
LINEAGES = ["B cells", "Endothelial", "Epithelial", "Mesenchymal", "Myeloid", "Neuronal", "Plasma cells", "T cells"]

# Ancestry-matched 1000 Genomes Phase 3 reference files (Zenodo 10515792)
ANC = {
    "EAS": {
        "bfile": REF / "1000G_Phase3_EAS_plinkfiles/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.{chr}",
        "frq": REF / "1000G_Phase3_EAS_plinkfiles/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.{chr}",
        "baseline": REF / "1000G_Phase3_EAS_baseline_v1.2_ldscores/baseline.{chr}",
        "baselineLD": REF / "1000G_Phase3_EAS_baselineLD_v2.2_ldscores/baselineLD.{chr}",
        "weights": REF / "1000G_Phase3_EAS_weights_hm3_no_MHC/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC.{chr}",
    },
    "EUR": {
        "bfile": REF / "1000G_Phase3_plinkfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}",
        "frq": REF / "1000G_Phase3_frq/1000G_Phase3_frq/1000G.EUR.QC.{chr}",
        "baseline": REF / "1000G_Phase3_baseline_v1.2_ldscores/baseline_v1.2/baseline.{chr}",
        "baselineLD": REF / "1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.{chr}",
        "weights": REF / "1000G_Phase3_weights_hm3_no_MHC/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.{chr}",
    },
}


def slug(lineage):
    return lineage.replace(" ", "_")
