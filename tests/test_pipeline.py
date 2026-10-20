"""Correctness and leakage tests.

Leakage in this design means: gene sets or annotations built with any GWAS data, one ancestry's reference data used
for the other, or the real GWAS outcome columns read before the primary run. Run from the repo root:
    .venv/Scripts/python -m unittest tests/test_pipeline.py -v
"""
import hashlib
import importlib.util
import pathlib
import subprocess
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
sys.path.insert(0, str(CODE))
PY = sys.executable


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", "").lstrip("0123456789_"), CODE / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


prep = load("07_prepare_sumstats.py")
test = load("09_primary_test.py")


class Audit:
    """Record every file opened and every subprocess launched while active (sys.addaudithook is permanent, so the
    hook checks a flag)."""
    active, opened, spawned = False, [], []

    @classmethod
    def hook(cls, event, args):
        if not cls.active:
            return
        if event == "open" and isinstance(args[0], (str, pathlib.Path)):
            cls.opened.append(str(args[0]).replace("\\", "/"))
        elif event == "subprocess.Popen":
            cmd = args[1] if isinstance(args[1], str) else " ".join(map(str, args[1]))  # Windows passes one string
            cls.spawned.append(cmd.replace("\\", "/"))


sys.addaudithook(Audit.hook)


class TestArithmetic(unittest.TestCase):
    def test_eas_effective_n(self):
        w = [4 / (1 / a + 1 / b) for a, b in [(2552, 6279), (1611, 3724), (1619, 4419), (1590, 1034)]]
        got = prep.eas_neff_from_direction(["+-+-", "?+??", "????", "+?+-"])
        np.testing.assert_allclose(got, [sum(w), w[1], 0, w[0] + w[2] + w[3]])

    def test_allele_match_and_ambiguous(self):
        s = lambda *v: pd.Series(v)
        got = prep.allele_match(s("a", "a", "c", "a"), s("g", "g", "t", "c"), s("A", "G", "G", "T"),
                                s("G", "A", "A", "T"))
        self.assertEqual(list(got), [1, -1, 1, 0])  # same; swapped; strand flip (c/t -> G/A); inconsistent
        self.assertEqual(list(prep.is_ambiguous(s("A", "C", "A"), s("T", "G", "G"))), [True, True, False])

    def test_prepare_fin_on_synthetic_finngen_rows(self):
        # rs62224609 is C/T and rs587646183 is C/T in 1000G EUR chr22 (bim A1/A2)
        import tempfile
        rows = pd.DataFrame({"#chrom": ["22"] * 4, "pos": [1, 2, 3, 4],
                             "ref": ["T", "C", "A", "A"], "alt": ["C", "T", "G", "G"],
                             "rsids": ["rs62224609", "rs587646183", "rsNOTREAL", "rsNOTREAL2,rsALSO_NOT"],
                             "beta": [0.1, 0.2, 0.3, 0.4], "sebeta": [0.05, 0.1, 0.1, 0.1]})
        with tempfile.TemporaryDirectory() as d:
            f = pathlib.Path(d) / "fin.tsv.gz"
            rows.to_csv(f, sep="\t", index=False)
            out = prep.prepare_fin(f, {"rs62224609", "rs587646183", "rsNOTREAL"}).set_index("SNP")
        self.assertEqual(sorted(out.index), ["rs587646183", "rs62224609"])
        # alt is the effect allele: alt=C equals bim A1=C -> A1 C; alt=T equals bim A2 -> reported as A1 T
        self.assertEqual((out.loc["rs62224609", "A1"], out.loc["rs62224609", "A2"]), ("C", "T"))
        self.assertEqual((out.loc["rs587646183", "A1"], out.loc["rs587646183", "A2"]), ("T", "C"))
        np.testing.assert_allclose(out.loc[["rs62224609", "rs587646183"], "Z"], [2.0, 2.0])
        self.assertTrue((out["N"].round() == 10378).all())

    def test_jackknife_cov_matches_ldsc_formula(self):
        rng = np.random.default_rng(0)
        est, dele = rng.normal(size=3), rng.normal(size=(200, 3))
        pseudo = 200 * est - 199 * dele
        np.testing.assert_allclose(test.jackknife_cov(est, dele), np.cov(pseudo.T, ddof=1) / 200)

    def test_wald_equals_sum_of_z2_when_covariance_is_diagonal(self):
        # Build delete values whose pseudovalue deviations are orthogonal across lineages, so each ancestry's
        # jackknife covariance is exactly diagonal; then W must equal the sum of squared per-lineage z computed
        # independently from the per-lineage variances.
        B = 4

        def deletes(est, c1, c2):
            dev = np.array([[c1, 0], [-c1, 0], [0, c2], [0, -c2]], float)  # orthogonal columns, mean zero
            pseudo = np.asarray(est)[None, :] + dev
            return ((B * np.asarray(est)[None, :] - pseudo) / (B - 1)).tolist()

        a = {"names": ["x", "y"], "tau_star": [0.30, -0.10], "delete": deletes([0.30, -0.10], 0.4, 0.2)}
        b = {"names": ["x", "y"], "tau_star": [0.10, 0.00], "delete": deletes([0.10, 0.00], 0.3, 0.6)}
        var = lambda c: 2 * c ** 2 / (B - 1) / B  # sample variance of [c, -c, 0, 0] over B
        z1 = (0.30 - 0.10) / np.sqrt(var(0.4) + var(0.3))
        z2 = (-0.10 - 0.00) / np.sqrt(var(0.2) + var(0.6))
        r = test.wald(a, b)
        self.assertAlmostEqual(r["W"], z1 ** 2 + z2 ** 2)
        np.testing.assert_allclose(r["z"], [z1, z2])

    def test_tau_star_arithmetic(self):
        fm = load("08_fit_models.py")
        # 3 annotations: base (all 1000 common SNPs), control (400), lineage (100); tau per SNP-unit
        tau, M = np.array([1e-4, 2e-4, 5e-4]), np.array([1000.0, 400.0, 100.0])
        part = np.tile(tau, (5, 1)) * np.array([[1.0], [1.1], [0.9], [1.0], [1.2]])
        tot = part @ M
        r = fm.tau_star_from(tau, M, part, tot)
        h2 = 1e-4 * 1000 + 2e-4 * 400 + 5e-4 * 100  # = 0.23
        sd = np.sqrt(0.1 * 0.9)
        self.assertAlmostEqual(r["h2"], h2)
        self.assertAlmostEqual(r["tau_star"], 5e-4 * sd * 1000 / h2)
        # scaling every block's coefficients by the same factor leaves per-block tau* unchanged
        np.testing.assert_allclose(r["delete"], [5e-4 * sd * 1000 / h2] * 5)

    def test_holm(self):
        np.testing.assert_allclose(test.holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])


class TestBlindingAndLeakage(unittest.TestCase):
    def test_real_gwas_blocked_without_unblind(self):
        for src, f in [("EAS", "data/raw/gwas/GCST90446793.tsv"),
                       ("EUR", "data/raw/gwas/cd_build37_40266_20161107.txt.gz")]:
            r = subprocess.run([PY, str(CODE / "07_prepare_sumstats.py"), src, f, str(ROOT / "results/_x")],
                               cwd=ROOT, capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("Blinded", r.stderr + r.stdout)

    def test_step_a_reads_no_gwas_and_is_reproducible(self):
        # Rerun Step A into a temporary folder (production files untouched) and compare with the production manifest.
        import os
        import tempfile
        manifest = ROOT / "data/processed/genesets/MANIFEST.sha256"
        step_a = load("05_lineage_genesets.py")
        with tempfile.TemporaryDirectory() as tmp:
            step_a.OUT, step_a.RES = pathlib.Path(tmp) / "genesets", pathlib.Path(tmp) / "res"
            Audit.opened.clear()
            Audit.active = True
            cwd = os.getcwd()
            try:
                os.chdir(ROOT)
                step_a.main()
            finally:
                Audit.active = False
                os.chdir(cwd)
            rerun = (step_a.OUT / "MANIFEST.sha256").read_text()
        self.assertFalse([p for p in Audit.opened if "/gwas/" in p or "finngen" in p.lower()])
        self.assertEqual(manifest.read_text(), rerun, "gene sets changed when Step A was rerun")

    def test_step_b_uses_only_own_ancestry_and_no_gwas(self):
        # Rerun Step B for EAS chr22 into a temporary folder and compare one LD-score file with production.
        import tempfile
        step_b = load("06_annotations_ldscores.py")
        sets = step_b.gene_sets()
        prod = step_b.ANNOT
        with tempfile.TemporaryDirectory() as tmp:
            step_b.ANNOT = pathlib.Path(tmp)
            Audit.opened.clear()
            Audit.spawned.clear()
            Audit.active = True
            try:
                step_b.run_chrom(("EAS", 22, sets))
            finally:
                Audit.active = False
            name = "L_Myeloid_top10"
            new = pd.read_csv(pathlib.Path(tmp) / "EAS" / name / f"{name}.22.l2.ldscore.gz", sep="\t")
        old = pd.read_csv(prod / "EAS" / name / f"{name}.22.l2.ldscore.gz", sep="\t")
        pd.testing.assert_frame_equal(new, old)
        touched = Audit.opened + Audit.spawned
        self.assertFalse([p for p in touched if "/gwas/" in p or "finngen" in p.lower()])
        self.assertFalse([p for p in touched if "EUR" in p or "1000G_Phase3_plinkfiles" in p
                          or "1000G_Phase3_frq" in p or "annot/EUR" in p])
        self.assertTrue(any("1000G.EAS.QC.22" in p for p in Audit.spawned))

    def test_genesets_match_manifest(self):
        d = ROOT / "data/processed/genesets"
        for line in (d / "MANIFEST.sha256").read_text().splitlines():
            h, name = line.split()
            self.assertEqual(hashlib.sha256((d / name).read_bytes()).hexdigest(), h, name)


if __name__ == "__main__":
    unittest.main()
