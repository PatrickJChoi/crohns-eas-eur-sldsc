"""Step E of the registered plan: joint Wald test of ancestry differences in tau* across lineages.

For each ancestry, the K x K covariance of tau* comes from block-jackknife pseudovalues of the K separate model fits
(identical SNP blocks within an ancestry): pseudo_b = B * est - (B - 1) * delete_b, cov = cov(pseudo) / B (as in LDSC).
D = tau*_A - tau*_B; Cov(D) = Cov_A + Cov_B (independent people); W = D' Cov(D)^-1 D ~ chi2(K).
Secondary: per-lineage z = D_L / se(D_L), Holm-adjusted.

Library use: wald(res_a, res_b) with dicts from code/08_fit_models.fit_set.
CLI: .venv/Scripts/python code/09_primary_test.py A.json B.json
"""
import json
import sys

import numpy as np
from scipy import stats


def jackknife_cov(est, delete):
    est, delete = np.asarray(est, float), np.asarray(delete, float)
    B = delete.shape[0]
    pseudo = B * est[None, :] - (B - 1) * delete
    return np.atleast_2d(np.cov(pseudo.T, ddof=1)) / B


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(p) - rank) * p[i]))
        adj[i] = running
    return adj


def wald(res_a, res_b):
    if res_a["names"] != res_b["names"]:
        raise ValueError("gene sets differ between the two fits")
    ta, tb = np.asarray(res_a["tau_star"]), np.asarray(res_b["tau_star"])
    ca, cb = jackknife_cov(ta, res_a["delete"]), jackknife_cov(tb, res_b["delete"])
    D, C = ta - tb, ca + cb
    W = float(D @ np.linalg.solve(C, D))
    K = len(D)
    se = np.sqrt(np.diag(C))
    z = D / se
    p_l = 2 * stats.norm.sf(np.abs(z))
    return {"names": res_a["names"], "K": K, "W": W, "p": float(stats.chi2.sf(W, K)),
            "D": D.tolist(), "se_D": se.tolist(), "z": z.tolist(), "p_lineage": p_l.tolist(),
            "p_holm": holm(p_l).tolist(), "se_a": np.sqrt(np.diag(ca)).tolist(), "se_b": np.sqrt(np.diag(cb)).tolist()}


if __name__ == "__main__":
    a, b = (json.loads(open(f).read()) for f in sys.argv[1:3])
    print(json.dumps(wald(a, b), indent=1))
