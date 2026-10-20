# Software environment

Python 3.11.5 (Windows 11). Seeds are set inside each script that uses randomness.

## Environments
- System Python: numpy 2.1.3, pandas 2.2.3, scipy 1.17.1, certifi 2023.07.22. Used by code/01–02.
- `.venv` (analysis), pins in `requirements-analysis.txt`: numpy 2.1.3, pandas 2.2.3, scipy 1.17.1, h5py 3.16.0,
  pyliftover 0.4.1, certifi 2026.7.22. Used by code/03 onward.
- `.venv-ldsc` (LDSC only), pins in `requirements-ldsc.txt`: numpy 1.23.5, pandas 1.5.3, scipy 1.9.3, bitarray 2.9.3.

## Tools
- LDSC: https://github.com/CBIIT/ldsc at commit 1f09cf0c5e61a38a281016fefcc1e0583cc78402 (2026-08-10), cloned to
  `tools/ldsc/` (not tracked in git). Rebuild: `git clone https://github.com/CBIIT/ldsc.git tools/ldsc` then
  `git -C tools/ldsc checkout 1f09cf0`.

Recreate an environment: `python -m venv .venv && .venv/Scripts/python -m pip install -r requirements-analysis.txt`.
