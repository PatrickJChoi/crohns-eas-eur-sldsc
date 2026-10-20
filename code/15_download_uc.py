"""Download the ulcerative colitis GWAS files for registered exploratory analysis 3.

GCST90446795  East Asian UC (Liu et al. 2023), GRCh38, MD5 from the Catalog's md5sum.txt
GCST004133    European UC (de Lange et al. 2017), GRCh37, no MD5 published
Rows appended to data/raw/PROVENANCE.tsv (same format as code/02).
Run from the repo root:  .venv/Scripts/python code/15_download_uc.py
"""
import datetime as dt
import hashlib
import os
import pathlib
import ssl
import stat
import urllib.request

import certifi

CAT = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
FILES = [("GCST90446795.tsv", f"{CAT}/GCST90446001-GCST90447000/GCST90446795"),
         ("uc_build37_45975_20161107.txt.gz", f"{CAT}/GCST004001-GCST005000/GCST004133")]
OUT = pathlib.Path("data/raw/gwas")
PROV = pathlib.Path("data/raw/PROVENANCE.tsv")
CTX = ssl.create_default_context(cafile=certifi.where())


def file_md5(path):
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            md5.update(chunk)
    return md5.hexdigest()


def fetch_resumable(url, tmp, tries=20):
    """Download with a 60 s socket timeout, resuming from the partial file (HTTP Range) after a stall or error.
    The MD5 is computed over the finished file, so resumed downloads are still verified."""
    for attempt in range(tries):
        have = tmp.stat().st_size if tmp.exists() else 0
        req = urllib.request.Request(url, headers={"Range": f"bytes={have}-"} if have else {})
        try:
            with urllib.request.urlopen(req, context=CTX, timeout=60) as r, open(tmp, "ab" if have else "wb") as f:
                if have and r.status != 206:
                    raise RuntimeError("server ignored the Range request")
                while chunk := r.read(1 << 20):
                    f.write(chunk)
            return
        except (OSError, TimeoutError) as e:
            print(f"attempt {attempt + 1}: {e!r}; resuming from {tmp.stat().st_size if tmp.exists() else 0:,} bytes",
                  flush=True)
    raise SystemExit(f"download failed after {tries} attempts: {url}")


def main():
    sums = {}
    with urllib.request.urlopen(f"{FILES[0][1]}/md5sum.txt", context=CTX) as r:
        for line in r.read().decode().splitlines():
            h, n = line.split()
            sums[n] = h
    for name, base in FILES:
        dest = OUT / name
        if dest.exists():
            print("exists:", dest)
            continue
        tmp = dest.with_name(name + ".part")
        fetch_resumable(f"{base}/{name}", tmp)
        got, pub = file_md5(tmp), sums.get(name)
        if pub and got != pub:
            raise SystemExit(f"MD5 mismatch for {name}")
        tmp.rename(dest)
        os.chmod(dest, stat.S_IREAD)
        when = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        os.chmod(PROV, stat.S_IREAD | stat.S_IWRITE)
        with open(PROV, "a") as f:
            f.write(f"gwas\t{name}\t{base}/{name}\t{when}\t{dest.stat().st_size}\t{got}\t{pub or ''}\t"
                    f"{'True' if pub else 'n/a'}\n")
        os.chmod(PROV, stat.S_IREAD)
        print("OK", name, dest.stat().st_size, "md5 match:", pub == got if pub else "n/a", flush=True)


if __name__ == "__main__":
    main()
