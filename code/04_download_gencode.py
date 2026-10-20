"""Download GENCODE release 19 (GRCh37/hg19) gene annotation (gene coordinates).

Source: EBI GENCODE FTP; MD5 checked against the release's MD5SUMS file. Provenance row
appended to data/raw/PROVENANCE.tsv (same format as code/02_download_data.py).
Run from the repo root:  .venv/Scripts/python code/04_download_gencode.py
"""
import datetime as dt
import hashlib
import os
import pathlib
import ssl
import stat
import urllib.request

import certifi

BASE = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19"
NAME = "gencode.v19.annotation.gtf.gz"
OUT = pathlib.Path("data/raw/genes")
PROV = pathlib.Path("data/raw/PROVENANCE.tsv")
CTX = ssl.create_default_context(cafile=certifi.where())


def main():
    dest = OUT / NAME
    if dest.exists():
        print("exists:", dest)
        return
    OUT.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(f"{BASE}/MD5SUMS", context=CTX) as r:
        sums = dict(line.split()[::-1] for line in r.read().decode().splitlines() if line.strip())
    with urllib.request.urlopen(f"{BASE}/{NAME}", context=CTX) as r:
        data = r.read()
    got = hashlib.md5(data).hexdigest()
    if got != sums[NAME]:
        raise SystemExit(f"MD5 mismatch for {NAME}")
    dest.write_bytes(data)
    os.chmod(dest, stat.S_IREAD)
    when = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    # PROVENANCE.tsv is read-only once written; reopen for this one append, then restore
    os.chmod(PROV, stat.S_IREAD | stat.S_IWRITE)
    with open(PROV, "a") as f:
        f.write(f"genes\t{NAME}\t{BASE}/{NAME}\t{when}\t{len(data)}\t{got}\t{sums[NAME]}\tTrue\n")
    os.chmod(PROV, stat.S_IREAD)
    print("OK", NAME, len(data))


if __name__ == "__main__":
    main()
