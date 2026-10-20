"""Archive the IIBDGC readme that documents the East Asian meta-analysis files (cohort order in the Direction column).

Source: Google Drive folder "IBD_EAS/sumstats" linked from https://www.ibdgenetics.org/ ("East Asian and European meta
analysis GWAS"), file readme.sst.txt (Drive id 1taTigC3PM9TwK-RBx8ZSprci59yIW9qO). Saved as
data/raw/gwas/IIBDGC_IBD_EAS_sumstats_readme.sst.txt with MD5 and SHA-256 recorded in data/raw/PROVENANCE.tsv
(Drive publishes no checksum).
Run from the repo root:  .venv/Scripts/python code/22_archive_iibdgc_readme.py
"""
import datetime as dt
import hashlib
import os
import pathlib
import ssl
import stat
import urllib.request

import certifi

URL = "https://drive.google.com/uc?export=download&id=1taTigC3PM9TwK-RBx8ZSprci59yIW9qO"
DEST = pathlib.Path("data/raw/gwas/IIBDGC_IBD_EAS_sumstats_readme.sst.txt")
PROV = pathlib.Path("data/raw/PROVENANCE.tsv")


def main():
    if DEST.exists():
        print("exists:", DEST)
        return
    ctx = ssl.create_default_context(cafile=certifi.where())
    with urllib.request.urlopen(URL, context=ctx, timeout=60) as r:
        data = r.read()
    if data.lstrip()[:1] == b"<":
        raise SystemExit("got an HTML page instead of the text file")
    DEST.write_bytes(data)
    os.chmod(DEST, stat.S_IREAD)
    md5, sha = hashlib.md5(data).hexdigest(), hashlib.sha256(data).hexdigest()
    when = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    os.chmod(PROV, stat.S_IREAD | stat.S_IWRITE)
    with open(PROV, "a") as f:
        f.write(f"gwas\t{DEST.name}\t{URL}\t{when}\t{len(data)}\t{md5}\t\tn/a (sha256 {sha})\n")
    os.chmod(PROV, stat.S_IREAD)
    print("OK", DEST, len(data), "sha256", sha)
    print(data.decode("utf-8", "replace"))


if __name__ == "__main__":
    main()
