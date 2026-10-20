"""Download FinnGen R13 GWAS summary statistics (K11_CD_STRICT2 is the European validation dataset).

Source: public bucket gs://finngen-public-data-r13 (access instructions emailed by
finngen-servicedesk@helsinki.fi on 2026-09-07). Each file is checked against the
MD5 that Google Cloud Storage reports, and a provenance row (URL, UTC download
time, bytes, MD5) is appended to data/raw/finngen_r13/PROVENANCE.tsv.

Run from the repo root:  python code/01_download_finngen.py
"""
import base64
import datetime as dt
import hashlib
import json
import pathlib
import urllib.request

ENDPOINTS = ["CHRONSMALL", "CHRONLARGE", "K11_UC_STRICT2", "K11_CD_STRICT2"]
BUCKET = "finngen-public-data-r13"
OUT = pathlib.Path("data/raw/finngen_r13")
MANIFEST_URL = f"https://storage.googleapis.com/{BUCKET}/summary_stats/finngen_R13_manifest.tsv"


def remote_md5(obj):
    meta_url = f"https://storage.googleapis.com/storage/v1/b/{BUCKET}/o/{urllib.request.quote(obj, safe='')}"
    with urllib.request.urlopen(meta_url) as r:
        meta = json.load(r)
    return base64.b64decode(meta["md5Hash"]).hex(), int(meta["size"])


def local_md5(path):
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            md5.update(chunk)
    return md5.hexdigest()


def fetch(url, dest):
    md5 = hashlib.md5()
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as r, open(tmp, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
            md5.update(chunk)
    tmp.rename(dest)
    return md5.hexdigest(), dest.stat().st_size


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    prov = OUT / "PROVENANCE.tsv"
    if not prov.exists():
        prov.write_text("file\turl\tdownloaded_utc\tbytes\tmd5\tmd5_matches_gcs\n")
    objs = ["summary_stats/finngen_R13_manifest.tsv"]
    for ep in ENDPOINTS:
        objs += [f"summary_stats/finngen_R13_{ep}.gz", f"summary_stats/finngen_R13_{ep}.gz.tbi"]
    for obj in objs:
        dest = OUT / pathlib.Path(obj).name
        url = f"https://storage.googleapis.com/{BUCKET}/{obj}"
        want_md5, want_size = remote_md5(obj)
        if dest.exists():
            # A file already on disk is kept only if it matches the bucket exactly.
            if dest.name in {ln.split("\t")[0] for ln in prov.read_text().splitlines()}:
                print("exists, already recorded:", dest)
                continue
            got_md5, got_size = local_md5(dest), dest.stat().st_size
            if (got_md5, got_size) != (want_md5, want_size):
                raise SystemExit(f"{dest} exists but does not match the bucket; move it aside and rerun")
        else:
            got_md5, got_size = fetch(url, dest)
        ok = got_md5 == want_md5 and got_size == want_size
        when = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        with open(prov, "a") as f:
            f.write(f"{dest.name}\t{url}\t{when}\t{got_size}\t{got_md5}\t{ok}\n")
        print(("OK  " if ok else "MISMATCH  ") + dest.name, got_size)
        if not ok:
            raise SystemExit(f"checksum mismatch for {dest}")


if __name__ == "__main__":
    main()
