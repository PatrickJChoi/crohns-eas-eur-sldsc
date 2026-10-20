"""Download the public data (Crohn's disease heritability by cell type,
East Asian vs European ancestry).

- GWAS summary statistics (GWAS Catalog FTP):
    GCST90446793  East Asian CD (Liu et al. 2023), GRCh38, MD5 published by the Catalog
    GCST004132    European CD (de Lange et al. 2017), GRCh37, no MD5 published
- S-LDSC reference files (Zenodo record 10515792; MD5 published by Zenodo):
    1000 Genomes Phase 3 EAS and EUR: baselineLD v2.2 and baseline v1.2 LD scores,
    PLINK files, regression weights, EUR allele frequencies, HapMap3 SNP list
- Gut Cell Atlas full object, raw counts (Elmentaite et al. 2021; cellgeni.cog.sanger.ac.uk)
- UCSC hg38 -> hg19 liftOver chain (MD5 from UCSC's md5sum.txt)

Every file's MD5 is computed locally and compared with the published value where one
exists; rows go to data/raw/PROVENANCE.tsv.
Run from the repo root:  python code/02_download_data.py
"""
import datetime as dt
import hashlib
import json
import pathlib
import ssl
import urllib.request

import certifi

OUT = pathlib.Path("data/raw")
CAT = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
ZEN = "https://zenodo.org/records/10515792/files"
ZEN_FILES = [
    "1000G_Phase3_EAS_baselineLD_v2.2_ldscores.tgz", "1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz",
    "1000G_Phase3_EAS_plinkfiles.tgz", "1000G_Phase3_EAS_weights_hm3_no_MHC.tgz",
    "1000G_Phase3_baselineLD_v2.2_ldscores.tgz", "1000G_Phase3_baseline_v1.2_ldscores.tgz",
    "1000G_Phase3_plinkfiles.tgz", "1000G_Phase3_weights_hm3_no_MHC.tgz", "1000G_Phase3_frq.tgz",
    "hm3_no_MHC.list.txt",
]
UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver"
# Certificates are verified against certifi's CA bundle; the Windows default store here lacks
# an intermediate certificate that UCSC's server does not send.
CTX = ssl.create_default_context(cafile=certifi.where())


def text(url):
    with urllib.request.urlopen(url, context=CTX) as r:
        return r.read().decode()


def published_md5():
    """Map file name -> published MD5 (None when the source publishes none)."""
    md5 = {}
    for line in text(f"{CAT}/GCST90446001-GCST90447000/GCST90446793/md5sum.txt").splitlines():
        h, name = line.split()
        md5[name] = h
    # The Catalog's listed MD5 for the small metadata yaml is out of date (the yaml was edited on
    # 2025-01-05; three fresh downloads agree with each other but not with the list). The data file's
    # MD5 does match, so only the yaml is treated as having no usable published checksum.
    md5.pop("GCST90446793.tsv-meta.yaml", None)
    zen = json.loads(text("https://zenodo.org/api/records/10515792"))
    for f in zen["files"]:
        md5[f["key"]] = f["checksum"].split(":")[1]
    for line in text(f"{UCSC}/md5sum.txt").splitlines():
        h, name = line.split()
        md5[name] = h
    return md5


def plan():
    yield "gwas", "GCST90446793.tsv", f"{CAT}/GCST90446001-GCST90447000/GCST90446793/GCST90446793.tsv"
    yield "gwas", "GCST90446793.tsv-meta.yaml", f"{CAT}/GCST90446001-GCST90447000/GCST90446793/GCST90446793.tsv-meta.yaml"
    yield "gwas", "cd_build37_40266_20161107.txt.gz", f"{CAT}/GCST004001-GCST005000/GCST004132/cd_build37_40266_20161107.txt.gz"
    for f in ZEN_FILES:
        yield "ldsc_ref", f, f"{ZEN}/{f}?download=1"
    yield "atlas", "Full_obj_raw_counts_nosoupx_v2.h5ad", \
        "https://cellgeni.cog.sanger.ac.uk/gutcellatlas/Full_obj_raw_counts_nosoupx_v2.h5ad"
    yield "liftover", "hg38ToHg19.over.chain.gz", f"{UCSC}/hg38ToHg19.over.chain.gz"


def fetch(url, dest):
    md5 = hashlib.md5()
    tmp = dest.with_name(dest.name + ".part")
    with urllib.request.urlopen(url, context=CTX) as r, open(tmp, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
            md5.update(chunk)
    tmp.rename(dest)
    return md5.hexdigest(), dest.stat().st_size


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    prov = OUT / "PROVENANCE.tsv"
    if not prov.exists():
        prov.write_text("subdir\tfile\turl\tdownloaded_utc\tbytes\tmd5_local\tmd5_published\tmatch\n")
    done = {tuple(ln.split("\t")[:2]) for ln in prov.read_text().splitlines()}
    want = published_md5()
    for sub, name, url in plan():
        dest = OUT / sub / name
        dest.parent.mkdir(exist_ok=True)
        if (sub, name) in done and dest.exists():
            print("exists, already recorded:", sub, name)
            continue
        got, size = fetch(url, dest)
        pub = want.get(name)
        match = "n/a" if pub is None else str(got == pub)
        when = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        with open(prov, "a") as f:
            f.write(f"{sub}\t{name}\t{url}\t{when}\t{size}\t{got}\t{pub or ''}\t{match}\n")
        print(f"{'OK' if match != 'False' else 'MISMATCH'}  {sub}/{name}  {size:,} bytes  md5 match: {match}")
        if match == "False":
            raise SystemExit(f"checksum mismatch: {dest}")


if __name__ == "__main__":
    main()
