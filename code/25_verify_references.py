"""Verify every manuscript reference by DOI on Crossref.

For each DOI: fetch https://api.crossref.org/works/<doi>, compare the first author's family name, year and a title
fragment with the citation as written in the manuscript. Writes results/manuscript/reference_check.tsv.
Run from the repo root:  .venv/Scripts/python code/25_verify_references.py
"""
import json
import pathlib
import ssl
import time
import urllib.request

import certifi

# key: (doi, expected first-author family name, expected year, expected title fragment, as cited)
REFS = {
    "Liu2023": ("10.1038/s41588-023-01384-0", "Liu", 2023, "Genetic architecture of the inflammatory bowel diseases"),
    "deLange2017": ("10.1038/ng.3760", "de Lange", 2017, "integrin genes in inflammatory bowel disease"),
    "Kurki2023": ("10.1038/s41586-022-05473-8", "Kurki", 2023, "FinnGen provides genetic insights"),
    "Elmentaite2021": ("10.1038/s41586-021-03852-1", "Elmentaite", 2021, "Cells of the human intestinal tract"),
    "BulikSullivan2015": ("10.1038/ng.3211", "Bulik-Sullivan", 2015, "LD Score regression distinguishes confounding"),
    "Finucane2015": ("10.1038/ng.3404", "Finucane", 2015, "Partitioning heritability by functional annotation"),
    "Finucane2018": ("10.1038/s41588-018-0081-4", "Finucane", 2018, "Heritability enrichment of specifically expressed"),
    "Gazal2017": ("10.1038/ng.3954", "Gazal", 2017, "Linkage disequilibrium"),
    "Skene2016": ("10.3389/fnins.2016.00016", "Skene", 2016, "Identification of Vulnerable Cell Types"),
    "Bryois2020": ("10.1038/s41588-020-0610-9", "Bryois", 2020, "Genetic identification of cell types"),
    "Jagadeesh2022": ("10.1038/s41588-022-01187-9", "Jagadeesh", 2022, "Identifying disease-critical cell types"),
    "Squair2021": ("10.1038/s41467-021-25960-2", "Squair", 2021, "Confronting false discoveries"),
    "Brown2016": ("10.1016/j.ajhg.2016.05.001", "Brown", 2016, "Transethnic Genetic-Correlation Estimates"),
    "Shi2021": ("10.1038/s41467-021-21286-1", "Shi", 2021, "Population-specific causal disease effect sizes"),
    "Amariuta2020": ("10.1038/s41588-020-00740-8", "Amariuta", 2020, "trans-ancestry portability of polygenic risk"),
    "Krzak2026": ("10.1038/s41588-026-02634-7", "Krzak", 2026, "Single-cell RNA sequencing of terminal ileal biopsies"),
    "Zhao2026": ("10.1038/s41588-026-02755-z", "Zhao", 2026, "Single-cell analysis of chromatin accessibility"),
    "TGP2015": ("10.1038/nature15393", None, 2015, "A global reference for human genetic variation"),
    "Harrow2012": ("10.1101/gr.135350.111", "Harrow", 2012, "GENCODE"),
    "Breeze2025": ("10.64898/2025.12.19.695639", "Breeze", 2025, "LDscore"),
}
CTX = ssl.create_default_context(cafile=certifi.where())


def fetch(doi, tries=5):
    for i in range(tries):
        try:
            req = urllib.request.Request(f"https://api.crossref.org/works/{doi}",
                                         headers={"User-Agent": "crohns-reference-check (mailto:none)"})
            with urllib.request.urlopen(req, context=CTX, timeout=30) as r:
                return json.load(r)["message"]
        except Exception as e:  # rate limits / transient errors
            time.sleep(3 * (i + 1))
            last = e
    raise RuntimeError(f"{doi}: {last}")


def main():
    out = pathlib.Path("results/manuscript")
    out.mkdir(parents=True, exist_ok=True)
    rows = ["key\tdoi\tok\tfirst_author\tyear\ttitle\tjournal"]
    for key, (doi, author, year, frag) in REFS.items():
        m = fetch(doi)
        title = (m.get("title") or [""])[0]
        people = [a for a in (m.get("author") or []) if a.get("family")]  # skip consortium entries listed first
        fam = people[0]["family"] if people else ""
        yr = m["issued"]["date-parts"][0][0]
        ok = (frag.lower() in title.lower()) and (yr == year) and (author is None or author.lower() in fam.lower())
        rows.append(f"{key}\t{doi}\t{ok}\t{fam}\t{yr}\t{title}\t{(m.get('container-title') or [''])[0]}")
        print(("OK  " if ok else "FAIL") + f" {key}: {fam} {yr} {title[:70]}", flush=True)
        time.sleep(1)
    (out / "reference_check.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
