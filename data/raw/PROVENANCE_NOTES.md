# Provenance notes, data/raw (2026-09-07)

- GCST90446793.tsv-meta.yaml: the MD5 listed in the GWAS Catalog's md5sum.txt (8dfa8092ce55ea22b6729e531fabebb6) does
  not match the file served; three separate downloads (http and https) all gave 0e5f4e24f3ab20bd8aca05d7ccd6d010. The
  yaml says it was last modified 2025-01-05, so the listed checksum appears stale. The summary-statistics file itself
  (GCST90446793.tsv) matches its listed MD5. The yaml is recorded with no published checksum.
- GCST90446793 is the East Asian fixed-effect meta-analysis of four cohorts. The cohort order in the Direction column
  is SHA1, ICH1, KOR1, JPN1, per readme.sst.txt (version 12202022) in the IIBDGC Google Drive folder "IBD_EAS/sumstats"
  linked from https://www.ibdgenetics.org/ (viewed 2026-09-07). ICH1 was genotyped on ImmunoChip.
- IIBDGC_IBD_EAS_sumstats_readme.sst.txt: archived 2026-09-09 from the consortium's Google Drive (Drive id 1taTigC3PM9TwK-RBx8ZSprci59yIW9qO); Drive publishes no checksum, so MD5 and SHA-256 were computed locally (PROVENANCE.tsv).
