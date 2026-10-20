# Provenance notes (2026-09-07)

- Source: public bucket gs://finngen-public-data-r13/summary_stats/ (FinnGen data freeze 13), accessed via the
  instructions sent by finngen-servicedesk@helsinki.fi on 2026-09-07 after registration.
- All files were fetched by code/01_download_finngen.py and each MD5 was checked against the MD5 that Google Cloud
  Storage reports (PROVENANCE.tsv, column md5_matches_gcs).
- finngen_R13_CHRONSMALL.gz was already on disk from an earlier, interrupted download (2026-09-07 16:38 local time);
  it was kept because its MD5 matches the bucket exactly. Its PROVENANCE.tsv timestamp is the time it was verified.
- A partial finngen_R13_CHRONLARGE.gz (505,262,080 of 797,702,339 bytes) from that interrupted download was discarded
  and a fresh copy downloaded and verified.
