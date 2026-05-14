---
# GEO Data Download Rules — v2
# Migrated from ~/claude-skills/shared/geo_download.md
# No project-specific content in v1; carried forward as-is.
---

# GEO Data Download Rules

---

## Priority Order for Downloading Public Datasets

Always attempt in this order — stop at the first success:

1. **Processed count matrix** — look for `*_matrix.h5`, `*_counts.h5`, `matrix.mtx.gz` (MEX format),
   or any processed count file posted directly to the GEO accession page. Download these.
2. **Per-sample processed files** — if no merged matrix exists, look for individual sample matrix
   files and download all of them.
3. **Supplementary files** — check GSE supplementary file listings on GEO for any count-level data
   before assuming realignment is needed.
4. **FASTQ / raw reads** — NEVER download FASTQs without explicit user approval. If only FASTQs
   are available, stop and report this to the user before proceeding.

---

## Efficiency Rules

- Before downloading anything, fetch the GEO accession page and summarize the available files
  to the user (file names, sizes, formats). Wait for approval before downloading.
- Prefer merged/combined matrices over per-sample files when both exist.
- Only collect **single-cell RNA sequencing** data. Skip bulk RNA-seq, ATAC-seq, or other
  modalities even if posted under the same GSE.
- If total download size exceeds 20 GB, flag this to the user before proceeding.
- Never download FASTQs. If realignment is the only option, ask the user explicitly.
