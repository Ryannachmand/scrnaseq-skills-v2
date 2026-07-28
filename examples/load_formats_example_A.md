---
# Load Formats — worked example (generic)
#
# This example shows how to invoke the load_formats module on a typical
# multi-source dataset with mixed input formats. Placeholders in angle brackets
# and values marked with "# REPLACE:" must be customized for your project.
#
# For the module reference itself, see: modules/load_formats.md
# Status: full-build
---

# Load Formats — worked example (generic)

## Instantiates
- @modules/load_formats.md

## Brief block

```yaml
inputs:
  datasets:
    - name: "dataset_1"
      type: mex                     # 10x Genomics MEX format
      path: "<YOUR_DATA_DIR>/dataset_1/"
      notes: "REPLACE: source publication and accession pattern for dataset_1"

    - name: "dataset_2"
      type: mex                     # MEX format with nested directory structure
      path: "<YOUR_DATA_DIR>/dataset_2/"
      notes: |
        REPLACE: dataset_2 source and any directory layout notes.
        If data has a nested structure, document the subdirectory pattern here.
        Use recursive directory traversal to locate barcodes/features/matrix files.

    - name: "dataset_3"
      type: rds                     # pre-processed RDS
      path: "<YOUR_DATA_DIR>/dataset_3/"
      notes: |
        REPLACE: dataset_3 source and any format notes.
        If from velocyto: spliced/unspliced format — load the spliced (RNA) assay.
        If UMAP reduction has a non-standard name, see note below on dynamic detection.
```

## File pattern conventions

```r
# REPLACE: update these patterns to match your actual file naming conventions

# dataset_1 — standard h5 or MEX pattern
dataset1_pattern <- "_filtered_feature_bc_matrix\\.h5$"   # REPLACE with actual pattern
# Usage: list.files(data_dir, pattern = dataset1_pattern, recursive = TRUE)

# dataset_2 — nested MEX directory example
# REPLACE: substitute your actual accession prefix and subdirectory structure
dataset2_accession <- "GSMXXXXXXX_sample_name"   # REPLACE
# Read with: Read10X(data.dir = file.path(data_dir, dataset2_accession, "filtered_feature_bc_matrix"))

# dataset_3 — multiple accessions merged into one object (if applicable)
dataset3_accessions <- c(
  # REPLACE: list actual GEO accessions
  # Pattern: paste0(accession, "_barcodes.tsv.gz")
)
```

## Non-standard UMAP reduction name

Some pre-processed RDS objects use a non-standard UMAP reduction name (e.g., "UMAP_dim30"
instead of "umap"). Always use dynamic detection:

```r
# Use dynamic UMAP detection for any dataset that may have a non-standard reduction name
umap_key <- names(so@reductions)[grepl("umap", names(so@reductions), ignore.case = TRUE)][1]
# Standard datasets: umap_key will be "umap"
# Non-standard example: umap_key might be "UMAP_dim30" or "wnn.umap"
# See @primitives/seurat_v5_rules.md Rule 5
```

## Validation notes

- Example: three source datasets representing three anatomical sites, each in a different
  input format (MEX standard, MEX nested, pre-processed RDS)
- Each dataset loaded separately before label harmonization and integration
- JoinLayers called after merge (Stage 5), before integration (Stage 6) — per Rule 1
- If any source dataset contains a non-standard UMAP reduction, use dynamic UMAP detection
  (Rule 5) rather than hardcoding "umap"

## Known issues / quirks

- Nested MEX directory structures require `recursive = TRUE` or explicit subdirectory paths;
  inspect the directory tree before writing the loader
- Velocyto-processed objects contain spliced/unspliced assays; load the spliced assay
  (or the RNA assay if pre-processed) as the count matrix
- REPLACE: add a marker gene QC check appropriate for your tissue (e.g., confirm cell type
  identity with canonical markers for your biological context)
- Ensembl ID to gene symbol conversion (if needed): use `RenameFeatures()` per
  @modules/load_formats.md P2-9 fix (NOT direct slot assignment to `rownames(counts)`)
