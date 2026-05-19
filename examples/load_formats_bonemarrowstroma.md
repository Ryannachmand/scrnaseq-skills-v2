---
# Example: BoneMarrowStroma Load Formats
# Status: full-build
# Validated: BoneMarrowStroma run
---

# Example: BoneMarrowStroma Load Formats

## Instantiates
- @modules/load_formats.md

## Project context
- Project: BoneMarrowStroma
- Validated: BoneMarrowStroma run
- Status: full-build

## Brief block

```yaml
inputs:
  datasets:
    - name: "Vertebrae"
      type: mex                     # 10x Genomics MEX format (Leimkuler dataset)
      path: "data/Vertebrae/"       # TODO: confirm actual path
      notes: "Leimkuler et al. (Vertebrae); accession pattern: _raw_gene_bc_matrices_h5.h5"

    - name: "IliacCrest"
      type: mex                     # Wang dataset — nested directory structure
      path: "data/IliacCrest/"      # TODO: confirm actual path
      notes: |
        Wang dataset (Iliac Crest); nested directory layout:
        GSM4423510_MNC_op/2.2.filtered_feature_bc_matrix/
        Each GEO accession has its own subdirectory containing the MEX files.
        Use recursive directory traversal to locate barcodes/features/matrix files.

    - name: "FemoralHead"
      type: rds                     # Li dataset — pre-processed RDS
      path: "data/FemoralHead/"     # TODO: confirm actual path
      notes: |
        Li et al. (Femoral Head); spliced/unspliced format from velocyto.
        Accession pattern: paste0(accession, "_barcodes.tsv.gz") for barcodes,
        paste0(accession, "_genes.tsv.gz") for features.
        Merged from multiple GEO accessions (GSE190965).
        IMPORTANT: stores UMAP as "UMAP_dim30" — non-standard reduction name.
```

## GSM file pattern conventions

```r
# Leimkuler (Vertebrae) — h5 format with project-specific file pattern
leimkuler_pattern <- "_raw_gene_bc_matrices_h5\\.h5$"
# Usage: list.files(data_dir, pattern = leimkuler_pattern, recursive = TRUE)

# Wang (Iliac Crest) — nested MEX directory
# Nested directory: GSM4423510_MNC_op/2.2.filtered_feature_bc_matrix/
wang_example <- "GSM4423510_MNC_op"   # GSM accession prefix
# Read with: Read10X(data.dir = file.path(data_dir, wang_example, "2.2.filtered_feature_bc_matrix"))

# Li (Femoral Head) — spliced/unspliced MEX
# File pattern: {accession}_barcodes.tsv.gz, {accession}_genes.tsv.gz, {accession}_matrix.mtx.gz
# GSE190965: multiple accessions merged into a single Seurat object
li_accessions <- c(
  # TODO: fill in actual GSM accessions from GSE190965
  # Pattern: paste0(accession, "_barcodes.tsv.gz")
)
```

## FemoralHead UMAP reduction name

The FemoralHead (Li) RDS object uses a non-standard UMAP reduction name:

```r
# FemoralHead stores UMAP as "UMAP_dim30" — always use dynamic detection:
umap_key <- names(so@reductions)[grepl("umap", names(so@reductions), ignore.case = TRUE)][1]
# For FemoralHead: umap_key will be "UMAP_dim30"
# For all other datasets: umap_key will be "umap" (standard)
# See @primitives/seurat_v5_rules.md Rule 5
```

## Validation notes

- Three source datasets representing three anatomical sites: Vertebrae, Iliac Crest, Femoral Head
- Each dataset loaded separately before label harmonization and integration
- JoinLayers called after merge (Stage 5), before integration (Stage 6) — per Rule 1
- FemoralHead pre-processed RDS contains the non-standard UMAP reduction; use dynamic
  UMAP detection (Rule 5) rather than hardcoding "umap"
- The `isolatedstroma` object reference in v1 was BoneMarrowStroma-specific;
  in v2 this is represented as a subset operation using target_celltypes from the brief

## Known issues / quirks

- Wang dataset (IliacCrest) has deeply nested directory structure; the loader must use
  recursive = TRUE or explicit subdirectory paths to find MEX files
- Li dataset (FemoralHead) has spliced/unspliced assays from velocyto; load the spliced
  assay (or the RNA assay if pre-processed) as the count matrix
- CLEC4G + LYVE1 > 100 FPKM threshold confirms liver sinusoidal ECs — this is an example
  quality control check pattern; substitute appropriate markers for bone marrow context
- Ensembl ID to gene symbol conversion (if needed): use `RenameFeatures()` per
  @modules/load_formats.md P2-9 fix (NOT direct slot assignment to `rownames(counts)`)
