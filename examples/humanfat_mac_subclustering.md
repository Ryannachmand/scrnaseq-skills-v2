---
# Example: HumanFat Macrophage/Monocyte Subclustering
# Status: full-build
# Validated: 2026-03-24
---

# Example: HumanFat Macrophage/Monocyte Subclustering

## Instantiates
- @modules/celltype_subclustering.md

## Project context
- Project: HumanFat
- Validated: 2026-03-24
- Status: full-build

## Brief block

```yaml
downstream_analyses:
  subclustering:
    target_celltypes:
      - "Macrophage"
      - "Monocyte"              # also matches "Naive Monocyte", "Polarized Monocyte",
                                # "Prolif. Monocyte" — use partial match or list all
    subset_name: "mac_mono"
    pipeline_params:
      batch_correction_var: source_file   # HumanFat-specific: use source_file, NOT sample_id
      n_variable_features: 2850           # lab default for subsets
      n_pcs: 40                           # lab default for subsets
      clustering_resolution: 0.39         # lab default for subsets

context_overrides:
  palettes:
    group_colors:
      "Subcutaneous Fat": "#F4A460"
      "Visceral Fat":     "#CD853F"
      "Lipoma":           "#DEB887"
      "Liposuction Fat":  "#D2691E"
      "Myelolipoma":      "#8B4513"
      "Breast Fat":       "#FFB6C1"
      "Orbital Fat":      "#20B2AA"
```

## Cluster-to-label mapping (validated 2026-03-24)

After running the subclustering pipeline, the user provides the cluster→label mapping.
The validated HumanFat macrophage/monocyte mapping from the 2026-03-24 run:

```r
# Cluster → label mapping (set after reviewing FindAllMarkers output)
# These values are project-specific and run-specific; document here for reproducibility
cluster_labels <- c(
  # REPLACE with actual validated cluster→label mapping from the 2026-03-24 run
  # Pattern from that run:
  # cluster 0 → "Resident Macrophage"
  # cluster 1 → "Infiltrating Macrophage"
  # cluster 2 → "Naive Monocyte"
  # cluster 3 → "Polarized Macrophage"
  # cluster 4 → "Prolif. Macrophage"
  # cluster 5 → "Dendritic Cell"  # cross-contamination; verify markers
  # TODO: retrieve exact mapping from 2026-03-24 run CLAUDE.md or project notes
)
```

## Validation notes

- Validated on HumanFat whole-object (Macrophage + Monocyte subset, two-phase subclustering)
- Two-phase: (1) subset and re-cluster, (2) user provides cluster→label mapping
- Batch correction variable: `source_file` (HumanFat-specific; not the lab default `sample_id`)
- SAMPLE_COL defaults to BATCH_CORRECTION_VAR = source_file (both are the same column here)
- GROUP_COL = tissue_type (adipose tissue type column)
- Outputs written to `output/subclustering/mac_mono/`

## Known issues / quirks

- HumanFat uses `source_file` as both the batch variable and the sample column — this
  is non-standard; the module config must set `BATCH_CORRECTION_VAR <- "source_file"` and
  `SAMPLE_COL <- "source_file"` explicitly
- Macrophage and Monocyte labels in `celltype` (whole-object column) include:
  "Macrophage", "Naive Monocyte", "Polarized Monocyte", "Prolif. Monocyte"
  Subsetting by CELLTYPE_LABELS should include all four labels
- The subclustering module assigns a new `subtype_label` column (not `mylabel`) to the
  re-clustered object; if this sub-object is merged back into the main object, reconcile
  the column names
- HumanFat UMAP is pre-computed for the EC subset; subclustering runs its own UMAP on the
  mac/mono subset — these are independent reductions, not the same UMAP space
