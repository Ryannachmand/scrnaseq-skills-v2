---
# Cohort Plots — worked example (generic)
#
# This example shows how to invoke the cohort_plots module on a typical
# multi-source harmonized dataset. Placeholders in angle brackets and values
# marked with "# REPLACE:" must be customized for your project.
#
# For the module reference itself, see: modules/cohort_plots.md
# Status: full-build
---

# Cohort Plots — worked example (generic)

## Instantiates
- @modules/cohort_plots.md

## Brief block

```yaml
cohort_plots: true

context_overrides:
  palettes:
    dataset_colors:
      # REPLACE: assign colors for your dataset/site labels; see context/color_palettes.md
      "site_1":  "#REPLACE_HEX"
      "site_2":  "#REPLACE_HEX"
      "site_3":  "#REPLACE_HEX"
    ct_colors:
      # REPLACE: cell type palette keyed by your unified_label values (stromal subtypes)
      # see context/color_palettes.md
      null
```

## Annotation strip configuration

Bone marrow stroma example uses four metadata columns in the annotation strip:

```r
# REPLACE: list the metadata columns relevant to your cohort annotation strip
ann_track_cols <- c("Dataset", "Organ", "Condition", "Data.Type")
# These map to standard IntegratePublicData metadata columns:
# Dataset   → which source dataset (e.g., site_1, site_2, site_3 — REPLACE)
# Organ     → anatomical site (REPLACE with your site labels)
# Condition → healthy / disease (or sort strategy for in vitro)
# Data.Type → fresh isolate vs cultured vs FACS-sorted
```

## Validated figure dimensions

Use the dynamic formula to size figures for your dataset:

```r
# Dynamic formula — REPLACE n_samples and n_celltypes with your dataset counts
fig_w <- max(n_samples * 0.3 + 4, 8)
fig_h <- max(n_celltypes * 0.35 + 2, 5)
# Validated example: n_samples ≈ 30, n_celltypes ≈ 14 → fig_w ≈ 13, fig_h ≈ 7
# Adjust multipliers if labels are wider or taller than typical
```

## ann_colors structure

The `ann_colors` argument to `pheatmap` is a named list matching `ann_col`:

```r
# REPLACE: structure for your annotation strip columns
ann_colors <- list(
  Dataset   = c(
    "site_1" = "#REPLACE_HEX",   # REPLACE: from context/color_palettes.md
    "site_2" = "#REPLACE_HEX",
    "site_3" = "#REPLACE_HEX"
  ),
  Organ     = c(
    "site_1" = "#REPLACE_HEX",   # REPLACE: anatomical site colors
    "site_2" = "#REPLACE_HEX",
    "site_3" = "#REPLACE_HEX"
  ),
  Condition = c(
    "GroupA" = "#REPLACE_HEX",   # REPLACE: your condition labels and colors
    "GroupB" = "#REPLACE_HEX"
  ),
  Data.Type = c(
    "TypeA" = "#REPLACE_HEX",    # REPLACE: your data type labels and colors
    "TypeB" = "#REPLACE_HEX"
  )
)
# Note: ann_colors keys must match the factor levels of each ann_col column exactly
# (case-sensitive, including spaces/periods)
```

## Validation notes

- Example validated on a harmonized UMAP from multiple source datasets (bone marrow stroma)
- Four annotation strips: Dataset, Organ, Condition, Data.Type
- Outputs written to `<YOUR_OUTPUT_DIR>/cohort_plots/`
- Three-part output: correlation heatmap (A), chord diagram (B), proportion heatmap (C)

## Known issues / quirks

- `pdf()` / `dev.off()` is required for correlation heatmap with styled title overlay
  (pheatmap + grid.text exception — CONVENTIONS.md §4 exception #2); this is correct behavior
- `circlize::chordDiagram()` requires `pdf()` / `dev.off()` (base R graphics; CONVENTIONS.md §4 exception #3)
- All other plots use `ggsave()` per CONVENTIONS.md rules
- JoinLayers should be called after merge (Stage 5), before integration (Stage 6); by the
  time cohort_plots runs, JoinLayers should already have been applied to the merged object
- `AverageExpression()` output requires `as.matrix()` conversion — see @primitives/seurat_v5_rules.md Rule 6
