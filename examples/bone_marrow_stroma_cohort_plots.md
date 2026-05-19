---
# Example: BoneMarrowStroma Cohort Plots
# Status: full-build
# Validated: BoneMarrowStroma run (date TODO — not confirmed in registry)
---

# Example: BoneMarrowStroma Cohort Plots

## Instantiates
- @modules/cohort_plots.md

## Project context
- Project: BoneMarrowStroma
- Validated: BoneMarrowStroma run (see validated_examples.yaml)
- Status: full-build

## Brief block

```yaml
cohort_plots: true

context_overrides:
  palettes:
    dataset_colors:
      "Vertebrae":    null    # TODO: confirm hex values from BoneMarrowStroma project records
      "IliacCrest":   null    # TODO: confirm hex values
      "FemoralHead":  null    # TODO: confirm hex values
    ct_colors:
      null                    # TODO: cell type palette for BoneMarrowStroma labels
                              # keyed by unified_label values (stromal subtypes)
```

## Annotation strip configuration

BoneMarrowStroma uses four metadata columns in the annotation strip:

```r
ann_track_cols <- c("Dataset", "Organ", "Condition", "Data.Type")
# These map to standard IntegratePublicData metadata columns:
# Dataset   → which source dataset (Leimkuler/Vertebrae, Wang/IliacCrest, Li/FemoralHead)
# Organ     → anatomical site
# Condition → healthy / disease (or sort strategy for in vitro)
# Data.Type → fresh isolate vs cultured vs FACS-sorted
```

## Validated figure dimensions

Validated figure output dimensions from BoneMarrowStroma run:

```r
# Correlation heatmap (Part A)
fig_w <- 13.25   # validated width in inches
fig_h <- 7.00    # validated height in inches
# These values correspond to the BoneMarrowStroma dataset size:
# ~3 datasets × multiple samples; 13.25 × 7 was the appropriate figure footprint
# for the sample-level correlation matrix annotation

# Dynamic formula for other datasets:
# fig_w <- max(n_samples * 0.3 + 4, 8)
# fig_h <- max(n_celltypes * 0.35 + 2, 5)
# BoneMarrowStroma case validates: n_samples ≈ 31, n_celltypes ≈ 14
```

## ann_colors structure

The `ann_colors` argument to `pheatmap` is a named list matching `ann_col`:

```r
# Structure for BoneMarrowStroma annotation strip (4 columns)
ann_colors <- list(
  Dataset   = c(
    "Vertebrae"   = "#project_specific",   # TODO: confirm from project records
    "IliacCrest"  = "#project_specific",   # TODO: confirm from project records
    "FemoralHead" = "#project_specific"    # TODO: confirm from project records
  ),
  Organ     = c(
    "Vertebrae"   = "#project_specific",
    "Iliac Crest" = "#project_specific",
    "Femoral Head"= "#project_specific"
  ),
  Condition = c(
    "project_specific" = "#project_specific"   # TODO: fill in condition labels
  ),
  Data.Type = c(
    "project_specific" = "#project_specific"   # TODO: fill in data type labels
  )
)
# Note: ann_colors keys must match the factor levels of each ann_col column exactly
# (case-sensitive, including spaces/periods)
```

## Validation notes

- Validated on BoneMarrowStroma harmonized UMAP (multi-source: Leimkuler, Wang, Li)
- Four annotation strips: Dataset, Organ, Condition, Data.Type
- Validated figure size 13.25 × 7 in (see above)
- Outputs written to `output/cohort_plots/`
- Three-part output: correlation heatmap (A), chord diagram (B), proportion heatmap (C)

## Known issues / quirks

- `pdf()` / `dev.off()` is required for correlation heatmap with styled title overlay
  (pheatmap + grid.text exception — CONVENTIONS.md §4 exception #2); this is correct behavior
- `circlize::chordDiagram()` requires `pdf()` / `dev.off()` (base R graphics; CONVENTIONS.md §4 exception #3)
- All other plots use `ggsave()` per CONVENTIONS.md rules
- BoneMarrowStroma JoinLayers timing: called in Stage 5 (after merge), before Stage 6 integration;
  by the time cohort_plots runs, JoinLayers has already been applied to the merged object
- `AverageExpression()` output requires `as.matrix()` conversion — see @primitives/seurat_v5_rules.md Rule 6
