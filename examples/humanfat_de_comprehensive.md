---
# Example: HumanFat Comprehensive DE Statistics
# Status: full-build
# Validated: 2026-04-23
---

# Example: HumanFat Comprehensive DE Statistics

## Instantiates
- @modules/de_comprehensive_csv.md

## Project context
- Project: HumanFat
- Validated: 2026-04-23
- Status: full-build

## Brief block

```yaml
downstream_analyses:
  de_comprehensive_csv:
    group_col: tissue_type          # REQUIRED: tissue comparison column
                                    # values: Subcutaneous Fat, Visceral Fat, Orbital Fat,
                                    #         Liposuction Fat, Breast Fat, Lipoma
                                    # (Myelolipoma excluded — see note below)
    label_col: mylabel              # optional: EC subtype column for top_subtype annotation
    group_order:                    # optional: display ordering for tissue groups
      - "Subcutaneous Fat"
      - "Visceral Fat"
      - "Liposuction Fat"
      - "Breast Fat"
      - "Lipoma"
      - "Orbital Fat"
```

## Pre-processing note

```r
# Exclude Myelolipoma before running de_comprehensive_csv
# Rationale: Myelolipoma is a rare benign tumor; excluded to focus on
# adipose EC biology across the 6 canonical tissue groups
ec_subset <- ec[, ec$tissue_type != "Myelolipoma"]
```

## Validation notes

- Dataset: ~37,692 EC cells after Myelolipoma exclusion
- Group variable: tissue_type (6 tissue groups)
- Label variable: mylabel (6 EC subtypes: AEC, CapEC, CapEC2, VenEC1, VenEC2, VenEC3)
  Note: RibHighEC is present in mylabel but treated as a group in the Kruskal-Wallis test;
  however, RibHighEC genes are caught by `is_confound()` (ambient/quality markers) and
  filtered from the final output
- Result: ~253 significant genes after `is_confound()` filter (padj < 0.05, vectorized KW)
- `is_confound()` is loaded from @primitives/differential_expression.md — it extends
  `is_ambient()` with sex-linked gene exclusion, HLA class II, histones, and unannotated IDs
- Outputs written to `output/de_comprehensive_csv/comprehensive_DE_stats.csv`
- Cache file: `output/de_comprehensive_csv/kw_results_cache.Rds` (preserves intermediate
  Kruskal-Wallis results to avoid re-running the expensive KW step)
- top_subtype column identifies which EC subtype drives the significant variation per gene

## Known issues / quirks

- HumanFat uses `source_file` as batch variable — this does NOT affect de_comprehensive_csv
  directly (the module runs KW, not FindMarkers), but JoinLayers must be called first
- tissue_type has an additional_notes column for higher-resolution location; do NOT use
  additional_notes as GROUP_COL — it contains empty strings and inconsistent labeling
- RibHighEC will appear in the top_subtype column for some genes; these entries are valid
  (the KW test runs on all cells including RibHighEC) but should be interpreted with caution
- GROUP_ORDER omits Myelolipoma since it is excluded in pre-processing; if Myelolipoma is
  accidentally included, the KW test will run on 7 groups (not 6) — verify row count before running
