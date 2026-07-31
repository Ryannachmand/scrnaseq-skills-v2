---
# Comprehensive DE Statistics — worked example (generic)
#
# This example shows how to invoke the de_comprehensive_csv module on a typical
# dataset. Placeholders in angle brackets and values marked with
# "# REPLACE:" must be customized for your project.
#
# For the module reference itself, see: modules/de_comprehensive_csv.md
# Status: full-build
---

# Comprehensive DE Statistics — worked example (generic)

## Instantiates
- @modules/de_comprehensive_csv.md

## Brief block

```yaml
downstream_analyses:
  de_comprehensive_csv:
    group_col: "project_specific"   # REQUIRED: comparison column (REPLACE with your group column)
                                    # REPLACE: list your actual group values in a comment here
    label_col: "project_specific"   # optional: cell subtype column for top_subtype annotation
    group_order:                    # optional: display ordering for groups (REPLACE with your levels)
      - "GroupA"
      - "GroupB"
      - "GroupC"
```

## Pre-processing note

```r
# REPLACE: exclude any groups that should not be part of the analysis
# Example: exclude a rare or outlier group before running de_comprehensive_csv
so_subset <- so[, so$project_specific != "ExcludedGroup"]
# Rationale: document why this group is excluded for your project
```

## Validation notes

- Dataset: N cells (project-specific; update with your actual count after subsetting)
- Group variable: tissue_type (REPLACE: N groups — list them)
- Label variable: mylabel (REPLACE: N cell subtypes — list them)
  Note: any low-quality or ambient clusters in the label column will appear in top_subtype
  output but are caught by `is_confound()` filtering
- Result: N significant genes after `is_confound()` filter (padj < 0.05, vectorized KW)
- `is_confound()` is loaded from @primitives/differential_expression.md — it extends
  `is_ambient()` with sex-linked gene exclusion, HLA class II, histones, and unannotated IDs
- Outputs written to `<YOUR_OUTPUT_DIR>/de_comprehensive_csv/comprehensive_DE_stats.csv`
- Cache file: `<YOUR_OUTPUT_DIR>/de_comprehensive_csv/kw_results_cache.Rds` (preserves
  intermediate Kruskal-Wallis results to avoid re-running the expensive KW step)
- top_subtype column identifies which cell subtype drives the significant variation per gene

## Known issues / quirks

- If using a batch variable in your project, this does NOT affect de_comprehensive_csv
  directly (the module runs KW, not FindMarkers), but JoinLayers must be called first
- Verify the group_col does not contain empty strings or inconsistent labels before running;
  additional metadata columns with inconsistent formatting should not be used as GROUP_COL
- Any low-quality clusters will appear in the top_subtype column for some genes; these
  entries are valid (the KW test runs on all cells) but should be interpreted with caution
- GROUP_ORDER should match exactly the groups present after subsetting; verify row count
  matches expected N before running
