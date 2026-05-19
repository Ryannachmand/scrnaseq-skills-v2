---
# Example: CoCulture Feature UMAP Plot (Skeleton)
# Status: skeleton
# Validated: TODO
---

# Example: CoCulture Feature UMAP Plot

## Instantiates
- @modules/feature_umap_plot.md

## Project context
- Project: CoCulture
- Validated: TODO — date not confirmed in registry
- Status: skeleton

## Brief block

```yaml
downstream_analyses:
  feature_umap_plot:
    enabled: true
    genes:
      # TODO: fill in specific gene list from CoCulture run
      # CoCulture is the validation project for feature_umap_plot.md;
      # the specific genes used in the validated run are not documented in library records
      - "GENE1"    # TODO: replace with actual CoCulture gene list
      - "GENE2"
      - "GENE3"
    output_dir: "output"    # TODO: confirm output directory used in validated run

# TODO: fill in CoCulture-specific context_defaults once known:
# - batch_correction_var: null (TODO: confirm)
# - label_col: null (TODO: confirm cell type label column)
# - input RDS path (not documented in registry)
```

## Validation notes

CoCulture is the validation project for `@modules/feature_umap_plot.md`. The module was
validated on 19,348 cells in the CoCulture run (R 4.5 / Seurat v5, log-normalized RNA layer).

The seven critical gotchas documented in the module (AutoPointSize formula, oob=squish,
gradientn vs gradient, q95 limit compression, etc.) were all debugged against CoCulture
output vs native Seurat FeaturePlot side-by-side.

Specific details still needed:
- **TODO**: Gene list used in the validated CoCulture feature UMAP run
- **TODO**: CoCulture dataset path / input RDS
- **TODO**: Validated output dimensions and file names
- **TODO**: Validation date for the CoCulture run
- **TODO**: Cell type label column name for CoCulture

## Known issues / quirks

- CoCulture context is sparse in the project registry (validated_examples.yaml);
  the validated_examples.yaml entry has `context_defaults: null` for most fields
- Biology context for CoCulture is not documented in library records; the registry
  notes only that CoCulture validates feature_umap_plot.md
- The `SUBSET_COL` / `SUBSET_VAL` optional subsetting feature was added to the module
  in Phase 2A; confirm whether the CoCulture run used subsetting or plotted all cells

## TODO for future wrap-up agent

1. Retrieve CoCulture gene list from project CLAUDE.md or run records
2. Retrieve input RDS path
3. Confirm validation date and output dimensions
4. Add CoCulture to validated_examples.yaml with complete context_defaults
5. Mark this example as full-build once the above are resolved
