---
# Trajectory Analysis (Monocle3) — worked example (generic)
#
# This example shows how to invoke the trajectory_monocle3 module on a typical
# dataset. Placeholders in angle brackets and values marked with
# "# REPLACE:" must be customized for your project.
#
# For the module reference itself, see: modules/trajectory_monocle3.md
# Status: full-build
---

# Trajectory Analysis (Monocle3) — worked example (generic)

## Instantiates
- @modules/trajectory_monocle3.md

## Brief block

```yaml
downstream_analyses:
  trajectory:
    enabled: true
    label_col: mylabel              # REPLACE: label column containing your cell type annotations
    exclude_pattern: "^AEC"        # REPLACE: regex pattern to exclude cell types from trajectory
                                    # choose a biologically plausible exclusion for your data
                                    # (e.g., exclude a terminal/endpoint cell type that would
                                    # collapse branching structure into a single terminus)
    start_cluster: "CapEC2"        # REPLACE: root cluster for pseudotime ordering
                                    # choose a biologically plausible root cell type for your data
                                    # (e.g., a progenitor, transitional, or immature phenotype)
    end_cluster: null               # Monocle3 infers endpoints from graph structure

context_overrides:
  palettes:
    subtype_colors:
      # REPLACE: assign colors for your cell types; see context/color_palettes.md
      "AEC":    "#REPLACE_HEX"
      "CapEC":  "#REPLACE_HEX"
      "CapEC2": "#REPLACE_HEX"
      "VenEC1": "#REPLACE_HEX"
      "VenEC2": "#REPLACE_HEX"
      "VenEC3": "#REPLACE_HEX"
```

## Validation notes

- Dataset: N cells with M subtypes (project-specific; update with your actual count)
- Label column: `mylabel` — verify all expected subtypes are populated before running
- After applying exclude_pattern and any low-quality cluster filters, confirm the remaining
  subtypes make biological sense as trajectory states
- START_CLUSTER should be the cluster that sits at the origin of the principal graph;
  pseudotime should increase toward differentiated or endpoint subtypes
- Outputs written to `<YOUR_OUTPUT_DIR>/trajectory/`
- Module correctly uses SeuratWrappers bypass (`new_cell_data_set()` manually) because
  SeuratWrappers is not available in r-env
- `preprocess_cds()` skipped (causes `as_cholmod_sparse` ABI error) — this is correct behavior

## Known issues / quirks

- Low-quality or high-ribosomal clusters should be excluded before trajectory analysis;
  they are not biologically meaningful as trajectory states
- If a UMAP is pre-computed in the Seurat object, do NOT rerun UMAP; use the existing reduction
- UMAP name detection via `grep("umap", names(so@reductions), ...)` handles non-standard
  reduction names — use dynamic detection rather than hardcoding "umap"
- Assay switching: switch to RNA + JoinLayers before `GetAssayData()` calls in trajectory;
  switch back to SCT (if used for visualization) after
- Always use the post-annotation label column (not any intermediate annotation columns) for
  the trajectory input
