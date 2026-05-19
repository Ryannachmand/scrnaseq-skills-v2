---
# Example: HumanFat Trajectory Analysis
# Status: full-build
# Validated: 2026-03-05
---

# Example: HumanFat Trajectory Analysis

## Instantiates
- @modules/trajectory_monocle3.md

## Project context
- Project: HumanFat
- Validated: 2026-03-05 (HumanFat_Yang run)
- Status: full-build

## Brief block

```yaml
downstream_analyses:
  trajectory:
    enabled: true
    label_col: mylabel            # EC subset label column
    exclude_pattern: "^AEC"      # exclude Arterial ECs from trajectory
                                  # rationale: AECs are the arterial endpoint, not the
                                  # transition state; including them collapses the capillary
                                  # branching structure into a single arterial terminus
    start_cluster: "CapEC2"      # root cluster for pseudotime ordering
                                  # CapEC2 = immature/transitional capillary EC phenotype;
                                  # validated as the biologically appropriate root
    end_cluster: null            # Monocle3 infers endpoints from graph structure

context_overrides:
  palettes:
    subtype_colors:
      "AEC":    "#F4433C"
      "CapEC":  "#FF9800"
      "CapEC2": "#FFEB3B"
      "VenEC1": "#4CAF50"
      "VenEC2": "#2196F3"
      "VenEC3": "#9C27B0"
```

## Validation notes

- Dataset: 53,487 EC cells, 13 subtypes (HumanFat_Yang run, 2026-03-05)
- Label column: `mylabel` contains AEC, CapEC, CapEC2, VenEC1, VenEC2, VenEC3, RibHighEC
- After `^AEC` exclusion and RibHighEC filter: trajectory runs on CapEC, CapEC2, VenEC1-3 subtypes
- START_CLUSTER = "CapEC2" validated as root — CapEC2 cells cluster at the origin of the
  principal graph; pseudotime increases toward VenEC endpoints
- Outputs written to `output/trajectory/`
- Module correctly uses SeuratWrappers bypass (`new_cell_data_set()` manually) because
  SeuratWrappers is not available in r-env
- `preprocess_cds()` skipped (causes `as_cholmod_sparse` ABI error) — this is correct behavior

## Known issues / quirks

- RibHighEC must be excluded before trajectory analysis (low-quality/high-ribosomal cluster;
  not biologically meaningful as a trajectory state)
- HumanFat UMAP is pre-computed — do NOT rerun UMAP; use the existing reduction
- UMAP name detection via `grep("umap", names(so@reductions), ...)` handles non-standard
  reduction names (FemoralHead pattern not relevant here, but the dynamic detection is robust)
- Assay switching: switch to RNA + JoinLayers before `GetAssayData()` calls in trajectory;
  switch back to SCT (if used for visualization) after
- The `ec_subtype` column (pre-annotation) vs `mylabel` (post-annotation): always use `mylabel`
  for the trajectory input in HumanFat
