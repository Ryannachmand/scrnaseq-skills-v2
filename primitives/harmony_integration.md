---
# Harmony Integration — Canonical Recipe — v2
---

# Harmony Integration — Canonical Recipe

The standard Seurat + Harmony processing sequence. Referenced by both LargeDataset (Stage 3)
and IntegratePublicData (Stage 6). Parameters come from the brief or lab defaults.

---

## Configuration

```r
# ── Integration parameters — set from brief ────────────────────────────────────
N_VARIABLE_FEATURES <- 4000    # lab default for whole objects; 2850 for subset re-clustering
N_PCS               <- 30      # lab default for whole objects; 40 for subset re-clustering
CLUSTERING_RES      <- 0.5     # lab default for whole objects; 0.39 for subsets
BATCH_COL           <- "sample_id"   # REQUIRED: metadata column passed to RunHarmony
                                     # lab default is sample_id — override in brief if different
ASSAY               <- "RNA"   # always RNA; never SCT for integration
```

---

## Processing Sequence

```r
# 1. Normalize and find variable features
so <- NormalizeData(so)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = N_VARIABLE_FEATURES)

# 2. Scale and PCA (variable features only — see seurat_v5_rules.md Rule 2)
so <- ScaleData(so, features = VariableFeatures(so))
gc()
so <- RunPCA(so, features = VariableFeatures(so), npcs = N_PCS)

# 3. Harmony batch correction
#    REQUIRED: BATCH_COL must exist in so@meta.data before this step
so <- harmony::RunHarmony(so, group.by.vars = BATCH_COL, assay.use = ASSAY)
gc()

# 4. Build neighbor graph on Harmony embeddings
so <- FindNeighbors(so, reduction = "harmony", dims = 1:N_PCS)
so <- FindClusters(so, resolution = CLUSTERING_RES)

# 5. UMAP on Harmony embeddings
so <- RunUMAP(so, reduction = "harmony", dims = 1:N_PCS)
```

---

## Critical Notes

- **Always JoinLayers first:** Before this sequence on a merged multi-sample object,
  call `so <- JoinLayers(so, assay = "RNA")` (see @primitives/seurat_v5_rules.md Rule 1)
- **Never SCT here:** All integration uses the RNA assay. SCT can be run for visualization
  after Harmony if needed, but Harmony and FindMarkers always use RNA.
- **BATCH_COL must be the actual column name in metadata.** Lab default is `sample_id`.
  Confirm this column exists before running: `table(so@meta.data[[BATCH_COL]])`
- **gc() after RunPCA and RunHarmony:** Both steps are memory-intensive on large objects.

---

## Parameter Defaults by Context

| Context | n_variable_features | n_pcs | clustering_resolution | batch_correction_var |
|---|---|---|---|---|
| Whole object (LargeDataset) | 4000 | 30 | 0.5 | sample_id |
| Subset re-clustering | 2850 | 40 | 0.39 | sample_id |
| IntegratePublicData exploratory | project_specific | 30 | 0.3 | sample_id |

All defaults come from lab_context.md. Override in the brief's `pipeline_params` block.

---

## Multi-Variable Harmony Correction

For two-factor correction (e.g., batch AND dataset), use:
```r
so <- harmony::RunHarmony(so, group.by.vars = c("dataset", "patient_id"), assay.use = ASSAY)
```
This is undocumented in the lab default but follows the Harmony API. Use when the brief
specifies multiple batch variables.
