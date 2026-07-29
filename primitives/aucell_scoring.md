---
# AUCell Scoring — Reusable Helpers — v2
---

# AUCell Scoring — Reusable Helpers

Rank-based gene set enrichment scoring using AUCell. These two functions handle
the core computation; gene set definitions are always caller-provided.

**Important distinction:**
- `run_aucell()` + `add_auc_to_seurat()` are general-purpose primitives (this file)
- Specific gene sets (metabolic pathways, PPARG signatures, etc.) belong in examples/

---

## Core Functions

```r
run_aucell <- function(seurat_obj, gene_sets, assay = "RNA") {
  # seurat_obj: Seurat object
  # gene_sets:  named list of character vectors — e.g. list(Glycolysis = c("HK1","HK2",...))
  #             Gene sets are CALLER-PROVIDED — never hardcode here.
  # assay:      assay to use for count matrix (always "RNA" for AUCell)
  #
  # Returns: data frame, cells × gene_set_names, with AUC scores

  counts_mat    <- GetAssayData(seurat_obj, assay = assay, layer = "counts")
  gene_sets_f   <- lapply(gene_sets, function(g) g[g %in% rownames(counts_mat)])
  cell_rankings <- AUCell_buildRankings(counts_mat, plotStats = FALSE, verbose = FALSE)
  auc_scores    <- AUCell_calcAUC(gene_sets_f, cell_rankings, verbose = FALSE)
  as.data.frame(t(getAUC(auc_scores)))
}

add_auc_to_seurat <- function(seurat_obj, auc_df) {
  # seurat_obj: Seurat object to add scores to
  # auc_df:     data frame returned by run_aucell() — cells in rows, gene sets in columns
  #
  # Adds AUC score columns to Seurat metadata via AddMetaData.
  # Handles cell barcode intersection gracefully (common cells only).

  common_cells <- intersect(rownames(auc_df), colnames(seurat_obj))
  for (col in colnames(auc_df)) {
    vals <- setNames(auc_df[[col]], rownames(auc_df))
    seurat_obj <- AddMetaData(seurat_obj, metadata = vals[common_cells], col.name = col)
  }
  seurat_obj
}
```

---

## Memory Management for Large Objects

For whole objects > 100k cells, subsample before AUCell to avoid out-of-memory:

```r
# Stratified 2500/group gives near-identical group means and reduces memory ~5x
keep_cells <- so@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(!!sym(label_col)) %>%
  slice_sample(n = 2500) %>%
  pull(cell)

obj_sub <- subset(so, cells = keep_cells)
auc_df  <- run_aucell(obj_sub, gene_sets)
```

Note: AUCell `buildRankings` holds the full count matrix rank-transformed in memory.
For >200k cells, even with subsampling, memory may be tight. Monitor with `gc()` calls.

---

## Usage Pattern

```r
# 1. Define gene sets (caller-provided — project-specific, not hardcoded here)
gene_sets <- list(
  PathwayA = c("GENE1", "GENE2", "GENE3"),   # project_specific
  PathwayB = c("GENE4", "GENE5", "GENE6")    # project_specific
)

# 2. Score cells
auc_df <- run_aucell(so, gene_sets, assay = "RNA")

# 3. Add scores to Seurat object
so <- add_auc_to_seurat(so, auc_df)
# Scores now available as so$PathwayA, so$PathwayB in metadata

# 4. Visualize
FeaturePlot(so, features = names(gene_sets))
VlnPlot(so, features = names(gene_sets), group.by = label_col)
```

---

## Gene Set Sources

Gene sets for specific biological contexts are maintained in examples/.
For common pathway gene sets, consider:
- MSigDB Hallmark or KEGG collections via `msigdbr` package
- Lab-curated sets from analysis briefs (see context/validated_examples.yaml per-project entries)

Do NOT hardcode pathway gene sets in this primitive.
