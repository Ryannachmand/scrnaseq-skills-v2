---
# Seurat v5 Rules — v2
# Migrated from ~/claude-skills/shared/seurat_v5_rules.md
# Permanent corrections for Seurat v5 behavior. Apply in every generated R script.
# Last validated: TFEC Expression Atlas run, 2026-04-10
---

# Seurat v5 Rules

Permanent bug fixes for Seurat v5 behavior. Apply in every generated R script
that uses Seurat — no exceptions.

---

## 1. JoinLayers Before ScaleData

After `merge()` on multi-sample objects, always call:
```r
so <- JoinLayers(so, assay = "RNA")
```
before `FindVariableFeatures()` or `ScaleData()`.

**Why:** Seurat v5 stores each sample as a separate layer after `merge()`.
Without joining, `ScaleData()` silently fails or hits memory issues.

---

## 2. Scale Only Variable Features

Always:
```r
so <- ScaleData(so, features = VariableFeatures(so))
```
Never `features = rownames(so)`.

**Why:** Scaling all genes (~40k+) causes OOM kill (exit code 137) on
large datasets. Only variable features are needed for PCA.

---

## 3. Metadata Assignment from Named Vectors

When assigning cluster-derived labels (e.g., celltype from a named vector
keyed by cluster ID), always use direct slot assignment:
```r
so@meta.data$celltype <- named_vector[as.character(so@meta.data$seurat_clusters)]
```
Never: `so$celltype <- named_vector`

**Why:** In Seurat v5, `$<-` with a named vector attempts cell-name matching
on the vector names (cluster IDs), producing "No cell overlap" error.

---

## 4. External Metadata CSV Merge

When a Seurat RDS has outdated labels and a newer annotation CSV exists separately,
merge the CSV metadata into the object before any downstream work.

```r
meta <- read.csv(meta_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
shared <- intersect(colnames(so), rownames(meta))
cat("Shared cells:", length(shared), "/", ncol(so), "\n")

# Add/overwrite columns from CSV using direct slot assignment (safe in v5)
for (col in colnames(meta)) {
  so@meta.data[[col]] <- NA
  so@meta.data[shared, col] <- meta[shared, col]
}
```

**Why:** `so$col <- meta$col` in Seurat v5 attempts cell-name matching and will silently
produce NA for all cells if the vector is not named by barcode. Direct `@meta.data` assignment
is always safe. Always report the shared cell count — a low overlap (< 95%) warrants investigation.

---

## 5. Dynamic UMAP Reduction Name Detection

Never hardcode `reduction = "umap"`. Objects from other labs frequently use non-standard names
(`UMAP_dim30`, `wnn.umap`, `ref.umap`, etc.).

```r
red_name <- grep("umap", names(so@reductions), ignore.case = TRUE, value = TRUE)[1]
if (is.na(red_name)) stop("No UMAP found. Available: ",
                           paste(names(so@reductions), collapse = ", "))
cat("Using reduction:", red_name, "\n")
coords <- as.data.frame(Embeddings(so, reduction = red_name))
```

---

## 6. Differential Abundance — Matrix Conversion

When building a log2 odds ratio matrix for heatmaps, always use base R:
```r
lor_mat <- as.matrix(lor_wide[, -1])
rownames(lor_mat) <- lor_wide$subtype
```
Never `tibble::column_to_rownames()`.

**Why:** `tibble` package not reliably available in r-env; base R is more robust.
