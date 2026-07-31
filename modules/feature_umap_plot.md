---
# Module: Feature UMAP Plot (Gene Expression Overlay)
# Pipeline: LargeDataset (compatible with any Seurat-based pipeline)
# Migrated from ~/claude-skills/pipelines/LargeDataset/methods/feature_umap_plot.md
# Phase 2 Group A — gene list and paths parameterized; primitives referenced
# Validated on: multiple single-cell RNA-seq datasets
requires_context:
  palettes: []
  metadata_columns:
    required: []
    optional:
      - label_col     # optional: column to use for subsetting before plotting
  brief_keys:
    required:
      - output_dir
    optional:
      - downstream_analyses.feature_umap_plot.genes    # list of gene symbols to plot
references:
  - "@primitives/visualization.md"
  - "@primitives/aesthetics.md"
  - "@primitives/seurat_v5_rules.md"
---

# Module: Feature UMAP Plot (Gene Expression Overlay)

Plots gene expression on the UMAP embedding, one panel per gene, assembled into a
single-row multi-panel PDF. Matches native Seurat FeaturePlot behavior exactly
(same color scale, same point size) while applying project aesthetic standards.

**Key references:**
- Color scale rules: @primitives/aesthetics.md (UMAP Feature Plots section)
- Typography standards: @primitives/aesthetics.md
- Dynamic UMAP reduction name: @primitives/seurat_v5_rules.md Rule 5

---

## Critical Gotchas — Do Not Repeat These Mistakes

These were all debugged against real project output vs native FeaturePlot side-by-side.
Each looked correct in isolation but produced visibly wrong output in comparison.

| Wrong | Correct | Why |
|---|---|---|
| `gene_max <- max(expr_vals, 0.001)` | `gene_max <- max(expr_vals)` with a separate `if (gene_max == 0)` guard | The `0.001` floor compresses all-zero genes to a near-zero pseudo-gradient instead of flat grey |
| `limits = c(0, gene_max)` with hardcoded 0 | `gene_min <- 0` (explicit, not `min(expr_vals)`) | For log-normalized data the minimum is always 0 and should be anchored there |
| `scale_color_gradient(low, high)` | `scale_color_gradientn(colors = c("lightgrey", "blue"))` | Seurat uses `gradientn` internally; gradient vs gradientn produce different intermediate stops |
| No `oob` argument | `oob = scales::squish` | Without squish, values marginally above `gene_max` due to floating point render as NA/transparent — indistinguishable from zero-expressors |
| `quantile(expr_vals[expr_vals > 0], 0.95)` as gene_max | `max(expr_vals)` | q95 of non-zero cells is too aggressive for sparse markers; scale compresses by ~2x |
| `geom_point(size = 0.8)` hardcoded | `Seurat:::AutoPointSize(data = data.frame(x = seq_len(ncol(so))))` | 0.8 is ~10x Seurat's computed size for ~19k cells; stacked large points produce over-dark saturated color |
| `AutoPointSize(so)` (Seurat object as input) | `AutoPointSize(data.frame(x = seq_len(ncol(so))))` | `nrow(so)` returns features (~38k), not cells (~19k); wrong input gives wrong pt.size |
| `df <- df[order(df$expr), ]` unconditionally | Wrap in `if (order_by_expression)` flag | Seurat default is `order = FALSE`; sorting changes visual appearance and should be opt-in |

---

## Configuration

```r
# ── Inputs ────────────────────────────────────────────────────────────────────
RDS_PATH  <- "project_specific"   # REPLACE: path to annotated Seurat RDS

# Genes to plot — from brief downstream_analyses block or specified directly
# Stage 8 dispatch: reads from downstream_analyses.feature_umap_plot.genes
GENES <- brief$downstream_analyses$feature_umap_plot$genes %||%
         brief$genes %||%    # fallback: top-level brief$genes (legacy compatibility)
         c("GENE1", "GENE2", "GENE3")  # REPLACE with actual gene symbols

# ── Behavioral options ────────────────────────────────────────────────────────
ORDER_BY_EXPRESSION <- FALSE  # FALSE matches Seurat default (order = FALSE)
                               # set TRUE to sort high-expressors on top

# ── Optional: subset before plotting ─────────────────────────────────────────
SUBSET_COL <- NULL   # e.g., "cell_type" — set NULL to plot all cells
SUBSET_VAL <- NULL   # e.g., c("Endothelial") — set NULL to plot all cells

# ── Output ────────────────────────────────────────────────────────────────────
OUTPUT_DIR <- file.path(
  brief$downstream_analyses$feature_umap_plot$output_dir %||%
  brief$output_dir %||%    # fallback: top-level brief$output_dir (legacy compatibility)
  "output",
  "feature_umap"
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
```

---

## Validated R Template

```r
library(Seurat)
library(ggplot2)
library(patchwork)

so <- readRDS(RDS_PATH)

# Optional subset
if (!is.null(SUBSET_COL) && !is.null(SUBSET_VAL)) {
  cells_keep <- colnames(so)[so[[SUBSET_COL]] %in% SUBSET_VAL]
  cat("Subsetting to", length(cells_keep), "of", ncol(so), "cells\n")
  so <- so[, cells_keep]
}

# Filter to genes present in the object
genes_present <- GENES[GENES %in% rownames(so)]
if (length(setdiff(GENES, genes_present)) > 0)
  cat("Genes not found:", paste(setdiff(GENES, genes_present), collapse = ", "), "\n")

# Dynamic UMAP reduction name — never hardcode "umap"
# See @primitives/seurat_v5_rules.md Rule 5
red_name <- grep("umap", names(so@reductions), ignore.case = TRUE, value = TRUE)[1]
if (is.na(red_name)) stop("No UMAP. Available: ", paste(names(so@reductions), collapse = ", "))
cat("Using reduction:", red_name, "\n")

coords <- as.data.frame(Embeddings(so, reduction = red_name))
colnames(coords) <- c("UMAP_1", "UMAP_2")

expr_mat <- GetAssayData(so, assay = "RNA", layer = "data")[genes_present, , drop = FALSE]

# Match Seurat AutoPointSize — must use ncol(so) (cells), NOT nrow(so) (features)
# See Gotcha 6 above
pt_size <- Seurat:::AutoPointSize(data = data.frame(x = seq_len(ncol(so))))
cat("pt_size:", pt_size, "\n")

plot_list <- lapply(genes_present, function(gene) {
  expr_vals <- as.numeric(expr_mat[gene, ])
  gene_min  <- 0
  gene_max  <- max(expr_vals)
  if (gene_max == 0) gene_max <- 0.001  # guard for all-zero genes (see Gotcha 1)

  df      <- coords
  df$expr <- expr_vals
  if (ORDER_BY_EXPRESSION) df <- df[order(df$expr), ]

  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = expr)) +
    geom_point(size = pt_size) +
    scale_color_gradientn(
      colors = c("lightgrey", "blue"),   # matches Seurat FeaturePlot default
      limits = c(gene_min, gene_max),
      oob    = scales::squish,           # required — see Gotcha 4
      name   = "Expression"
    ) +
    labs(title = gene, x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
    theme_classic() +
    theme(
      plot.title   = element_text(size = 16, face = "italic"),
      axis.title   = element_text(size = 15),
      axis.text    = element_text(size = 13, color = "grey20"),
      legend.text  = element_text(size = 12.5),
      legend.title = element_text(size = 12.5),
      axis.line    = element_line(color = "grey60", linewidth = 0.3),
      axis.ticks   = element_line(color = "grey60", linewidth = 0.3)
    )
})

p_final <- wrap_plots(plot_list, nrow = 1)

fig_w <- 5.5 * length(genes_present) + 1
fig_h <- 5.5

ggsave(file.path(OUTPUT_DIR, "feature_umap.pdf"),
       plot = p_final,
       width = fig_w, height = fig_h, units = "in",
       device = "pdf", useDingbats = FALSE)
```

---

## How to Verify Against Native Seurat

Run before finalizing to confirm ranges and pt.size match:

```r
fp <- FeaturePlot(so, genes_present[1], order = TRUE)

# Range should match
cat("Seurat range:", range(fp$data[[genes_present[1]]]), "\n")
cat("Our range:   ", range(GetAssayData(so, assay = "RNA", layer = "data")[genes_present[1], ]), "\n")

# pt.size should match
cat("Seurat pt.size:", fp$layers[[1]]$aes_params$size, "\n")
cat("Our pt_size:   ", Seurat:::AutoPointSize(data = data.frame(x = seq_len(ncol(so)))), "\n")
```

---

## Figure Sizing

| n_genes | Width (in) | Height (in) |
|---|---|---|
| 1 | 6.5 | 5.5 |
| 5 | 28.5 | 5.5 |
| Formula | `5.5 * n_genes + 1` | `5.5` |

Height is fixed at 5.5 in for all panel counts (single-row layout).
For multi-row layouts, use `nrow = n_rows` in `wrap_plots()` and multiply height accordingly.
