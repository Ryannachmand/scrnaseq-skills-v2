---
# Module: Cohort Summary Plots
# Pipeline: IntegratePublicData
# Migrated from ~/claude-skills/pipelines/IntegratePublicData/methods/cohort_plots.md
# Phase 2 Group A — column names and colors parameterized; hardcoded size removed
requires_context:
  palettes:
    - dataset_colors    # required — named vector: dataset value → hex color
    - ct_colors         # required — named vector: cell type label → hex color
  metadata_columns:
    required:
      - dataset_col     # column identifying dataset (default: "dataset")
      - label_col       # column holding harmonized cell type labels (default: "unified_label")
      - sample_col      # column identifying samples (default: "sample_id")
    optional:
      - condition_col   # column holding condition (default: "condition")
      - organ_col       # column holding organ (default: "organ")
      - data_type_col   # column holding data type (default: "data_type")
  brief_keys:
    required:
      - project_name
      - output_dir
    optional:
      - n_hvg           # number of HVGs for correlation heatmap (default: 2250)
references:
  - "@primitives/visualization.md"
  - "@primitives/aesthetics.md"
  - "@context/color_palettes.md"
  - "@primitives/seurat_v5_rules.md"
---

# Module: Cohort Summary Plots

Three cohort-level plot types for the IntegratePublicData pipeline:
A. Correlation heatmap (transcriptional similarity across cell types)
B. Chord diagram (dataset × cell type composition)
C. Proportion heatmap (cell type percentages per sample with metadata strips)

Set these variables at the top of the script before calling any section:

```r
# ── Configuration ────────────────────────────────────────────────────────────
DATASET_COL    <- "dataset"       # brief: metadata.dataset_col
LABEL_COL      <- "unified_label" # brief: metadata.label_col
SAMPLE_COL     <- "sample_id"     # brief: metadata.sample_col
CONDITION_COL  <- "condition"     # brief: metadata.condition_col
ORGAN_COL      <- "organ"         # brief: metadata.organ_col
DATA_TYPE_COL  <- "data_type"     # brief: metadata.data_type_col
N_HVG          <- 2250            # brief: n_hvg (default 2250)
OUTPUT_DIR     <- "output/cohort_plots"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Palettes — must be named vectors keyed by actual data values
dataset_colors <- project_specific  # REPLACE: named vector of dataset → hex color
ct_colors      <- project_specific  # REPLACE: named vector of cell type → hex color
```

---

## A. Correlation Heatmap

Cell type × cell type Pearson correlation using highly variable genes.
Uses z-scored within-dataset averages before cross-dataset averaging.
This corrects for batch/depth effects without limma or embedding-space operations.

```r
library(pheatmap)
library(Matrix)

# Get HVGs
hvg <- head(VariableFeatures(so), N_HVG)

# Compute scaled average per cell type per dataset
datasets <- sort(unique(so[[DATASET_COL]]))
all_ct   <- sort(unique(so[[LABEL_COL]]))

avg_list <- list()
for (ds in datasets) {
  so_ds <- subset(so, cells = colnames(so)[so[[DATASET_COL]] == ds])
  # IMPORTANT: AverageExpression returns dgCMatrix in Seurat v5 — always wrap as.matrix()
  # See @primitives/seurat_v5_rules.md Rule 6
  avg <- as.matrix(AverageExpression(so_ds, assay = "RNA", layer = "data",
                                     group.by = LABEL_COL,
                                     features = hvg)$RNA)
  avg_z          <- scale(avg)  # z-score within dataset (scales columns = cell types)
  avg_list[[ds]] <- avg_z
}

# Average z-scores across datasets for cell types present in multiple datasets
all_ct_all <- unique(unlist(lapply(avg_list, colnames)))
z_mat <- matrix(NA, nrow = length(hvg), ncol = length(all_ct_all),
                dimnames = list(hvg, all_ct_all))
for (ct in all_ct_all) {
  ct_zscores <- lapply(avg_list, function(m) if (ct %in% colnames(m)) m[, ct] else NULL)
  ct_zscores <- Filter(Negate(is.null), ct_zscores)
  if (length(ct_zscores) > 0)
    z_mat[, ct] <- rowMeans(do.call(cbind, ct_zscores), na.rm = TRUE)
}

# Pearson correlation
cor_mat <- cor(z_mat, use = "pairwise.complete.obs")
```

### Plotting

```r
# Diverging blue → white → red palette (6-stop; see @context/color_palettes.md)
pal <- colorRampPalette(c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#D6604D","#B2182B"))(100)

ph <- pheatmap(cor_mat,
               color             = pal,
               clustering_method = "ward.D2",
               display_numbers   = TRUE,
               number_format     = "%.2f",
               fontsize_number   = 6.5,
               fontsize_row      = 11,
               fontsize_col      = 11,
               border_color      = "white",
               main              = "")   # leave blank when using grid title overlay below
```

**Simple output (no styled subtitle) — use ggsave:**
```r
n     <- nrow(cor_mat)
sz    <- 13
fig_h <- n * (sz / 72) + 3.5
fig_w <- fig_h + 2   # extra width to prevent left-side clipping

ggsave(file.path(OUTPUT_DIR, "correlation_heatmap.pdf"),
       plot = ph$gtable,
       width = fig_w, height = fig_h, units = "in",
       device = "pdf", useDingbats = FALSE)
```

**With styled title + subtitle — intentional pdf()/dev.off() exception:**
pheatmap + grid.draw requires an active device for the title overlay.
This is a documented exception analogous to ComplexHeatmap.
```r
title    <- project_specific  # REPLACE: e.g., paste(project_name, "Cell Type Correlation")
subtitle <- project_specific  # REPLACE: e.g., "Pearson correlation of scaled HVG expression"

pdf(file.path(OUTPUT_DIR, "correlation_heatmap.pdf"), width = fig_w, height = fig_h)
grid::grid.newpage()
grid::pushViewport(grid::viewport(y = 0, height = 0.93, just = "bottom"))
grid::grid.draw(ph$gtable)
grid::popViewport()
grid::grid.text(title,    x = 0.5, y = 0.985, gp = gpar(fontsize = 15, fontface = "bold"))
grid::grid.text(subtitle, x = 0.5, y = 0.955, gp = gpar(fontsize = 10, col = "grey55"))
dev.off()
```

### Factor level mismatch fix (annotation_colors)

```r
# Always fill unmapped factor levels to avoid pheatmap error:
# "Factor levels on variable X do not match with annotation_colors"
unmapped <- setdiff(levels(ann$`Cell Type`), names(ct_cols))
if (length(unmapped) > 0) {
  new_cols <- setNames(scales::hue_pal()(length(unmapped)), unmapped)
  ct_cols  <- c(ct_cols, new_cols)
}
```

---

## B. Chord Diagram (circlize)

Visualizes dataset ↔ cell type composition — how many cells of each cell type
come from each dataset.

```r
library(circlize)
library(dplyr)

# Build count matrix — CRITICAL: use a for loop, NOT pivot_wider
# pivot_wider produces NA for absent combinations; chordDiagramFromMatrix
# fails with "missing value where TRUE/FALSE needed"
meta_df <- so@meta.data
all_ds  <- sort(unique(meta_df[[DATASET_COL]]))
all_ct  <- sort(unique(meta_df[[LABEL_COL]]))

chord_mat <- matrix(0, nrow = length(all_ds), ncol = length(all_ct),
                    dimnames = list(all_ds, all_ct))
for (i in seq_len(nrow(meta_df))) {
  d <- meta_df[[DATASET_COL]][i]
  c <- meta_df[[LABEL_COL]][i]
  chord_mat[d, c] <- chord_mat[d, c] + 1
}

# Color matrix — chord color by source dataset
chord_col_mat <- matrix(dataset_colors[all_ds], nrow = length(all_ds),
                        ncol = length(all_ct), dimnames = dimnames(chord_mat))

# Apply transparency — CRITICAL: adjustcolor on matrix loses dimensions
# Must restore dimensions and dimnames after adjustcolor
chord_col_mat <- matrix(
  adjustcolor(chord_col_mat, alpha.f = 0.55),
  nrow = nrow(chord_col_mat), ncol = ncol(chord_col_mat),
  dimnames = dimnames(chord_col_mat)
)

# Plot
pdf(file.path(OUTPUT_DIR, "chord_diagram.pdf"),
    width = project_specific, height = project_specific)  # REPLACE: typically 7–9 inches square

circos.clear()
circos.par(gap.after = c(rep(2, length(all_ds) - 1), 10,
                          rep(2, length(all_ct) - 1), 10))
chordDiagramFromMatrix(
  chord_mat,
  col             = chord_col_mat,
  grid.col        = c(dataset_colors, ct_colors),
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.08)
)

# Sector labels — use circos.track only, NOT circos.trackPlotRegion
# Using both causes duplicate labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.1,
              CELL_META$sector.index, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 0.75)
}, bg.border = NA)

title("Dataset x Cell Type", cex.main = 1.2)
circos.clear()
dev.off()
```

**Note:** circlize uses base R graphics and has no ggsave path.
The `pdf()/dev.off()` pattern here is a documented exception for base-graphics output.

---

## C. Proportion Heatmap

Percentage of each cell type per sample, with configurable metadata annotation strips.

```r
library(pheatmap)
library(dplyr)

# Build proportion matrix: rows = cell types, cols = samples
prop_mat <- prop.table(table(so[[LABEL_COL]], so[[SAMPLE_COL]]), margin = 2) * 100
prop_mat  <- as.matrix(prop_mat)

# Order columns: by dataset, then condition, then sample
# Adjust sort columns to match your brief's metadata structure
sort_cols <- c(DATASET_COL, CONDITION_COL, SAMPLE_COL)
sort_cols_present <- intersect(sort_cols, colnames(so@meta.data))

col_order <- so@meta.data %>%
  select(all_of(sort_cols_present)) %>%
  distinct() %>%
  arrange(across(all_of(sort_cols_present)))

prop_mat <- prop_mat[, col_order[[SAMPLE_COL]], drop = FALSE]

# Annotation tracks — specify which metadata columns to include as strips
# ann_track_cols: character vector of metadata column names to show as annotation
# Adjust to match your brief; omit columns with sparse or uninformative data
ann_track_cols <- c(DATASET_COL, ORGAN_COL, CONDITION_COL, DATA_TYPE_COL)
ann_track_cols <- intersect(ann_track_cols, colnames(so@meta.data))

ann_col <- as.data.frame(
  col_order[, ann_track_cols, drop = FALSE]
)
rownames(ann_col) <- col_order[[SAMPLE_COL]]

# Annotation color list — must cover all levels of all annotation columns
# Caller provides ann_colors as a named list (one named vector per column)
# Example structure:
# ann_colors <- list(
#   dataset   = dataset_colors,      # named vector
#   condition = condition_colors,    # named vector — project_specific
#   organ     = organ_colors,        # named vector — project_specific
#   data_type = data_type_colors     # named vector — project_specific
# )
ann_colors <- project_specific  # REPLACE with list of named color vectors from brief

# Fill palette
fill_pal <- colorRampPalette(c("#F5F5F5","#BBDEFB","#1565C0"))(100)

# Dataset separator positions
dataset_gap_positions <- which(diff(as.integer(factor(col_order[[DATASET_COL]]))) != 0)

ph <- pheatmap(prop_mat,
               color             = fill_pal,
               annotation_col    = ann_col,
               annotation_colors = ann_colors,
               cluster_cols      = FALSE,
               cluster_rows      = TRUE,
               display_numbers   = FALSE,
               fontsize          = 16,
               fontsize_row      = 15,
               fontsize_col      = 12,
               border_color      = "white",
               gaps_col          = dataset_gap_positions,
               main              = "")
```

### Figure sizing and output

```r
# Proportion heatmaps scale poorly with dynamic formulas —
# calibrate empirically for your number of samples and cell types.
# Starting formula: fig_w = n_samples * 0.3 + 4, fig_h = n_celltypes * 0.35 + 2
n_samples   <- ncol(prop_mat)
n_celltypes <- nrow(prop_mat)
fig_w <- max(n_samples * 0.3 + 4, 8)   # REVIEW: adjust after initial render
fig_h <- max(n_celltypes * 0.35 + 2, 5) # REVIEW: adjust after initial render

ggsave(file.path(OUTPUT_DIR, "proportion_heatmap.pdf"),
       plot   = ph$gtable,
       width  = fig_w, height = fig_h,
       units  = "in", device = "pdf", useDingbats = FALSE)
```

### Metadata strip notes
- Include strips that have meaningful data — omit sparse ones
- Continuous variables (e.g., Age): use gradient color scale
- Categorical variables: use consistent colors matching other plots in project
- Dataset strip colors: use `dataset_colors` from brief
