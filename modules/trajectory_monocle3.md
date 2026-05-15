---
requires_context:
  palettes:
    - subtype_colors   # optional — named vector for cell subtype colors in celltype UMAP
  metadata_columns:
    required:
      - label_col      # cell subtype label column used for cluster injection into CDS
    optional:
      - group_col      # grouping column for pseudotime violin (tissue, condition, etc.)
  brief_keys:
    required:
      - output_dir
      - downstream_analyses.trajectory.label_col
      - downstream_analyses.trajectory.start_cluster
    optional:
      - downstream_analyses.trajectory.exclude_pattern
      - downstream_analyses.trajectory.end_cluster
      - downstream_analyses.trajectory.supplementary_plots
      - downstream_analyses.trajectory.group_col
references:
  - "@primitives/visualization.md"
  - "@primitives/seurat_v5_rules.md"
  - "@context/color_palettes.md"
---

# Module: Trajectory Analysis — Monocle3

Pseudotime trajectory analysis on a Seurat subset object using Monocle3.
Bypasses SeuratWrappers and preprocess_cds entirely — both are unavailable or
cause ABI errors in r-env. Injects Seurat Harmony embeddings and UMAP coordinates
directly into the CDS object.

---

## Critical Constraints

| ❌ Don't use | ✅ Use instead | Why |
|---|---|---|
| `SeuratWrappers::as.cell_data_set()` | `new_cell_data_set()` manually (Step 2) | SeuratWrappers not installed in r-env |
| `preprocess_cds()` | Skip entirely | Triggers `as_cholmod_sparse` ABI error with irlba/Matrix |
| `learn_graph(use_partition=TRUE)` | `use_partition=FALSE` | Partitions fragment trajectory on subset objects |
| `pdf()` / `dev.off()` for ggplot outputs | `ggsave(..., device = "pdf", useDingbats = FALSE)` | Convention — see CONVENTIONS.md §4 |

---

## R Packages Required

```r
library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)
```

---

## Configuration Block

```r
# ── CONFIG ────────────────────────────────────────────────────────────────────
LABEL_COL       <- "project_specific"   # REPLACE: metadata column holding cell subtype labels
                                          # (brief: downstream_analyses.trajectory.label_col)
GROUP_COL       <- NULL                  # optional: grouping column for violin faceting
                                          # (brief: downstream_analyses.trajectory.group_col)
EXCLUDE_PATTERN <- NULL                  # optional: regex to exclude cell types before trajectory
                                          # (brief: downstream_analyses.trajectory.exclude_pattern)
                                          # Example: "^AEC" catches AEC1, AEC2, etc.
START_CLUSTER   <- "project_specific"   # REPLACE: subtype label value to use as trajectory root
                                          # (brief: downstream_analyses.trajectory.start_cluster)
SUBSET_NAME     <- "celltype_traj"       # short name used in output filenames
OUTPUT_DIR      <- "output/trajectory"
SUPPLEMENTARY_PLOTS <- TRUE
```

---

## Step 1: Exclude Non-Trajectory Cell Types (Optional)

Use pattern matching, not exact string matching — subtype names may be numbered:

```r
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

obj_traj <- readRDS("project_specific")   # REPLACE: path to Seurat subset object

if (!is.null(EXCLUDE_PATTERN)) {
  keep_cells <- colnames(obj_traj)[
    !grepl(EXCLUDE_PATTERN, obj_traj@meta.data[[LABEL_COL]])
  ]
  obj_traj <- obj_traj[, keep_cells]
  message(sprintf("After exclusion: %d cells remaining", ncol(obj_traj)))
}
```

---

## Step 2: Build CDS Manually

```r
# Extract expression matrix
expr_mat <- GetAssayData(obj_traj, assay = "RNA", layer = "counts")

# Cell metadata and gene metadata
cell_meta <- obj_traj@meta.data
gene_meta  <- data.frame(gene_short_name = rownames(expr_mat),
                          row.names       = rownames(expr_mat))

# Create CDS — no preprocess_cds() call
cds <- new_cell_data_set(
  expression_data = expr_mat,
  cell_metadata   = cell_meta,
  gene_metadata   = gene_meta
)
```

---

## Step 3: Inject Seurat Embeddings

```r
# Inject Harmony embeddings as PCA — see @primitives/seurat_v5_rules.md Rule 1
# (JoinLayers must have been called before GetAssayData above if multi-batch)
reducedDims(cds)[["PCA"]] <- Embeddings(obj_traj, "harmony")

# Detect UMAP reduction name — see @primitives/seurat_v5_rules.md Rule 5
umap_key <- names(obj_traj@reductions)[
  grepl("umap", names(obj_traj@reductions), ignore.case = TRUE)][1]
umap_coords <- Embeddings(obj_traj, umap_key)
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
reducedDims(cds)[["UMAP"]] <- umap_coords

# Inject cluster assignments from the label column
cluster_vec   <- setNames(factor(obj_traj@meta.data[[LABEL_COL]]), colnames(obj_traj))
partition_vec <- setNames(factor(rep(1, ncol(obj_traj))),          colnames(obj_traj))
cds@clusters[["UMAP"]] <- list(
  clusters   = cluster_vec,
  partitions = partition_vec
)
```

---

## Step 4: Learn Graph and Order Cells

```r
cds <- learn_graph(cds, use_partition = FALSE)

# Find root cells from start_cluster
root_cells <- colnames(cds)[cds@clusters[["UMAP"]]$clusters == START_CLUSTER]
if (length(root_cells) == 0) {
  stop(sprintf("START_CLUSTER '%s' not found in %s. Available: %s",
               START_CLUSTER, LABEL_COL,
               paste(unique(obj_traj@meta.data[[LABEL_COL]]), collapse = ", ")))
}

cds <- order_cells(cds, root_cells = root_cells)
```

---

## Step 5: Save Pseudotime Back to Seurat Object

```r
obj_traj$pseudotime <- pseudotime(cds)
saveRDS(obj_traj, file.path(OUTPUT_DIR, paste0(SUBSET_NAME, "_with_pseudotime.Rds")))
```

---

## Output Plots

### Main Trajectory Plots

```r
# Pseudotime colored UMAP — magma palette (NEVER inferno; magma is the validated convention)
# See @context/color_palettes.md PALETTE_PSEUDOTIME
p_pseudo <- plot_cells(cds,
  color_cells_by      = "pseudotime",
  cell_size           = 0.6,
  label_groups_by_cluster = FALSE,
  label_leaves        = TRUE,
  label_branch_points = TRUE
) +
  theme_classic(base_size = 16) +
  scale_color_viridis_c(option = "magma")   # magma, NOT inferno

ggsave(file.path(OUTPUT_DIR, paste0(SUBSET_NAME, "_trajectory_pseudotime.pdf")),
       p_pseudo, width = 7, height = 6, device = "pdf", useDingbats = FALSE)

# Cell type colored UMAP with trajectory overlaid
p_celltype <- plot_cells(cds,
  color_cells_by          = LABEL_COL,
  cell_size               = 0.6,
  label_groups_by_cluster = TRUE
) + theme_classic(base_size = 16)

# Apply subtype_colors if provided (from brief context_overrides.palettes.subtype_colors)
# subtype_colors <- c("project_specific" = "#project_specific")
# Uncomment and populate:
# p_celltype <- p_celltype + scale_color_manual(values = subtype_colors)

ggsave(file.path(OUTPUT_DIR, paste0(SUBSET_NAME, "_trajectory_celltypes.pdf")),
       p_celltype, width = 7, height = 6, device = "pdf", useDingbats = FALSE)
```

---

### Supplementary Plots

```r
if (SUPPLEMENTARY_PLOTS) {
  pseudo_df <- data.frame(
    pseudotime = pseudotime(cds),
    subtype    = cds@clusters[["UMAP"]]$clusters
  ) %>% filter(is.finite(pseudotime))

  # 1. Pseudotime density by cell subtype
  p_density <- ggplot(pseudo_df, aes(x = pseudotime, fill = subtype, color = subtype)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    theme_classic(base_size = 14) +
    labs(x = "Pseudotime", y = "Density",
         title = "Pseudotime Distribution by Cell Subtype")

  ggsave(file.path(OUTPUT_DIR, paste0(SUBSET_NAME, "_pseudotime_density.pdf")),
         p_density, width = 8, height = 4, device = "pdf", useDingbats = FALSE)

  # 2. Pseudotime violin by group variable (if GROUP_COL is set)
  if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj_traj@meta.data)) {
    pseudo_df$group_var <- obj_traj@meta.data[rownames(pseudo_df), GROUP_COL]

    median_order <- pseudo_df %>%
      group_by(group_var) %>%
      summarise(med = median(pseudotime)) %>%
      arrange(med) %>%
      pull(group_var)

    p_violin <- ggplot(pseudo_df, aes(x = factor(group_var, levels = median_order),
                                       y = pseudotime, fill = group_var)) +
      geom_violin(width = 0.8, linewidth = 0.3) +
      geom_boxplot(width = 0.07, outlier.shape = NA, fill = "white", linewidth = 0.3) +
      stat_summary(fun = median, geom = "point",
                   size = 2.5, color = "#E53935") +   # red median dot — validated convention
      theme_classic(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(x = GROUP_COL, y = "Pseudotime")

    ggsave(file.path(OUTPUT_DIR, paste0(SUBSET_NAME, "_pseudotime_violin_", GROUP_COL, ".pdf")),
           p_violin, width = 8, height = 5, device = "pdf", useDingbats = FALSE)
  }
}
```

---

## Brief Configuration

```yaml
downstream_analyses:
  trajectory:
    enabled: true
    method: monocle3
    label_col: project_specific      # REPLACE: cell subtype label column
    exclude_pattern: null            # optional regex — e.g. "^AEC" catches AEC1, AEC2, etc.
    start_cluster: project_specific  # REPLACE: must match exact subtype label in annotated object
    end_cluster: project_specific    # informational — monocle3 infers endpoints from graph
    group_col: null                  # optional: grouping column for violin plot
    supplementary_plots: true
```

---

## Output File Summary

| File | Description |
|---|---|
| `{SUBSET_NAME}_with_pseudotime.Rds` | Seurat object with pseudotime added as metadata column |
| `{SUBSET_NAME}_trajectory_pseudotime.pdf` | Pseudotime UMAP (magma palette) |
| `{SUBSET_NAME}_trajectory_celltypes.pdf` | Cell type UMAP with trajectory overlaid |
| `{SUBSET_NAME}_pseudotime_density.pdf` | Density curves by subtype (supplementary) |
| `{SUBSET_NAME}_pseudotime_violin_{GROUP_COL}.pdf` | Pseudotime violin by group (supplementary; if GROUP_COL set) |

---

## Project-Specific Values (Stage for Phase 4 examples/)

- `examples/humanfat_trajectory.md`:
  - `LABEL_COL = "ec_subtype"`, `EXCLUDE_PATTERN = "^AEC"`, `START_CLUSTER = "CapEC2"`
  - HumanFat EC subtype colors
  - Validated on 53,487 EC cells, 13 subtypes (HumanFat_Yang run 2026-03-05)
