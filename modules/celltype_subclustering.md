---
requires_context:
  palettes:
    - group_colors     # optional — named vector keyed by group_col values (tissue/condition)
    - subtype_colors   # optional — named vector for final subtype labels (Phase 2 endpoint)
  metadata_columns:
    required:
      - label_col      # cell type label column in whole object (for subsetting)
    optional:
      - group_col      # group/tissue/condition column for proportion plots and exclusion
      - sample_col     # sample/batch column (default: batch_correction_var)
      - batch_col      # batch correction column (default: batch_correction_var from brief)
  brief_keys:
    required:
      - output_dir
      - pipeline_params.batch_correction_var
    optional:
      - downstream_analyses.subclustering.target_celltypes
      - downstream_analyses.subclustering.resolution
      - downstream_analyses.subclustering.n_pcs
      - downstream_analyses.subclustering.keep_clusters
      - downstream_analyses.subclustering.exclude_groups
references:
  - "@primitives/harmony_integration.md"
  - "@primitives/visualization.md"
  - "@primitives/seurat_v5_rules.md"
  - "@context/color_palettes.md"
---

# Module: Cell Type Subclustering, Sub-labeling, and Endpoint Plots

Takes a broad cell type from a pre-annotated whole object, subsets and re-clusters it,
explores cluster identity through diagnostic plots and DEG analysis, then applies
user-approved labels to produce a clean set of endpoint publication figures.

**Two-phase structure:**
- **Phase 1 — Exploration (`output/subclustering/`):** subset → re-cluster → diagnostic plots → DEG
  → user reviews and provides label mapping + clusters to remove
- **Phase 2 — Endpoint (`output/subclustering/endpoint/`):** apply labels, subset to
  keep-clusters, produce final UMAP + proportion plots + curated marker dotplot

---

## Critical Constraints

| ❌ Don't | ✅ Do instead | Why |
|---|---|---|
| `subset(obj, seurat_clusters %in% c(0,1,2))` | `obj[, colnames(obj)[as.character(obj$seurat_clusters) %in% as.character(c(0,1,2))]]` | `seurat_clusters` is a factor; `subset()` throws "No cell overlap" with integer comparisons |
| `obj$label <- label_map[as.character(obj$seurat_clusters)]` | `obj$label <- unname(label_map[as.character(obj$seurat_clusters)])` | Without `unname()`, the vector carries cluster numbers as names; Seurat matches against cell barcodes and throws "No cell overlap" |
| `Embeddings(obj, "umap")[,"UMAP_1"]` | `as.data.frame(Embeddings(obj, "umap"))$umap_1` | UMAP column names are `umap_1`/`umap_2` (lowercase underscore) in r-env |
| Pass gene vector directly to `DotPlot` | `unique(genes_vec)` before passing | Duplicate gene names cause a factor level duplication error |
| Use base R | Always run via `@primitives/r_environment.md` | `sp` package is broken in base R |
| Use `pdf()` / `dev.off()` for ggplot outputs | `ggsave(..., device = "pdf", useDingbats = FALSE)` | Convention — see CONVENTIONS.md §4 |

---

## Phase 1: Exploration

### Step 1a: Subset and Re-cluster

```r
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(patchwork)
set.seed(42)

# ── CONFIG ────────────────────────────────────────────────────────────────────
# Target cell types to subset from the whole object
CELLTYPE_LABELS <- c("project_specific")  # REPLACE: values from whole-object label_col
                    # Example: c("Macrophage", "Monocyte") — use blank template
                    # Structure: character vector of one or more cell type label strings

LABEL_COL       <- "project_specific"   # REPLACE: metadata column in whole object
                                         # holding the cell type label (brief: metadata.label_col)

BATCH_CORRECTION_VAR <- "project_specific"  # REPLACE: from brief pipeline_params.batch_correction_var
                                             # column used for Harmony batch correction

GROUP_COL       <- "project_specific"   # REPLACE: metadata column for tissue/condition grouping
                                         # used in proportion plots (brief: metadata.group_col)
                                         # Set to NULL to skip group-aware proportion plots

SAMPLE_COL      <- BATCH_CORRECTION_VAR  # sample/batch column for per-sample plots
                                          # defaults to batch_correction_var; override if different

RDS_IN          <- "project_specific"   # REPLACE: path to whole-object RDS
RDS_OUT         <- file.path(OUTPUT_DIR, "CellType_subclustered.Rds")
OUTPUT_DIR      <- "output/subclustering"
ENDPOINT_DIR    <- "output/subclustering/endpoint"

N_PCS           <- 50
NEIGHBOR_DIMS   <- 1:25
CLUSTER_RES     <- 0.4

# ── LOAD AND SUBSET ──────────────────────────────────────────────────────────
dir.create(OUTPUT_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(ENDPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

whole <- readRDS(RDS_IN)
obj   <- subset(whole, subset = !!sym(LABEL_COL) %in% CELLTYPE_LABELS)
rm(whole); gc()

# ── RE-CLUSTER (see @primitives/harmony_integration.md for full recipe) ──────
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = N_PCS)
obj <- RunHarmony(obj, group.by.vars = BATCH_CORRECTION_VAR, reduction.use = "pca",
                  reduction.save = "harmony")
obj <- FindNeighbors(obj, reduction = "harmony", dims = NEIGHBOR_DIMS)
obj <- FindClusters(obj, resolution = CLUSTER_RES)
obj <- RunUMAP(obj, reduction = "harmony", dims = NEIGHBOR_DIMS)

saveRDS(obj, RDS_OUT)
```

**If a subset of clusters needs further sub-clustering:**
```r
mac_clusters <- c("2", "4", "7")   # cluster IDs that are enriched for target subtype
obj2 <- obj[, as.character(obj$seurat_clusters) %in% mac_clusters]
# re-run NormalizeData → RunUMAP as above
```

**If a contaminating cluster must be removed:**
```r
obj3 <- obj2[, as.character(obj2$seurat_clusters) != "1"]
# re-run NormalizeData → RunUMAP
```

---

### Step 1b: Diagnostic Plots (Exploration)

**P2-7 fix:** All plots use `ggsave()` — never `pdf()` / `dev.off()` for ggplot outputs.

**4-panel UMAP overview (patchwork):**
```r
# Build group panel only if GROUP_COL is specified
p1 <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE, label.size = 5)
p2 <- DimPlot(obj, group.by = LABEL_COL,         label = TRUE, label.size = 4, repel = TRUE)

if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  p3 <- DimPlot(obj, group.by = GROUP_COL)           # group_colors added in Phase 2 once known
  p4 <- DimPlot(obj, group.by = SAMPLE_COL)
  panel4 <- (p1 | p2) / (p3 | p4)
} else {
  p3 <- DimPlot(obj, group.by = SAMPLE_COL)
  panel4 <- (p1 | p2) / (p3 | plot_spacer())
}

ggsave(file.path(OUTPUT_DIR, "CellType_umap_overview.pdf"),
       panel4, width = 16, height = 12, device = "pdf", useDingbats = FALSE)
```

**Cluster UMAP (individual):**
```r
ggsave(file.path(OUTPUT_DIR, "CellType_umap_clusters.pdf"),
       DimPlot(obj, group.by = "seurat_clusters", label = TRUE, label.size = 5),
       width = 8, height = 7, device = "pdf", useDingbats = FALSE)
```

**Canonical marker feature plots with cluster number overlaid:**
```r
umap_df   <- as.data.frame(Embeddings(obj, "umap"))   # columns: umap_1, umap_2
umap_df$cluster <- as.character(obj$seurat_clusters)
centroids <- umap_df %>%
  group_by(cluster) %>%
  summarise(x = median(umap_1), y = median(umap_2), .groups = "drop")

marker_genes <- unique(c("project_specific"))   # REPLACE: canonical markers for this cell type

plots <- lapply(marker_genes, function(g) {
  FeaturePlot(obj, features = g, reduction = "umap", order = TRUE,
              min.cutoff = "q10") +
    geom_text(data = centroids, aes(x = x, y = y, label = cluster),
              inherit.aes = FALSE, size = 3, fontface = "bold", color = "black")
})

ggsave(file.path(OUTPUT_DIR, "CellType_feature_grid.pdf"),
       wrap_plots(plots, ncol = 4),
       width = 16, height = 12, device = "pdf", useDingbats = FALSE)
```

**Canonical marker dot plot:**
```r
Idents(obj) <- "seurat_clusters"
ggsave(file.path(OUTPUT_DIR, "CellType_marker_dotplot.pdf"),
       DotPlot(obj, features = unique(marker_genes), group.by = "seurat_clusters",
               dot.scale = 6) + RotatedAxis(),
       width = 16, height = max(5, length(levels(obj$seurat_clusters)) * 0.5 + 3),
       device = "pdf", useDingbats = FALSE)
```

**Group proportion — clusters on x (only when GROUP_COL is set):**
```r
if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  prop_df <- obj@meta.data %>%
    group_by(seurat_clusters, !!sym(GROUP_COL)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(prop = n / sum(n))

  p_prop <- ggplot(prop_df, aes(x = seurat_clusters, y = prop,
                                 fill = !!sym(GROUP_COL))) +
    geom_col(position = "stack") +
    # NOTE: group_colors — supply from brief context_overrides.palettes.group_colors
    # or leave unset for automatic ggplot palette
    theme_classic(base_size = 13) +
    labs(x = "Cluster", y = "Proportion")

  ggsave(file.path(OUTPUT_DIR, "CellType_group_proportion.pdf"),
         p_prop, width = 10, height = 5, device = "pdf", useDingbats = FALSE)
}
```

**Per-sample proportion (faceted by group) — colorblind-safe Paul Tol palette for clusters:**
```r
# Paul Tol muted palette — safe for red-green color blindness, up to 10 clusters
tol_muted <- c("#332288","#88CCEE","#44AA99","#117733","#999933",
               "#DDCC77","#CC6677","#882255","#AA4499","#DDDDDD")
cluster_levels <- as.character(sort(unique(as.integer(as.character(obj$seurat_clusters)))))
cluster_colors <- setNames(tol_muted[seq_len(length(cluster_levels))], cluster_levels)

if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  prop_sample <- obj@meta.data %>%
    group_by(!!sym(SAMPLE_COL), !!sym(GROUP_COL), seurat_clusters) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(!!sym(SAMPLE_COL)) %>%
    mutate(prop = n / sum(n),
           seurat_clusters = factor(seurat_clusters, levels = cluster_levels))

  sample_order <- prop_sample %>%
    distinct(!!sym(SAMPLE_COL), !!sym(GROUP_COL)) %>%
    arrange(!!sym(GROUP_COL), !!sym(SAMPLE_COL)) %>%
    pull(!!sym(SAMPLE_COL))
  prop_sample[[SAMPLE_COL]] <- factor(prop_sample[[SAMPLE_COL]], levels = sample_order)

  ggsave(file.path(OUTPUT_DIR, "CellType_proportion_by_sample_clusters.pdf"),
         ggplot(prop_sample, aes(x = !!sym(SAMPLE_COL), y = prop,
                                  fill = seurat_clusters)) +
           geom_col(position = "stack", color = "white", linewidth = 0.3) +
           scale_fill_manual(values = cluster_colors, name = "Cluster") +
           facet_grid(~ !!sym(GROUP_COL), scales = "free_x", space = "free_x") +
           scale_x_discrete(guide = guide_axis(angle = 45)) +
           labs(x = "Sample", y = "Proportion") +
           theme_classic(base_size = 14) +
           theme(strip.background = element_rect(fill = "grey92", color = NA),
                 strip.text = element_text(size = 12, face = "bold")),
         width = 14, height = 5.5, device = "pdf", useDingbats = FALSE)
}
```

**FindAllMarkers + CSV handoff:**
```r
# See @primitives/seurat_v5_rules.md Rule 1 — JoinLayers required before FindAllMarkers
obj <- JoinLayers(obj)
all_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = 0.4, test.use = "wilcox")
write.csv(all_markers, file.path(OUTPUT_DIR, "CellType_cluster_markers.csv"), row.names = FALSE)

top5 <- all_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
write.csv(top5, file.path(OUTPUT_DIR, "CellType_top5_markers.csv"), row.names = FALSE)

top100 <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 100)
write.csv(top100, file.path(OUTPUT_DIR, "CellType_DEG_top100_significant.csv"), row.names = FALSE)
```

**DEG between specific ambiguous clusters (optional):**
```r
# Example: compare clusters 0, 1, 2 pairwise or 0 vs rest
obj_sub <- obj[, as.character(obj$seurat_clusters) %in% c("0","1","2")]
obj_sub <- JoinLayers(obj_sub)
Idents(obj_sub) <- "seurat_clusters"
deg_012 <- FindAllMarkers(obj_sub, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.3)
write.csv(deg_012, file.path(OUTPUT_DIR, "CellType_DEG_clusters012.csv"), row.names = FALSE)
```

---

### Step 1c: Handoff to User

After generating exploration plots, output to the user:
1. Path to cluster UMAP PDF
2. Path to feature grid PDF(s)
3. Path to `CellType_DEG_top100_significant.csv`
4. Summary table of cluster sizes and dominant group per cluster

Request from user:
- Cluster → label mapping (many-to-one: multiple clusters → same label)
- Which clusters to remove entirely
- (Optional) Which group values to exclude from endpoint plots

---

## Phase 2: Endpoint Plots

### Step 2a: Apply Labels and Remove Unwanted Clusters

```r
# ── USER-PROVIDED MAPPING ────────────────────────────────────────────────────
# Fill in cluster → subtype label mapping (blank template below)
label_map <- c(
  "0" = "project_specific",   # REPLACE with actual subtype names
  "1" = "project_specific",
  "2" = "project_specific"
  # Map multiple clusters to the same label when biology warrants it
  # Omit clusters to remove them (they will not appear in keep_clusters)
)
keep_clusters <- names(label_map)

# ── OPTIONAL: exclude specific group values ───────────────────────────────────
# Set to NULL to keep all; or c("GroupA", "GroupB") to remove
EXCLUDE_GROUPS <- NULL  # REPLACE from brief: downstream_analyses.subclustering.exclude_groups

# ── FILTER CELLS ─────────────────────────────────────────────────────────────
keep_cells <- colnames(obj)[as.character(obj$seurat_clusters) %in% keep_clusters]

if (!is.null(EXCLUDE_GROUPS) && !is.null(GROUP_COL) &&
    GROUP_COL %in% colnames(obj@meta.data)) {
  keep_cells <- intersect(
    keep_cells,
    colnames(obj)[!obj@meta.data[[GROUP_COL]] %in% EXCLUDE_GROUPS]
  )
}
obj <- obj[, keep_cells]

# ── ASSIGN LABELS ─────────────────────────────────────────────────────────────
# MUST use unname() — see Critical Constraints above
obj$subtype_label <- unname(label_map[as.character(obj$seurat_clusters)])
obj$subtype_label <- factor(obj$subtype_label, levels = unique(label_map))

Idents(obj) <- "subtype_label"
saveRDS(obj, file.path(ENDPOINT_DIR, "CellType_labeled.Rds"))
```

---

### Step 2b: Labeled UMAP

```r
# ── SUBTYPE COLORS ────────────────────────────────────────────────────────────
# Provide from brief context_overrides.palettes.subtype_colors, or define here
subtype_cols <- c(
  "project_specific" = "#project_specific"  # REPLACE: one entry per subtype label
)
# Example structure: c("SubtypeA" = "#E05C5C", "SubtypeB" = "#4A90D9")

p <- DimPlot(obj, group.by = "subtype_label", cols = subtype_cols,
             pt.size = 0.6, label = FALSE) +
  ggtitle("Cell Type Subtypes") +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.title = element_blank(),
    axis.title   = element_text(size = 12),
    axis.text    = element_blank(),
    axis.ticks   = element_blank()
  )

ggsave(file.path(ENDPOINT_DIR, "CellType_labeled_umap.pdf"), p,
       width = 9, height = 5.5, device = "pdf", useDingbats = FALSE)
ggsave(file.path(ENDPOINT_DIR, "CellType_labeled_umap.png"), p,
       width = 9, height = 5.5, dpi = 200)
```

---

### Step 2c: Proportion Plots (Group on x, Labels as fill)

```r
# ── GROUP COLORS ──────────────────────────────────────────────────────────────
# Provide from brief context_overrides.palettes.group_colors
group_cols <- c("project_specific" = "#project_specific")  # REPLACE

if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  # By group — x = group, fill = subtype label
  prop_by_group <- obj@meta.data %>%
    group_by(!!sym(GROUP_COL), subtype_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(!!sym(GROUP_COL)) %>%
    mutate(prop = n / sum(n),
           subtype_label = factor(subtype_label, levels = names(subtype_cols)))

  p_group <- ggplot(prop_by_group, aes(x = !!sym(GROUP_COL), y = prop,
                                        fill = subtype_label)) +
    geom_col(position = "stack", color = "white", linewidth = 0.3) +
    scale_fill_manual(values = subtype_cols, name = "Subtype") +
    scale_x_discrete(guide = guide_axis(angle = 35)) +
    labs(x = "Group", y = "Proportion") +
    theme_classic(base_size = 14) +
    theme(axis.text  = element_text(size = 13),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14))

  ggsave(file.path(ENDPOINT_DIR, "CellType_proportion_by_group.pdf"), p_group,
         width = 9, height = 5.5, device = "pdf", useDingbats = FALSE)
  ggsave(file.path(ENDPOINT_DIR, "CellType_proportion_by_group.png"), p_group,
         width = 9, height = 5.5, dpi = 200)

  # Inverse — x = subtype label, fill = group
  prop_by_label <- obj@meta.data %>%
    group_by(subtype_label, !!sym(GROUP_COL)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(subtype_label) %>%
    mutate(prop = n / sum(n))

  p_label <- ggplot(prop_by_label, aes(x = subtype_label, y = prop,
                                         fill = !!sym(GROUP_COL))) +
    geom_col(position = "stack", color = "white", linewidth = 0.3) +
    scale_fill_manual(values = group_cols, name = "Group") +
    scale_x_discrete(guide = guide_axis(angle = 35)) +
    labs(x = "Subtype", y = "Proportion") +
    theme_classic(base_size = 14) +
    theme(axis.text  = element_text(size = 13),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14))

  ggsave(file.path(ENDPOINT_DIR, "CellType_proportion_by_label.pdf"), p_label,
         width = 9, height = 5.5, device = "pdf", useDingbats = FALSE)
  ggsave(file.path(ENDPOINT_DIR, "CellType_proportion_by_label.png"), p_label,
         width = 9, height = 5.5, dpi = 200)

  # Per sample — x = sample, faceted by group, fill = subtype label
  prop_sample <- obj@meta.data %>%
    group_by(!!sym(SAMPLE_COL), !!sym(GROUP_COL), subtype_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(!!sym(SAMPLE_COL)) %>%
    mutate(prop = n / sum(n),
           subtype_label = factor(subtype_label, levels = names(subtype_cols)))

  sample_order <- prop_sample %>%
    distinct(!!sym(SAMPLE_COL), !!sym(GROUP_COL)) %>%
    arrange(!!sym(GROUP_COL), !!sym(SAMPLE_COL)) %>%
    pull(!!sym(SAMPLE_COL))
  prop_sample[[SAMPLE_COL]] <- factor(prop_sample[[SAMPLE_COL]], levels = sample_order)

  p_sample <- ggplot(prop_sample, aes(x = !!sym(SAMPLE_COL), y = prop,
                                       fill = subtype_label)) +
    geom_col(position = "stack", color = "white", linewidth = 0.3) +
    scale_fill_manual(values = subtype_cols, name = "Subtype") +
    facet_grid(~ !!sym(GROUP_COL), scales = "free_x", space = "free_x") +
    labs(x = "Sample", y = "Proportion") +
    theme_classic(base_size = 14) +
    theme(axis.text.x      = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y      = element_text(size = 13),
          axis.title       = element_text(size = 15),
          legend.text      = element_text(size = 13),
          legend.title     = element_text(size = 14),
          strip.text       = element_text(size = 12, face = "bold"),
          strip.background = element_rect(fill = "grey92", color = NA),
          panel.spacing    = unit(0.4, "lines"))

  ggsave(file.path(ENDPOINT_DIR, "CellType_proportion_by_sample.pdf"), p_sample,
         width = 14, height = 5.5, device = "pdf", useDingbats = FALSE)
  ggsave(file.path(ENDPOINT_DIR, "CellType_proportion_by_sample.png"), p_sample,
         width = 14, height = 5.5, dpi = 200)
}
```

---

### Step 2d: Curated Marker Dot Plot

User provides a curated gene list grouped by subtype. Use the order exactly as provided.

```r
# Build gene vector in user-specified order (grouped by subtype)
genes_ordered <- c(
  # Subtype A-specific
  "project_specific",   # REPLACE with actual curated gene list
  # Shared / pan-markers
  # Subtype B-specific
  # Subtype C-specific
)

genes_use <- unique(genes_ordered[genes_ordered %in% rownames(obj)])
cat("Missing genes:", paste(setdiff(genes_ordered, rownames(obj)), collapse = ", "), "\n")

Idents(obj) <- "subtype_label"

p_dot <- DotPlot(obj, features = genes_use, group.by = "subtype_label", dot.scale = 7) +
  RotatedAxis() +
  scale_color_gradient2(low = "#2166AC", mid = "lightyellow", high = "#B2182B",
                        midpoint = 0, name = "Avg. Expr.") +
  labs(x = NULL, y = NULL, size = "% Expressed") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1, face = "italic"),
    axis.text.y      = element_text(size = 13),
    legend.text      = element_text(size = 11),
    legend.title     = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave(file.path(ENDPOINT_DIR, "CellType_dotplot_curated.pdf"), p_dot,
       width = max(8, length(genes_use) * 0.45 + 3), height = 4.5,
       device = "pdf", useDingbats = FALSE)
ggsave(file.path(ENDPOINT_DIR, "CellType_dotplot_curated.png"), p_dot,
       width = max(8, length(genes_use) * 0.45 + 3), height = 4.5, dpi = 200)
```

---

## Output File Summary

### Phase 1 — Exploration (`output/subclustering/`)

| File | Description |
|---|---|
| `CellType_umap_overview.pdf` | 4-panel: clusters / original labels / group / batch |
| `CellType_umap_clusters.pdf` | Single cluster UMAP |
| `CellType_feature_grid.pdf` | Canonical marker FeaturePlots with cluster number overlaid |
| `CellType_marker_dotplot.pdf` | DotPlot of canonical markers × clusters |
| `CellType_group_proportion.pdf` | Proportion: clusters on x, group fill (if GROUP_COL set) |
| `CellType_proportion_by_sample_clusters.pdf` | Per-sample proportion faceted by group |
| `CellType_cluster_markers.csv` | Full FindAllMarkers output |
| `CellType_top5_markers.csv` | Top 5 per cluster |
| `CellType_DEG_top100_significant.csv` | Top 100 significant per cluster — for label decisions |

### Phase 2 — Endpoint (`output/subclustering/endpoint/`)

| File | Description |
|---|---|
| `CellType_labeled_umap.pdf/png` | UMAP colored by subtype label |
| `CellType_proportion_by_group.pdf/png` | Proportion: group on x, label fill |
| `CellType_proportion_by_label.pdf/png` | Proportion (inverse): label on x, group fill |
| `CellType_proportion_by_sample.pdf/png` | Per-sample proportion faceted by group, label fill |
| `CellType_dotplot_curated.pdf/png` | Curated marker dot plot in user-specified gene/group order |
| `CellType_labeled.Rds` | Labeled object with `subtype_label` column |

---

## Iteration Patterns

**Removing a cluster after initial clustering** (common — user sees one contaminating cluster):
- Load `output/subclustering/CellType_subclustered.Rds`
- Remove cluster via `obj[, as.character(obj$seurat_clusters) != "X"]`
- Re-run from `NormalizeData` forward — do NOT go back to the whole object

**Two-stage subclustering** (e.g. MonoMac → Mac):
- First pass: broad subset (Macrophages + Monocytes together) at moderate resolution
- Identify macrophage-enriched clusters by group pattern + markers
- Second pass: subset those clusters only, re-run pipeline at same or slightly higher resolution
- This avoids monocyte contamination distorting the macrophage embedding

**Metadata assignment safety:** see `@primitives/seurat_v5_rules.md` Rule 3 —
always use `so@meta.data$col <-` not `so$col <-`.

---

## Project-Specific Values (Stage for Phase 4 examples/)

The following values are project-specific and must not appear in this module:

- `examples/humanfat_mac_subclustering.md` — CELLTYPE_LABELS for Macrophage/Monocyte subset,
  adipose_type_colors palette, source_file as BATCH_CORRECTION_VAR and SAMPLE_COL,
  tissue_type as GROUP_COL, validated 2026-03-24 cluster → label mapping
