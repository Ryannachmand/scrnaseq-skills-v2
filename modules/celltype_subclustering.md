---
requires_context:
  palettes:
    - group_colors     # optional -- named vector keyed by group_col values (tissue/condition)
    - subtype_colors   # optional -- named vector for final subtype labels (endpoint)
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
      - downstream_analyses.subclustering.mode               # "auto" (default) | "user_handoff"
      - downstream_analyses.subclustering.require_user_labels # bool, default false
      - downstream_analyses.subclustering.broad_type_markers  # list of {broad_type: [genes]}
references:
  - "@primitives/aesthetics.md"
  - "@primitives/doublet_removal.md"
  - "@primitives/harmony_integration.md"
  - "@primitives/visualization.md"
  - "@primitives/seurat_v5_rules.md"
  - "@context/color_palettes.md"
---

# Module: Cell Type Subclustering, Sub-labeling, and Endpoint Plots

Takes a broad cell type from a pre-annotated whole object, subsets and re-clusters it,
explores cluster identity through diagnostic plots and DEG analysis, then applies
labels to produce a clean set of endpoint publication figures.

Supports two execution paths (see Execution modes below):
- **Mode A (default):** autonomous annotation via marker scoring, no blocking checkpoint.
- **Mode B (opt-in):** user-handoff checkpoint for explicit label control.

---

## Execution modes

**Mode A (default -- autonomous platform):**
Subcluster -> diagnostic plots -> FindAllMarkers -> autonomous annotation
(AverageExpression scoring + `which.max()` assignment) -> duplicate-label guard
(append numeric suffixes when multiple Seurat clusters receive the same broad-type label)
-> labeled endpoint outputs. NO blocking checkpoint.

Enable Mode A by omitting `subclustering.require_user_labels` from the brief, or setting
it to `false`. This is the default for all autonomous pipeline runs.

**Mode B (opt-in -- user-handoff):**
Subcluster -> diagnostic plots -> FindAllMarkers -> BLOCKING checkpoint awaiting
user-provided cluster-to-label mapping -> label application -> labeled endpoint outputs.

Enable Mode B by setting `downstream_analyses.subclustering.require_user_labels: true`
in the brief. Use when autonomous annotation quality is insufficient or when the user
wants explicit label control.

---

## Critical Constraints

| Do NOT | Do instead | Why |
|---|---|---|
| `subset(obj, seurat_clusters %in% c(0,1,2))` | `obj[, colnames(obj)[as.character(obj$seurat_clusters) %in% as.character(c(0,1,2))]]` | `seurat_clusters` is a factor; `subset()` throws "No cell overlap" with integer comparisons |
| `obj$label <- label_map[as.character(obj$seurat_clusters)]` | `obj$label <- unname(label_map[as.character(obj$seurat_clusters)])` | Without `unname()`, the vector carries cluster numbers as names; Seurat matches against cell barcodes and throws "No cell overlap" |
| `Embeddings(obj, "umap")[,"UMAP_1"]` | `as.data.frame(Embeddings(obj, "umap"))$umap_1` | UMAP column names are `umap_1`/`umap_2` (lowercase underscore) in r-env |
| Pass gene vector directly to `DotPlot` | `unique(genes_vec)` before passing | Duplicate gene names cause a factor level duplication error |
| Use base R | Always run via `@primitives/r_environment.md` | `sp` package is broken in base R |
| Use `pdf()` / `dev.off()` for ggplot outputs | `ggsave(..., device = "pdf", useDingbats = FALSE)` | Convention -- see CONVENTIONS.md section 4 |
| `DimPlot(label = TRUE)` | `make_labeled_umap()` from `@primitives/visualization.md` | visualization.md mandates `LabelClusters()`, not `label = TRUE` |
| `DotPlot() + RotatedAxis()` bare | `make_canonical_dotplot()` from `@primitives/visualization.md` | Ensures canonical expression palette; `scale_color_gradient2(lightyellow)` is wrong |
| Inline proportion ggplot | `make_proportion_plot()` from `@primitives/visualization.md` | Ensures numeric x-positions and consistent separator logic |

---

## Step 1: Subset and Re-cluster

```r
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(patchwork)
set.seed(42)

# -- CONFIG -------------------------------------------------------------------
# Target cell types to subset from the whole object
CELLTYPE_LABELS <- c("project_specific")  # REPLACE: values from whole-object label_col
                    # Example: c("Macrophage", "Monocyte") -- use blank template
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
OUTPUT_DIR      <- "output/celltype_subclustering"
ENDPOINT_DIR    <- "output/celltype_subclustering/endpoint"

N_PCS           <- 50
NEIGHBOR_DIMS   <- 1:25
CLUSTER_RES     <- 0.4

# Color vectors -- provide from brief context_overrides.palettes or define here.
# GROUP_COLORS must be defined before Step 2 (used in Phase 1 diagnostic proportion plots).
GROUP_COLORS    <- c("project_specific" = "#project_specific")  # REPLACE from brief
# SUBTYPE_COLORS is needed for endpoint plots -- fill in when labels are known.
SUBTYPE_COLORS  <- c("project_specific" = "#project_specific")  # REPLACE after labeling

# -- LOAD AND SUBSET ----------------------------------------------------------
dir.create(OUTPUT_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(ENDPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

whole <- readRDS(RDS_IN)
obj   <- subset(whole, subset = !!sym(LABEL_COL) %in% CELLTYPE_LABELS)
rm(whole); gc()

# -- DOUBLET REMOVAL ON SUBSET ------------------------------------------------
# scDblFinder requires raw counts; must run before normalization.
# Reuses run_doublet_removal() from @primitives/doublet_removal.md.
obj <- run_doublet_removal(obj, sample_col = SAMPLE_COL,
                           output_dir = file.path(OUTPUT_DIR, "doublets/"))

# -- RE-CLUSTER ON SUBSET -----------------------------------------------------
# After subsetting, obj retains parent-object pca/harmony/umap embeddings for
# these cells. ALL steps below OVERWRITE them with embeddings computed on just
# this subset. Do NOT skip or reorder any step -- parent geometry will persist.
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj, features = VariableFeatures(obj))
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
# obj2 inherits pca/harmony/umap from obj. Re-run the FULL embedding pipeline:
# NormalizeData -> FindVariableFeatures -> ScaleData(features=VariableFeatures) ->
# RunPCA -> RunHarmony -> FindNeighbors -> FindClusters -> RunUMAP
```

**If a contaminating cluster must be removed:**
```r
obj3 <- obj2[, as.character(obj2$seurat_clusters) != "1"]
# obj3 inherits embeddings from obj2. Re-run the full embedding pipeline (same steps as above).
```

---

## Step 2: Diagnostic Plots (Both Modes)

Re-read `@primitives/aesthetics.md` and `@primitives/visualization.md` before writing
any plot code -- mandatory per pipeline Universal Rules.

All plots use `ggsave()` -- never `pdf()` / `dev.off()` for ggplot outputs.

**4-panel UMAP overview (patchwork):**
```r
# make_labeled_umap() uses LabelClusters(), NOT label = TRUE (see Critical Constraints)
p1 <- make_labeled_umap(obj, label_col = "seurat_clusters")
p2 <- make_labeled_umap(obj, label_col = LABEL_COL)

if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  p3 <- make_labeled_umap(obj, label_col = GROUP_COL)
  p4 <- make_labeled_umap(obj, label_col = SAMPLE_COL)
  panel4 <- (p1 | p2) / (p3 | p4)
} else {
  p3 <- make_labeled_umap(obj, label_col = SAMPLE_COL)
  panel4 <- (p1 | p2) / (p3 | plot_spacer())
}

ggsave(file.path(OUTPUT_DIR, "CellType_umap_overview.pdf"),
       panel4, width = 16, height = 12, device = "pdf", useDingbats = FALSE)
```

**Cluster UMAP (individual):**
```r
ggsave(file.path(OUTPUT_DIR, "CellType_umap_clusters.pdf"),
       make_labeled_umap(obj, label_col = "seurat_clusters"),
       width = 8, height = 7, device = "pdf", useDingbats = FALSE)
```

**Canonical marker feature plots with cluster number overlaid:**
```r
# aesthetics.md spec: c("lightgrey","blue") scale, explicit limits anchored at zero,
# raster = FALSE (rasterized points render invisible at small pt.size)
umap_df   <- as.data.frame(Embeddings(obj, "umap"))   # columns: umap_1, umap_2
umap_df$cluster <- as.character(obj$seurat_clusters)
centroids <- umap_df %>%
  group_by(cluster) %>%
  summarise(x = median(umap_1), y = median(umap_2), .groups = "drop")

marker_genes <- unique(c("project_specific"))   # REPLACE: canonical markers for this cell type

plots <- lapply(marker_genes, function(g) {
  g_max <- max(GetAssayData(obj, assay = "RNA", layer = "data")[g, ])
  FeaturePlot(obj, features = g, reduction = "umap", order = TRUE,
              raster = FALSE, cols = c("lightgrey", "blue")) +
    scale_color_gradient(low = "lightgrey", high = "blue",
                         limits = c(0, g_max)) +
    geom_text(data = centroids, aes(x = x, y = y, label = cluster),
              inherit.aes = FALSE, size = 3, fontface = "bold", color = "black")
})

ggsave(file.path(OUTPUT_DIR, "CellType_feature_grid.pdf"),
       wrap_plots(plots, ncol = 4),
       width = 16, height = 12, device = "pdf", useDingbats = FALSE)
```

**Canonical marker dot plot:**
```r
# Extract dot data from Seurat, then render via make_canonical_dotplot()
genes_valid <- unique(marker_genes)[unique(marker_genes) %in% rownames(obj)]
Idents(obj) <- "seurat_clusters"
dp_data <- DotPlot(obj, features = genes_valid, group.by = "seurat_clusters")$data
dot_df <- dp_data %>%
  rename(gene = features.plot, group = id, pct_exp = pct.exp, avg_exp_scaled = avg.exp.scaled) %>%
  mutate(gene = factor(gene, levels = genes_valid))

make_canonical_dotplot(
  dot_df      = dot_df,
  gene_col    = "gene",
  group_col   = "group",
  pct_col     = "pct_exp",
  avg_col     = "avg_exp_scaled",
  output_file = file.path(OUTPUT_DIR, "CellType_marker_dotplot.pdf"),
  n_genes     = length(genes_valid),
  n_groups    = length(unique(dot_df$group))
)
```

**Group proportion -- clusters on x (only when GROUP_COL is set):**
```r
if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  prop_df <- obj@meta.data %>%
    group_by(seurat_clusters, !!sym(GROUP_COL)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(proportion = n / sum(n),
           seurat_clusters = as.character(seurat_clusters))

  make_proportion_plot(
    df          = prop_df,
    sample_col  = "seurat_clusters",
    group_col   = "seurat_clusters",
    subtype_col = GROUP_COL,
    group_colors = GROUP_COLORS,
    output_file = file.path(OUTPUT_DIR, "CellType_group_proportion.pdf")
  )
}
```

**Per-sample proportion (faceted by group):**
```r
# Paul Tol muted palette -- safe for red-green color blindness, up to 10 clusters
tol_muted <- c("#332288","#88CCEE","#44AA99","#117733","#999933",
               "#DDCC77","#CC6677","#882255","#AA4499","#DDDDDD")
cluster_levels <- as.character(sort(unique(as.integer(as.character(obj$seurat_clusters)))))
cluster_colors <- setNames(tol_muted[seq_len(length(cluster_levels))], cluster_levels)

if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  prop_sample <- obj@meta.data %>%
    group_by(!!sym(SAMPLE_COL), !!sym(GROUP_COL), seurat_clusters) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(!!sym(SAMPLE_COL)) %>%
    mutate(proportion = n / sum(n),
           seurat_clusters = factor(as.character(seurat_clusters), levels = cluster_levels))

  make_proportion_plot(
    df          = prop_sample,
    sample_col  = SAMPLE_COL,
    group_col   = GROUP_COL,
    subtype_col = "seurat_clusters",
    group_colors = cluster_colors,
    output_file = file.path(OUTPUT_DIR, "CellType_proportion_by_sample_clusters.pdf")
  )
}
```

---

## Step 3: FindAllMarkers (Both Modes)

```r
# See @primitives/seurat_v5_rules.md Rule 1 -- JoinLayers required before FindAllMarkers
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

## Mode A: Autonomous annotation (default)

Run this section when `downstream_analyses.subclustering.require_user_labels` is absent
or `false`. No blocking checkpoint. Labels are assigned programmatically from marker
gene scoring.

### Autonomous annotation algorithm

```r
# -- BROAD TYPE MARKERS -------------------------------------------------------
# broad_type_markers: named list of gene vectors (from brief or defined here)
# Each name is a broad cell subtype; each vector is its canonical marker genes.
# Example:
#   broad_type_markers <- list(
#     "Capillary"    = c("CA4", "PLVAP", "FCN3"),
#     "ArterialEC"   = c("HEY1", "GJA5", "EFNB2"),
#     "VenousEC"     = c("ACKR1", "NR2F2", "VWF"),
#     "LymphaticEC"  = c("LYVE1", "PROX1", "PDPN")
#   )
broad_type_markers <- list(
  "project_specific" = c("project_specific")  # REPLACE from brief
)

# -- SCORE EACH CLUSTER -------------------------------------------------------
all_marker_genes <- unique(unlist(broad_type_markers))
avg_expr <- AverageExpression(obj, assay = "RNA",
                               features = all_marker_genes[all_marker_genes %in% rownames(obj)],
                               group.by = "seurat_clusters")$RNA

score_mat <- do.call(rbind, lapply(names(broad_type_markers), function(bt) {
  genes <- intersect(broad_type_markers[[bt]], rownames(avg_expr))
  if (length(genes) == 0) return(rep(0, ncol(avg_expr)))
  colMeans(avg_expr[genes, , drop = FALSE])
}))
rownames(score_mat) <- names(broad_type_markers)

# -- ASSIGN BROAD TYPE --------------------------------------------------------
cluster_ids    <- colnames(avg_expr)
assigned_types <- rownames(score_mat)[apply(score_mat, 2, which.max)]
names(assigned_types) <- cluster_ids

# -- DUPLICATE-LABEL GUARD ---------------------------------------------------
# When multiple clusters receive the same broad-type label, append numeric suffixes
# in order of cluster ID. Single-cluster labels stay unsuffixed.
# Example: clusters 1, 2, 5 all assigned "Capillary" -> "Capillary_1", "Capillary_2", "Capillary_3"
label_counts <- table(assigned_types)
final_labels <- assigned_types
for (lbl in names(label_counts)[label_counts > 1]) {
  dup_clusters <- sort(names(assigned_types)[assigned_types == lbl])
  for (i in seq_along(dup_clusters)) {
    final_labels[dup_clusters[i]] <- paste0(lbl, "_", i)
  }
}

# -- WRITE ASSIGNMENT TABLE ---------------------------------------------------
score_for_cluster <- sapply(cluster_ids, function(cl) max(score_mat[, cl]))
assignment_df <- data.frame(
  seurat_cluster   = cluster_ids,
  broad_type       = assigned_types[cluster_ids],
  score            = score_for_cluster,
  final_label      = final_labels[cluster_ids],
  received_suffix  = label_counts[assigned_types[cluster_ids]] > 1,
  stringsAsFactors = FALSE
)
write.csv(assignment_df,
          file.path(OUTPUT_DIR, "autonomous_label_assignment.csv"),
          row.names = FALSE)
message(sprintf("Autonomous annotation: %d clusters -> %d final labels",
                length(cluster_ids), length(unique(final_labels))))

# -- APPLY LABELS AND SAVE ---------------------------------------------------
dir.create(ENDPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
label_map_auto <- setNames(final_labels[cluster_ids], cluster_ids)
obj$subtype_label <- unname(label_map_auto[as.character(obj$seurat_clusters)])
obj$subtype_label <- factor(obj$subtype_label, levels = unique(final_labels))
Idents(obj) <- "subtype_label"
saveRDS(obj, file.path(ENDPOINT_DIR, "CellType_labeled.Rds"))
```

### Mode A endpoint plots

```r
# Auto-generate subtype color palette (Tol bright 12-color, colorblind-safe)
tol_bright <- c("#4477AA","#EE6677","#228833","#CCBB44","#66CCEE","#AA3377",
                "#BBBBBB","#000000","#E69F00","#56B4E9","#009E73","#F0E442")
subtype_levels_auto <- levels(obj$subtype_label)
subtype_cols_auto <- setNames(
  tol_bright[seq_len(length(subtype_levels_auto))],
  subtype_levels_auto
)

# Labeled UMAP
p_umap_auto <- make_labeled_umap(obj, label_col = "subtype_label") +
  scale_color_manual(values = subtype_cols_auto)
ggsave(file.path(ENDPOINT_DIR, "CellType_labeled_umap.pdf"), p_umap_auto,
       width = 9, height = 5.5, device = "pdf", useDingbats = FALSE)
ggsave(file.path(ENDPOINT_DIR, "CellType_labeled_umap.png"), p_umap_auto,
       width = 9, height = 5.5, dpi = 200)

# Proportion plots (only when GROUP_COL is set)
if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  # By group: x = GROUP_COL, fill = subtype
  prop_by_group_a <- obj@meta.data %>%
    group_by(!!sym(GROUP_COL), subtype_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(!!sym(GROUP_COL)) %>%
    mutate(proportion = n / sum(n),
           subtype_label = factor(subtype_label, levels = subtype_levels_auto))

  p_grp_a <- make_proportion_plot(
    df = prop_by_group_a, sample_col = GROUP_COL, group_col = GROUP_COL,
    subtype_col = "subtype_label", group_colors = subtype_cols_auto,
    output_file = file.path(ENDPOINT_DIR, "CellType_proportion_by_group.pdf")
  )
  ggsave(file.path(ENDPOINT_DIR, "CellType_proportion_by_group.png"),
         plot = p_grp_a, width = 9, height = 5.5, dpi = 200)

  # By sample: x = SAMPLE_COL, groups = GROUP_COL, fill = subtype
  prop_sample_a <- obj@meta.data %>%
    group_by(!!sym(SAMPLE_COL), !!sym(GROUP_COL), subtype_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(!!sym(SAMPLE_COL)) %>%
    mutate(proportion = n / sum(n),
           subtype_label = factor(subtype_label, levels = subtype_levels_auto))

  p_samp_a <- make_proportion_plot(
    df = prop_sample_a, sample_col = SAMPLE_COL, group_col = GROUP_COL,
    subtype_col = "subtype_label", group_colors = subtype_cols_auto,
    output_file = file.path(ENDPOINT_DIR, "CellType_proportion_by_sample.pdf")
  )
  ggsave(file.path(ENDPOINT_DIR, "CellType_proportion_by_sample.png"),
         plot = p_samp_a, width = 14, height = 5.5, dpi = 200)
}

# Marker dot plot using broad_type_markers genes
genes_auto_valid <- unique(unlist(broad_type_markers))
genes_auto_valid <- genes_auto_valid[genes_auto_valid %in% rownames(obj)]
dp_data_auto <- DotPlot(obj, features = genes_auto_valid, group.by = "subtype_label")$data
dot_df_auto <- dp_data_auto %>%
  rename(gene = features.plot, group = id, pct_exp = pct.exp, avg_exp_scaled = avg.exp.scaled) %>%
  mutate(gene = factor(gene, levels = genes_auto_valid))

make_canonical_dotplot(
  dot_df      = dot_df_auto,
  gene_col    = "gene",
  group_col   = "group",
  pct_col     = "pct_exp",
  avg_col     = "avg_exp_scaled",
  group_colors = subtype_cols_auto,
  output_file = file.path(ENDPOINT_DIR, "CellType_dotplot_markers.pdf"),
  n_genes     = length(genes_auto_valid),
  n_groups    = length(unique(dot_df_auto$group))
)
```

---

## Mode B: User-handoff (opt-in)

Activated when `downstream_analyses.subclustering.require_user_labels: true` in the brief.
Steps 2 and 3 (diagnostic plots and FindAllMarkers) have already run. This section
handles the blocking checkpoint, label application, and endpoint plots.

### Step B1: Handoff to User

After generating exploration plots and FindAllMarkers output:
1. Path to cluster UMAP PDF
2. Path to feature grid PDF(s)
3. Path to `CellType_DEG_top100_significant.csv`
4. Summary table of cluster sizes and dominant group per cluster

Request from user:
- Cluster -> label mapping (many-to-one: multiple clusters -> same label)
- Which clusters to remove entirely
- (Optional) Which group values to exclude from endpoint plots

Write a checkpoint JSON to `output/celltype_subclustering/checkpoint_mode_b.json` with the
paths above, then exit. The pipeline resumes when the user provides the mapping.

### Step B2: Apply Labels and Remove Unwanted Clusters

```r
# -- USER-PROVIDED MAPPING ----------------------------------------------------
# Fill in cluster -> subtype label mapping (blank template below)
label_map <- c(
  "0" = "project_specific",   # REPLACE with actual subtype names
  "1" = "project_specific",
  "2" = "project_specific"
  # Map multiple clusters to the same label when biology warrants it
  # Omit clusters to remove them (they will not appear in keep_clusters)
)
keep_clusters <- names(label_map)

# -- OPTIONAL: exclude specific group values ----------------------------------
# Set to NULL to keep all; or c("GroupA", "GroupB") to remove
EXCLUDE_GROUPS <- NULL  # REPLACE from brief: downstream_analyses.subclustering.exclude_groups

# -- FILTER CELLS -------------------------------------------------------------
keep_cells <- colnames(obj)[as.character(obj$seurat_clusters) %in% keep_clusters]

if (!is.null(EXCLUDE_GROUPS) && !is.null(GROUP_COL) &&
    GROUP_COL %in% colnames(obj@meta.data)) {
  keep_cells <- intersect(
    keep_cells,
    colnames(obj)[!obj@meta.data[[GROUP_COL]] %in% EXCLUDE_GROUPS]
  )
}
obj <- obj[, keep_cells]

# -- ASSIGN LABELS ------------------------------------------------------------
# MUST use unname() -- see Critical Constraints above
obj$subtype_label <- unname(label_map[as.character(obj$seurat_clusters)])
obj$subtype_label <- factor(obj$subtype_label, levels = unique(label_map))

Idents(obj) <- "subtype_label"
saveRDS(obj, file.path(ENDPOINT_DIR, "CellType_labeled.Rds"))
```

### Step B3: Labeled UMAP

```r
# -- SUBTYPE COLORS -----------------------------------------------------------
# Provide from brief context_overrides.palettes.subtype_colors, or define here
subtype_cols <- c(
  "project_specific" = "#project_specific"  # REPLACE: one entry per subtype label
)
# Example structure: c("SubtypeA" = "#E05C5C", "SubtypeB" = "#4A90D9")

p <- make_labeled_umap(obj, label_col = "subtype_label", pt_size = 0.6) +
  scale_color_manual(values = subtype_cols)

ggsave(file.path(ENDPOINT_DIR, "CellType_labeled_umap.pdf"), p,
       width = 9, height = 5.5, device = "pdf", useDingbats = FALSE)
ggsave(file.path(ENDPOINT_DIR, "CellType_labeled_umap.png"), p,
       width = 9, height = 5.5, dpi = 200)
```

### Step B4: Proportion Plots (Group on x, Labels as fill)

```r
# -- GROUP COLORS -------------------------------------------------------------
# Provide from brief context_overrides.palettes.group_colors; overrides GROUP_COLORS if needed
group_cols <- c("project_specific" = "#project_specific")  # REPLACE

if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  # By group: x = GROUP_COL, fill = subtype label
  prop_by_group <- obj@meta.data %>%
    group_by(!!sym(GROUP_COL), subtype_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(!!sym(GROUP_COL)) %>%
    mutate(proportion = n / sum(n),
           subtype_label = factor(subtype_label, levels = names(subtype_cols)))

  p_group <- make_proportion_plot(
    df          = prop_by_group,
    sample_col  = GROUP_COL,
    group_col   = GROUP_COL,
    subtype_col = "subtype_label",
    group_colors = subtype_cols,
    output_file = file.path(ENDPOINT_DIR, "CellType_proportion_by_group.pdf")
  )
  ggsave(file.path(ENDPOINT_DIR, "CellType_proportion_by_group.png"),
         plot = p_group, width = 9, height = 5.5, dpi = 200)

  # Inverse: x = subtype label, fill = group
  prop_by_label <- obj@meta.data %>%
    group_by(subtype_label, !!sym(GROUP_COL)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(subtype_label) %>%
    mutate(proportion = n / sum(n))

  p_label <- make_proportion_plot(
    df          = prop_by_label,
    sample_col  = "subtype_label",
    group_col   = "subtype_label",
    subtype_col = GROUP_COL,
    group_colors = group_cols,
    output_file = file.path(ENDPOINT_DIR, "CellType_proportion_by_label.pdf")
  )
  ggsave(file.path(ENDPOINT_DIR, "CellType_proportion_by_label.png"),
         plot = p_label, width = 9, height = 5.5, dpi = 200)

  # Per sample: x = SAMPLE_COL, grouped by GROUP_COL, fill = subtype label
  prop_sample <- obj@meta.data %>%
    group_by(!!sym(SAMPLE_COL), !!sym(GROUP_COL), subtype_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(!!sym(SAMPLE_COL)) %>%
    mutate(proportion = n / sum(n),
           subtype_label = factor(subtype_label, levels = names(subtype_cols)))

  p_sample <- make_proportion_plot(
    df          = prop_sample,
    sample_col  = SAMPLE_COL,
    group_col   = GROUP_COL,
    subtype_col = "subtype_label",
    group_colors = subtype_cols,
    output_file = file.path(ENDPOINT_DIR, "CellType_proportion_by_sample.pdf")
  )
  ggsave(file.path(ENDPOINT_DIR, "CellType_proportion_by_sample.png"),
         plot = p_sample, width = 14, height = 5.5, dpi = 200)
}
```

### Step B5: Curated Marker Dot Plot

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

# Extract dot data from Seurat, render via make_canonical_dotplot()
dp_data_b <- DotPlot(obj, features = genes_use, group.by = "subtype_label")$data
dot_df_b <- dp_data_b %>%
  rename(gene = features.plot, group = id, pct_exp = pct.exp, avg_exp_scaled = avg.exp.scaled) %>%
  mutate(gene = factor(gene, levels = genes_use))

p_dot <- make_canonical_dotplot(
  dot_df      = dot_df_b,
  gene_col    = "gene",
  group_col   = "group",
  pct_col     = "pct_exp",
  avg_col     = "avg_exp_scaled",
  group_colors = subtype_cols,
  output_file = file.path(ENDPOINT_DIR, "CellType_dotplot_curated.pdf"),
  n_genes     = length(genes_use),
  n_groups    = length(unique(dot_df_b$group))
)
ggsave(file.path(ENDPOINT_DIR, "CellType_dotplot_curated.png"),
       plot = p_dot,
       width = max(8, length(genes_use) * 0.45 + 3), height = 4.5, dpi = 200)
```

---

## Output File Summary

### Shared (both modes) -- Exploration (`output/celltype_subclustering/`)

| File | Description |
|---|---|
| `CellType_umap_overview.pdf` | 4-panel: clusters / original labels / group / batch |
| `CellType_umap_clusters.pdf` | Single cluster UMAP |
| `CellType_feature_grid.pdf` | Canonical marker FeaturePlots with cluster number overlaid |
| `CellType_marker_dotplot.pdf` | Dot plot of canonical markers x clusters (make_canonical_dotplot) |
| `CellType_group_proportion.pdf` | Proportion: clusters on x, group fill (if GROUP_COL set) |
| `CellType_proportion_by_sample_clusters.pdf` | Per-sample proportion faceted by group |
| `CellType_cluster_markers.csv` | Full FindAllMarkers output |
| `CellType_top5_markers.csv` | Top 5 per cluster |
| `CellType_DEG_top100_significant.csv` | Top 100 significant per cluster -- for label decisions |

### Mode A only -- Exploration (`output/celltype_subclustering/`)

| File | Description |
|---|---|
| `autonomous_label_assignment.csv` | Table: seurat_cluster, broad_type, score, final_label, received_suffix |

### Endpoint -- both modes (`output/celltype_subclustering/endpoint/`)

| File | Description |
|---|---|
| `CellType_labeled_umap.pdf/png` | UMAP colored by subtype label |
| `CellType_proportion_by_group.pdf/png` | Proportion: group on x, label fill |
| `CellType_proportion_by_sample.pdf/png` | Per-sample proportion, label fill |
| `CellType_labeled.Rds` | Labeled object with `subtype_label` column |
| `CellType_dotplot_markers.pdf` | Mode A: dot plot of broad_type_markers genes |
| `CellType_dotplot_curated.pdf/png` | Mode B: curated marker dot plot in user-specified order |
| `CellType_proportion_by_label.pdf/png` | Mode B: inverse proportion (label on x, group fill) |

---

## Iteration Patterns

**Removing a cluster after initial clustering** (common -- user sees one contaminating cluster):
- Load `output/celltype_subclustering/CellType_subclustered.Rds`
- Remove cluster via `obj[, as.character(obj$seurat_clusters) != "X"]`
- Re-run from `NormalizeData` forward -- do NOT go back to the whole object

**Two-stage subclustering** (e.g. MonoMac -> Mac):
- First pass: broad subset (Macrophages + Monocytes together) at moderate resolution
- Identify macrophage-enriched clusters by group pattern + markers
- Second pass: subset those clusters only, re-run pipeline at same or slightly higher resolution
- This avoids monocyte contamination distorting the macrophage embedding

**Metadata assignment safety:** see `@primitives/seurat_v5_rules.md` Rule 3 --
always use `so@meta.data$col <-` not `so$col <-`.

---

## Project-Specific Values (Stage for Phase 4 examples/)

The following values are project-specific and must not appear in this module:

- `examples/humanfat_mac_subclustering.md` -- CELLTYPE_LABELS for Macrophage/Monocyte subset,
  adipose_type_colors palette, source_file as BATCH_CORRECTION_VAR and SAMPLE_COL,
  tissue_type as GROUP_COL, validated 2026-03-24 cluster -> label mapping
