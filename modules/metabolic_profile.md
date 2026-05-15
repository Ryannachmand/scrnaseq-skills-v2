---
requires_context:
  palettes:
    - subtype_colors   # optional — named vector for cell subtype colors in faceted plots
    - group_colors     # optional — named vector for group/tissue colors in proportion context
  metadata_columns:
    required:
      - label_col      # broad cell type label column (for whole-object faceting)
    optional:
      - subtype_col    # finer subtype label column (for EC/subset-level plots)
      - group_col      # tissue/condition grouping column (for group-level plots)
  brief_keys:
    required:
      - output_dir
      - downstream_analyses.metabolic_profile.gene_sets
    optional:
      - downstream_analyses.metabolic_profile.label_col
      - downstream_analyses.metabolic_profile.subtype_col
      - downstream_analyses.metabolic_profile.group_col
      - downstream_analyses.metabolic_profile.subsample_ceiling
references:
  - "@primitives/aucell_scoring.md"
  - "@primitives/seurat_v5_rules.md"
  - "@context/color_palettes.md"
---

# Module: Single-Cell Metabolic Profiling (AUCell-based)

Score every cell for metabolic pathway activity using AUCell (rank-based enrichment),
then visualize scores across cell types, subtypes, and tissue/condition groups. Applies
to any pre-processed Seurat object — no re-QC or re-clustering.

The module shows HOW to score AUCell against caller-provided gene sets.
The gene sets themselves (which specific genes belong to Glycolysis, OXPHOS, etc.)
are project-specific curation choices and belong in `examples/`. See the brief schema below.

**Note on palette:** This module uses `PALETTE_DIVERGING_6` (blue→white→red, 6-stop) from
`@context/color_palettes.md` for AUCell UMAPs. The v1 source used a distinct 3-stop
blue-yellow-red palette (`#4575B4`, `#FFFFBF`, `#D73027`). The 6-stop diverging palette
is the lab canonical choice (duplication_report.md §3.1). If a 3-stop diverging is needed
for visual simplicity, use `PALETTE_DIVERGING_3` from `@context/color_palettes.md`.

---

## Brief Schema

```yaml
downstream_analyses:
  metabolic_profile:
    enabled: true
    label_col: project_specific      # broad cell type label column (whole object)
    subtype_col: project_specific    # finer subtype label column (subset object, if available)
    group_col: project_specific      # tissue/condition column for group-level plots
    subsample_ceiling: 2500          # max cells per label group before AUCell; NULL = no subsampling
    gene_sets: project_specific      # REQUIRED: named list of gene vectors (see note below)
    # NOTE: gene_sets must be defined explicitly in your project CLAUDE.md or brief.
    # Gene set membership reflects curation choices — see examples/ for validated gene sets.
    # Example structure:
    #   gene_sets:
    #     Glycolysis:     ["HK1","HK2","PFKL","ALDOA","GAPDH","PKM","LDHA"]
    #     Beta_Oxidation: ["ACADM","HADHA","CPT1A","ACSL1"]
    #     TCA_Cycle:      ["CS","IDH1","OGDH","SDHA","MDH2"]
    #     OXPHOS:         ["NDUFA1","NDUFB1","SDHB","UQCRB","COX5A","ATP5B"]
    #     FA_Synthesis:   ["FASN","ACACA","ACLY","SCD","FADS1","DGAT1"]
    #   # Add additional gene sets as needed (e.g. PPARG_Targets for a TF OE project)
```

---

## Critical Constraints

| ❌ Don't | ✅ Do instead | Why |
|---|---|---|
| Define gene sets in this module | Pass as `gene_sets` argument from brief | Gene set membership is curation-choice, not universal biology |
| Use invalid range syntax `"NDUFA1"-"NDUFA10"` | Enumerate genes explicitly as character vector | R range operator does not work on strings — will throw an error |
| Apply project-specific sample exclusion filters | Accept the object as-is | Exclusions (e.g. tissue type filters) belong in the project CLAUDE.md |
| Inline `run_aucell` or `add_auc_to_seurat` | Use functions from `@primitives/aucell_scoring.md` | Canonical implementations live in the primitive |
| Mix `scale_color_manual` + `scale_color_identity` | Pre-compute hex strings into the data frame | Two scale calls silently conflict; the second replaces the first |

---

## R Packages Required

```r
library(Seurat)
library(AUCell)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(Matrix)
library(ggplotify)
library(cowplot)
```

---

## Configuration Block

```r
# ── CONFIG ───────────────────────────────────────────────────────────────────
LABEL_COL    <- "project_specific"   # REPLACE: broad cell type label column
SUBTYPE_COL  <- "project_specific"   # REPLACE: finer subtype label column (or NULL)
GROUP_COL    <- "project_specific"   # REPLACE: tissue/condition column (or NULL)
RDS_IN       <- "project_specific"   # REPLACE: path to Seurat object
OUTPUT_DIR   <- "output/metabolic"
SUBSAMPLE_CEILING <- 2500            # max cells per label group before AUCell; NULL = no subsampling

# GENE SETS — load from brief or define explicitly in project CLAUDE.md
# See brief schema above for expected structure.
# DO NOT define gene sets here; they are caller-provided named lists.
gene_sets <- brief$downstream_analyses$metabolic_profile$gene_sets
# Structure: list(PathwayName = c("GENE1","GENE2",...), ...)
stopifnot(!is.null(gene_sets), is.list(gene_sets), length(gene_sets) > 0)

PATHWAY_ORDER <- names(gene_sets)   # controls row order in heatmaps

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
```

---

## Step 1: Load and Optionally Subsample

```r
obj <- readRDS(RDS_IN)

# Stratified subsampling before AUCell to avoid OOM on large objects
# (2500 cells/group gives near-identical group means, reduces memory ~5x)
if (!is.null(SUBSAMPLE_CEILING) && ncol(obj) > SUBSAMPLE_CEILING * length(unique(obj@meta.data[[LABEL_COL]]))) {
  meta <- obj@meta.data
  keep_cells <- meta %>%
    tibble::rownames_to_column("cell") %>%
    group_by(!!sym(LABEL_COL)) %>%
    slice_sample(n = SUBSAMPLE_CEILING) %>%
    pull(cell)
  obj_score <- obj[, keep_cells]
  message(sprintf("Subsampled to %d cells for AUCell scoring", ncol(obj_score)))
} else {
  obj_score <- obj
}
```

---

## Step 2: AUCell Scoring

Uses `run_aucell()` and `add_auc_to_seurat()` from `@primitives/aucell_scoring.md`.

```r
# Run AUCell — gene_sets is the caller-provided named list from brief
auc_df    <- run_aucell(obj_score, gene_sets = gene_sets, assay = "RNA")
obj_score <- add_auc_to_seurat(obj_score, auc_df)

# If subsampled: propagate scores back to full object via merge on cell barcodes
if (!identical(obj_score, obj)) {
  for (pw in names(gene_sets)) {
    vals <- setNames(obj_score@meta.data[[pw]], colnames(obj_score))
    obj  <- AddMetaData(obj, metadata = vals, col.name = pw)
  }
}

# Save scored object
saveRDS(obj, file.path(OUTPUT_DIR, "metabolic_scored.Rds"))

# Export per-cell scores to CSV
auc_export <- obj@meta.data[, c(LABEL_COL, names(gene_sets))]
write.csv(auc_export, file.path(OUTPUT_DIR, "auc_scores_per_cell.csv"), quote = FALSE)
```

---

## Step 3: Visualization

### 3.1 AUCell UMAP — one panel per pathway

```r
# PALETTE_DIVERGING_6 from @context/color_palettes.md
PALETTE_DIVERGING_6 <- c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#D6604D","#B2182B")

# Detect UMAP reduction name (see @primitives/seurat_v5_rules.md Rule 5)
umap_key <- names(obj@reductions)[grepl("umap", names(obj@reductions), ignore.case = TRUE)][1]
umap_df  <- as.data.frame(Embeddings(obj, umap_key))
colnames(umap_df) <- c("umap_1", "umap_2")

umap_plots <- lapply(PATHWAY_ORDER, function(pw) {
  score_col <- umap_df
  score_col$score <- obj@meta.data[[pw]]
  lims <- quantile(score_col$score, c(0.02, 0.98), na.rm = TRUE)
  ggplot(score_col, aes(x = umap_1, y = umap_2, colour = score)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_colour_gradientn(
      colours = PALETTE_DIVERGING_6,
      limits  = lims,
      oob     = scales::squish,
      name    = "AUC"
    ) +
    ggtitle(pw) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 18, hjust = 0.5))
})

ggsave(file.path(OUTPUT_DIR, "umap_metabolic_scores.pdf"),
       wrap_plots(umap_plots, ncol = min(length(PATHWAY_ORDER), 3)),
       width = 7 * min(length(PATHWAY_ORDER), 3),
       height = 6 * ceiling(length(PATHWAY_ORDER) / 3),
       device = "pdf", useDingbats = FALSE)
```

---

### 3.2 AUCell Violin by Cell Type Label

```r
vln_data <- obj@meta.data %>%
  select(all_of(c(LABEL_COL, PATHWAY_ORDER))) %>%
  tidyr::pivot_longer(cols = all_of(PATHWAY_ORDER),
                      names_to = "pathway", values_to = "score") %>%
  mutate(pathway = factor(pathway, levels = PATHWAY_ORDER))

p_vln <- ggplot(vln_data, aes(x = !!sym(LABEL_COL), y = score, fill = !!sym(LABEL_COL))) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.85, linewidth = 0.3, color = "grey30") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7, linewidth = 0.3) +
  facet_wrap(~ pathway, scales = "free_y", ncol = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(x = NULL, y = "AUC Score") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "violin_metabolic_by_celltype.pdf"),
       p_vln,
       width = 14,
       height = ceiling(length(PATHWAY_ORDER) / 2) * 4,
       device = "pdf", useDingbats = FALSE)
```

---

### 3.3 AUC Heatmap (Groups × Pathways) — ComplexHeatmap exception

Z-scored scaled mean AUC across groups. Uses ComplexHeatmap — pdf()/dev.off() required.
See CONVENTIONS.md §4 exception #1.

```r
# Compute mean AUC per label group
mean_auc <- sapply(PATHWAY_ORDER, function(pw) {
  tapply(obj@meta.data[[pw]], obj@meta.data[[LABEL_COL]], mean, na.rm = TRUE)
})   # rows = label groups, cols = pathways

# Z-score across groups per pathway
z_mat <- scale(mean_auc)

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D6604D"))
ht <- Heatmap(
  z_mat,
  name             = "Z-score\nAUC",
  col              = col_fun,
  cluster_rows     = TRUE,
  cluster_columns  = FALSE,
  column_order     = PATHWAY_ORDER,
  cell_fun         = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", z_mat[i, j]), x, y, gp = gpar(fontsize = 9))
  },
  row_names_gp     = gpar(fontsize = 12),
  column_names_gp  = gpar(fontsize = 12, fontface = "bold"),
  column_names_rot = 45
)

n_groups <- nrow(z_mat)
hm_h <- n_groups * 0.45 + 2
hm_w <- length(PATHWAY_ORDER) * 0.7 + 2.5

# ComplexHeatmap exception — pdf()/dev.off() required (see CONVENTIONS.md §4)
pdf(file.path(OUTPUT_DIR, "heatmap_celltype_x_pathway.pdf"),
    width = hm_w, height = hm_h)
draw(ht)
dev.off()
```

---

### 3.4 Evidence Dot Plot (Pathway × Group)

Direct gene-expression evidence; does NOT use AUC scores.

```r
make_evidence_dotplot <- function(seurat_obj, gene_sets, group_col, output_path) {
  group_levels <- unique(seurat_obj@meta.data[[group_col]])
  expr_mat <- GetAssayData(seurat_obj, layer = "data")

  rows <- lapply(names(gene_sets), function(pw) {
    genes_present <- intersect(gene_sets[[pw]], rownames(expr_mat))
    lapply(group_levels, function(g) {
      cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[group_col]] == g]
      if (length(cells) < 5 || length(genes_present) == 0) return(NULL)
      sub_mat <- as.matrix(expr_mat[genes_present, cells, drop = FALSE])
      data.frame(
        pathway     = pw,
        group       = g,
        pct_detect  = mean(colMeans(sub_mat > 0.1)),
        mean_expr   = mean(sub_mat)
      )
    })
  })

  dot_df <- bind_rows(unlist(rows, recursive = FALSE)) %>%
    mutate(pathway = factor(pathway, levels = rev(names(gene_sets))),
           group   = factor(group, levels = group_levels))

  PALETTE_EXPRESSION <- c("#F5F5F5","#FFF9C4","#FFB300","#E53935")

  p <- ggplot(dot_df, aes(x = group, y = pathway, size = pct_detect, fill = mean_expr)) +
    geom_point(shape = 21, color = "grey30", stroke = 0.32) +
    scale_fill_gradientn(colors = PALETTE_EXPRESSION, name = "Mean expr.") +
    scale_size_continuous(range = c(1, 8), name = "Pct genes > 0.1") +
    theme_classic(base_size = 13) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y  = element_text(size = 12, face = "bold"),
          axis.title   = element_blank(),
          legend.text  = element_text(size = 10))

  ggsave(output_path, p,
         width  = length(group_levels) * 0.75 + 3.0,
         height = length(gene_sets) * 0.52 + 2.5,
         device = "pdf", useDingbats = FALSE)
}

make_evidence_dotplot(
  obj, gene_sets, LABEL_COL,
  file.path(OUTPUT_DIR, "dotplot_evidence_celltype.pdf")
)

if (!is.null(SUBTYPE_COL) && SUBTYPE_COL %in% colnames(obj@meta.data)) {
  make_evidence_dotplot(
    obj, gene_sets, SUBTYPE_COL,
    file.path(OUTPUT_DIR, "dotplot_evidence_subtype.pdf")
  )
}
```

---

### 3.5 Per-Gene Expression Heatmap (ComplexHeatmap exception)

All pathway genes × groups (subtypes or cell types). Uses ComplexHeatmap.
See CONVENTIONS.md §4 exception #1.

```r
make_gene_heatmap <- function(seurat_obj, gene_sets, group_col, output_path) {
  expr_mat    <- GetAssayData(seurat_obj, layer = "data")
  group_levels <- unique(seurat_obj@meta.data[[group_col]])

  all_genes <- unlist(lapply(names(gene_sets), function(pw)
    intersect(gene_sets[[pw]], rownames(expr_mat))))
  pw_vec <- rep(names(gene_sets),
                sapply(names(gene_sets), function(pw)
                  length(intersect(gene_sets[[pw]], rownames(expr_mat)))))

  mat <- do.call(cbind, lapply(group_levels, function(g) {
    cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[group_col]] == g]
    rowMeans(as.matrix(expr_mat[all_genes, cells, drop = FALSE]))
  }))
  colnames(mat) <- group_levels

  col_fun <- colorRamp2(
    c(0, quantile(mat, 0.5), quantile(mat, 0.95)),
    c("#F5F5F5", "#FFB300", "#E53935")
  )

  ht <- Heatmap(
    mat,
    name              = "Mean expr.",
    col               = col_fun,
    row_split         = pw_vec,
    cluster_rows      = FALSE,
    cluster_columns   = TRUE,
    row_names_gp      = gpar(fontsize = 9, fontface = "italic"),
    column_names_gp   = gpar(fontsize = 12),
    column_names_rot  = 45
  )

  n_genes  <- length(all_genes)
  n_groups <- length(group_levels)
  hm_h     <- n_genes * 0.18 + 3
  hm_w     <- n_groups * 0.5 + 3

  # ComplexHeatmap exception — pdf()/dev.off() required (see CONVENTIONS.md §4)
  pdf(output_path, width = hm_w, height = hm_h)
  draw(ht)
  dev.off()
}

make_gene_heatmap(obj, gene_sets, LABEL_COL,
                  file.path(OUTPUT_DIR, "gene_heatmap_celltype.pdf"))

if (!is.null(SUBTYPE_COL) && SUBTYPE_COL %in% colnames(obj@meta.data)) {
  make_gene_heatmap(obj, gene_sets, SUBTYPE_COL,
                    file.path(OUTPUT_DIR, "gene_heatmap_subtype.pdf"))
}
```

---

## Output File Summary

| File | Description |
|---|---|
| `metabolic_scored.Rds` | Seurat object with AUC scores added as metadata columns |
| `auc_scores_per_cell.csv` | Per-cell AUC scores + label column |
| `umap_metabolic_scores.pdf` | One UMAP panel per pathway, colored by AUC |
| `violin_metabolic_by_celltype.pdf` | Violin + boxplot per pathway, x = label group |
| `heatmap_celltype_x_pathway.pdf` | Z-scored mean AUC heatmap (ComplexHeatmap) |
| `dotplot_evidence_celltype.pdf` | Evidence dot plot: pathway × label group |
| `dotplot_evidence_subtype.pdf` | Evidence dot plot: pathway × subtype (if SUBTYPE_COL set) |
| `gene_heatmap_celltype.pdf` | Per-gene expression heatmap (ComplexHeatmap) |
| `gene_heatmap_subtype.pdf` | Per-gene expression heatmap × subtype (if SUBTYPE_COL set) |

---

## Key Pitfalls and Fixes

### 1. Seurat v5 layer vs slot
Use `layer = "counts"` / `layer = "data"` in `GetAssayData()`, not `slot = "counts"`.

### 2. scale_color_identity conflict
Never use `scale_color_manual()` and `scale_color_identity()` in the same ggplot.
**Fix:** Pre-compute actual hex strings into the data frame column, use only `scale_color_identity()`.

### 3. AUCell memory
For objects >200k cells, subsample before AUCell (`SUBSAMPLE_CEILING = 2500` cells/group).
The `run_aucell()` primitive builds full cell rankings in memory. See Step 1.

### 4. Invalid range syntax
Gene sets with pseudo-notation (e.g. `"NDUFA1"-"NDUFA10"`) are NOT valid R. Always
enumerate gene vectors explicitly. See examples/ for validated explicit gene lists.

---

## Project-Specific Values (Stage for Phase 4 examples/)

- `examples/humanfat_metabolic_profile.md`:
  - Five metabolic pathway gene sets (Glycolysis, Beta_Oxidation, TCA_Cycle, OXPHOS, FA_Synthesis)
    with explicit gene vectors replacing the invalid range syntax from v1
  - PPARG_Targets gene set (biology-specific to PPARG OE project)
  - Myelolipoma exclusion rationale and filter
  - ec_subtype_colors from HumanFat validated_examples.yaml
  - Output subdirectory structure from MetabolicProfile run
