---
requires_context:
  palettes:
    - group_colors    # optional: named vector for column background highlighting
  metadata_columns:
    required:
      - source_col    # metadata column identifying dataset source
      - subtype_col   # cell type column
    optional:
      - col_group_col # metadata column for column background grouping
  brief_keys:
    required:
      - output_dir
      - downstream_analyses.cross_dataset_dotplot.source_col
      - downstream_analyses.cross_dataset_dotplot.subtype_col
      - downstream_analyses.cross_dataset_dotplot.sources_order
      - downstream_analyses.cross_dataset_dotplot.marker_genes
      - downstream_analyses.cross_dataset_dotplot.col_group_col
    optional:
      - downstream_analyses.cross_dataset_dotplot.n_variable
      - downstream_analyses.cross_dataset_dotplot.z_cap
      - downstream_analyses.cross_dataset_dotplot.group_colors
      - downstream_analyses.cross_dataset_dotplot.sec1_force
references:
  - "@primitives/visualization.md"
  - "@primitives/aesthetics.md"
  - "@context/color_palettes.md"
  - "@primitives/seurat_v5_rules.md"
---

# Module: Cross-Dataset Comparison Dotplot

Compares cell subtypes across two or more datasets with depth correction via within-dataset
z-scoring. Produces a multi-section dotplot and a centroid correlation heatmap in Harmony
space.

**Known-Atlas Convention (Phase 3 concern):** When a brief names a recognized atlas (e.g.
"Tabula Sapiens"), the deployment agent can auto-populate `source_col`, `sources_order`,
and column label conventions from the context registry (`context/known_atlases` — proposed,
not yet implemented). See PHASE2C_REPORT.md for the design. This module is atlas-agnostic;
all dataset identity flows through parameters.

---

## Depth Correction Assumption — Read Before Using

**This module z-scores expression within each dataset independently before combining.**
This corrects sequencing depth differences between datasets but rests on the assumption
that **both datasets have been processed with comparable normalization and the same gene
panel**. If the datasets use fundamentally different protocols (e.g., very different
transcript detection rates or partially overlapping gene panels), within-dataset z-scoring
may shift but not fully correct the batch effect — and the combined dotplot color scale
will still be skewed toward the dataset with higher absolute expression values.

Verify this assumption by checking the per-gene z-score distributions after combining.
Any gene where one dataset consistently hits the z_cap in one direction is a candidate
depth-correction failure — investigate before including in figures.

---

## `sec1_force` — Explicit Analysis Decision Required

`sec1_force` allows specific genes to be forced into Section 1 regardless of whether they
pass the normal marker selection criteria. This is a **project-specific analysis decision**
that must be justified per project — it is not a general-purpose feature. Document the
rationale explicitly in the examples/ file. Use it only when there is a biological reason
to include genes that appear in one dataset's markers but not the other's.

---

## Brief Schema

```yaml
downstream_analyses:
  cross_dataset_dotplot:
    enabled: true
    source_col: project_specific        # REQUIRED: metadata column identifying dataset source
    subtype_col: project_specific       # REQUIRED: cell type column
    sources_order: project_specific     # REQUIRED: character vector of source values in column display order
    marker_genes: project_specific      # REQUIRED: named list of gene section vectors
                                        # e.g. list("Section 1" = c(...), "Section 2" = c(...))
    col_group_col: project_specific     # REQUIRED: metadata column for column background highlight grouping
    n_variable: 50                      # top variable genes for centroid correlation
    z_cap: 2                            # z-score cap applied after within-dataset scaling
    group_colors: null                  # optional: named vector for column background highlight colors
    sec1_force: null                    # optional: genes forced into Section 1
                                        # requires per-project justification — see module header
```

---

## Critical Constraints

| ❌ Don't | ✅ Do instead | Why |
|---|---|---|
| Combine avg_exp across datasets without z-scoring | Z-score per gene within each dataset first | Depth differences inflate atlas values relative to in-house |
| Use `expm1(mean(expression))` as color scale | Use `mean(expression)` then z-score | `expm1` amplifies depth differences |
| Use shared color scale across datasets pre-correction | Diverging scale centered at 0 after z-scoring | Z-scores are on the same scale only after the correction |
| Hardcode column order | Build from `sources_order` parameter | Column order changes with different source values |
| Read correlation matrix from CSV without fixing dots | Apply `gsub("\\.", " ", ...)` on read | `write.csv` converts spaces to dots in column names |
| Use `annotate("segment")` for section dividers on discrete y | `geom_hline(yintercept = ...)` | `annotate("segment")` with numeric y fails on discrete scales |
| Use `useDingbats = FALSE` with `cairo_pdf` | Omit `useDingbats` when device is `cairo_pdf` | `useDingbats` is only valid for the base `pdf()` device |
| Include cell type names with special characters without handling | Convert spaces consistently before CSV round-trip | Whitespace-to-dot conversion in `write.csv` breaks column matching |
| Use `select()` without namespace | `dplyr::select()` when Seurat loaded | Seurat masks dplyr's `select` |

---

## R Packages Required

```r
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)    # centroid correlation heatmap
library(circlize)          # colorRamp2
```

---

## Configuration Block

```r
# ── CONFIG ────────────────────────────────────────────────────────────────────
SOURCE_COL    <- "project_specific"    # REPLACE: metadata column identifying dataset source
SUBTYPE_COL   <- "project_specific"   # REPLACE: cell type column
COL_GROUP_COL <- "project_specific"   # REPLACE: column for background group highlighting

SOURCES_ORDER <- c("project_specific")  # REPLACE: source values in display order, left to right
                                         # in-house sources first, atlas sources after

# Marker genes: named list of gene sections for the dotplot
# Section names become section divider labels in the figure
marker_genes <- list(
  "Section 1" = c("project_specific"),  # REPLACE: e.g. shared markers
  "Section 2" = c("project_specific"),  # REPLACE: e.g. in-house-specific markers
  "Section 3" = c("project_specific")   # REPLACE: e.g. atlas-specific markers
)

# Optional: force specific genes into Section 1 regardless of marker selection criteria
# Requires explicit per-project justification — document rationale in examples/ file
# Leave as NULL for standard marker selection
FORCE_GENES_SEC1 <- NULL   # REPLACE: e.g. c("GENE1", "GENE2") — with justification

# Background highlight colors for column groups
# e.g. warm yellow for in-house columns, light red for atlas columns of interest
GROUP_COLORS <- c("project_specific" = "#project_specific")  # REPLACE

N_VARIABLE <- 50    # top variable genes for centroid correlation Harmony space
Z_CAP      <- 2     # z-score cap: expression clamped to [-Z_CAP, Z_CAP] after scaling

RDS_IN     <- "project_specific"   # REPLACE: integrated Seurat object with Harmony reduction
OUTPUT_DIR <- file.path("output", "cross_dataset_dotplot")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
```

---

## Part A: Centroid Correlation Heatmap in Harmony Space

Compute per-group centroids in Harmony embedding, then Pearson correlation.
This avoids gene-space batch effects that distort pseudobulk comparisons.

```r
so <- readRDS(RDS_IN)
so <- JoinLayers(so)   # @primitives/seurat_v5_rules.md Rule 1

# Detect Harmony reduction — @primitives/seurat_v5_rules.md Rule 5
harm_key <- names(so@reductions)[grepl("harmony", names(so@reductions), ignore.case = TRUE)][1]
if (is.na(harm_key)) stop("No Harmony reduction found in object.")

# Build combined group label: one label per (source, subtype) combination
# Convention: paste source prefix + subtype value — adapt separator to project naming
so@meta.data$group_label <- paste0(so@meta.data[[SOURCE_COL]], " ", so@meta.data[[SUBTYPE_COL]])

harm_emb <- Embeddings(so, harm_key)
group_labels <- sort(unique(so@meta.data$group_label))

centroids <- t(sapply(group_labels, function(g) {
  cells <- rownames(so@meta.data)[so@meta.data$group_label == g]
  colMeans(harm_emb[cells, , drop = FALSE])
}))

# Use top N_VARIABLE variable genes to weight the centroid space
var_genes <- head(VariableFeatures(so), N_VARIABLE)
cor_mat <- cor(t(centroids))

# Save correlation matrix — fix space→dot conversion on read-back
cor_csv <- file.path(OUTPUT_DIR, "centroid_correlation.csv")
write.csv(cor_mat, cor_csv)

# Read back and restore spaces
cor_mat_read <- read.csv(cor_csv, row.names = 1, check.names = FALSE)
colnames(cor_mat_read) <- gsub("\\.", " ", colnames(cor_mat_read))
rownames(cor_mat_read) <- gsub("\\.", " ", rownames(cor_mat_read))
cor_mat <- as.matrix(cor_mat_read)

# Order: rows/cols follow SOURCES_ORDER, then alphabetical within each source
source_order_idx <- order(match(
  sub(" .*", "", group_labels),   # extract source prefix
  SOURCES_ORDER
))
group_order <- group_labels[source_order_idx]
group_order <- group_order[group_order %in% rownames(cor_mat)]

col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))

ht <- ComplexHeatmap::Heatmap(
  cor_mat[group_order, group_order],
  col              = col_fun,
  row_order        = seq_along(group_order),
  column_order     = seq_along(group_order),
  cell_fun         = function(j, i, x, y, w, h, f)
    grid.text(sprintf("%.2f", cor_mat[group_order[i], group_order[j]]),
              x, y, gp = gpar(fontsize = 7)),
  row_names_gp    = gpar(fontsize = 11),
  column_names_gp = gpar(fontsize = 11),
  heatmap_legend_param = list(title = "Pearson r", at = c(-1, -0.5, 0, 0.5, 1))
)

# ComplexHeatmap exception: pdf()/draw()/dev.off() required (CONVENTIONS.md §4 exception #1)
pdf(file.path(OUTPUT_DIR, "centroid_correlation_heatmap.pdf"), width = 10, height = 9)
ComplexHeatmap::draw(ht)
dev.off()
```

---

## Part B: Depth-Corrected Multi-Section Dotplot

Within-dataset z-scoring corrects sequencing depth differences. See depth correction
assumption in the module header — verify this assumption for your datasets before
interpreting color differences between datasets as biological.

```r
DefaultAssay(so) <- "RNA"
meta <- so@meta.data

# Build dot statistics per (source, subtype, gene) — mean log-normalized expression
all_genes <- unique(unlist(marker_genes))
if (!is.null(FORCE_GENES_SEC1)) all_genes <- unique(c(all_genes, FORCE_GENES_SEC1))
all_genes <- intersect(all_genes, rownames(so))

expr <- as.data.frame(t(as.matrix(
  GetAssayData(so, assay = "RNA", layer = "data")[all_genes, ]
)))
expr[[SOURCE_COL]]  <- meta[[SOURCE_COL]]
expr[[SUBTYPE_COL]] <- meta[[SUBTYPE_COL]]

dot_long <- tidyr::pivot_longer(expr, cols = all_of(all_genes),
                                  names_to = "gene", values_to = "expression") %>%
  group_by(!!rlang::sym(SOURCE_COL), !!rlang::sym(SUBTYPE_COL), gene) %>%
  summarise(
    avg_exp = mean(expression),   # mean log-normalized, NOT expm1 — see Critical Constraints
    pct_exp = mean(expression > 0) * 100,
    .groups = "drop"
  )

# Depth correction: z-score each gene within each source independently
dot_inhouse <- dot_long %>%
  filter(!!rlang::sym(SOURCE_COL) != SOURCES_ORDER[length(SOURCES_ORDER)]) %>%   # adjust slice if >2 sources
  group_by(gene) %>%
  mutate(avg_exp_z = as.numeric(scale(avg_exp))) %>%
  ungroup()

dot_atlas <- dot_long %>%
  filter(!!rlang::sym(SOURCE_COL) == SOURCES_ORDER[length(SOURCES_ORDER)]) %>%
  group_by(gene) %>%
  mutate(avg_exp_z = as.numeric(scale(avg_exp))) %>%
  ungroup()

# For > 2 sources: apply group_by(source, gene) %>% mutate(avg_exp_z = scale(avg_exp))
# to each source independently, then bind_rows.

dot_all <- bind_rows(dot_inhouse, dot_atlas) %>%
  mutate(avg_exp_z = pmax(pmin(avg_exp_z, Z_CAP), -Z_CAP))   # cap at ±Z_CAP

# Column order: sources in SOURCES_ORDER, subtypes alphabetical within each source.
# NOTE: Intentionally not purely alphabetical — SOURCES_ORDER governs source grouping.
# This is the documented exception to the global alphabetical-column rule in visualization.md.
# Handle whitespace in subtype names: column_order uses group_label (source + subtype)
dot_all$col_id <- paste0(dot_all[[SOURCE_COL]], " ", dot_all[[SUBTYPE_COL]])

col_labels <- sort(unique(dot_all$col_id))
col_order  <- col_labels[order(match(sub(" .*", "", col_labels), SOURCES_ORDER))]
# If col_order still has whitespace issues after CSV round-trip, apply gsub("\\.", " ", ...)

# Gene (row) order: section membership from marker_genes, apply FORCE_GENES_SEC1 if set.
# NOTE: Gene order is caller-specified via marker_genes — diagonal ordering does not apply here.
section_df <- bind_rows(lapply(names(marker_genes), function(s) {
  genes <- if (s == names(marker_genes)[1] && !is.null(FORCE_GENES_SEC1)) {
    unique(c(FORCE_GENES_SEC1, marker_genes[[s]]))   # forced genes prepended to Section 1
  } else {
    marker_genes[[s]]
  }
  data.frame(gene = genes, section = s, stringsAsFactors = FALSE)
})) %>%
  filter(gene %in% all_genes) %>%
  filter(!duplicated(gene))

gene_order   <- section_df$gene
section_vec  <- section_df$section
n_sections   <- length(names(marker_genes))

dot_plot_df <- dot_all %>%
  filter(gene %in% gene_order) %>%
  mutate(
    gene   = factor(gene, levels = rev(gene_order)),
    col_id = factor(col_id, levels = col_order)
  ) %>%
  filter(!is.na(col_id))

n_cols   <- length(col_order)
n_genes  <- length(gene_order)

# Section divider y-positions (for geom_hline on discrete y-axis)
section_sizes  <- rev(table(factor(section_vec, levels = names(marker_genes))))
divider_ys     <- cumsum(section_sizes)[-length(section_sizes)] + 0.5

p <- ggplot(dot_plot_df, aes(x = col_id, y = gene, size = pct_exp, fill = avg_exp_z)) +
  geom_point(shape = 21, color = "grey30", stroke = 0.32) +
  geom_hline(yintercept = divider_ys, linetype = "dotted", color = "grey50", linewidth = 0.3) +
  scale_fill_gradientn(
    colors = c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B"),
    values = scales::rescale(c(-2, -0.5, 0, 0.5, 1, 2)),
    limits = c(-Z_CAP, Z_CAP),
    name   = "Z-score\n(depth-corrected)"
  ) +
  scale_size_continuous(range = c(0.3, 6), limits = c(0, 100), name = "% Expressing") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x.top = element_text(size = 13, angle = 45, hjust = 0, vjust = 0, face = "bold"),
    axis.text.y     = element_text(size = 14, face = "italic"),
    panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.4),
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10)
  )

# Optional: column background highlighting for source groups
# Requires COL_GROUP_COL values to match col_id source prefix pattern
if (!is.null(GROUP_COLORS) && length(GROUP_COLORS) > 0) {
  bg_df <- data.frame(
    col_id    = factor(col_order, levels = col_order),
    col_group = sub(" .*", "", col_order)   # extract source prefix
  ) %>%
    group_by(col_group) %>%
    mutate(
      x_min = as.numeric(col_id) - 0.5,
      x_max = as.numeric(col_id) + 0.5
    ) %>%
    summarise(x_min = min(x_min), x_max = max(x_max), .groups = "drop")

  for (grp in bg_df$col_group) {
    if (grp %in% names(GROUP_COLORS)) {
      row_grp <- bg_df[bg_df$col_group == grp, ]
      p <- p + annotate("rect",
                         xmin = row_grp$x_min, xmax = row_grp$x_max,
                         ymin = -Inf, ymax = Inf,
                         fill = GROUP_COLORS[grp], alpha = 0.15)
    }
  }
}

fig_w <- n_cols  * 0.62 + 4.5
fig_h <- n_genes * 0.38 + 4.0

# Use cairo_pdf for dotplots — omit useDingbats (incompatible with cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "cross_dataset_dotplot.pdf"),
       p, width = fig_w, height = fig_h, units = "in", device = cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "cross_dataset_dotplot.png"),
       p, width = fig_w, height = fig_h, units = "in", dpi = 200, bg = "white")
```

---

## Part C: TF Diamond Variant

Same structure as Part B but:
- `shape = 23` (diamond) instead of `21` (circle)
- Gene universe filtered to transcription factors from a caller-provided TF list
- `sec1_force` may differ from the standard dotplot (document separately if so)

```r
# TF list — caller-provided, biology-specific
# For EC biology: Lambert et al. 2018 ETS/KLF/SOX/FOX list; see examples/
# For other cell types: provide curated TF list appropriate to the biology
human_tfs <- c("project_specific")   # REPLACE: curated TF list

tf_genes   <- intersect(human_tfs, all_genes)
tf_in_data <- intersect(tf_genes, unique(dot_plot_df$gene))

if (length(tf_in_data) == 0) {
  message("No TF genes found in dotplot gene list — skipping diamond variant.")
} else {
  p_tf <- p %+%
    filter(dot_plot_df, gene %in% tf_in_data) +
    geom_point(shape = 23, color = "grey30", stroke = 0.32)   # override to diamond

  n_tfs <- length(tf_in_data)
  fig_h_tf <- n_tfs * 0.38 + 4.0

  ggsave(file.path(OUTPUT_DIR, "cross_dataset_tf_diamond.pdf"),
         p_tf, width = fig_w, height = fig_h_tf, units = "in", device = cairo_pdf)
}
```

---

## Figure Sizing Reference

| Plot | Width | Height |
|---|---|---|
| Dotplot | `n_cols * 0.62 + 4.5` | `n_genes * 0.38 + 4.0` |
| TF diamond variant | same width | `n_tfs * 0.38 + 4.0` |
| Centroid heatmap | 10 | 9 |

