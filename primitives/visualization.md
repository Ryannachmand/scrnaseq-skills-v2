---
# Visualization — Parameterized R Plot Recipes — v2
# Migrated from ~/claude-skills/pipelines/LargeDataset/methods/aesthetics.md
# F5 FIX applied throughout:
#   - Proportion plot: patient_id / cell_fraction / source_file → function arguments
#   - ec_colors / tissue_colors / adipose_type_colors → group_colors / subtype_colors args
#   - 'EC Functional Dot Plot' → make_canonical_dotplot() with no EC-specific assumptions
#   - EC/HumanFat-specific size comments removed from stacked violin
# Last validated: HumanFat_Yang run (source), 2026-03-05
---

# Visualization — Parameterized R Plot Recipes

Concrete R implementations of the visual principles in @primitives/aesthetics.md.
Every function is parameterized — no hardcoded column names, color vectors, or project assumptions.

---

## UMAP — Labeled Cell Type Plot

```r
# Always use LabelClusters — never DimPlot(label=TRUE)
make_labeled_umap <- function(so, label_col, reduction = NULL, pt_size = 0.3) {
  if (is.null(reduction)) {
    reduction <- grep("umap", names(so@reductions), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(reduction)) stop("No UMAP found in object")
  }
  plot <- DimPlot(so, group.by = label_col, reduction = reduction,
                  pt.size = pt_size, raster = FALSE)  # raster=FALSE: raster renders invisible at small pt.size
  LabelClusters(plot = plot, id = label_col,
                box = TRUE, repel = TRUE, size = 4.65) +
    theme(legend.text = element_text(size = 14),
          title       = element_text(size = 16)) +
    xlab("UMAP Dimension 1") +
    ylab("UMAP Dimension 2")
}
```

---

## Stacked Violin Plot

```r
make_stacked_violin <- function(df, label_col, gene_col = "gene",
                                value_col = "value", fill_col = NULL,
                                output_file, n_genes = NULL, n_ctypes = NULL) {
  # df: tidy data frame with columns: label_col, gene_col, value_col
  # fill_col: column to use for fill; defaults to label_col
  fill_col <- if (is.null(fill_col)) label_col else fill_col
  n_genes  <- if (is.null(n_genes))  length(unique(df[[gene_col]]))   else n_genes
  n_ctypes <- if (is.null(n_ctypes)) length(unique(df[[label_col]])) else n_ctypes

  p <- ggplot(df, aes_string(x = label_col, y = value_col, fill = fill_col)) +
    geom_violin(scale = "width", trim = TRUE,
                linewidth = 0.3, color = "grey30", adjust = 1.2) +
    scale_y_continuous(
      breaks = function(x) c(0, floor(max(x))),
      labels = function(x) c("0", as.character(floor(x[2]))),
      expand = expansion(mult = c(0.05, 0.15))
    ) +
    labs(y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      legend.position  = "none",
      axis.text.y      = element_text(size = 9.85, color = "grey30"),
      axis.title.y     = element_text(size = 13, face = "italic",
                                      angle = 0, hjust = 1, vjust = 0.5),
      axis.line.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.line.y      = element_line(color = "grey60", linewidth = 0.3),
      axis.ticks.y     = element_line(color = "grey60", linewidth = 0.3),
      plot.margin      = margin(t = 1, b = 1, l = 2, r = 2)
    )

  fig_h <- n_genes  * 0.51 + 1.5
  fig_w <- n_ctypes * 0.35 + 1.5
  ggsave(output_file, plot = p, width = fig_w, height = fig_h,
         units = "in", device = "pdf", useDingbats = FALSE)
  p
}
# Bottom panel only: add axis.text.x = element_text(size=13.5, angle=45, hjust=1, color="grey20")
# Stacking: wrap_plots(plot_list, ncol=1) — NOT patchwork / operator (fails >~10 panels)
```

---

## Cell Type Proportion Plot

X-axis uses numeric positions, not discrete sample labels.
Samples grouped by `group_col`, labeled by `subtype_col`.

```r
make_proportion_plot <- function(df,
                                 sample_col,      # column identifying individual samples
                                 group_col,       # column grouping samples (e.g. patient_id)
                                 subtype_col,     # column holding cell type proportions
                                 group_colors,    # named vector: subtype → hex color
                                 output_file,
                                 meta_tracks = list()) {
  # df: one row per sample × subtype combination with a 'proportion' column
  # meta_tracks: optional named list of additional metadata data frames for tile tracks
  #   each element: list(df=..., fill_scale=..., label_col=..., label_size=1.8)

  # Numeric x positions — samples within each group are adjacent
  df <- df %>%
    arrange(!!sym(group_col), !!sym(sample_col)) %>%
    mutate(x_pos = as.numeric(factor(!!sym(sample_col),
                                     levels = unique(!!sym(sample_col)))))

  separator_x <- df %>%
    group_by(!!sym(group_col)) %>%
    summarise(sep = max(x_pos) + 0.5, .groups = "drop") %>%
    filter(sep < max(df$x_pos)) %>%
    pull(sep)

  p_main <- ggplot(df, aes(x = x_pos, y = proportion, fill = !!sym(subtype_col))) +
    geom_bar(stat = "identity", width = 0.85, color = NA) +
    geom_vline(xintercept = separator_x, color = "white", linewidth = 0.8) +
    scale_fill_manual(values = group_colors) +
    scale_x_continuous(expand = expansion(add = 0.5)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = 0)) +
    theme_classic(base_size = 14) +
    theme(legend.position = "right")

  n_samples <- length(unique(df[[sample_col]]))
  fig_h <- 5
  fig_w <- max(10 / 1.65, n_samples * 0.38 / 1.65 + 3.5)

  # If meta_tracks provided, assemble with patchwork
  if (length(meta_tracks) > 0) {
    track_plots <- lapply(meta_tracks, function(tr) {
      ggplot(tr$df, aes(x = x_pos, y = 1, fill = !!sym(tr$label_col))) +
        geom_tile() + tr$fill_scale +
        theme_void() + theme(legend.position = "none")
    })
    heights <- c(10, rep(0.6, length(track_plots)))
    p_final <- wrap_plots(c(list(p_main), track_plots), ncol = 1) +
      plot_layout(heights = heights, guides = "collect")
  } else {
    p_final <- p_main
  }

  ggsave(output_file, plot = p_final, width = fig_w, height = fig_h,
         units = "in", device = "pdf", useDingbats = FALSE)
  p_final
}
```

---

## Canonical Dot Plot (`make_canonical_dotplot`)

General-purpose dot plot with shape=21 (filled circle), size=pct_exp, fill=avg_exp_scaled.
Replaces the v1 "EC Functional Dot Plot" — no EC-specific assumptions.

```r
make_canonical_dotplot <- function(dot_df,
                                   gene_col     = "gene",
                                   group_col    = "group",
                                   pct_col      = "pct_exp",
                                   avg_col      = "avg_exp_scaled",
                                   group_colors = NULL,  # named vector for x-axis label colors
                                   output_file,
                                   n_genes = NULL, n_groups = NULL,
                                   section_dividers = NULL,   # numeric yintercepts
                                   section_labels   = NULL,   # data.frame(y, label, pin_to_group)
                                   x_position = "top") {
  # dot_df: tidy data frame with gene_col, group_col, pct_col, avg_col
  # section_dividers: optional numeric vector of yintercepts for dashed divider lines
  # section_labels: optional df for geom_text labels pinned to one facet (avoids annotate() in facets)
  n_genes  <- if (is.null(n_genes))  length(unique(dot_df[[gene_col]]))   else n_genes
  n_groups <- if (is.null(n_groups)) length(unique(dot_df[[group_col]])) else n_groups

  p <- ggplot(dot_df, aes_string(x = group_col, y = gene_col,
                                  size = pct_col, fill = avg_col)) +
    geom_point(shape = 21, color = "grey30", stroke = 0.32) +
    scale_fill_gradientn(colors = c("#F5F5F5", "#FFF9C4", "#FFB300", "#E53935"),
                         name = "Scaled\nExpression") +
    scale_size_continuous(range = c(0.3, 6), limits = c(0, 100),
                          breaks = c(10, 25, 50, 75),
                          name = "% Expressing") +
    scale_x_discrete(position = x_position) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x.top  = element_text(size = 14.5, angle = 45, hjust = 0, vjust = 0,
                                       face = "bold"),
      axis.text.x      = element_text(size = 14.5, angle = 45, hjust = 1),
      axis.text.y      = element_text(size = 12.5, color = "grey10", face = "italic"),
      panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      plot.margin      = margin(t = 10, r = 90, b = 10, l = 10),
      legend.position  = "none"
    ) +
    coord_cartesian(clip = "off")

  # Color x-axis labels by cell type if group_colors provided
  if (!is.null(group_colors)) {
    group_order <- levels(factor(dot_df[[group_col]]))
    label_colors <- group_colors[group_order]
    p <- p + theme(axis.text.x.top = element_text(color = label_colors))
  }

  # Section dividers
  if (!is.null(section_dividers)) {
    p <- p + geom_hline(yintercept = section_dividers, color = "grey40",
                        linewidth = 0.5, linetype = "dashed")
  }

  # Section labels — pinned to one group to avoid rendering in all facet panels
  # section_labels df must include: y (numeric), label (character), <group_col> (factor pinned to one group)
  if (!is.null(section_labels)) {
    p <- p + geom_text(data = section_labels,
                       aes_string(x = as.name(group_col), y = "y", label = "label"),
                       inherit.aes = FALSE, size = 4.5, hjust = 0,
                       fontface = "bold.italic", color = "grey35")
  }

  fig_h <- n_genes  * 0.28 + 2
  fig_w <- n_groups * 0.75 + 3.0
  ggsave(output_file, plot = p, width = fig_w, height = fig_h,
         units = "in", device = "pdf", useDingbats = FALSE)
  p
}
```

**Legend handling:** Save legend separately to avoid cramped main plot:
```r
leg <- cowplot::get_legend(p + theme(legend.position = "right"))
ggsave(legend_file, plot = cowplot::plot_grid(leg), width = 2.5, height = 4,
       units = "in", device = "pdf", useDingbats = FALSE)
```

---

## TF Diamond Dot Plot (differences from canonical dot plot)

```r
make_tf_diamond_plot <- function(dot_df, gene_col = "gene", group_col = "group",
                                  pct_col = "pct_exp", lfc_col = "avg_log2FC",
                                  output_file, n_genes = NULL, n_groups = NULL) {
  # Same structure as make_canonical_dotplot but:
  # - shape = 23 (diamond) with stroke = 0.4, color = "grey25"
  # - fill = diverging blue-white-red (TF fold change, not expression)
  # - size range = c(0.8, 5.5)
  # - section dividers: solid, not dashed
  # - gene names italic on y-axis
  n_genes  <- if (is.null(n_genes))  length(unique(dot_df[[gene_col]]))   else n_genes
  n_groups <- if (is.null(n_groups)) length(unique(dot_df[[group_col]])) else n_groups

  p <- ggplot(dot_df, aes_string(x = group_col, y = gene_col,
                                  size = pct_col, fill = lfc_col)) +
    geom_point(shape = 23, color = "grey25", stroke = 0.4) +
    scale_fill_gradientn(
      colors = c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#D6604D","#B2182B"),
      values = scales::rescale(c(0, 0.2, 0.5, 0.7, 0.85, 1)),
      name = "log2FC"
    ) +
    scale_size_continuous(range = c(0.8, 5.5), name = "% Expressing") +
    scale_x_discrete(position = "top") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x.top = element_text(size = 14.5, angle = 45, hjust = 0, vjust = 0, face = "bold"),
      axis.text.y     = element_text(size = 12.5, color = "grey10", face = "italic"),
      panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      plot.margin     = margin(t = 10, r = 90, b = 10, l = 10),
      legend.position = "none"
    ) +
    coord_cartesian(clip = "off")

  fig_h <- n_genes  * 0.29 + 2
  fig_w <- n_groups * 0.75 + 3.0
  ggsave(output_file, plot = p, width = fig_w, height = fig_h,
         units = "in", device = "pdf", useDingbats = FALSE)
  p
}
```

---

## Differential Abundance Heatmap

```r
make_diff_abundance_heatmap <- function(results_df,
                                         x_col = "group",    # columns (e.g. tissue types)
                                         y_col = "cluster",  # rows (e.g. cell types)
                                         fill_col = "log2OR",
                                         sig_col = "p_adj",
                                         sig_threshold = 0.01,
                                         or_threshold = 1.5,
                                         output_file) {
  df <- results_df %>%
    mutate(label = ifelse(!!sym(sig_col) < sig_threshold & abs(!!sym(fill_col)) > log2(or_threshold),
                          "*", ""))

  p <- ggplot(df, aes_string(x = x_col, y = y_col, fill = fill_col)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 7, vjust = 0.75) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
    theme_classic() +
    theme(
      axis.text.x   = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y   = element_text(size = 14),
      axis.title.x  = element_text(size = 15),
      axis.title.y  = element_text(size = 15),
      plot.title    = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12.5),
      legend.text   = element_text(size = 12.5),
      legend.title  = element_text(size = 12.5)
    )

  ggsave(output_file, plot = p, units = "in", device = "pdf", useDingbats = FALSE)
  p
}
```

---

## Trajectory / Pseudotime Plots

```r
make_pseudotime_umap <- function(cds, output_file) {
  p <- plot_cells(cds, color_cells_by = "pseudotime", cell_size = 0.6) +
    theme_classic(base_size = 16) +
    scale_color_viridis_c(option = "magma")  # magma — NOT inferno
  ggsave(output_file, plot = p, units = "in", device = "pdf", useDingbats = FALSE)
  p
}

make_pseudotime_violin <- function(df, celltype_col, pseudotime_col,
                                    group_colors = NULL, output_file) {
  # Sorted by median pseudotime per cell type
  ct_order <- df %>%
    group_by(!!sym(celltype_col)) %>%
    summarise(med = median(!!sym(pseudotime_col)), .groups = "drop") %>%
    arrange(med) %>% pull(!!sym(celltype_col))

  df[[celltype_col]] <- factor(df[[celltype_col]], levels = ct_order)

  p <- ggplot(df, aes_string(x = celltype_col, y = pseudotime_col,
                               fill = if (!is.null(group_colors)) celltype_col else NULL)) +
    geom_violin(trim = TRUE, scale = "width") +
    stat_summary(fun = median, geom = "point", size = 2.5, color = "#E53935") +
    geom_boxplot(width = 0.07, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (!is.null(group_colors)) p <- p + scale_fill_manual(values = group_colors)
  ggsave(output_file, plot = p, units = "in", device = "pdf", useDingbats = FALSE)
  p
}
```
