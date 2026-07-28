---
requires_context:
  palettes:
    - cell_colors        # required — named vector mapping every cell type in analysis to a hex color
    - source_colors      # optional — subset of cell_colors for source cell types (EC subtypes);
                         #   used for x-axis label coloring in bubble and bar plots;
                         #   auto-derived as cell_colors[source_types] if not provided separately
  metadata_columns:
    required:
      - label_col        # cell type label column in Seurat metadata
    optional:
      - group_col        # condition/tissue column for Script 6 comparison mode (null = skip Script 6)
  brief_keys:
    required:
      - output_dir
      - downstream_analyses.cellchat.label_col
      - downstream_analyses.cellchat.organism
    optional:
      - downstream_analyses.cellchat.group_col
      - downstream_analyses.cellchat.n_workers
      - downstream_analyses.cellchat.signaling_pathways
      - downstream_analyses.cellchat.source_celltypes
      - downstream_analyses.cellchat.target_celltypes
references:
  - "@primitives/seurat_v5_rules.md"
  - "@primitives/r_environment.md"
  - "@primitives/visualization.md"
  - "@primitives/aesthetics.md"
  - "@context/color_palettes.md"
---

# Module: CellChat v2 Interactome Analysis

Cell-cell communication analysis using CellChat v2 on a pre-annotated Seurat object.
Six-script architecture: inference → basic plots → stacked bubbles → bar plots →
circos diagrams → tissue comparison (conditional).

**Validated on:** two independent in-house projects (2026-Q1).
Validated against up to 19 cell types, 104 pathways, and 16k+ LR interactions.

**Never re-run inference for plotting iterations.** Script 1 produces the saved
CellChat object; Scripts 2–6 load that object and generate all downstream plots.

---

## Brief Schema

```yaml
downstream_analyses:
  cellchat:
    enabled: true
    label_col: project_specific         # cell type label column in Seurat metadata
    group_col: null                     # condition/tissue column; null = skip Script 6
    organism: "human"                   # "human" or "mouse" (selects PPI database)
    n_workers: 1                        # MUST remain 1 — see Critical Constraints
    signaling_pathways: null            # optional list of pathways to focus on; null = top by net weight
    source_celltypes: project_specific  # cell types to use as senders/sources
    target_celltypes: project_specific  # cell types to use as receivers/targets
    output_dir: project_specific        # defaults to output/cellchat/
```

---

## Critical Constraints

| Do not use | Use instead | Why |
|---|---|---|
| `projectData(cellchat, PPI.human)` | *(remove this line entirely)* | **Function removed in CellChat v2.** Throws `could not find function "projectData"`. The v1 library erroneously included this call — do not copy it. |
| `netVisual_aggregate(signaling=NULL)` | `netVisual_circle()` on `@net$count` | CellChat v2 crashes when `signaling=NULL` |
| `netVisual_heatmap()` with pathway subset | `plot_heatmap()` (defined below) | Internal annotation dimension bug when filtering pathways |
| `pairLR.use=` parameter in `netVisual_chord_gene` | `net=` parameter with pre-built df from `get_filtered_net()` | Precedence bug — validation always triggers an error regardless of input |
| `top=` in `netVisual_chord_gene` | `get_filtered_net()` helper | Parameter does not exist; throws "unused argument" |
| `plan("multisession")` or default future | `plan("sequential")` at the very top of every script | Parallelization crashes chord diagram rendering in CellChat v2 |
| `cell.level = TRUE` in `mergeCellChat()` | *(remove this argument)* | Parameter removed in CellChat v2; throws "unused argument" |
| `netAnalysis_computeCentrality()` on merged object | Run on each individual object before merging | Throws `rowSums` error on merged objects |
| Load all tissue subsets simultaneously | Load one at a time, extract matrix, `rm(obj); gc()` | OOM / exit 137 on large objects |

**`plan("sequential")` is mandatory.** Place it at the very top of every CellChat script —
before any library loads, not just before chord diagram calls. CellChat v2 interacts with
the future package in ways that crash chord diagram rendering if any parallelization plan
was previously set.

---

## Module-Level Helper Functions

These functions are module-level helpers defined here because they exist to work around
specific CellChat v2 bugs. They are NOT promoted to primitives (self-contained workarounds).
Include them at the top of every plotting script.

### `subset_net()` — slice pathway probability array for category plots

```r
subset_net <- function(cellchat, pathway_names) {
  avail <- pathway_names[pathway_names %in% cellchat@netP$pathways]
  if (length(avail) == 0) return(NULL)
  prob_arr   <- cellchat@netP$prob
  idx        <- which(cellchat@netP$pathways %in% avail)
  sub_arr    <- prob_arr[, , idx, drop = FALSE]
  count_mat  <- apply(sub_arr, c(1, 2), function(x) sum(x > 0))
  weight_mat <- apply(sub_arr, c(1, 2), sum)
  rownames(count_mat)  <- colnames(count_mat)  <- rownames(prob_arr)
  rownames(weight_mat) <- colnames(weight_mat) <- rownames(prob_arr)
  list(count = count_mat, weight = weight_mat, paths = avail)
}
```

### `get_filtered_net()` — top-N LR pairs per pathway for chord diagrams

Bypasses the `pairLR.use` validation bug by building a pre-filtered `net`
data frame and passing it via the `net=` parameter (not `pairLR.use=`).

```r
get_filtered_net <- function(cellchat, pathways, sources, targets, top_n = 2) {
  prob <- slot(cellchat, "net")$prob
  pval <- slot(cellchat, "net")$pval
  prob[pval > 0.05] <- 0
  net <- reshape2::melt(prob, value.name = "prob")
  colnames(net)[1:3] <- c("source", "target", "interaction_name")

  pairLR <- cellchat@LR$LRsig
  cols <- intersect(c("interaction_name_2", "pathway_name", "ligand",
                      "receptor", "annotation", "evidence"), colnames(pairLR))
  idx  <- match(net$interaction_name, rownames(pairLR))
  net  <- cbind(net, pairLR[idx, cols, drop = FALSE])

  net  <- net[net$pathway_name %in% pathways & net$prob > 0, ]
  net  <- net[net$source %in% sources & net$target %in% targets, ]

  lr_totals <- aggregate(prob ~ interaction_name + pathway_name, net, sum)
  keep <- do.call(rbind, lapply(unique(lr_totals$pathway_name), function(p) {
    sub <- lr_totals[lr_totals$pathway_name == p, ]
    head(sub[order(-sub$prob), ], top_n)
  }))
  net[net$interaction_name %in% keep$interaction_name, ]
}
```

**top_n guidance:**
- Category-level chord diagrams: `top_n = 2`
- Spotlight / single-cell-type chord diagrams: `top_n = 1` (cleaner gene labels)

### `plot_heatmap()` — custom ggplot heatmap (replaces buggy `netVisual_heatmap`)

`netVisual_heatmap()` has an annotation dimension bug when filtering to a pathway subset.
Use this ggplot replacement for all per-category heatmaps. (The merged-object comparison
heatmap in Script 6 uses `netVisual_heatmap()` + `draw()` — that case is safe because
no pathway subset is applied.)

```r
plot_heatmap <- function(mat, title, grad_high, ct_order, cell_colors) {
  df <- reshape2::melt(mat, varnames = c("source", "target"), value.name = "count")
  df$source <- factor(df$source, levels = rev(ct_order))
  df$target <- factor(df$target, levels = ct_order)

  ggplot(df, aes(x = target, y = source, fill = count)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(count > 0, count, "")),
              size = 3.5, color = "grey20") +
    scale_fill_gradientn(
      colors = c("#F5F5F5", grad_high),
      name   = "Interactions",
      limits = c(0, max(mat))
    ) +
    scale_x_discrete(position = "top") +
    labs(title = title, x = "Target", y = "Source") +
    theme_classic(base_size = 13) +
    theme(
      plot.title   = element_text(size = 14, face = "bold"),
      axis.text.x  = element_text(size = 11, angle = 45, hjust = 0, vjust = 0,
                                  face = "bold", color = cell_colors[ct_order]),
      axis.text.y  = element_text(size = 11, face = "bold",
                                  color = cell_colors[rev(ct_order)]),
      axis.title   = element_text(size = 12),
      legend.text  = element_text(size = 11),
      legend.title = element_text(size = 11),
      axis.line    = element_blank(),
      axis.ticks   = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5)
    )
}
```

### `reorder_bubble_by_ec()` — sort LR pairs by dominant source cell type

```r
reorder_bubble_by_ec <- function(p, source_types, direction = "source") {
  d <- p$data
  if (direction == "source" && "source" %in% colnames(d)) {
    ec_col_vals <- as.character(d$source)
  } else if (direction == "target" && "target" %in% colnames(d)) {
    ec_col_vals <- as.character(d$target)
  } else {
    parts <- strsplit(as.character(d$source.target), " --> | -> | \\| ")
    ec_col_vals <- if (direction == "source")
      trimws(sapply(parts, `[`, 1)) else
      trimws(sapply(parts, function(x) x[length(x)]))
  }
  d$ec_col <- ec_col_vals
  d_ec <- d[d$ec_col %in% source_types & !is.na(d$prob) & d$prob > 0, ]
  if (nrow(d_ec) == 0) return(p)
  lr_names <- unique(as.character(d_ec$interaction_name_2))
  best_df <- do.call(rbind, lapply(lr_names, function(lr) {
    sub_d    <- d_ec[d_ec$interaction_name_2 == lr, ]
    best_row <- sub_d[which.max(sub_d$prob), ]
    data.frame(lr = lr, best_ec = as.character(best_row$ec_col),
               best_prob = best_row$prob, stringsAsFactors = FALSE)
  }))
  best_df$ec_rank <- match(best_df$best_ec, source_types)
  best_df <- best_df[order(best_df$ec_rank, -best_df$best_prob), ]
  p$data$interaction_name_2 <- factor(p$data$interaction_name_2,
                                       levels = rev(best_df$lr))
  p
}
```

### `rebuild_bubble()` — custom ggplot from CellChat bubble data

Rebuilds from `p$data` rather than modifying CellChat's internal scales via `+`.
CellChat's `netVisual_bubble()` hard-codes color and size scales; adding a second
scale causes silent conflicts. Call `netVisual_bubble()` to get the data, then use
this function to build a clean plot.

```r
rebuild_bubble <- function(p_raw, source_types, source_colors,
                           direction, section_label, show_x = TRUE) {
  d <- p_raw$data
  x_levels <- levels(factor(d$source.target))
  x_colors <- sapply(x_levels, function(st) {
    parts   <- trimws(strsplit(st, " --> | -> | \\| | - ")[[1]])
    ec_part <- if (direction == "source") parts[1] else parts[length(parts)]
    if (ec_part %in% names(source_colors)) source_colors[ec_part] else "grey40"
  })
  d$source.target      <- factor(d$source.target, levels = x_levels)
  d$interaction_name_2 <- factor(d$interaction_name_2,
                                  levels = levels(p_raw$data$interaction_name_2))
  raw_size <- -log10(d$pval + 1e-10)
  raw_size <- pmin(raw_size, quantile(raw_size, 0.95, na.rm = TRUE))
  d$size_val <- (raw_size - min(raw_size, na.rm = TRUE)) /
                (max(raw_size, na.rm = TRUE) - min(raw_size, na.rm = TRUE) + 1e-10)

  p <- ggplot(d, aes(x = source.target, y = interaction_name_2)) +
    geom_point(aes(size = size_val, color = prob), alpha = 1) +
    scale_color_gradientn(
      colors = c("#BDBDBD", "#EF9A9A", "#E53935", "#B71C1C"),
      name   = "Communication\nProbability",
      guide  = guide_colorbar(barheight = 4, barwidth = 0.8,
                              title.position = "top", title.hjust = 0.5)
    ) +
    scale_size_continuous(
      name   = "-log10\n(p-value)",
      range  = c(2.5, 8),
      guide  = guide_legend(override.aes = list(color = "grey60"),
                            title.position = "top", title.hjust = 0.5)
    ) +
    scale_x_discrete(position = "bottom") +
    labs(title = section_label, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(size = 12, color = "grey10"),
      axis.title       = element_blank(),
      plot.title       = element_text(size = 15, face = "bold.italic",
                                      hjust = 0, color = "grey35",
                                      margin = margin(b = 4)),
      legend.title     = element_text(size = 12),
      legend.text      = element_text(size = 11)
    )

  if (show_x) {
    p <- p + theme(
      axis.text.x  = element_text(size = 11, angle = 40, hjust = 1, vjust = 1,
                                   color = x_colors),
      plot.margin  = margin(t = 5, r = 10, b = 45, l = 10)
    ) + coord_cartesian(clip = "off")
  } else {
    p <- p + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(t = 5, r = 10, b = 5, l = 10)
    )
  }
  p
}
```

### `assign_category()`, `make_bar_data()`, `make_bar_plot()`, `save_both()` — bar plot helpers

```r
assign_category <- function(pathway, pathway_categories) {
  for (cat in names(pathway_categories)) {
    if (pathway %in% pathway_categories[[cat]]$paths) return(cat)
  }
  return("Other")
}

make_bar_data <- function(df, source_types, partner_types,
                          source_as = "sender",  # "sender" or "receiver"
                          pathway_categories) {
  if (source_as == "sender") {
    d <- df[df$source %in% source_types & df$target %in% partner_types, ]
    d$source_ct <- d$source
    d$partner   <- d$target
  } else {
    d <- df[df$target %in% source_types & df$source %in% partner_types, ]
    d$source_ct <- d$target
    d$partner   <- d$source
  }
  d$category <- vapply(d$pathway_name, assign_category,
                       FUN.VALUE = character(1),
                       pathway_categories = pathway_categories)
  d %>%
    dplyr::group_by(source_ct, partner, category) %>%
    dplyr::summarise(n = dplyr::n(), strength = sum(prob), .groups = "drop") %>%
    dplyr::mutate(source_ct = factor(source_ct, levels = source_types))
}

make_bar_plot <- function(bar_data, title_label, cat_colors,
                          source_colors, metric = "n",
                          show_legend = TRUE) {
  y_var  <- if (metric == "n") "n" else "strength"
  y_lab  <- if (metric == "n") "Interaction count" else "Communication strength (prob sum)"
  x_cols <- source_colors[levels(bar_data$source_ct)]
  x_cols[is.na(x_cols)] <- "grey40"

  ggplot(bar_data, aes(x = source_ct, y = .data[[y_var]], fill = category)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = cat_colors, name = "Pathway\nCategory", drop = FALSE) +
    facet_wrap(~partner, nrow = 1, scales = "fixed") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = title_label, x = NULL, y = y_lab) +
    theme_minimal(base_size = 13) +
    theme(
      panel.background   = element_rect(fill = "white", color = NA),
      panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey96", color = "grey80", linewidth = 0.4),
      strip.text         = element_text(size = 12, face = "bold", color = "grey25"),
      axis.text.x        = element_text(size = 11, angle = 40, hjust = 1, vjust = 1,
                                        color = x_cols),
      axis.text.y        = element_text(size = 11, color = "grey30"),
      axis.title.y       = element_text(size = 12, color = "grey30"),
      plot.title         = element_text(size = 14, face = "bold.italic",
                                        hjust = 0, color = "grey30",
                                        margin = margin(b = 6)),
      legend.title       = element_text(size = 12, face = "bold"),
      legend.text        = element_text(size = 11),
      legend.key.size    = unit(0.9, "lines"),
      plot.margin        = margin(t = 5, r = 10, b = 40, l = 10)
    ) +
    coord_cartesian(clip = "off") +
    if (!show_legend) theme(legend.position = "none") else NULL
}

save_both <- function(send_data, recv_data, base_name, output_dir, width,
                      cat_colors, source_colors) {
  for (metric in c("n", "strength")) {
    suffix <- if (metric == "n") "count" else "strength"
    p <- make_bar_plot(send_data, "Source cell types as Senders",
                       cat_colors, source_colors, metric, show_legend = FALSE) /
         make_bar_plot(recv_data, "Source cell types as Receivers",
                       cat_colors, source_colors, metric, show_legend = TRUE)
    ggsave(file.path(output_dir, paste0(base_name, "_", suffix, ".pdf")),
           plot = p, width = width, height = 10,
           device = "pdf", useDingbats = FALSE)
  }
}
```

### `make_circos()` — circos chord diagram (exception #3 from CONVENTIONS.md §4)

Uses `cairo_pdf()` (a variant of the pdf()/dev.off() exception #3 pattern) because
cairo embeds text as proper text objects editable in Adobe Illustrator. Standard
`pdf()` bakes sector labels into paths.

```r
make_circos <- function(df, sources, targets,
                        cell_colors,
                        all_cat_paths,
                        output_dir,
                        file_name,
                        title        = "",
                        top_n        = 10,
                        min_prob     = 1e-4,
                        pval_thresh  = 0.05,
                        source_order = NULL,
                        target_order = NULL,
                        link_alpha   = 0.80) {
  qual_pal <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#A65628", "#F781BF", "#4DBEEE", "#BCBD22", "#17BECF",
    "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5"
  )

  d <- df[df$source %in% sources & df$target %in% targets &
          df$pathway_name %in% all_cat_paths &
          df$prob > min_prob & df$pval <= pval_thresh, ]
  if (nrow(d) == 0) {
    cat("No interactions for:", file_name, "\n")
    return(invisible(NULL))
  }

  path_counts <- sort(table(d$pathway_name), decreasing = TRUE)
  top_paths   <- names(path_counts)[seq_len(min(top_n, length(path_counts)))]
  d <- d[d$pathway_name %in% top_paths, ]

  path_cols <- stats::setNames(qual_pal[seq_along(top_paths)], top_paths)

  agg <- stats::aggregate(prob ~ source + target + pathway_name, data = d, FUN = sum)
  agg$link_col <- grDevices::adjustcolor(path_cols[agg$pathway_name], alpha.f = link_alpha)
  agg <- agg[order(match(agg$pathway_name, top_paths)), ]

  src_present  <- intersect(if (!is.null(source_order)) source_order else sources,
                             unique(agg$source))
  tgt_present  <- intersect(if (!is.null(target_order)) target_order else targets,
                             unique(agg$target))
  sector_order <- c(src_present, tgt_present)
  n_src <- length(src_present)
  n_tgt <- length(tgt_present)

  all_sectors <- unique(c(agg$source, agg$target))
  grid_col <- stats::setNames(
    sapply(all_sectors, function(s) {
      if (s %in% names(cell_colors)) cell_colors[s] else "grey70"
    }),
    all_sectors)

  legend_labels <- paste0(top_paths, "  (n=", path_counts[top_paths], ")")
  legend_cols   <- grDevices::adjustcolor(path_cols[top_paths], alpha.f = link_alpha)

  path_out <- file.path(output_dir, file_name)
  grDevices::cairo_pdf(path_out, width = 11, height = 13)  # cairo for Illustrator text

  graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(5.5, 1))

  graphics::par(mar = c(1, 2, 3, 2))
  circlize::circos.clear()
  circlize::circos.par(
    canvas.xlim              = c(-1.35, 1.35),
    canvas.ylim              = c(-1.35, 1.35),
    gap.after                = c(rep(2, n_src - 1), 12, rep(2, n_tgt - 1), 12),
    start.degree             = 90,
    clock.wise               = TRUE,
    track.margin             = c(0.005, 0.005),
    points.overflow.warning  = FALSE
  )
  circlize::chordDiagram(
    agg[, c("source", "target", "prob")],
    order           = sector_order,
    grid.col        = grid_col,
    col             = agg$link_col,
    directional     = 1,
    direction.type  = c("diffHeight", "arrows"),
    link.arr.type   = "big.arrow",
    link.arr.length = 0.12,
    link.arr.width  = 0.08,
    link.sort       = TRUE,
    annotationTrack = "grid",
    preAllocateTracks = list(list(track.height = circlize::mm_h(5)))
  )
  circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim   <- circlize::get.cell.meta.data("xlim")
    ylim   <- circlize::get.cell.meta.data("ylim")
    sector <- circlize::get.cell.meta.data("sector.index")
    circlize::circos.text(mean(xlim), ylim[1] + 0.5, sector,
                          facing = "clockwise", niceFacing = TRUE,
                          adj = c(0, 0.5), cex = 0.78, col = "grey10")
  }, bg.border = NA)
  if (nchar(title) > 0)
    graphics::title(title, cex.main = 1.25, col.main = "grey10", font.main = 1, line = 0.5)
  circlize::circos.clear()

  graphics::par(mar = c(0, 1, 0, 1))
  graphics::plot.new()
  graphics::legend("center",
    legend = legend_labels, fill = legend_cols, border = "grey60",
    bty = "n", cex = 0.88, ncol = ceiling(length(top_paths) / 2),
    title = "Pathway  (# sig. interactions)", title.font = 2,
    x.intersp = 0.6, y.intersp = 1.0)

  grDevices::dev.off()
  cat("Saved:", path_out, "\n  top pathways:", paste(top_paths, collapse = ", "), "\n")
}
```

---

## Pathway Categorization Strategy

Organize all significant pathways into 2–4 biological categories before writing Scripts 2–5.
This is the single most important structural decision — it makes outputs interpretable.

**How to select pathways:**

1. Load `CellChat_all_interactions.csv` (written by Script 1)
2. Compute `total_prob` per pathway: `aggregate(prob ~ pathway_name, df, sum)`
3. Drop pathways with `total_prob < ~0.0005` AND `≤ 10 LR pairs` (noise)
4. Drop single-LR-pair pathways (uninformative for chord diagrams)
5. Group remaining pathways by biological theme (ECM, adhesion, growth factor, etc.)
6. For chord diagrams: if one pathway's total_prob is >3–4× the next strongest in its
   category, **exclude it from chord only** (keep in bubble/heatmap/circle). This prevents
   one dominant pathway from crushing width scaling of all others.

**CHECKPOINT after Script 1:** Report total pathways and interaction count; top 15 pathways
by total_prob. Ask user to review and confirm category assignments before writing Script 2.

---

## Script 1: Core Inference (`01_cellchat_inference.R`)

```r
library(future)
plan("sequential")  # MUST be first — before any CellChat or parallel library loads
library(CellChat)
library(Seurat)

# ── CONFIG ────────────────────────────────────────────────────────────────────
LABEL_COL  <- "project_specific"  # REPLACE: cell type column name in Seurat metadata
ORGANISM   <- "human"             # "human" or "mouse"
INPUT_RDS  <- "project_specific"  # REPLACE: path to annotated Seurat object
OUTPUT_DIR <- file.path("output", "cellchat")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
# ─────────────────────────────────────────────────────────────────────────────

so <- readRDS(INPUT_RDS)

# Seurat v5: join layers before extracting counts (seurat_v5_rules Rule 1)
DefaultAssay(so) <- "RNA"
so <- JoinLayers(so)

data_input <- GetAssayData(so, assay = "RNA", layer = "data")  # log-normalized
meta       <- so@meta.data[, LABEL_COL, drop = FALSE]

cellchat <- createCellChat(object = data_input, meta = meta, group.by = LABEL_COL)

# Select database: CellChatDB.human or CellChatDB.mouse
CellChatDB <- if (ORGANISM == "human") CellChatDB.human else CellChatDB.mouse
cellchat@DB <- CellChatDB

# NOTE: projectData(cellchat, PPI.human) is NOT called here.
# This function was removed in CellChat v2. The v1 library erroneously included
# this call in its inference skeleton. Any v1 script that included projectData()
# must have that line removed before running on a v2 install.

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file.path(OUTPUT_DIR, "cellchat_object.Rds"))

# Write interaction table for pathway selection decisions
df_net <- subsetCommunication(cellchat)
write.csv(df_net, file.path(OUTPUT_DIR, "CellChat_all_interactions.csv"),
          row.names = FALSE)

cat("Total pathways:", length(cellchat@netP$pathways), "\n")
cat("Top 15 by prob:\n")
pw_total <- aggregate(prob ~ pathway_name, df_net, sum)
print(head(pw_total[order(-pw_total$prob), ], 15))
```

---

## Script 2: Basic Plots (`02_cellchat_plots.R`)

```r
library(future)
plan("sequential")  # MUST be first
library(CellChat)
library(ggplot2)
library(patchwork)
library(reshape2)

# ── CONFIG ────────────────────────────────────────────────────────────────────
OUTPUT_DIR <- file.path("output", "cellchat")

# Pathway categories (define after reviewing Script 1 checkpoint output)
PATHWAY_CATEGORIES <- list(
  "CategoryA" = list(
    label       = "project_specific",  # display name, e.g. "Angiocrine & Growth Factor"
    grad_high   = "project_specific",  # hex color for heatmap gradient high end
    paths       = c("project_specific"),  # pathway names in this category
    chord_paths = c("project_specific")   # subset for chord (exclude dominant outliers)
  )
  # Add more categories as needed
)

# Cell type vectors
CT_ORDER     <- "project_specific"  # REPLACE: ordered character vector of all cell type names
SOURCE_TYPES <- "project_specific"  # REPLACE: cell types used as sources (e.g., EC subtypes)
TARGET_TYPES <- "project_specific"  # REPLACE: cell types used as targets

# Color vectors (must cover all cell types in CT_ORDER)
CELL_COLORS   <- c("project_specific" = "#000000")  # REPLACE with full palette
SOURCE_COLORS <- CELL_COLORS[SOURCE_TYPES]  # auto-derived; override if needed
# ─────────────────────────────────────────────────────────────────────────────

# Paste module-level helpers here (subset_net, get_filtered_net, plot_heatmap)

cellchat   <- readRDS(file.path(OUTPUT_DIR, "cellchat_object.Rds"))
cellchat   <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
group_size <- as.numeric(table(cellchat@idents))

# ── Global plots ──────────────────────────────────────────────────────────────

# All-pathway circle: DO NOT use netVisual_aggregate(signaling=NULL) — it crashes.
# Use netVisual_circle() on the pre-aggregated count matrix instead.
pdf(file.path(OUTPUT_DIR, "CellChat_global_circle.pdf"),
    width = 8, height = 8, useDingbats = FALSE)
par(mar = c(0, 0, 2, 0))
netVisual_circle(cellchat@net$count, vertex.weight = group_size,
                 weight.scale = TRUE, label.edge = FALSE,
                 color.use = CELL_COLORS, title.name = "All Interactions (count)")
dev.off()

# Pathway strength ranking
p_rank <- rankNet(cellchat, mode = "single", stacked = FALSE, do.stat = TRUE)
ggsave(file.path(OUTPUT_DIR, "CellChat_pathway_ranking.pdf"), plot = p_rank,
       width = 8, height = max(6, length(cellchat@netP$pathways) * 0.18 + 2),
       device = "pdf", useDingbats = FALSE)

# Signaling role scatter
p_scatter <- netAnalysis_signalingRole_scatter(cellchat, color.use = CELL_COLORS,
                                               title = "Signaling Role")
ggsave(file.path(OUTPUT_DIR, "CellChat_signaling_role_scatter.pdf"), plot = p_scatter,
       width = 7, height = 6, device = "pdf", useDingbats = FALSE)

# Outgoing / incoming heatmaps (cell type x pathway)
p_out <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",
                                           color.use = CELL_COLORS,
                                           title = "Outgoing Signaling")
p_in  <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                           color.use = CELL_COLORS,
                                           title = "Incoming Signaling")
ggsave(file.path(OUTPUT_DIR, "CellChat_signaling_outgoing.pdf"), plot = p_out,
       width = 10, height = 7, device = "pdf", useDingbats = FALSE)
ggsave(file.path(OUTPUT_DIR, "CellChat_signaling_incoming.pdf"), plot = p_in,
       width = 10, height = 7, device = "pdf", useDingbats = FALSE)

# ── Per-category plots ────────────────────────────────────────────────────────

for (cat_id in names(PATHWAY_CATEGORIES)) {
  ci         <- PATHWAY_CATEGORIES[[cat_id]]
  cat_label  <- ci$label
  cat_paths  <- ci$paths
  chord_use  <- ci$chord_paths
  prefix     <- file.path(OUTPUT_DIR, paste0("CellChat_", cat_id))

  # 1. Bubble — source to target and back
  p_bub_out <- netVisual_bubble(cellchat,
    sources.use = SOURCE_TYPES, targets.use = TARGET_TYPES,
    signaling = cat_paths, remove.isolate = TRUE,
    title.name = paste0(cat_label, ": Sources to Targets"))
  bub_h <- max(6, length(cat_paths) * 0.55 + 3)
  ggsave(paste0(prefix, "_bubble_out.pdf"), plot = p_bub_out,
         width = 10, height = bub_h, device = "pdf", useDingbats = FALSE)

  p_bub_in <- netVisual_bubble(cellchat,
    sources.use = TARGET_TYPES, targets.use = SOURCE_TYPES,
    signaling = cat_paths, remove.isolate = TRUE,
    title.name = paste0(cat_label, ": Targets to Sources"))
  ggsave(paste0(prefix, "_bubble_in.pdf"), plot = p_bub_in,
         width = 10, height = bub_h, device = "pdf", useDingbats = FALSE)

  # 2. Circle — interaction count, all cell types
  net <- subset_net(cellchat, cat_paths)
  if (!is.null(net)) {
    pdf(paste0(prefix, "_circle.pdf"), width = 8, height = 8, useDingbats = FALSE)
    par(mar = c(0, 0, 2, 0))
    netVisual_circle(net$count, vertex.weight = group_size,
                     weight.scale = TRUE, label.edge = FALSE,
                     color.use = CELL_COLORS,
                     title.name = paste0(cat_label, " (interactions)"))
    dev.off()

    # 3. Chord — gene-level, source senders only
    net_chord <- get_filtered_net(cellchat, chord_use,
                                   sources = SOURCE_TYPES,
                                   targets = TARGET_TYPES, top_n = 2)
    if (nrow(net_chord) > 0) {
      pdf(paste0(prefix, "_chord.pdf"), width = 10, height = 10, useDingbats = FALSE)
      par(mar = c(1, 1, 3, 1))
      netVisual_chord_gene(cellchat,
        net         = net_chord,
        sources.use = SOURCE_TYPES,
        targets.use = TARGET_TYPES,
        color.use   = CELL_COLORS,
        lab.cex     = 0.9,
        big.gap     = 20,
        small.gap   = 3,
        title.name  = paste0(cat_label, ": Source Outgoing"))
      dev.off()
    }

    # 4. Heatmap (custom ggplot — replaces buggy netVisual_heatmap with pathway subset)
    p_ht <- plot_heatmap(net$count,
                         title       = paste0(cat_label, " - Interaction Count"),
                         grad_high   = ci$grad_high,
                         ct_order    = CT_ORDER,
                         cell_colors = CELL_COLORS)
    ggsave(paste0(prefix, "_heatmap.pdf"), plot = p_ht,
           width = 7, height = 6, device = "pdf", useDingbats = FALSE)
  }
}
```

---

## Script 3: Stacked Composite Bubble Plots (`03_cellchat_bubbles.R`)

Produces a multi-panel stacked figure with one panel per signaling category.
Rebuilds from `netVisual_bubble()`'s `p$data` to avoid CellChat's hard-coded
internal scales conflicting with custom color/size aesthetics.

```r
library(future)
plan("sequential")  # MUST be first
library(CellChat)
library(ggplot2)
library(patchwork)

# ── CONFIG ────────────────────────────────────────────────────────────────────
OUTPUT_DIR <- file.path("output", "cellchat")

# Pathway categories (same as Script 2)
PATHWAY_CATEGORIES <- list(...)  # project_specific — copy from Script 2 CONFIG

CT_ORDER     <- "project_specific"  # REPLACE
SOURCE_TYPES <- "project_specific"  # REPLACE
TARGET_TYPES <- "project_specific"  # REPLACE
CELL_COLORS  <- c("project_specific" = "#000000")  # REPLACE
SOURCE_COLORS <- CELL_COLORS[SOURCE_TYPES]

# Categories to include in the stacked figure (exclude if dominated by one pathway)
STACKED_CATS <- names(PATHWAY_CATEGORIES)  # or subset, e.g. c("CategoryA", "CategoryB")
# ─────────────────────────────────────────────────────────────────────────────

# Paste module-level helpers here (reorder_bubble_by_ec, rebuild_bubble)

cellchat <- readRDS(file.path(OUTPUT_DIR, "cellchat_object.Rds"))

build_stacked <- function(direction) {
  p_list  <- list()
  heights <- c()

  for (i in seq_along(STACKED_CATS)) {
    cat_id    <- STACKED_CATS[i]
    ci        <- PATHWAY_CATEGORIES[[cat_id]]
    show_x    <- (i == length(STACKED_CATS))

    p_raw <- netVisual_bubble(cellchat,
      sources.use    = if (direction == "source") SOURCE_TYPES else TARGET_TYPES,
      targets.use    = if (direction == "source") TARGET_TYPES else SOURCE_TYPES,
      signaling      = ci$paths,
      remove.isolate = TRUE)

    if (is.null(p_raw) || nrow(p_raw$data) == 0) next

    p_raw  <- reorder_bubble_by_ec(p_raw, source_types = SOURCE_TYPES, direction)
    n_lr   <- length(unique(p_raw$data$interaction_name_2))
    p_list <- c(p_list, list(
      rebuild_bubble(p_raw, SOURCE_TYPES, SOURCE_COLORS, direction, ci$label, show_x)
    ))
    heights <- c(heights, n_lr)
  }

  combined <- Reduce(function(a, b) a / b, p_list) +
    plot_layout(heights = heights, guides = "collect") &
    theme(legend.position = "right")

  total_h <- max(12, sum(heights) * 0.3 + length(p_list) * 2.5)
  width   <- max(10, length(SOURCE_TYPES) * length(TARGET_TYPES) * 0.9)
  label   <- if (direction == "source") "sources_to_targets" else "targets_to_sources"
  ggsave(file.path(OUTPUT_DIR, paste0("CellChat_stacked_", label, ".pdf")),
         plot = combined, width = width, height = total_h,
         device = "pdf", useDingbats = FALSE)
}

build_stacked("source")
build_stacked("target")
```

**Note:** Do NOT use `& coord_cartesian(clip = "off")` on the patchwork when x-axis
is at the bottom — it causes a patchwork theme warning. Only use `clip = "off"` for
individual panels with `show_x = TRUE` (handled inside `rebuild_bubble()`).

### Receiver x-axis column ordering

When the target population is biologically mixed (e.g., multiple partner types),
sort x-axis columns by source cell type rank so all columns from the same source
appear together:

```r
if (direction == "target") {
  ec_parts <- sapply(x_levels, function(st) {
    parts <- trimws(strsplit(st, " --> | -> | \\| | - ")[[1]])
    parts[length(parts)]
  })
  ec_rank <- match(ec_parts, SOURCE_TYPES)
  ec_rank[is.na(ec_rank)] <- 999
  x_levels <- x_levels[order(ec_rank)]
}
```

### Trimming global plots to a probability cutoff

To trim `rankNet` or `netAnalysis_signalingRole_heatmap` to specific pathways:

```r
keep_pws <- c("project_specific")  # REPLACE with pathway names to keep
keep_idx <- cellchat@netP$pathways %in% keep_pws
cc_sub <- cellchat
cc_sub@netP$pathways <- cellchat@netP$pathways[keep_idx]
cc_sub@netP$prob     <- cellchat@netP$prob[, , keep_idx, drop = FALSE]
rankNet(cc_sub, mode = "single", stacked = FALSE, do.stat = TRUE)
```

---

## Script 4: Pathway Category Bar Plots (`04_cellchat_bars.R`)

Stacked bar plots showing interaction count and total communication probability
per source cell type × partner cell type, broken down by pathway category.
Fixed y-axis scale across all panels (`scales = "fixed"`).

```r
library(future)
plan("sequential")  # MUST be first
library(ggplot2)
library(dplyr)
library(patchwork)

# ── CONFIG ────────────────────────────────────────────────────────────────────
OUTPUT_DIR <- file.path("output", "cellchat")

PATHWAY_CATEGORIES <- list(...)  # project_specific — copy from Script 2 CONFIG

SOURCE_TYPES <- "project_specific"  # REPLACE
TARGET_TYPES <- "project_specific"  # REPLACE
CELL_COLORS  <- c("project_specific" = "#000000")  # REPLACE
SOURCE_COLORS <- CELL_COLORS[SOURCE_TYPES]

# Category colors (project_specific — choose colors that match your biological categories)
CAT_COLORS <- c(
  "CategoryA" = "project_specific",  # REPLACE with one hex per category
  "Other"     = "#BDBDBD"
)
# ─────────────────────────────────────────────────────────────────────────────

# Paste module-level helpers here (assign_category, make_bar_data, make_bar_plot, save_both)

df <- read.csv(file.path(OUTPUT_DIR, "CellChat_all_interactions.csv"),
               stringsAsFactors = FALSE)
df <- df[df$prob > 0 & df$pval < 0.05, ]
df$category <- vapply(df$pathway_name, assign_category,
                      FUN.VALUE = character(1),
                      pathway_categories = PATHWAY_CATEGORIES)
df$category <- factor(df$category, levels = names(CAT_COLORS))

send_data <- make_bar_data(df, SOURCE_TYPES, TARGET_TYPES, "sender", PATHWAY_CATEGORIES)
recv_data <- make_bar_data(df, SOURCE_TYPES, TARGET_TYPES, "receiver", PATHWAY_CATEGORIES)

bar_width <- max(8, length(SOURCE_TYPES) * length(TARGET_TYPES) * 1.2)
save_both(send_data, recv_data, "pathway_dist_all",
          OUTPUT_DIR, bar_width, CAT_COLORS, SOURCE_COLORS)
```

---

## Script 5: Pathway Circos Plots (`05_cellchat_circos.R`)

**CONVENTIONS.md §4 Exception #3** applies: `circlize::chordDiagram()` uses base R
graphics with no ggplot/ggsave path. This script uses `cairo_pdf()` + `dev.off()`.
`cairo_pdf()` is used instead of `pdf()` to produce Illustrator-editable text objects.

```r
library(future)
plan("sequential")  # MUST be first
library(circlize)
library(dplyr)

# ── CONFIG ────────────────────────────────────────────────────────────────────
OUTPUT_DIR <- file.path("output", "cellchat")

PATHWAY_CATEGORIES <- list(...)  # project_specific — copy from Script 2 CONFIG

SOURCE_TYPES <- "project_specific"  # REPLACE
TARGET_TYPES <- "project_specific"  # REPLACE
CELL_COLORS  <- c("project_specific" = "#000000")  # REPLACE

# All pathway names across all categories (passed to make_circos for curated filtering)
ALL_CAT_PATHS <- unlist(lapply(PATHWAY_CATEGORIES, `[[`, "paths"))

# Preferred sector ordering (project_specific — usually match UMAP or biological order)
SOURCE_ORDER <- "project_specific"  # REPLACE: ordered character vector of source types
TARGET_ORDER <- "project_specific"  # REPLACE: ordered character vector of target types
# ─────────────────────────────────────────────────────────────────────────────

# Paste make_circos() helper here

df <- read.csv(file.path(OUTPUT_DIR, "CellChat_all_interactions.csv"),
               stringsAsFactors = FALSE)

# Sources-to-targets
make_circos(df, SOURCE_TYPES, TARGET_TYPES, CELL_COLORS, ALL_CAT_PATHS,
            OUTPUT_DIR, "circos_sources_to_targets_bypathway.pdf",
            title = "Sources -> Targets",
            source_order = SOURCE_ORDER, target_order = TARGET_ORDER)

# Targets-to-sources
make_circos(df, TARGET_TYPES, SOURCE_TYPES, CELL_COLORS, ALL_CAT_PATHS,
            OUTPUT_DIR, "circos_targets_to_sources_bypathway.pdf",
            title = "Targets -> Sources",
            source_order = TARGET_ORDER, target_order = SOURCE_ORDER)
```

### Common circos adjustments

| Situation | Fix |
|---|---|
| Legend still crowded | Reduce `top_n` to 8; or increase Panel 2 height in `layout()` |
| Important pathway not shown | Lower `top_n` OR manually add it to `top_paths` after ranking |
| Sector labels still clipping | Widen `canvas.xlim/ylim` to `c(-1.4, 1.4)` |
| Chords too faint | Increase `link_alpha` toward 1.0 |
| Too many small chords cluttering | Raise `min_prob` to `1e-3` to suppress weak interactions |

---

## Script 6: Tissue-Type Comparison (`06_cellchat_comparison.R`) [Conditional]

**Run only when `group_col` is set in the brief.** Runs separate CellChat inference per
condition/tissue, then compares objects using CellChat's differential framework.

### Memory-efficient tissue extraction helper

```r
extract_tissue <- function(obj_path, group_val, group_col, label_col, prefix) {
  cat(sprintf("Loading %s for %s = %s ...\n", obj_path, group_col, group_val))
  obj   <- readRDS(obj_path)
  cells <- colnames(obj)[obj[[group_col, drop = TRUE]] == group_val]
  cat(sprintf("  %d cells\n", length(cells)))
  sub_obj <- subset(obj, cells = cells)
  sub_obj <- JoinLayers(sub_obj)  # Seurat v5: join before GetAssayData
  mat  <- GetAssayData(sub_obj, assay = "RNA", layer = "data")
  meta <- sub_obj@meta.data[, label_col, drop = FALSE]
  rm(obj, sub_obj); gc()
  colnames(mat)  <- paste0(prefix, "_", colnames(mat))
  rownames(meta) <- paste0(prefix, "_", rownames(meta))
  list(mat = mat, meta = meta)
}
```

### Per-group inference with caching

```r
run_cellchat_from_parts <- function(mat, meta, label_col, organism, out_rds) {
  if (file.exists(out_rds)) {
    cat("Loading cached:", out_rds, "\n")
    return(readRDS(out_rds))
  }
  db <- if (organism == "human") CellChatDB.human else CellChatDB.mouse
  cc <- createCellChat(object = mat, meta = meta, group.by = label_col)
  cc@DB <- subsetDB(db)
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  saveRDS(cc, out_rds)
  cat("Saved:", out_rds, "\n")
  cc
}
```

### Comparison workflow

**Order of operations is strict** (violating it causes known CellChat v2 bugs):

```r
compare_pair <- function(cc_list, names_vec, col_use = c("#2166AC", "#B2182B")) {
  # 1. Find common cell types — objects may differ after filtering
  common_cts <- Reduce(intersect, lapply(cc_list, function(cc) levels(cc@idents)))

  # 2. Subset EACH object to common cell types
  cc_list <- lapply(cc_list, function(cc) {
    cc <- subsetCellChat(cc, idents.use = common_cts)
    # 3. Compute centrality on INDIVIDUAL objects — NEVER on the merged object
    #    (running on merged throws a rowSums error — known CellChat v2 bug)
    cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
    cc
  })

  # 4. Merge — do NOT pass cell.level = TRUE (removed in CellChat v2)
  mergeCellChat(cc_list, add.names = names_vec)
}
```

### Saving comparison plots

Different CellChat comparison functions return different object types — each requires
a different save pattern:

```r
# compareInteractions → ggplot → ggsave
g_int <- compareInteractions(cc_merged, show.legend = FALSE, group = c(1, 2),
                              color.use = COL_USE)
ggsave(file.path(OUTPUT_DIR, "compare_interactions.pdf"),
       g_int, width = 7, height = 5, device = cairo_pdf)

# netVisual_heatmap → ComplexHeatmap → draw() inside pdf()/dev.off()
# CONVENTIONS.md §4 Exception #1: ComplexHeatmap::draw() sends output to active device
# Note: here we use the merged-object netVisual_heatmap (no pathway subset) — the
# annotation bug does NOT apply to this merged comparison context.
ht1 <- netVisual_heatmap(cc_merged,
  color.heatmap = c("#2166AC", "white", "#B2182B"))
ht2 <- netVisual_heatmap(cc_merged, measure = "weight",
  color.heatmap = c("#2166AC", "white", "#B2182B"))
pdf(file.path(OUTPUT_DIR, "compare_heatmap.pdf"), width = 14, height = 6,
    useDingbats = FALSE)
ComplexHeatmap::draw(ht1 + ht2)  # draw() is mandatory — auto-print does NOT work
dev.off()

# rankNet → ggplot → ggsave
p_rank <- rankNet(cc_merged, mode = "comparison", stacked = TRUE,
                  do.stat = TRUE, color.use = COL_USE)
ggsave(file.path(OUTPUT_DIR, "compare_pathways.pdf"),
       p_rank, width = 10, height = 9, device = cairo_pdf)

# Signaling role scatter → ggplot per individual object → ggsave
g_sc <- lapply(seq_along(cc_list), function(i) {
  netAnalysis_signalingRole_scatter(cc_list[[i]]) + ggtitle(names_vec[i])
})
ggsave(file.path(OUTPUT_DIR, "compare_scatter.pdf"),
       Reduce(`+`, g_sc), width = 13, height = 6, device = cairo_pdf)
```

### Top-differential pathway filtering

Pre-filter `rankNet` to the N most differentially expressed pathways:

```r
pathway_stats <- function(cc) {
  pathways <- cc@netP$pathways
  flow <- sapply(pathways, function(pw) {
    if (pw %in% dimnames(cc@netP$prob)[[3]])
      sum(cc@netP$prob[, , pw], na.rm = TRUE)
    else 0
  })
  n_int <- sapply(pathways, function(pw) {
    df <- subsetCommunication(cc, slot.name = "net")
    sum(df$pathway_name == pw, na.rm = TRUE)
  })
  data.frame(pathway = pathways, flow = flow, n_int = n_int, row.names = pathways)
}

top_diff_pathways <- function(cc1, cc2, top_n = 25, min_int = 3) {
  s1 <- pathway_stats(cc1);  s2 <- pathway_stats(cc2)
  all_pw <- union(rownames(s1), rownames(s2))
  flow1  <- stats::setNames(s1$flow[match(all_pw, rownames(s1))], all_pw)
  flow2  <- stats::setNames(s2$flow[match(all_pw, rownames(s2))], all_pw)
  nint1  <- stats::setNames(s1$n_int[match(all_pw, rownames(s1))], all_pw)
  nint2  <- stats::setNames(s2$n_int[match(all_pw, rownames(s2))], all_pw)
  flow1[is.na(flow1)] <- 0;  flow2[is.na(flow2)] <- 0
  nint1[is.na(nint1)] <- 0;  nint2[is.na(nint2)] <- 0
  data.frame(pathway = all_pw, flow1 = flow1, flow2 = flow2,
             abs_diff = abs(flow1 - flow2), total_int = nint1 + nint2) %>%
    dplyr::filter(total_int >= min_int) %>%
    dplyr::arrange(dplyr::desc(abs_diff)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(pathway)
}
```

### Y-axis label coloring by dominant condition

```r
top_pw    <- top_diff_pathways(cc_list[[1]], cc_list[[2]])
cc_merged <- mergeCellChat(cc_list, add.names = names_vec)

s1 <- pathway_stats(cc_list[[1]]);  s2 <- pathway_stats(cc_list[[2]])
flow1 <- stats::setNames(s1$flow[match(top_pw, rownames(s1))], top_pw)
flow2 <- stats::setNames(s2$flow[match(top_pw, rownames(s2))], top_pw)
flow1[is.na(flow1)] <- 0;  flow2[is.na(flow2)] <- 0
dominant <- ifelse(flow1 >= flow2, COL_USE[1], COL_USE[2])

p_base <- rankNet(cc_merged, mode = "comparison", stacked = TRUE,
                  do.stat = TRUE, color.use = COL_USE, signaling = top_pw)

pb       <- ggplot_build(p_base)
y_labels <- pb$layout$panel_params[[1]]$y$get_labels()
label_cols <- dominant[y_labels]
label_cols[is.na(label_cols)] <- "grey30"

p_final <- p_base +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13, color = label_cols))
```

### Same-axis scatter for N-condition comparison

```r
condition_list <- list(
  "project_specific" = cc_a,  # REPLACE keys with actual condition names
  "project_specific" = cc_b
)

plots <- lapply(names(condition_list), function(nm) {
  netAnalysis_signalingRole_scatter(condition_list[[nm]],
    label.size = 5, font.size = 12, font.size.title = 13) +
    ggplot2::ggtitle(nm) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))
})

get_range <- function(p) {
  d <- ggplot_build(p)$data[[1]]
  list(x = range(d$x, na.rm = TRUE), y = range(d$y, na.rm = TRUE))
}
ranges  <- lapply(plots, get_range)
x_range <- range(sapply(ranges, function(r) r$x))
y_range <- range(sapply(ranges, function(r) r$y))
xpad <- diff(x_range) * 0.12
ypad <- diff(y_range) * 0.12

plots_fixed <- lapply(plots, function(p) {
  p + ggplot2::coord_cartesian(
    xlim = c(x_range[1] - xpad, x_range[2] + xpad),
    ylim = c(y_range[1] - ypad, y_range[2] + ypad))
})

ggsave(file.path(OUTPUT_DIR, "compare_scatter_sameaxis.pdf"),
       patchwork::wrap_plots(plots_fixed, nrow = 2),
       width = 14, height = 11, device = cairo_pdf)
```

### Script 6 common pitfalls

| Error / Symptom | Cause | Fix |
|---|---|---|
| OOM / exit 137 | Loading all Seurat subsets simultaneously | `extract_tissue()` one at a time; `rm(obj); gc()` after extraction |
| "subscript out of bounds" on sparse matrix | Seurat v5 layers not joined | `JoinLayers(sub_obj)` before `GetAssayData()` |
| "unused argument (cell.level=TRUE)" | CellChat v2 removed parameter | Remove `cell.level = TRUE` from `mergeCellChat()` |
| "non-conformable arrays" in heatmap | Cell types differ between objects | `subsetCellChat(cc, idents.use = common_cts)` before merging |
| "rowSums" error in computeCentrality | Called on merged object | Run `netAnalysis_computeCentrality` on individual objects only |
| Empty / 3.6 KB PDF for comparison heatmap | ComplexHeatmap not printed | `draw(ht1 + ht2)` inside `pdf()` / `dev.off()` |

---

## Aesthetic Consistency

- Use the **same color palette as the UMAP** across all CellChat plots
- Apply `color.use = CELL_COLORS` in every `netVisual_circle()` and `netVisual_chord_gene()` call
- Color axis labels in heatmaps and dot grids with the matching palette vectors
- Compute per-axis color vectors **before** the `ggplot()` call (not inside `theme()`)
- Use `useDingbats = FALSE` in every `pdf()` and `ggsave()` call
- ASCII arrows in titles only — use `->` not `→`. Unicode arrows cause encoding warnings.
- Size the bubble height: `max(6, length(cat_paths) * 0.55 + 3)` inches

---

## Output Directory Convention

```
output/cellchat/
  CellChat_all_interactions.csv           -- full LR table from inference (Script 1)
  cellchat_object.Rds                     -- saved inference object (Script 1)
  CellChat_{CATEGORY}_bubble_out.pdf      -- source-to-target bubbles (Script 2)
  CellChat_{CATEGORY}_bubble_in.pdf       -- target-to-source bubbles (Script 2)
  CellChat_{CATEGORY}_circle.pdf          -- circle interaction count (Script 2)
  CellChat_{CATEGORY}_chord.pdf           -- gene-level chord (Script 2)
  CellChat_{CATEGORY}_heatmap.pdf         -- interaction count heatmap (Script 2)
  CellChat_global_circle.pdf              -- all-pathway circle (Script 2)
  CellChat_pathway_ranking.pdf            -- rankNet (Script 2)
  CellChat_signaling_role_scatter.pdf     -- scatter (Script 2)
  CellChat_signaling_outgoing.pdf         -- role heatmap (Script 2)
  CellChat_signaling_incoming.pdf         -- role heatmap (Script 2)
  CellChat_stacked_sources_to_targets.pdf -- stacked bubbles (Script 3)
  CellChat_stacked_targets_to_sources.pdf -- stacked bubbles (Script 3)
  pathway_dist_all_count.pdf              -- bar plots (Script 4)
  pathway_dist_all_strength.pdf           -- bar plots (Script 4)
  circos_sources_to_targets_bypathway.pdf -- circos (Script 5)
  circos_targets_to_sources_bypathway.pdf -- circos (Script 5)
  compare_interactions.pdf                -- Script 6 (conditional)
  compare_heatmap.pdf                     -- Script 6 (conditional)
  compare_pathways.pdf                    -- Script 6 (conditional)
  compare_scatter.pdf                     -- Script 6 (conditional)
  compare_scatter_sameaxis.pdf            -- Script 6, N conditions (conditional)
```

---

## pdf()/dev.off() Exceptions Used

| Script | Usage | CONVENTIONS.md §4 exception |
|---|---|---|
| Scripts 2, 3, 5 | `pdf(...); netVisual_circle(...); dev.off()` and `pdf(...); netVisual_chord_gene(...); dev.off()` | Exception #3 — circlize/base R graphics (CellChat chord and circle internally use base R graphics devices) |
| Script 5 (`make_circos`) | `cairo_pdf(...); chordDiagram(...); dev.off()` | Exception #3 — circlize/base R graphics; cairo_pdf variant used for Illustrator-editable text |
| Script 6 (comparison) | `pdf(...); ComplexHeatmap::draw(ht1 + ht2); dev.off()` | Exception #1 — ComplexHeatmap::draw() sends output to the active graphics device |

---

## Phase 4 Staging

The following project-specific values appeared in the v1 source. They must be moved to
`examples/` files in Phase 4. Do NOT use these values in any generated R script without
first confirming they match the current project state.

### examples/cellchat_example_A.md (stage for Phase 4)

Validated 2026-03-30 (HumanFat Interactions run).

```
LABEL_COL  <- "mylabel"   # EC subset label column; whole-object uses "celltype"
ORGANISM   <- "human"
INPUT_RDS  <- "project_specific"    # TODO: fill in from HumanFat project records
GROUP_COL  <- "tissue_type"         # Script 6: tissue comparison column

SOURCE_TYPES: AEC, CapEC, CapEC2, VenEC1, VenEC2, VenEC3
  (exclude RibHighEC — low-quality/high-ribosomal cluster)
TARGET_TYPES: Mac types + ASPC types (project_specific — confirm from project records)

Script 6 comparison tissue pairs:
  LiposuctionFat vs BreastFat
  (four-tissue variant also validated: project_specific — confirm tissue names)

EC subtype colors (from validated_examples.yaml HumanFat entry):
  AEC:    "#F4433C"
  CapEC:  "#FF9800"
  CapEC2: "#FFEB3B"
  VenEC1: "#4CAF50"
  VenEC2: "#2196F3"
  VenEC3: "#9C27B0"

Tissue comparison color convention (HumanFat):
  Orbital Fat:   "#20B2AA"
  Lipoma:        "#DEB887"
  Visceral Fat:  "#CD853F"
  Subcutaneous Fat: "#F4A460"

Script 3 stacked bubble validated figure:
  3 panels: Angiocrine (top), Chemokine & Immune (middle), Adhesion (bottom)
  ECM category intentionally excluded (dominated by COLLAGEN/LAMININ)
  Separate plots for: EC-to-Mac, EC-to-ASPC-merged, EC-to-ASPC16

Script 5 circos validated:
  4 standard plots: EC-to-ASPC16, ASPC16-to-EC, EC-to-Mac, Mac-to-EC
  Two inference objects: collapsed ASPC labels (v3) + uncollapsed ASPC16 labels
  (EC-to-ASPC16 requires separate inference run with ASPC subtypes uncollapsed)

PATHWAY_CATEGORIES (HumanFat — project_specific, confirm from project records):
  Similar biological groupings to KidneyNew but with adipose-specific pathway focus.
  TODO: extract from HumanFat project CLAUDE.md.

CAT_COLORS (bar plots, validated HumanFat):
  Angiocrine:         "#FF9800"
  Metabolic:          "#4CAF50"
  Chemokine/Immune:   "#2196F3"
  Adhesion:           "#9C27B0"
  ECM:                "#795548"
  Guidance:           "#F06292"
  Other:              "#BDBDBD"
```
