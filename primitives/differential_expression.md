---
# Differential Expression — Core Functions — v2
# Migrated from ~/claude-skills/pipelines/LargeDataset/methods/differential_expression.md
# F2 FIX: is_ambient() unified to full dual-filter (AMBIENT_PATTERNS + AMBIENT_EXPLICIT 16 genes).
#         The v1 anatomical_de_analysis.md version (patterns-only) is wrong and not carried forward.
# F3 FIX: is_confound() authored new — superset of is_ambient() adding sex-linked genes,
#         HLA class II, histones, and unannotated identifiers. Phantom in v1 (called but never defined).
# F4 FIX: functional_gene_sets removed from this primitive. It is a placeholder only.
#         EC-specific gene sets from v1 belong in examples/ (Phase 4).
---

# Differential Expression — Core Functions

Compares expression between two categorical groups within a subset object.
Designed to run across multiple comparison pairs × multiple cell-type scopes.

**Filter hierarchy:**
- `is_ambient()` — narrower filter for standard volcano/DE labeling (ambient contamination only)
- `is_confound()` — broader filter for comprehensive DE sweeps (ambient + sex-linked + HLA + histones + unannotated)
  Use `is_confound()` in de_comprehensive_csv.md; use `is_ambient()` for standard visualization functions.

---

## Critical Constraints

| Do NOT | Do Instead | Why |
|---|---|---|
| Use `!!sym(col)` in `group_by` inside pivot functions | Rename the column to a fixed name before pivoting | `!!sym()` fails in certain dplyr/tidyr contexts — cryptic type error |
| Use `annotate("text")` with `facet_wrap` | Use `geom_text(data = standalone_df, inherit.aes = FALSE)` pinned to one facet | `annotate()` expands vectors by n_panels, causing `check_aesthetics` error |
| Place text outside a facet panel without pinning the facet variable | Pin the label data frame to one facet level | Without pinning, the label renders in EVERY panel |
| Draw direction squares as a geom inside the main ggplot | Use `patchwork` to attach a separate narrow direction-strip panel | Numeric x in a `scale_x_discrete` plot renders inside panel regardless of x value |
| Use em-dash (—, U+2014) in PDF plot titles | Use plain ASCII equivalent | grid.Call.graphics rendering error |
| Use base R system Rscript | Always run via `conda run --no-capture-output -n r-env Rscript` | sp package is broken in base R |

---

## Configuration Block

```r
# ── Thresholds ────────────────────────────────────────────────────────────────
PADJ_CUT        <- 0.05
LFC_CUT         <- 0.5        # |log2FC| minimum
N_CELLS_HM      <- 300        # max cells per group sampled into heatmaps
N_TOP_GENES_HM  <- 30         # top N genes per direction for the overall heatmap
N_LABEL_VOLCANO <- 25         # genes labelled per direction on volcano
set.seed(42)

# ── Comparisons ───────────────────────────────────────────────────────────────
# ident1 = shown on the RIGHT of heatmaps and volcano (positive log2FC)
# ident2 = shown on the LEFT  of heatmaps and volcano (negative log2FC)
# col    = the metadata column that holds the group labels (project_specific)
comparisons <- list(
  comp_a = list(
    label  = "GroupA_vs_GroupB",      # used in filenames
    ident1 = "project_specific",      # REPLACE: right side, positive log2FC
    ident2 = "project_specific",      # REPLACE: left side, negative log2FC
    col    = "project_specific"       # REPLACE: metadata column name
  )
)

# ── Cell-type scopes ──────────────────────────────────────────────────────────
# NULL = use the whole subset object as-is; character vector = subset to these labels
cell_subsets <- list(
  AllCells = NULL,
  SubtypeA = "project_specific"     # REPLACE: character vector of subtype labels
)

# ── Label column ──────────────────────────────────────────────────────────────
LABEL_COL <- "project_specific"   # REPLACE: metadata column holding cell type labels

# ── Color references ─────────────────────────────────────────────────────────
# group_colors and subtype_colors are defined in the project CLAUDE.md context block.
# Do not hardcode here. See context/validated_examples.yaml for per-project palettes.
# group_colors   <- project_specific  # named vector: group value → hex color
# subtype_colors <- project_specific  # named vector: subtype label → hex color

# ── Functional gene sets ─────────────────────────────────────────────────────
# functional_gene_sets is a brief-supplied named list of gene vectors.
# It is NOT defined here — the caller must provide it.
# EC-specific gene sets from v1 belong in examples/ (Phase 4), not this primitive.
# Example structure:
#   functional_gene_sets <- list(
#     "Section A" = c("GENE1", "GENE2", ...),   # project_specific
#     "Section B" = c("GENE3", "GENE4", ...)    # project_specific
#   )
# If functional_gene_sets is not provided, make_functional_heatmap() and
# make_functional_dotplot() will error — this is intentional.
```

---

## Ambient RNA Exclusion

**Critical:** Ambient RNA causes immunoglobulins, haemoglobins, ribosomal, and
other contaminant genes to appear as DE. Exclude from all visualization layers
while keeping raw DE CSV tables unfiltered.

### `is_ambient()` — Standard Filter

Use for: volcano labels, heatmap gene selection, standard dot plots.

```r
AMBIENT_PATTERNS <- c(
  "^IGK", "^IGL", "^IGH",           # immunoglobulin chains
  "^HBA", "^HBB", "^HBD", "^HBE",  # haemoglobin alpha/beta/delta/epsilon
  "^HBG", "^HBM", "^HBQ", "^HBZ",  # haemoglobin gamma/mu/theta/zeta
  "^MT-",                            # mitochondrial (technical noise)
  "^RPS", "^RPL"                     # ribosomal proteins small/large subunit
)

AMBIENT_EXPLICIT <- c(
  "ALB", "FGA", "FGB", "FGG", "FGL1",          # plasma/liver proteins
  "PF4", "PPBP", "GP9", "GP1BA", "GP1BB", "ITGA2B",  # platelet markers
  "MALAT1", "NEAT1",                             # abundant nuclear lncRNAs
  "TPSAB1", "TPSB2", "CPA3"                     # mast cell tryptases
)

is_ambient <- function(genes) {
  pattern_hit <- Reduce(`|`, lapply(AMBIENT_PATTERNS, function(p) grepl(p, genes)))
  genes %in% AMBIENT_EXPLICIT | pattern_hit
}
```

**Note:** Extend `AMBIENT_EXPLICIT` with project-specific contaminants found during review.

---

### `is_confound()` — Comprehensive Filter

Use for: comprehensive DE CSV generation (de_comprehensive_csv module).
This is a superset of `is_ambient()` — every gene flagged by `is_ambient()` is
also flagged by `is_confound()`, plus additional categories.

```r
CONFOUND_PATTERNS <- c(
  # All ambient patterns
  "^IGK", "^IGL", "^IGH",
  "^HBA", "^HBB", "^HBD", "^HBE",
  "^HBG", "^HBM", "^HBQ", "^HBZ",
  "^MT-",
  "^RPS", "^RPL",
  # HLA class II — context-dependent antigen presentation, often confounds DE
  "^HLA-D[RPQ]",
  # Histones — cell cycle / chromatin state, rarely of interest in expression DE
  "^HIST[0-9]",
  # Unannotated identifiers — ENSEMBL IDs, chromosomal loci, lincRNAs without symbols
  "^ENSG[0-9]+", "^ENSMUSG[0-9]+",
  "^AC[0-9]+\\.[0-9]+",
  "^AL[0-9]+\\.[0-9]+",
  "^AP[0-9]+\\.[0-9]+",
  "^LINC[0-9]+"
)

CONFOUND_EXPLICIT <- c(
  # All ambient explicit genes
  "ALB", "FGA", "FGB", "FGG", "FGL1",
  "PF4", "PPBP", "GP9", "GP1BA", "GP1BB", "ITGA2B",
  "MALAT1", "NEAT1",
  "TPSAB1", "TPSB2", "CPA3",
  # Sex-linked genes — chromosome X/Y markers that confound disease/condition DE
  "XIST",                                   # X-inactive specific transcript (female)
  "RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D",
  "NLGN4Y", "RPS4Y2", "USP9Y", "UTY", "ZFY"  # Y-chromosome markers (male)
)

is_confound <- function(genes) {
  pattern_hit <- Reduce(`|`, lapply(CONFOUND_PATTERNS, function(p) grepl(p, genes)))
  genes %in% CONFOUND_EXPLICIT | pattern_hit
}
```

---

## Step 1: FindMarkers

```r
run_findmarkers <- function(so_obj, comp) {
  # so_obj: Seurat object with comp$col in metadata
  # comp:   list with ident1, ident2, col, label
  cells_keep  <- so_obj@meta.data[[comp$col]] %in% c(comp$ident1, comp$ident2)
  so_comp     <- so_obj[, cells_keep]
  Idents(so_comp) <- comp$col

  n1 <- sum(Idents(so_comp) == comp$ident1)
  n2 <- sum(Idents(so_comp) == comp$ident2)

  if (n1 < 10 || n2 < 10) {
    message(sprintf("SKIPPED — too few cells (%s: %d | %s: %d)",
                    comp$ident1, n1, comp$ident2, n2))
    return(NULL)
  }

  FindMarkers(
    so_comp,
    ident.1 = comp$ident1,
    ident.2 = comp$ident2,
    assay    = "RNA",
    test.use = "wilcox",
    logfc.threshold = 0,    # return all genes — filter later for flexibility
    min.pct = 0.05
  ) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(p_val_adj)
}
# Save full unfiltered results: DE_full_{subset}_{comp_label}.csv
```

---

## Step 2: Volcano Plot

```r
make_volcano <- function(markers, comp, subset_name, n_label = N_LABEL_VOLCANO,
                          functional_gene_sets = NULL) {
  # functional_gene_sets: caller-provided named list; used for gene prioritization in labels.
  #   Pass NULL to skip biological prioritization (all significant genes weighted equally).

  df <- markers %>%
    mutate(
      sig = p_val_adj < PADJ_CUT & abs(avg_log2FC) > LFC_CUT,
      significance = case_when(
        sig & avg_log2FC >  0 ~ "up_ident1",
        sig & avg_log2FC <= 0 ~ "up_ident2",
        TRUE                  ~ "ns"
      )
    )

  # Prioritise biologically known genes in top-N label selection
  all_bio_genes <- if (!is.null(functional_gene_sets)) unique(unlist(functional_gene_sets)) else character(0)
  df <- df %>% mutate(is_bio = gene %in% all_bio_genes)

  get_labels <- function(direction) {
    sub_df  <- df %>% filter(significance == direction, !is_ambient(gene))
    bio_top <- sub_df %>% filter(is_bio)  %>% arrange(p_val_adj) %>%
                head(ceiling(n_label * 0.6)) %>% pull(gene)
    rest    <- sub_df %>% filter(!is_bio) %>% arrange(p_val_adj) %>%
                head(n_label - length(bio_top)) %>% pull(gene)
    c(bio_top, rest)
  }

  label_genes <- c(get_labels("up_ident1"), get_labels("up_ident2"))
  df <- df %>% mutate(label = ifelse(gene %in% label_genes, gene, ""))

  colors <- c(up_ident1 = "#B2182B", up_ident2 = "#2166AC", ns = "grey78")

  ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj),
                 color = significance, size = significance)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", color = "grey50") +
    ggrepel::geom_text_repel(
      aes(label = label), size = 3.5, fontface = "italic",
      max.overlaps = 30, box.padding = 0.4, segment.size = 0.3,
      segment.color = "grey60"
    ) +
    scale_color_manual(values = colors, guide = "none") +
    scale_size_manual(values = c(up_ident1 = 1.3, up_ident2 = 1.3, ns = 0.7), guide = "none") +
    labs(
      title    = sprintf("%s | %s", subset_name, gsub("_", " ", comp$label)),
      subtitle = sprintf("Up in %s: %d  |  Up in %s: %d\n(padj < %g, |log2FC| > %g)",
                         comp$ident1, sum(df$significance == "up_ident1"),
                         comp$ident2, sum(df$significance == "up_ident2"),
                         PADJ_CUT, LFC_CUT),
      x = "log2 Fold Change", y = expression(-log[10](p[adj]))
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title    = element_text(hjust = 0.5, size = 17, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey35"),
      axis.title    = element_text(size = 15),
      axis.text     = element_text(size = 13),
      panel.grid.major = element_blank()
    )
}
```

---

## Step 3: Overall Cell-Level Heatmap (ComplexHeatmap)

Samples up to `N_CELLS_HM` cells per group, orders top N DE genes per direction,
draws a cell-level z-score heatmap with cells clustered within each group.

```r
make_overall_heatmap <- function(so_obj, markers, comp, subset_name,
                                  label_col,         # metadata column for subtype annotation
                                  subtype_colors,    # named vector: subtype → hex color
                                  group_colors,      # named vector: group → hex color
                                  n_cells = N_CELLS_HM, n_genes = N_TOP_GENES_HM) {
  sig_df <- markers %>%
    filter(p_val_adj < PADJ_CUT, abs(avg_log2FC) > LFC_CUT, !is_ambient(gene))

  top_genes <- c(
    sig_df %>% filter(avg_log2FC >  0) %>% arrange(p_val_adj) %>% head(n_genes) %>% pull(gene),
    sig_df %>% filter(avg_log2FC <= 0) %>% arrange(p_val_adj) %>% head(n_genes) %>% pull(gene)
  )
  if (length(top_genes) == 0) return(NULL)

  meta   <- so_obj@meta.data
  c1_all <- rownames(meta)[meta[[comp$col]] == comp$ident1]
  c2_all <- rownames(meta)[meta[[comp$col]] == comp$ident2]
  c1_use <- if (length(c1_all) > n_cells) sample(c1_all, n_cells) else c1_all
  c2_use <- if (length(c2_all) > n_cells) sample(c2_all, n_cells) else c2_all

  cells_use  <- c(c1_use, c2_use)
  so_sub     <- so_obj[, cells_use]
  expr_mat   <- as.matrix(GetAssayData(so_sub, assay = "RNA", layer = "data")[top_genes, ])
  expr_z     <- t(scale(t(expr_mat)))
  expr_z     <- pmin(pmax(expr_z, -3), 3)

  col_split <- factor(so_sub@meta.data[[comp$col]],
                      levels = c(comp$ident1, comp$ident2))

  gene_dir <- sapply(top_genes, function(g) {
    lfc <- markers$avg_log2FC[markers$gene == g]
    if (length(lfc) == 0 || is.na(lfc)) return("Unknown")
    if (lfc > 0) sprintf("Up: %s", comp$ident1) else sprintf("Up: %s", comp$ident2)
  })
  dir_colors <- setNames(
    c(group_colors[comp$ident1], group_colors[comp$ident2]),
    c(sprintf("Up: %s", comp$ident1), sprintf("Up: %s", comp$ident2))
  )

  col_ann <- HeatmapAnnotation(
    Group   = so_sub@meta.data[[comp$col]],
    Subtype = so_sub@meta.data[[label_col]],
    col = list(
      Group   = group_colors[names(group_colors) %in% c(comp$ident1, comp$ident2)],
      Subtype = subtype_colors
    )
  )

  row_ann <- rowAnnotation(
    Direction = gene_dir,
    col       = list(Direction = dir_colors),
    width     = unit(4, "mm"),
    annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
    show_legend = FALSE
  )

  ComplexHeatmap::Heatmap(
    expr_z,
    name                  = "Z-score",
    col                   = colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "#F7F7F7", "#B2182B")),
    top_annotation        = col_ann,
    right_annotation      = row_ann,
    cluster_rows          = FALSE,
    cluster_columns       = TRUE,
    cluster_column_slices = FALSE,
    column_split          = col_split,
    show_column_names     = FALSE,
    row_names_gp          = gpar(fontsize = 10, fontface = "italic"),
    column_title_gp       = gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_param  = list(
      title_gp  = gpar(fontsize = 11, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    )
  )
}
# After building: call draw(ht) inside a pdf() block for ComplexHeatmap output
```

---

## Step 4: Functional Heatmap

Same structure as the overall heatmap but genes come from `functional_gene_sets`
rather than top-N DE. functional_gene_sets MUST be caller-provided.

```r
make_functional_heatmap <- function(so_obj, markers, comp, subset_name,
                                     label_col, subtype_colors, group_colors,
                                     functional_gene_sets,   # REQUIRED: caller-provided named list
                                     n_cells = N_CELLS_HM) {
  # functional_gene_sets: named list of gene vectors (project-specific biology)
  # Genes present in functional_gene_sets that are also significant DE appear in heatmap.
  all_fg_genes <- unique(unlist(functional_gene_sets))

  sig_df <- markers %>%
    filter(p_val_adj < PADJ_CUT, abs(avg_log2FC) > LFC_CUT,
           !is_ambient(gene), gene %in% all_fg_genes)

  if (nrow(sig_df) == 0) { message("No significant functional genes found."); return(NULL) }

  # Gene ordering: ident1-up then ident2-up, within each section sorted by p_adj
  genes_up1 <- sig_df %>% filter(avg_log2FC > 0) %>% arrange(p_val_adj) %>% pull(gene)
  genes_up2 <- sig_df %>% filter(avg_log2FC <= 0) %>% arrange(p_val_adj) %>% pull(gene)
  genes_ordered <- c(genes_up1, genes_up2)

  # Build section factor for row_split
  section_df <- lapply(names(functional_gene_sets), function(s) {
    data.frame(gene = functional_gene_sets[[s]], section = s)
  }) %>% do.call(rbind, .)

  row_split <- factor(
    section_df$section[match(genes_ordered, section_df$gene)],
    levels = names(functional_gene_sets)
  )

  meta   <- so_obj@meta.data
  c1_use <- sample(rownames(meta)[meta[[comp$col]] == comp$ident1],
                   min(n_cells, sum(meta[[comp$col]] == comp$ident1)))
  c2_use <- sample(rownames(meta)[meta[[comp$col]] == comp$ident2],
                   min(n_cells, sum(meta[[comp$col]] == comp$ident2)))
  so_sub <- so_obj[, c(c1_use, c2_use)]

  expr_mat <- as.matrix(GetAssayData(so_sub, assay = "RNA", layer = "data")[genes_ordered, ])
  expr_z   <- pmin(pmax(t(scale(t(expr_mat))), -3), 3)

  col_split <- factor(so_sub@meta.data[[comp$col]], levels = c(comp$ident1, comp$ident2))

  col_ann <- HeatmapAnnotation(
    Group   = so_sub@meta.data[[comp$col]],
    Subtype = so_sub@meta.data[[label_col]],
    col     = list(Group = group_colors, Subtype = subtype_colors)
  )

  ComplexHeatmap::Heatmap(
    expr_z,
    name             = "Z-score",
    col              = colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "#F7F7F7", "#B2182B")),
    top_annotation   = col_ann,
    cluster_rows     = FALSE,
    cluster_columns  = TRUE,
    cluster_column_slices = FALSE,
    column_split     = col_split,
    row_split        = row_split,
    row_title_gp     = gpar(fontsize = 11, fontface = "bold.italic"),
    row_gap          = unit(3, "mm"),
    show_column_names = FALSE,
    row_names_gp     = gpar(fontsize = 10, fontface = "italic")
  )
}
```

---

## Step 5: Top-Gene Dot Plot

```r
make_topgene_dotplot <- function(so_obj, markers, comp, subset_name,
                                  label_col,       # metadata column for subtype labels
                                  group_col,       # metadata column for group (facet variable)
                                  subtype_colors,  # named vector for x-axis label colors
                                  n_each = 12,     # top N genes per direction
                                  output_file) {
  sig_df <- markers %>%
    filter(p_val_adj < PADJ_CUT, abs(avg_log2FC) > LFC_CUT, !is_ambient(gene))

  top_up1 <- sig_df %>% filter(avg_log2FC >  0) %>% arrange(p_val_adj) %>% head(n_each) %>% pull(gene)
  top_up2 <- sig_df %>% filter(avg_log2FC <= 0) %>% arrange(p_val_adj) %>% head(n_each) %>% pull(gene)
  gene_order <- c(top_up1, top_up2)

  # Use fixed column names to avoid !!sym() pivot failures (Critical Constraints)
  expr <- as.data.frame(t(as.matrix(GetAssayData(so_obj, assay = "RNA", layer = "data")[gene_order, ])))
  expr$group_var  <- so_obj@meta.data[[group_col]]   # fixed column name
  expr$subtype_var <- so_obj@meta.data[[label_col]]  # fixed column name

  dot_df <- tidyr::pivot_longer(expr, cols = all_of(gene_order),
                                 names_to = "gene", values_to = "expression") %>%
    group_by(group_var, subtype_var, gene) %>%
    summarise(
      avg_exp = mean(expm1(expression)),
      pct_exp = mean(expression > 0) * 100,
      .groups = "drop"
    ) %>%
    mutate(
      avg_exp_scaled = scale(avg_exp)[, 1],
      gene = factor(gene, levels = rev(gene_order))
    )

  n_subtypes <- length(unique(dot_df$subtype_var))
  n_genes    <- length(gene_order)

  p <- ggplot(dot_df, aes(x = subtype_var, y = gene, size = pct_exp, fill = avg_exp_scaled)) +
    geom_point(shape = 21, color = "grey30", stroke = 0.32) +
    facet_wrap(~ group_var, ncol = length(unique(dot_df$group_var)), scales = "free_x") +
    scale_fill_gradientn(colors = c("#F5F5F5", "#FFF9C4", "#FFB300", "#E53935")) +
    scale_size_continuous(range = c(0.3, 6), limits = c(0, 100)) +
    scale_x_discrete(position = "top") +
    geom_hline(yintercept = n_each + 0.5, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x.top = element_text(size = 13.5, angle = 45, hjust = 0, vjust = 0,
                                      face = "bold", color = subtype_colors[levels(factor(dot_df$subtype_var))]),
      axis.text.y     = element_text(size = 11, color = "grey10", face = "italic"),
      strip.text      = element_text(size = 13, face = "bold"),
      panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      legend.position = "right"
    )

  fig_h <- n_genes  * 0.28 + 2
  fig_w <- n_subtypes * 0.75 * length(unique(dot_df$group_var)) + 3
  ggsave(output_file, plot = p, width = fig_w, height = fig_h,
         units = "in", device = "pdf", useDingbats = FALSE)
  p
}
```

---

## Step 6: Functional Dot Plot with Direction Strip

```r
make_functional_dotplot <- function(so_obj, markers, comp, subset_name,
                                     label_col, group_col,
                                     functional_gene_sets,  # REQUIRED: caller-provided
                                     subtype_colors,
                                     show_direction = FALSE,
                                     output_file) {
  # Build dot_df from functional_gene_sets genes that are significant
  all_fg_genes <- unique(unlist(functional_gene_sets))
  sig_df <- markers %>%
    filter(p_val_adj < PADJ_CUT, abs(avg_log2FC) > LFC_CUT,
           !is_ambient(gene), gene %in% all_fg_genes) %>%
    arrange(p_val_adj)

  gene_order <- sig_df$gene
  if (length(gene_order) == 0) { message("No significant functional genes."); return(NULL) }

  expr <- as.data.frame(t(as.matrix(GetAssayData(so_obj, assay = "RNA", layer = "data")[gene_order, ])))
  expr$group_var   <- so_obj@meta.data[[group_col]]
  expr$subtype_var <- so_obj@meta.data[[label_col]]

  dot_df <- tidyr::pivot_longer(expr, cols = all_of(gene_order),
                                 names_to = "gene", values_to = "expression") %>%
    group_by(group_var, subtype_var, gene) %>%
    summarise(avg_exp = mean(expm1(expression)),
              pct_exp = mean(expression > 0) * 100, .groups = "drop") %>%
    mutate(avg_exp_scaled = scale(avg_exp)[, 1],
           gene = factor(gene, levels = rev(gene_order)))

  # Section dividers and pinned labels
  section_rows <- lapply(names(functional_gene_sets), function(s) {
    g <- intersect(functional_gene_sets[[s]], gene_order)
    if (length(g) == 0) return(NULL)
    data.frame(section = s, gene = g)
  }) %>% Filter(Negate(is.null), .) %>% do.call(rbind, .)

  section_label_df <- section_rows %>%
    group_by(section) %>%
    summarise(y = mean(as.numeric(factor(gene, levels = rev(gene_order)))), .groups = "drop") %>%
    mutate(
      label = section,
      x     = length(unique(dot_df$subtype_var)) + 0.65,
      # Pin to last group level to avoid rendering in all facet panels
      group_var = factor(tail(levels(factor(dot_df$group_var)), 1),
                         levels = levels(factor(dot_df$group_var)))
    )

  n_subtypes <- length(unique(dot_df$subtype_var))
  n_genes    <- length(gene_order)
  n_groups   <- length(unique(dot_df$group_var))

  p <- ggplot(dot_df, aes(x = subtype_var, y = gene, size = pct_exp, fill = avg_exp_scaled)) +
    geom_point(shape = 21, color = "grey30", stroke = 0.32) +
    facet_wrap(~ group_var, ncol = n_groups, scales = "free_x") +
    scale_fill_gradientn(colors = c("#F5F5F5", "#FFF9C4", "#FFB300", "#E53935")) +
    scale_size_continuous(range = c(0.3, 6), limits = c(0, 100)) +
    scale_x_discrete(position = "top") +
    geom_text(data = section_label_df,
              aes(x = x, y = y, label = label, group_var = group_var),
              inherit.aes = FALSE, size = 4.5, hjust = 0,
              fontface = "bold.italic", color = "grey30") +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x.top = element_text(size = 13.5, angle = 45, hjust = 0, vjust = 0, face = "bold",
                                      color = subtype_colors[levels(factor(dot_df$subtype_var))]),
      axis.text.y     = element_text(size = 11, color = "grey10", face = "italic"),
      strip.text      = element_text(size = 13, face = "bold"),
      panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      plot.margin     = margin(t = 8, r = 140, b = 8, l = 8),
      legend.position = "none"
    ) +
    coord_cartesian(clip = "off")

  if (show_direction) {
    dir_lookup <- setNames(
      ifelse(markers$avg_log2FC[match(gene_order, markers$gene)] > 0, "up_ident1", "up_ident2"),
      gene_order
    )
    dir_colors_dot <- c(up_ident1 = "#B2182B", up_ident2 = "#2166AC")

    dir_strip_df <- data.frame(
      gene      = factor(rev(gene_order), levels = rev(gene_order)),
      x         = 1L,
      direction = dir_lookup[rev(gene_order)]
    )

    p_strip <- ggplot(dir_strip_df, aes(x = x, y = gene, fill = direction)) +
      geom_tile(width = 0.7, height = 0.85, color = "white", linewidth = 0.2) +
      scale_fill_manual(values = dir_colors_dot, name = "DE Direction") +
      scale_y_discrete(limits = rev(gene_order)) +
      scale_x_continuous(limits = c(0.35, 1.65), expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(legend.position = "right",
            plot.margin = margin(t = 26, r = 5, b = 8, l = 5))  # 26pt top pad aligns with facet strip

    p <- p + theme(plot.margin = margin(t = 8, r = 5, b = 8, l = 8))
    p <- patchwork::wrap_plots(p, p_strip, widths = c(n_subtypes * n_groups, 1))
  }

  fig_h <- n_genes  * 0.28 + 2
  fig_w <- n_subtypes * 0.75 * n_groups + if (show_direction) 1.5 else 0 + 3
  ggsave(output_file, plot = p, width = fig_w, height = fig_h,
         units = "in", device = "pdf", useDingbats = FALSE)
  p
}
```

---

## Step 7: Pathway Bar Plot (GO:BP ORA)

```r
make_pathway_barplot <- function(markers, comp, subset_name, universe_genes,
                                  n_pathways = 15, output_file) {
  # universe_genes: character vector of all genes tested (e.g. rownames of so after filtering)
  sig_df <- markers %>%
    filter(p_val_adj < PADJ_CUT, abs(avg_log2FC) > LFC_CUT, !is_confound(gene))

  run_ora <- function(genes, direction_label) {
    if (length(genes) < 5) return(NULL)
    res <- clusterProfiler::enrichGO(
      gene          = genes,
      universe      = universe_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    res <- clusterProfiler::simplify(res, cutoff = 0.7)
    data.frame(res@result) %>%
      head(n_pathways) %>%
      mutate(direction = direction_label,
             neg_log10_p = -log10(p.adjust))
  }

  genes_up1 <- sig_df %>% filter(avg_log2FC >  0, !is_ambient(gene)) %>% pull(gene)
  genes_up2 <- sig_df %>% filter(avg_log2FC <= 0, !is_ambient(gene)) %>% pull(gene)

  bar_df <- bind_rows(
    run_ora(genes_up1, sprintf("Up: %s", comp$ident1)),
    run_ora(genes_up2, sprintf("Up: %s", comp$ident2))
  )

  if (is.null(bar_df) || nrow(bar_df) == 0) {
    message("No enriched pathways found."); return(NULL)
  }

  p <- ggplot(bar_df, aes(x = reorder(Description, neg_log10_p),
                           y = ifelse(direction == sprintf("Up: %s", comp$ident1),
                                      neg_log10_p, -neg_log10_p),
                           fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("#B2182B", "#2166AC")) +
    labs(title = sprintf("GO:BP -- %s | %s", subset_name, comp$label),
         x = NULL, y = expression(-log[10](p[adj]))) +
    theme_classic(base_size = 14)

  ggsave(output_file, plot = p, width = 10, height = max(6, nrow(bar_df) * 0.3 + 2),
         units = "in", device = "pdf", useDingbats = FALSE)
  p
}
```

---

## Main Analysis Loop

```r
# Requires: comparisons list, cell_subsets list, LABEL_COL, group_colors, subtype_colors,
#           functional_gene_sets (project_specific — must be defined in CLAUDE.md context)

for (comp in comparisons) {
  for (scope_name in names(cell_subsets)) {
    scope_labels <- cell_subsets[[scope_name]]
    so_scope <- if (is.null(scope_labels)) so_obj else
                  subset(so_obj, cells = rownames(so_obj@meta.data)[so_obj@meta.data[[LABEL_COL]] %in% scope_labels])

    markers <- run_findmarkers(so_scope, comp)
    if (is.null(markers)) next

    out_prefix <- file.path(output_dir, sprintf("%s_%s", scope_name, comp$label))
    dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

    # Save unfiltered CSV
    write.csv(markers, paste0(out_prefix, "_DE_full.csv"), row.names = FALSE)

    # Volcano
    p_vol <- make_volcano(markers, comp, scope_name, functional_gene_sets = functional_gene_sets)
    ggsave(paste0(out_prefix, "_volcano.pdf"), plot = p_vol, width = 9, height = 8,
           units = "in", device = "pdf", useDingbats = FALSE)

    # Heatmaps
    ht_overall <- make_overall_heatmap(so_scope, markers, comp, scope_name,
                                        label_col = LABEL_COL,
                                        subtype_colors = subtype_colors,
                                        group_colors = group_colors)
    if (!is.null(ht_overall)) {
      pdf(paste0(out_prefix, "_heatmap_overall.pdf"), width = 12, height = 10)
      ComplexHeatmap::draw(ht_overall)
      dev.off()
    }

    # Functional plots — only if functional_gene_sets is defined
    if (exists("functional_gene_sets") && !is.null(functional_gene_sets)) {
      ht_func <- make_functional_heatmap(so_scope, markers, comp, scope_name,
                                          label_col = LABEL_COL,
                                          subtype_colors = subtype_colors,
                                          group_colors = group_colors,
                                          functional_gene_sets = functional_gene_sets)
      if (!is.null(ht_func)) {
        pdf(paste0(out_prefix, "_heatmap_functional.pdf"), width = 12, height = 10)
        ComplexHeatmap::draw(ht_func)
        dev.off()
      }

      make_topgene_dotplot(so_scope, markers, comp, scope_name,
                            label_col = LABEL_COL, group_col = comp$col,
                            subtype_colors = subtype_colors, n_each = 12,
                            output_file = paste0(out_prefix, "_topgene_dotplot.pdf"))

      make_functional_dotplot(so_scope, markers, comp, scope_name,
                               label_col = LABEL_COL, group_col = comp$col,
                               functional_gene_sets = functional_gene_sets,
                               subtype_colors = subtype_colors, show_direction = FALSE,
                               output_file = paste0(out_prefix, "_functional_dotplot.pdf"))
    }

    # Pathway barplot
    universe <- rownames(so_scope)
    make_pathway_barplot(markers, comp, scope_name, universe_genes = universe,
                          output_file = paste0(out_prefix, "_pathways.pdf"))
  }
}
```
