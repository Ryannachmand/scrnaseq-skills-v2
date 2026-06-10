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

## Axis A: Cluster marker DE

Cluster marker DE (FindAllMarkers, condition-agnostic, one-vs-rest per cluster) is NOT
the responsibility of this primitive. It is handled by:

- `@modules/celltype_subclustering.md` Step 3 (FindAllMarkers within subclustering module)
- `@pipelines/large_dataset/pipeline.md` Stage 4 (whole-object annotation)
- `@pipelines/large_dataset/pipeline.md` Stage 6 (subset annotation)

This primitive (`differential_expression.md`) handles Axis B and Axis C only.

---

## Three-axis DE conventions

When both a condition column (`group_col`) and a cluster/subtype column (`label_col`)
exist on the object, the full DE pattern covers three axes:

**Axis A (cluster markers):** FindAllMarkers, condition-agnostic, one-vs-rest per cluster.
Handled by `@modules/celltype_subclustering.md`. NOT produced by this primitive.

**Axis B (condition global):** Two-group comparison across all cells (or a specified scope).
Standard use: `comparisons:` block with `per_cluster: false` (default). Output per scope.

**Axis C (condition within cluster):** Same two-group comparison run independently within
each unique value of `label_col`. Enabled by setting `per_cluster: true` on a comparison
entry and providing `label_col`. Gated by `MIN_CELLS_PER_CLUSTER` (default 100) per side
per cluster. Output filename: `DE_full_{comp_label}_within_{label_value}.csv`.

All three axes should run by default unless the brief explicitly suppresses one.

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
MIN_CELLS_PER_CLUSTER <- 100  # Axis C minimum: cells per side per cluster (condition_per_cluster)
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
  # Callers MUST NOT bypass the DE filter below via inline code.
  # See "Gene-set dotplots: filtered only" section for the policy.
  # If this function returns NULL, log the decision to output/decision_log.txt.

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
    mutate(avg_exp_scaled = scale(avg_exp)[, 1])

  # Build wide expression matrix: rows = genes, cols = group:subtype combinations
  clust_wide <- dot_df %>%
    filter(gene %in% gene_order) %>%
    mutate(col_id = paste(group_var, subtype_var, sep = ":")) %>%
    select(gene, col_id, avg_exp_scaled) %>%
    pivot_wider(names_from = col_id, values_from = avg_exp_scaled, values_fill = 0)
  clust_mat_m <- as.matrix(clust_wide[, -1])
  rownames(clust_mat_m) <- as.character(clust_wide$gene)

  # Rebuild section_df with clustered gene order within each section
  section_df <- bind_rows(lapply(names(functional_gene_sets), function(s) {
    g <- as.character(functional_gene_sets[[s]])
    g <- g[g %in% rownames(clust_mat_m)]
    if (length(g) > 2) {
      hc <- hclust(dist(clust_mat_m[g, , drop = FALSE], method = "euclidean"),
                   method = "ward.D2")
      g  <- g[hc$order]
    }
    data.frame(gene = g, section = s, stringsAsFactors = FALSE)
  })) %>% distinct(gene, .keep_all = TRUE)

  gene_order <- section_df$gene
  dot_df <- dot_df %>% mutate(gene = factor(gene, levels = rev(gene_order)))

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
                                 n_pathways = 15, output_file = NULL) {

  sig_df    <- markers %>%
    filter(p_val_adj < PADJ_CUT, abs(avg_log2FC) > LFC_CUT, !is_confound(gene))
  genes_up1 <- sig_df %>% filter(avg_log2FC >  0) %>% pull(gene)
  genes_up2 <- sig_df %>% filter(avg_log2FC <= 0) %>% pull(gene)

  run_ora <- function(genes, dir_label) {
    if (length(genes) < 5) return(NULL)
    ego <- enrichGO(
      gene          = genes,
      universe      = universe_genes,   # all genes expressed in the subset object
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      minGSSize     = 10,
      maxGSSize     = 500,
      readable      = TRUE
    )
    if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
    ego <- simplify(ego, cutoff = 0.6, by = "p.adjust", select_fun = min)
    as.data.frame(ego) %>%
      filter(p.adjust < 0.05) %>% arrange(p.adjust) %>%
      slice_head(n = n_pathways) %>%
      mutate(direction = dir_label)
  }

  res_all <- bind_rows(run_ora(genes_up1, comp$ident1),
                       run_ora(genes_up2, comp$ident2))
  if (nrow(res_all) == 0) return(NULL)

  # Signed x: positive = ident2 (RIGHT), negative = ident1 (LEFT)
  res_all <- res_all %>%
    mutate(
      neg_log_p  = -log10(p.adjust),
      x_val      = ifelse(direction == comp$ident2, neg_log_p, -neg_log_p),
      dot_x      = x_val * 1.03,
      label_desc = ifelse(nchar(Description) > 55,
                          paste0(substr(Description, 1, 52), "…"), Description),
      top_genes  = sapply(strsplit(geneID, "/"), function(g) paste(head(g,5), collapse=", ")),
      # Gene annotations: RIGHT bars go right of tip; LEFT bars go right of centre (x=0)
      # This prevents overlap with pathway name labels on the left margin.
      text_x     = ifelse(x_val >= 0,
                          dot_x + max(abs(x_val)) * 0.03,
                          max(abs(x_val)) * 0.03),
      text_hjust = 0,
      gene_text  = sprintf("n=%d  %s", Count, top_genes)
    ) %>%
    arrange(x_val) %>%
    mutate(label_desc = factor(label_desc, levels = unique(label_desc)))

  max_x <- max(abs(res_all$x_val))

  p <- ggplot(res_all, aes(x = x_val, y = label_desc)) +
    geom_col(aes(fill = direction), width = 0.55, alpha = 0.82) +
    geom_vline(xintercept = 0, color = "grey25", linewidth = 0.55) +
    geom_vline(xintercept = c(-log10(0.05), log10(0.05)),     # threshold lines
               color = "grey55", linewidth = 0.35, linetype = "dashed") +
    geom_point(
      aes(x = dot_x, fill = direction, size = Count),
      shape = 21, color = "grey20", stroke = 0.3, alpha = 0.95
    ) +
    geom_text(
      aes(x = text_x, label = gene_text, hjust = text_hjust),
      size = 3.5, fontface = "italic", color = "grey25", lineheight = 0.9
    ) +
    annotate("text", x = -max_x * 0.97, y = Inf,
             label = sprintf("<- Up in %s", comp$ident1),
             hjust = 0, vjust = 2, size = 3.8,
             color = dir_colors[comp$ident1], fontface = "bold") +
    annotate("text", x =  max_x * 0.97, y = Inf,
             label = sprintf("Up in %s ->", comp$ident2),
             hjust = 1, vjust = 2, size = 3.8,
             color = dir_colors[comp$ident2], fontface = "bold") +
    scale_fill_manual(values = dir_colors, name = "Enriched in") +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    scale_x_continuous(
      name   = expression(-log[10](p[adj])),
      labels = function(x) abs(round(x, 1)),
      expand = expansion(mult = c(0.02, 0.45))   # extra right room for gene text
    ) +
    labs(title = sprintf("%s | %s\nGO Biological Process Enrichment",
                         subset_name, gsub("_", " ", comp$label)),
         y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.y        = element_text(size = 13),
      axis.text.x        = element_text(size = 13),
      axis.title.x       = element_text(size = 14, face = "bold"),
      axis.line.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
      legend.position    = "right",
      plot.title         = element_text(hjust = 0.5, size = 15, face = "bold",
                                        lineheight = 1.2),
      plot.margin        = margin(t = 18, r = 10, b = 10, l = 10)
    ) +
    coord_cartesian(clip = "off")

  if (!is.null(output_file)) {
    n_bars <- nrow(res_all)
    ggsave(output_file, plot = p, width = 15, height = max(4, n_bars * 0.28 + 2),
           units = "in", device = "pdf", useDingbats = FALSE)
  }
  p
}
```

---

## make_pathway_barplot_msig — Hallmark + Reactome variant

Use this when the brief specifies Hallmark + Reactome (msigdbr) enrichment.
The canonical make_pathway_barplot() above uses GO:BP via enrichGO. If you
need Hallmark and/or Reactome, call make_pathway_barplot_msig() instead.
Do not write a third inline implementation.

```r
make_pathway_barplot_msig <- function(markers, comp, subset_name, universe_genes,
                                       n_pathways = N_PATHWAYS, output_file,
                                       flagged_keywords = NULL) {
  sig_df    <- markers %>%
    filter(p_val_adj < PADJ_CUT, abs(avg_log2FC) > LFC_CUT, !is_ambient(gene))
  genes_up1 <- sig_df %>% filter(avg_log2FC >  0) %>% pull(gene)
  genes_up2 <- sig_df %>% filter(avg_log2FC <= 0) %>% pull(gene)

  # Build TERM2GENE from msigdbr (Hallmark H + Reactome C2)
  msig_h  <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol)
  msig_r  <- msigdbr(species = "Homo sapiens", category = "C2",
                     subcategory = "CP:REACTOME") %>%
    dplyr::select(gs_name, gene_symbol)
  term2gene <- bind_rows(msig_h, msig_r) %>% distinct()

  run_ora <- function(genes, dir_label) {
    if (length(genes) < 5) return(NULL)
    res <- enricher(
      gene          = genes,
      universe      = universe_genes,
      TERM2GENE     = term2gene,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      minGSSize     = 10,
      maxGSSize     = 500
    )
    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    as.data.frame(res) %>%
      filter(p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      slice_head(n = n_pathways) %>%
      mutate(direction = dir_label)
  }

  res_all <- bind_rows(run_ora(genes_up1, comp$ident1),
                       run_ora(genes_up2, comp$ident2))

  if (is.null(res_all) || nrow(res_all) == 0) {
    message(sprintf("No significant pathways for %s -- barplot skipped", comp$label))
    return(NULL)
  }

  res_all <- res_all %>%
    mutate(
      neg_log_p  = -log10(p.adjust),
      x_val      = ifelse(direction == comp$ident2, neg_log_p, -neg_log_p),
      dot_x      = x_val * 1.03,
      label_desc = ifelse(nchar(Description) > 55,
                          paste0(substr(Description, 1, 52), "..."), Description),
      label_desc = gsub("_", " ", label_desc),
      top_genes  = sapply(strsplit(geneID, "/"), function(g) paste(head(g, 5), collapse = ", ")),
      text_x     = ifelse(x_val >= 0,
                          dot_x + max(abs(x_val)) * 0.03,
                          max(abs(x_val)) * 0.03),
      text_hjust = 0,
      gene_text  = sprintf("n=%d  %s", Count, top_genes),
      flagged    = if (!is.null(flagged_keywords))
                     grepl(flagged_keywords, tolower(Description))
                   else
                     FALSE
    ) %>%
    arrange(x_val) %>%
    mutate(label_desc = factor(label_desc, levels = unique(label_desc)))

  max_x <- max(abs(res_all$x_val))
  dir_colors <- setNames(c("#B2182B", "#2166AC"), c(comp$ident1, comp$ident2))

  p <- ggplot(res_all, aes(x = x_val, y = label_desc)) +
    geom_col(aes(fill = direction), width = 0.55, alpha = 0.82) +
    geom_vline(xintercept = 0, color = "grey25", linewidth = 0.55) +
    geom_vline(xintercept = c(-log10(0.05), log10(0.05)),
               color = "grey55", linewidth = 0.35, linetype = "dashed") +
    geom_point(aes(x = dot_x, fill = direction, size = Count),
               shape = 21, color = "grey20", stroke = 0.3, alpha = 0.95) +
    geom_text(aes(x = text_x, label = gene_text, hjust = text_hjust),
              size = 3.5, fontface = "italic", color = "grey25", lineheight = 0.9) +
    geom_text(data = res_all %>% filter(flagged),
              aes(x = 0, y = label_desc, label = "*"),
              inherit.aes = FALSE, size = 5, color = "darkred", hjust = 0.5) +
    scale_fill_manual(values = dir_colors, name = "Enriched in") +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    scale_x_continuous(
      name   = expression(-log[10](p[adj])),
      labels = function(x) abs(round(x, 1)),
      expand = expansion(mult = c(0.02, 0.45))
    ) +
    labs(title = sprintf("%s | %s\nHallmark + Reactome Pathway Enrichment",
                         subset_name, gsub("_", " ", comp$label)),
         y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.y        = element_text(size = 12),
      axis.text.x        = element_text(size = 13),
      axis.title.x       = element_text(size = 14, face = "bold"),
      axis.line.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
      legend.position    = "right",
      plot.title         = element_text(hjust = 0.5, size = 15, face = "bold",
                                        lineheight = 1.2),
      plot.margin        = margin(t = 18, r = 10, b = 10, l = 10)
    ) +
    coord_cartesian(clip = "off")

  n_bars <- nrow(res_all)
  ggsave(output_file, plot = p, width = 15, height = max(4, n_bars * 0.28 + 2),
         units = "in", device = "pdf", useDingbats = FALSE)
  message(sprintf("Pathway barplot: %d pathways, file: %s", n_bars, basename(output_file)))
  invisible(p)
}
```

Function signature mirrors make_pathway_barplot(): markers data frame with
columns gene/avg_log2FC/p_val_adj, comp list with label/ident1/ident2,
subset_name, universe_genes, n_pathways, output_file. The flagged_keywords
argument accepts a regex of pathway-name patterns to highlight with the
asterisk (e.g. "angiogenesis|immune|chemokine|ECM|remodeling|vascular|
cytokine|integrin|hypoxia").

---

## make_module_score_violin — per-subcluster module score violin

Canonical recipe for AddModuleScore visualization. Use this when showing
module score distributions across cell subtypes (e.g. EC subclusters)
with comparison between conditions or groups.

CRITICAL: x-axis is the subtype/subcluster column, fill is the comparison
column (Condition or Group). Do NOT put Condition or Group on the x-axis
-- that collapses the per-subtype dimension and loses the comparison
the plot is meant to show.

```r
# Module score violin: per-subtype breakdown, fill by condition or group.
#
# Args:
#   meta_df: data frame with at minimum subtype_col, fill_col, and score_col
#   score_col: name of the column containing the module score (numeric)
#   subtype_col: name of the column containing the subtype/subcluster labels
#                (this becomes the x-axis)
#   fill_col: name of the column containing the comparison variable
#             (this becomes the fill aesthetic -- typically Condition or Group)
#   fill_colors: named character vector of colors keyed by fill_col levels
#   output_file: PDF path
#   title: plot title (optional)
make_module_score_violin <- function(meta_df, score_col, subtype_col, fill_col,
                                     fill_colors, output_file,
                                     title = NULL) {
  df <- meta_df[, c(subtype_col, fill_col, score_col)]
  colnames(df) <- c("subtype", "fill_var", "score")

  n_sub <- length(unique(df$subtype))
  fig_w <- max(n_sub * 0.75 + 3, 6)

  p <- ggplot(df, aes(x = subtype, y = score, fill = fill_var)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.3, alpha = 0.8) +
    geom_boxplot(width = 0.05, outlier.shape = NA, fill = "white", linewidth = 0.2) +
    scale_fill_manual(values = fill_colors) +
    labs(title = title, x = NULL, y = "Module Score", fill = NULL) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11))

  ggsave(output_file, plot = p,
         width = fig_w, height = 5,
         device = "pdf", useDingbats = FALSE)
  invisible(p)
}
```

When calling with Condition as the fill, provide a 2-level color vector
(e.g. c(Normal = "#4477AA", Tumor = "#EE6677")). When calling with Group
as the fill, provide a 5-level color vector (Normal, PTC-FV, PTC-CT,
HCC, ATC). Project-specific color palettes should live in
context/lab_context.md or context/color_palettes.md as named vectors.

Always produce a separate plot per comparison axis: one violin grid
faceted on subtype with Condition fill, one with Group fill. Do not
combine them into a single plot.

---

## Adapter: make_pathway_barplot() for one-directional (signature) gene lists

When the input is a pre-filtered gene list rather than a bidirectional DE result,
construct a synthetic markers data frame and a direction-fixed comp object:

```r
# sig_genes: character vector of signature genes (already filtered, all positive direction)
# universe_genes: rownames(so) -- all expressed genes

# Build synthetic markers data frame
markers_synthetic <- data.frame(
  gene         = sig_genes,
  avg_log2FC   =  1.0,         # all positive -- will pass LFC_CUT filter
  p_val_adj    =  0.001,       # below PADJ_CUT -- all genes pass
  stringsAsFactors = FALSE
)

# Build synthetic comp object -- ident2 = NULL signals one-directional mode
comp_synthetic <- list(
  label  = "ThyroidEC_Signature",
  ident1 = "inhouse_thyroid",      # REPLACE: the enriched direction label
  ident2 = "other_organ",          # REPLACE: reference label (bars will be positive only)
  col    = "sig_group"
)

# dir_colors must be defined before calling make_pathway_barplot()
dir_colors <- c(
  "inhouse_thyroid" = "#4393C3",   # REPLACE: color for the enriched direction
  "other_organ"     = "#878787"    # REPLACE: color for reference (no bars, used for labels only)
)

make_pathway_barplot(
  markers       = markers_synthetic,
  comp          = comp_synthetic,
  subset_name   = "Thyroid EC",
  universe_genes = universe_genes,
  n_pathways    = 15,
  output_file   = "output/ec_organ_analysis/pathway/thyroid_EC_signature_pathway_barplot.pdf"
)
```

This adapter ensures the full visual layer specification (gene annotations, bubbles,
threshold lines, directional labels, coord_cartesian clip-off) is applied even for
one-directional use cases. The ident2 bars will be absent (no negative genes in the
synthetic markers) but the ident1 bars will carry all the patched visual layers.

Do NOT write an inline ggplot pathway barplot. Route to one of:
- make_pathway_barplot() for GO:BP (enrichGO-based; gene ontology biological process)
- make_pathway_barplot_msig() for Hallmark and/or Reactome (msigdbr-based)
- the one-directional adapter below if the input is a pre-filtered signature
  gene list rather than a bidirectional DE result
If none of these signatures fit, use the adapter rather than reinventing.

---

## Gene-set dotplots: filtered only

The sanctioned function for gene-set dotplots is `make_functional_dotplot()`. It filters
by `p_val_adj < PADJ_CUT` AND `|avg_log2FC| > LFC_CUT` (thresholds from the config block,
defaulting to 0.05 and 0.5) before plotting. Only genes that are BOTH in the gene set AND
significantly DE appear in the output.

**Rules:**
- Do NOT invent parallel inline functions (e.g., `make_geneset_dotplot()`) that skip the filter.
- If a gene set has zero genes passing the filter for a given comparison, `make_functional_dotplot()`
  returns NULL. The caller must log this to `output/decision_log.txt`:
  ```r
  p_fd <- make_functional_dotplot(...)
  if (is.null(p_fd)) {
    cat(sprintf("%s: %s -- no significant genes in gene set, plot skipped\n",
                Sys.time(), paste(comp$label, scope_name, sep = "_")),
        file = "output/decision_log.txt", append = TRUE)
  }
  ```
- The filter is not negotiable without explicit brief approval and deployment agent sign-off.

---

## Dotplot views: by Condition vs by Group

Both views (the "by Condition" 2-level view and the "by Group" multi-level view) filter
the gene set using the SAME DE result: the primary condition-vs-condition comparison
(e.g., Tumor vs Normal).

- **"by Condition" view:** columns = ident1, ident2 (2 levels); rows = significantly DE genes
  from the primary comparison filtered by `make_functional_dotplot()`.
- **"by Group" view:** columns = all unique group values (N levels); rows = the SAME filtered
  genes as the "by Condition" view. The group view shows how those genes distribute across
  all groups, not which genes differ across groups.

**Do NOT** filter the "by Group" view by a union of pairwise group comparisons or by a
multi-group test. Use the primary condition comparison filter for both views. This ensures:
- The same gene set appears in both views, enabling direct cross-view comparison.
- No additional DE computation is required for the expanded view.
- The interpretation is coherent: "among genes that differ in Tumor vs Normal, how does
  each group distribute?"

---

## Main Analysis Loop

```r
# Requires: comparisons list, cell_subsets list, LABEL_COL, group_colors, subtype_colors,
#           functional_gene_sets (project_specific -- must be defined in CLAUDE.md context)
# Axis C requires: comp$per_cluster = TRUE, comp$label_col set to cluster/subtype column

for (comp in comparisons) {
  if (isTRUE(comp$per_cluster)) {
    # Axis C: condition within each cluster
    stopifnot(!is.null(comp$label_col), comp$label_col %in% colnames(obj@meta.data))
    label_vals <- sort(unique(obj@meta.data[[comp$label_col]]))
    for (lv in label_vals) {
      so_lv <- subset(obj, cells = which(obj@meta.data[[comp$label_col]] == lv))
      n1 <- sum(so_lv@meta.data[[comp$col]] == comp$ident1)
      n2 <- sum(so_lv@meta.data[[comp$col]] == comp$ident2)
      if (n1 < MIN_CELLS_PER_CLUSTER || n2 < MIN_CELLS_PER_CLUSTER) {
        message(sprintf("SKIPPED Axis C %s within %s: too few cells (%d / %d, min %d)",
                        comp$label, lv, n1, n2, MIN_CELLS_PER_CLUSTER))
        next
      }
      markers <- run_findmarkers(so_lv, comp)
      if (is.null(markers)) next

      out_prefix <- file.path(output_dir, sprintf("%s_within_%s", comp$label, lv))
      dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
      write.csv(markers, paste0(out_prefix, "_DE_full.csv"), row.names = FALSE)
      message(sprintf("Axis C: %s within %s -- %d sig genes (padj<%.2f, |lfc|>%.1f)",
                      comp$label, lv,
                      sum(markers$p_val_adj < PADJ_CUT & abs(markers$avg_log2FC) > LFC_CUT),
                      PADJ_CUT, LFC_CUT))
    }
  } else {
    # Axis B: condition global (scope-based loop)
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

      # Functional plots -- only if functional_gene_sets is defined
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

        p_fd <- make_functional_dotplot(so_scope, markers, comp, scope_name,
                                 label_col = LABEL_COL, group_col = comp$col,
                                 functional_gene_sets = functional_gene_sets,
                                 subtype_colors = subtype_colors, show_direction = FALSE,
                                 output_file = paste0(out_prefix, "_functional_dotplot.pdf"))
        if (is.null(p_fd)) {
          cat(sprintf("%s: %s -- no significant genes in gene set, functional dotplot skipped\n",
                      Sys.time(), paste(scope_name, comp$label, sep = "_")),
              file = "output/decision_log.txt", append = TRUE)
        }
      }

      # Pathway barplot
      universe <- rownames(so_scope)
      make_pathway_barplot(markers, comp, scope_name, universe_genes = universe,
                            output_file = paste0(out_prefix, "_pathways.pdf"))
    }
  }
}
```

---

## Figure Sizing Reference

| Plot | Width | Height |
|---|---|---|
| Volcano | 9 | 7 |
| Overall heatmap | 12 | 10 |
| Functional heatmap | 12 | `max(8, n_genes * 0.28 + 3)` |
| Top-gene dot plot | `max(8, n_ctypes * 1.2 + 3)` | `max(6, 20 * 0.32 + 2.5)` |
| Functional dot plot | `n_ctypes * 0.55 * 2 + 7` | `max(6, n_fg * 0.28 + 3)` |
| Functional dot plot (with dir strip) | above + 1 | same |
| Pathway bar plot | 15 | `max(4, n_bars * 0.28 + 2)` |

---

## Output File Naming Convention

```
DE_full_{subset}_{CompLabel}.csv
volcano_{subset}_{CompLabel}.pdf
heatmap_overall_{subset}_{CompLabel}.pdf
heatmap_functional_{subset}_{CompLabel}.pdf
dotplot_topgenes_{subset}_{CompLabel}.pdf
dotplot_functional_{subset}_{CompLabel}.pdf
dotplot_functional_dir_{subset}_{CompLabel}.pdf
pathway_barplot_{subset}_{CompLabel}.pdf
{project}_proposed_gene_lists.yaml
```

`{CompLabel}` is `comp$label`, e.g. `ViscFat_vs_SubQFat`. Underscores only — no
spaces or special characters in filenames.

---

## Brief Configuration

```yaml
optional_analyses:
  deg:
    enabled: true
    axes:
      condition_global: true        # Axis B: condition comparison across all cells / scopes
      condition_per_cluster: true   # Axis C: condition comparison within each cluster value
    label_col: "subtype_label"      # required for Axis C; column holding cluster/subtype labels
    comparisons:
      - label: "GroupA_vs_GroupB"
        ident1: "Group A"           # shown RIGHT in heatmap, positive log2FC direction
        ident2: "Group B"           # shown LEFT  in heatmap, negative log2FC direction
        col: "tissue_type"          # metadata column holding group labels
      - label: "Tumor_vs_Normal"
        ident1: "Tumor"
        ident2: "Normal"
        col: "Condition"
        per_cluster: true           # opts this comparison into Axis C iteration
    scopes:
      - name: AllCells
        cell_types: null            # null = use whole subset object
      - name: SubtypeA
        cell_types: [SubtypeA1, SubtypeA2]
      - name: SubtypeB
        cell_types: [SubtypeB1, SubtypeB2, SubtypeB3]
    thresholds:
      padj: 0.05
      lfc:  0.5
      n_cells_heatmap: 300
      n_top_genes_heatmap: 30
      n_label_volcano: 25
      min_cells: 10                 # Axis B gate: minimum cells per side for global comparison
      min_cells_per_cluster: 100    # Axis C gate: minimum cells per side per cluster
    functional_gene_sets: project_specific   # define inline in CLAUDE.md
    pathway_analysis: true
    n_pathways: 15
```
