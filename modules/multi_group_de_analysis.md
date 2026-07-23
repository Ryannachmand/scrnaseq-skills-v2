---
requires_context:
  palettes:
    - group_colors    # required: named vector keyed by group values (for direction strip)
    - label_colors    # optional: named vector for x-axis cell type label coloring
  metadata_columns:
    required:
      - label_col     # cell type label column in Seurat metadata
      - group_col     # comparison group column (anatomical site, condition, time point, etc.)
    optional:
      - label_order   # ordered character vector of cell type names (default: alphabetical)
  brief_keys:
    required:
      - output_dir
      - downstream_analyses.multi_group_de.groups
      - downstream_analyses.multi_group_de.group_col
      - downstream_analyses.multi_group_de.label_col
      - downstream_analyses.multi_group_de.comparisons
    optional:
      - downstream_analyses.multi_group_de.functional_gene_sets  # null = skip functional plots
      - downstream_analyses.multi_group_de.label_order
      - downstream_analyses.multi_group_de.label_colors
      - downstream_analyses.multi_group_de.n_each
      - downstream_analyses.multi_group_de.padj_cut
      - downstream_analyses.multi_group_de.lfc_cut
references:
  - "@primitives/differential_expression.md"
  - "@primitives/visualization.md"
  - "@primitives/aesthetics.md"
  - "@primitives/seurat_v5_rules.md"
---

# Module: Multi-Group DE Analysis — Cell-Type-Resolved Dot Plots + GO Functional Plots

## Scope: N-Group Objects with Pairwise Comparisons

This module operates on a Seurat object containing cells from **N ≥ 2 discrete groups**
(anatomical sites, conditions, time points, or any categorical grouping). Comparisons are
always **pairwise** — each comparison is a `list(label, ident1, ident2)` specifying two
groups. Visualizations show expression per (cell_type × group) across all N panels
simultaneously, disentangling DE composition effects (different subtype abundances) from
expression effects (a subtype expressing a gene at different levels).

**This module does NOT perform N-way differential testing.** For a single ranked gene list
across all groups simultaneously, use `@modules/de_comprehensive_csv.md`.

---

## Function Name Collision Resolution

`make_topgene_dotplot()` is defined in `@primitives/differential_expression.md` with a
parameterized N-group signature. This module defines **`make_anatomical_dotplot()`** — a
separate function that extends the primitive's pattern with `groups` (explicit ordered
factor levels for the facet variable), `label_order` (x-axis cell type ordering), and
`label_colors` (per-cell-type axis label coloring). The primitive's `make_topgene_dotplot`
handles generic 2+-group comparisons; `make_anatomical_dotplot` handles the case where you
want canonical group ordering enforced and cell type axis labels colored.

`make_functional_dotplot()` from `@primitives/differential_expression.md` already supports
N groups dynamically (pins section labels to the last group level via
`tail(levels(factor(dot_df$group_var)), 1)`). It is referenced here without modification.
Set `GROUP_COL` factor levels to `GROUPS` on the Seurat object before calling — this
propagates the canonical ordering into the facet layout.

**`is_ambient()` must be loaded from `@primitives/differential_expression.md`** — the full
dual-filter version with both `AMBIENT_PATTERNS` and `AMBIENT_EXPLICIT`. Do not redefine
it here. The v1 source file used a patterns-only version that missed MALAT1 and platelet
markers; v2 resolves this by injecting the primitive.

---

## When to Use This Module

- N ≥ 2 discrete groups; want cell-type-resolved DE visualization per group
- Groups differ in composition and you want to distinguish "more Type-X cells" from "Type-X expresses differently"
- Need multi-panel dot plots with all N groups visible side by side
- Want data-driven GO:BP sections (make_go_functional_dotplot) derived from union DE

---

## Brief Schema

```yaml
downstream_analyses:
  multi_group_de:
    enabled: true
    group_col: project_specific           # REQUIRED: metadata column holding group values
    label_col: project_specific           # REQUIRED: metadata column holding cell type labels
    groups: project_specific              # REQUIRED: character vector of group levels, display order, N >= 2
    comparisons:                          # REQUIRED: list of pairwise comparison specs
      - label: "GroupA_vs_GroupB"         # filename-safe ASCII label — no Unicode, no em-dash
        ident1: project_specific          # group value for positive log2FC (shown right/top)
        ident2: project_specific          # group value for negative log2FC (shown left/bottom)
    functional_gene_sets: null              # OPTIONAL: named list of gene vectors, biology-specific
                                            # if null: make_functional_dotplot and make_functional_heatmap
                                            # steps are skipped entirely; base DE + dotplot still runs
    n_each: 12                            # top N genes per direction in make_anatomical_dotplot
    padj_cut: 0.05
    lfc_cut: 0.5
    label_order: null                     # optional: ordered cell type vector for x-axis ordering
                                          # if null: alphabetical
    label_colors: null                    # optional: named color vector for x-axis label coloring
                                          # if null: all labels grey30
```

---

## Critical Constraints

| ❌ Don't | ✅ Do instead | Why |
|---|---|---|
| Subset object to ident1/ident2 before make_anatomical_dotplot | Keep ALL cells in scope object | Filtering to the comparison pair drops the remaining group panels |
| Hardcode group factor levels | Pass `groups` argument; set factor on metadata | Hardcoded levels break for any other project |
| Use em-dash in comp$label | Use ASCII dash or underscore | `mbcsToSbcs` rendering error in grid PDF renderer |
| Call `is_ambient()` defined locally | Load from @primitives/differential_expression.md | v1 patterns-only version misses MALAT1 and platelet markers |
| Use `deframe()` from tibble | `{ setNames(.[[2]], .[[1]]) }` in a pipe block | `deframe` not always available |
| Omit `limitsize = FALSE` on ggsave | Add `limitsize = FALSE` | ggplot2 errors if height > 50 in |
| Use `select()` without namespace | Always `dplyr::select()` when Seurat loaded | Seurat masks dplyr's `select` |
| Use `annotate("text")` with `facet_wrap` | `geom_text(data = df, inherit.aes = FALSE)` pinned to one panel | `annotate()` renders in every panel |

---

## R Packages Required

```r
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)    # make_overall_heatmap, make_functional_heatmap
library(circlize)          # colorRamp2
library(clusterProfiler)   # make_go_functional_dotplot
library(org.Hs.eg.db)      # or org.Mm.eg.db for mouse
library(ggrepel)           # make_volcano
```

---

## Configuration Block

```r
# ── CONFIG ────────────────────────────────────────────────────────────────────
GROUPS     <- c("project_specific")  # REPLACE: character vector of group levels in display order
                                      # e.g. c("GroupA", "GroupB", "GroupC") — N >= 2
GROUP_COL  <- "project_specific"     # REPLACE: metadata column holding group values
LABEL_COL  <- "project_specific"     # REPLACE: metadata column holding cell type labels

LABEL_ORDER <- NULL   # REPLACE: optional ordered character vector of cell type names
                       # e.g. c("TypeA", "TypeB", "TypeC")
                       # if NULL: cell types appear in alphabetical order on x-axis
                       # NOTE: any project-specific ordered label vector (e.g. STROMA_ORDER,
                       #       CELLTYPE_ORDER) must be defined in that project's examples/ file
                       #       and passed in — the module does not hardcode label orderings.

LABEL_COLORS <- NULL  # REPLACE: optional named color vector for x-axis label coloring
                       # e.g. c("TypeA" = "#E41A1C", "TypeB" = "#377EB8")
                       # if NULL: all labels grey30

GROUP_COLORS <- c("project_specific" = "#project_specific")  # REPLACE: named vector, one entry per group value

# ── Pairwise comparisons within the N-group object ────────────────────────────
comparisons <- list(
  comp1 = list(
    label  = "GroupA_vs_GroupB",       # filename-safe ASCII label
    ident1 = "project_specific",       # REPLACE: positive log2FC group
    ident2 = "project_specific"        # REPLACE: negative log2FC group
  )
  # Add additional list() entries for additional pairwise comparisons
)

# ── Cell-type scopes ──────────────────────────────────────────────────────────
# NULL = all cells in object; character vector = subset to these label values
cell_subsets <- list(
  AllCells = NULL
  # FocusedSubset = c("project_specific")
)

# ── Functional gene sets — caller-provided, biology-specific ──────────────────
# Required for make_functional_dotplot and make_functional_heatmap.
# If not defined, those steps are skipped. EC biology sets belong in examples/.
functional_gene_sets <- list(
  "Pathway A" = c("project_specific"),  # REPLACE: gene vectors for your biology
  "Pathway B" = c("project_specific")
)

# ── Thresholds ────────────────────────────────────────────────────────────────
PADJ_CUT <- 0.05
LFC_CUT  <- 0.5
N_EACH   <- 12     # top genes per direction in make_anatomical_dotplot

# ── GO functional dotplot (make_go_functional_dotplot) ───────────────────────
N_GO_TERMS       <- 8    # top GO:BP terms as sections after simplify
N_GENES_PER_SECT <- 8    # max genes per GO section (top by mean |log2FC| across comparisons)

# ── Paths ─────────────────────────────────────────────────────────────────────
RDS_IN     <- "project_specific"   # REPLACE: path to Seurat subset object
OUTPUT_DIR <- file.path("output", "multi_group_de")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
```

---

## Load Object

```r
so <- readRDS(RDS_IN)
so <- JoinLayers(so)          # @primitives/seurat_v5_rules.md Rule 1
DefaultAssay(so) <- "RNA"

# Set group factor levels — required before calling plot functions so facets appear in
# the canonical order defined by GROUPS
so@meta.data[[GROUP_COL]] <- factor(so@meta.data[[GROUP_COL]], levels = GROUPS)

# Resolve subtype_colors for primitive function calls
subtype_colors_use <- if (!is.null(LABEL_COLORS)) LABEL_COLORS else {
  ct_levels <- sort(unique(so@meta.data[[LABEL_COL]]))
  setNames(rep("grey30", length(ct_levels)), ct_levels)
}
```

---

## Function: make_anatomical_dotplot

Extends the primitive's dotplot pattern with `groups` (ordered factor levels for facets),
`label_order` (x-axis cell type ordering and filtering), and `label_colors` (per-type
axis label coloring). Fetches from ALL cells in the scope object so all N group panels
have data — do NOT pre-filter to the comparison pair.

`is_ambient()` must be loaded from `@primitives/differential_expression.md`.

```r
make_anatomical_dotplot <- function(so_obj, markers, comp, subset_name,
                                     label_col,
                                     group_col,
                                     groups,             # all N group levels, display order
                                     label_order  = NULL,
                                     label_colors = NULL,
                                     n_each       = 12,
                                     output_file) {
  sig_df <- markers %>%
    filter(p_val_adj < PADJ_CUT, abs(avg_log2FC) > LFC_CUT, !is_ambient(gene))

  top_up1 <- sig_df %>% filter(avg_log2FC >  0) %>% arrange(p_val_adj) %>% head(n_each) %>% pull(gene)
  top_up2 <- sig_df %>% filter(avg_log2FC <= 0) %>% arrange(p_val_adj) %>% head(n_each) %>% pull(gene)
  gene_order <- c(top_up1, top_up2)
  if (length(gene_order) == 0) {
    message(sprintf("No significant genes for %s — skipping anatomical dotplot.", comp$label))
    return(NULL)
  }

  # Fetch from ALL cells in scope object — critical: do NOT pre-filter to comparison pair
  expr <- as.data.frame(t(as.matrix(
    GetAssayData(so_obj, assay = "RNA", layer = "data")[gene_order, ]
  )))
  expr$group_var   <- so_obj@meta.data[[group_col]]
  expr$subtype_var <- so_obj@meta.data[[label_col]]

  dot_df <- tidyr::pivot_longer(expr, cols = all_of(gene_order),
                                 names_to = "gene", values_to = "expression") %>%
    group_by(group_var, subtype_var, gene) %>%
    summarise(
      avg_exp = mean(expm1(expression)),
      pct_exp = mean(expression > 0) * 100,
      .groups = "drop"
    ) %>%
    group_by(gene) %>%
    mutate(avg_exp_scaled = pmax(pmin(scale(avg_exp)[, 1], 2.5), -2.5)) %>%
    ungroup()

  # Enforce group ordering — all N panels appear in canonical order
  dot_df$group_var <- factor(dot_df$group_var, levels = groups)

  # Apply label_order to filter and order cell types on x-axis (alphabetical when NULL)
  ct_present <- if (!is.null(label_order)) {
    intersect(label_order, unique(dot_df$subtype_var))
  } else {
    sort(unique(dot_df$subtype_var))
  }
  dot_df$subtype_var <- factor(dot_df$subtype_var, levels = ct_present)
  dot_df <- dot_df %>% filter(!is.na(subtype_var))

  # Diagonal gene ordering within each DE-direction group (uses x-axis column order for peak)
  avg_wide <- tapply(dot_df$avg_exp_scaled,
                     list(as.character(dot_df$gene), as.character(dot_df$subtype_var)),
                     mean)
  avg_wide <- avg_wide[, ct_present[ct_present %in% colnames(avg_wide)], drop = FALSE]
  avg_wide[is.na(avg_wide)] <- 0
  ct_ref <- colnames(avg_wide)
  up1_in <- top_up1[top_up1 %in% rownames(avg_wide)]
  if (length(up1_in) > 0) {
    m <- avg_wide[up1_in, , drop = FALSE]
    pg <- ct_ref[max.col(m, ties.method = "first")]
    pv <- m[cbind(seq_len(nrow(m)), max.col(m, ties.method = "first"))]
    up1_in <- up1_in[order(match(pg, ct_ref), -pv, up1_in)]
  }
  up2_in <- top_up2[top_up2 %in% rownames(avg_wide)]
  if (length(up2_in) > 0) {
    m <- avg_wide[up2_in, , drop = FALSE]
    pg <- ct_ref[max.col(m, ties.method = "first")]
    pv <- m[cbind(seq_len(nrow(m)), max.col(m, ties.method = "first"))]
    up2_in <- up2_in[order(match(pg, ct_ref), -pv, up2_in)]
  }
  gene_order_new <- c(up1_in, up2_in)
  dot_df$gene <- factor(dot_df$gene, levels = rev(gene_order_new))

  axis_colors <- if (!is.null(label_colors)) {
    cols <- label_colors[ct_present]
    cols[is.na(cols)] <- "grey30"
    cols
  } else {
    rep("grey30", length(ct_present))
  }

  n_groups   <- length(groups)
  n_subtypes <- length(ct_present)
  n_genes    <- length(gene_order_new)
  divider_y  <- length(up2_in) + 0.5

  p <- ggplot(dot_df, aes(x = subtype_var, y = gene, size = pct_exp, fill = avg_exp_scaled)) +
    geom_point(shape = 21, color = "grey30", stroke = 0.32) +
    geom_hline(yintercept = divider_y, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    facet_wrap(~ group_var, ncol = n_groups, scales = "free_x") +
    scale_fill_gradientn(colors = c("#F5F5F5", "#FFF9C4", "#FFB300", "#E53935"),
                         limits = c(-2.5, 2.5), oob = scales::squish) +
    scale_size_continuous(range = c(0.3, 6), limits = c(0, 100)) +
    scale_x_discrete(position = "top") +
    labs(title = sprintf("%s | %s", subset_name, gsub("_", " ", comp$label)),
         x = NULL, y = NULL) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x.top = element_text(size = 11, angle = 45, hjust = 0, vjust = 0,
                                      face = "bold", color = axis_colors),
      axis.text.y     = element_text(size = 10, face = "italic"),
      strip.text      = element_text(size = 13, face = "bold"),
      panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      legend.position = "right"
    )

  fig_h <- max(6, n_genes * 0.25 + 3)
  fig_w <- n_subtypes * 0.55 * n_groups + 7
  ggsave(output_file, plot = p, width = fig_w, height = fig_h,
         units = "in", device = "pdf", useDingbats = FALSE, limitsize = FALSE)
  p
}
```

---

## Main Analysis Loop

References `run_findmarkers`, `make_volcano`, `make_overall_heatmap`,
`make_functional_heatmap`, `make_functional_dotplot`, `make_pathway_barplot` from
`@primitives/differential_expression.md`. Calls `make_anatomical_dotplot` defined above.
`all_markers` accumulates results for `make_go_functional_dotplot` downstream.

```r
all_markers <- list()

for (scope_name in names(cell_subsets)) {
  scope_labels <- cell_subsets[[scope_name]]
  so_scope <- if (is.null(scope_labels)) so else {
    subset(so, cells = rownames(so@meta.data)[so@meta.data[[LABEL_COL]] %in% scope_labels])
  }
  # Propagate group factor levels to the scoped object
  so_scope@meta.data[[GROUP_COL]] <- factor(so_scope@meta.data[[GROUP_COL]], levels = GROUPS)

  for (comp in comparisons) {
    markers <- run_findmarkers(so_scope, comp)   # @primitives/differential_expression.md
    if (is.null(markers)) next

    out_pfx <- file.path(OUTPUT_DIR, sprintf("%s_%s", scope_name, comp$label))

    write.csv(markers, paste0(out_pfx, "_DE_full.csv"), row.names = FALSE)
    all_markers[[paste(scope_name, comp$label, sep = "_")]] <- markers

    # Volcano
    p_vol <- make_volcano(markers, comp, scope_name,
                           functional_gene_sets = if (exists("functional_gene_sets")) functional_gene_sets else NULL)
    ggsave(paste0(out_pfx, "_volcano.pdf"), plot = p_vol, width = 9, height = 7,
           units = "in", device = "pdf", useDingbats = FALSE)

    # Overall heatmap — ComplexHeatmap: pdf()/draw()/dev.off() required (CONVENTIONS.md §4 exception #1)
    ht_overall <- make_overall_heatmap(so_scope, markers, comp, scope_name,
                                        label_col      = LABEL_COL,
                                        subtype_colors = subtype_colors_use,
                                        group_colors   = GROUP_COLORS)
    if (!is.null(ht_overall)) {
      pdf(paste0(out_pfx, "_heatmap_overall.pdf"), width = 12, height = 10)
      ComplexHeatmap::draw(ht_overall)
      dev.off()
    }

    # Cell-type-resolved dotplot (this module's function)
    make_anatomical_dotplot(so_scope, markers, comp, scope_name,
                             label_col    = LABEL_COL,
                             group_col    = GROUP_COL,
                             groups       = GROUPS,
                             label_order  = LABEL_ORDER,
                             label_colors = LABEL_COLORS,
                             n_each       = N_EACH,
                             output_file  = paste0(out_pfx, "_anatomical_dotplot.pdf"))

    # Functional plots — only when functional_gene_sets is defined
    if (exists("functional_gene_sets") && !is.null(functional_gene_sets)) {
      ht_func <- make_functional_heatmap(so_scope, markers, comp, scope_name,
                                          label_col            = LABEL_COL,
                                          subtype_colors       = subtype_colors_use,
                                          group_colors         = GROUP_COLORS,
                                          functional_gene_sets = functional_gene_sets)
      if (!is.null(ht_func)) {
        pdf(paste0(out_pfx, "_heatmap_functional.pdf"), width = 12, height = 10)
        ComplexHeatmap::draw(ht_func)
        dev.off()
      }

      # make_functional_dotplot from @primitives: N-group aware; section labels auto-pin
      # to last group level (tail(levels(factor(group_var)), 1)) which equals GROUPS[length(GROUPS)]
      # because factor levels are set to GROUPS above.
      make_functional_dotplot(so_scope, markers, comp, scope_name,
                               label_col            = LABEL_COL,
                               group_col            = GROUP_COL,
                               functional_gene_sets = functional_gene_sets,
                               subtype_colors       = subtype_colors_use,
                               show_direction       = FALSE,
                               output_file          = paste0(out_pfx, "_functional_dotplot.pdf"))

      make_functional_dotplot(so_scope, markers, comp, scope_name,
                               label_col            = LABEL_COL,
                               group_col            = GROUP_COL,
                               functional_gene_sets = functional_gene_sets,
                               subtype_colors       = subtype_colors_use,
                               show_direction       = TRUE,
                               output_file          = paste0(out_pfx, "_functional_dotplot_dir.pdf"))
    }

    make_pathway_barplot(markers, comp, scope_name,
                          universe_genes = rownames(so_scope),
                          output_file    = paste0(out_pfx, "_pathway_barplot.pdf"))
  }
}
```

---

## Data-Driven GO Functional Dot Plot — make_go_functional_dotplot

Uses the union of significant DE genes across all comparisons for this scope. GO:BP
sections are derived from ORA rather than pre-curated lists. Direction strip uses one
color per group (N-way), not binary up/down. Call after the main loop above.

```r
make_go_functional_dotplot <- function(so_obj, scope_markers, scope_name,
                                        label_col, group_col, groups,
                                        group_colors,       # named vector: group → hex color
                                        label_order  = NULL,
                                        label_colors = NULL,
                                        n_go_terms       = N_GO_TERMS,
                                        n_genes_per_sect = N_GENES_PER_SECT,
                                        output_dir) {
  if (length(scope_markers) == 0) { message("No markers for scope: ", scope_name); return(NULL) }

  union_sig <- bind_rows(lapply(scope_markers, function(m) {
    m %>% filter(p_val_adj < PADJ_CUT, abs(avg_log2FC) > LFC_CUT, !is_ambient(gene)) %>%
      dplyr::select(gene)
  })) %>% distinct() %>% pull(gene)

  if (length(union_sig) < 5) { message("Too few union DE genes for GO analysis."); return(NULL) }

  ego <- clusterProfiler::enrichGO(
    gene          = union_sig,
    universe      = rownames(so_obj),
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  if (is.null(ego) || nrow(ego@result) == 0) { message("No GO:BP terms enriched."); return(NULL) }
  ego <- clusterProfiler::simplify(ego, cutoff = 0.5)

  top_terms <- head(ego@result, n_go_terms)$Description

  gene_to_term <- bind_rows(lapply(top_terms, function(term) {
    data.frame(
      gene    = strsplit(ego@result$geneID[ego@result$Description == term], "/")[[1]],
      section = term,
      stringsAsFactors = FALSE
    )
  })) %>%
    filter(gene %in% union_sig) %>%
    filter(!duplicated(gene))   # first-occurrence wins; term rank controls priority

  mean_lfc_df <- bind_rows(lapply(scope_markers, function(m) {
    m %>% filter(gene %in% gene_to_term$gene) %>% dplyr::select(gene, avg_log2FC)
  })) %>%
    group_by(gene) %>%
    summarise(mean_abs_lfc = mean(abs(avg_log2FC)), .groups = "drop")

  gene_order <- gene_to_term %>%
    left_join(mean_lfc_df, by = "gene") %>%
    group_by(section) %>%
    arrange(desc(mean_abs_lfc)) %>%
    slice_head(n = n_genes_per_sect) %>%
    ungroup() %>%
    mutate(section = factor(section, levels = top_terms)) %>%
    arrange(section) %>%
    pull(gene)

  if (length(gene_order) == 0) { message("No genes surviving GO section selection."); return(NULL) }

  # N-way direction: group with highest average expression per gene
  expr_avg <- as.data.frame(t(as.matrix(
    GetAssayData(so_obj, assay = "RNA", layer = "data")[gene_order, ]
  )))
  expr_avg$group_var <- so_obj@meta.data[[group_col]]

  dir_lookup <- tidyr::pivot_longer(expr_avg, cols = all_of(gene_order),
                                     names_to = "gene", values_to = "expression") %>%
    group_by(group_var, gene) %>%
    summarise(avg_exp = mean(expm1(expression)), .groups = "drop") %>%
    group_by(gene) %>%
    slice_max(avg_exp, n = 1, with_ties = FALSE) %>%
    dplyr::select(gene, group_var) %>%
    { setNames(.[[2]], .[[1]]) }

  ct_present <- if (!is.null(label_order)) {
    intersect(label_order, unique(so_obj@meta.data[[label_col]]))
  } else {
    sort(unique(so_obj@meta.data[[label_col]]))
  }

  expr2 <- as.data.frame(t(as.matrix(
    GetAssayData(so_obj, assay = "RNA", layer = "data")[gene_order, ]
  )))
  expr2$group_var   <- so_obj@meta.data[[group_col]]
  expr2$subtype_var <- so_obj@meta.data[[label_col]]

  dot_df <- tidyr::pivot_longer(expr2, cols = all_of(gene_order),
                                  names_to = "gene", values_to = "expression") %>%
    group_by(group_var, subtype_var, gene) %>%
    summarise(avg_exp = mean(expm1(expression)),
              pct_exp = mean(expression > 0) * 100, .groups = "drop") %>%
    group_by(gene) %>%
    mutate(avg_exp_scaled = pmax(pmin(scale(avg_exp)[, 1], 2.5), -2.5)) %>%
    ungroup() %>%
    mutate(group_var   = factor(group_var, levels = groups),
           subtype_var = factor(subtype_var, levels = ct_present)) %>%
    filter(!is.na(subtype_var))

  # Diagonal gene ordering within each GO section (preserves gene selection; changes within-section order)
  avg_by_subtype <- tapply(dot_df$avg_exp_scaled,
                           list(as.character(dot_df$gene), as.character(dot_df$subtype_var)),
                           mean)
  ct_ref <- ct_present[ct_present %in% colnames(avg_by_subtype)]
  avg_by_subtype <- avg_by_subtype[, ct_ref, drop = FALSE]
  avg_by_subtype[is.na(avg_by_subtype)] <- 0
  gene_order <- unlist(lapply(top_terms, function(s) {
    g <- gene_order[gene_order %in% gene_to_term$gene[gene_to_term$section == s]]
    if (length(g) > 0) {
      mat_s <- avg_by_subtype[g, , drop = FALSE]
      pg    <- ct_ref[max.col(mat_s, ties.method = "first")]
      pv    <- mat_s[cbind(seq_len(nrow(mat_s)), max.col(mat_s, ties.method = "first"))]
      g     <- g[order(match(pg, ct_ref), -pv, g)]
    }
    g
  }))

  dot_df <- dot_df %>% mutate(gene = factor(gene, levels = rev(gene_order)))

  n_groups   <- length(groups)
  n_subtypes <- length(ct_present)
  n_genes    <- length(gene_order)

  axis_colors <- if (!is.null(label_colors)) {
    cols <- label_colors[ct_present]; cols[is.na(cols)] <- "grey30"; cols
  } else rep("grey30", n_subtypes)

  section_label_df <- gene_to_term %>%
    filter(gene %in% gene_order) %>%
    mutate(gene = factor(gene, levels = rev(gene_order))) %>%
    group_by(section) %>%
    summarise(y = mean(as.numeric(gene)), .groups = "drop") %>%
    mutate(label     = section,
           x         = n_subtypes + 0.65,
           group_var = factor(groups[n_groups], levels = groups))   # rightmost panel

  p <- ggplot(dot_df, aes(x = subtype_var, y = gene, size = pct_exp, fill = avg_exp_scaled)) +
    geom_point(shape = 21, color = "grey30", stroke = 0.32) +
    facet_wrap(~ group_var, ncol = n_groups, scales = "free_x") +
    scale_fill_gradientn(colors = c("#F5F5F5", "#FFF9C4", "#FFB300", "#E53935"),
                         limits = c(-2.5, 2.5), oob = scales::squish) +
    scale_size_continuous(range = c(0.3, 6), limits = c(0, 100)) +
    scale_x_discrete(position = "top") +
    geom_text(data = section_label_df,
              aes(x = x, y = y, label = label, group = group_var),
              inherit.aes = FALSE, size = 4, hjust = 0,
              fontface = "bold.italic", color = "grey30") +
    labs(title = sprintf("%s | GO Functional", scope_name), x = NULL, y = NULL) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x.top = element_text(size = 11, angle = 45, hjust = 0, vjust = 0,
                                      face = "bold", color = axis_colors),
      axis.text.y     = element_text(size = 10, face = "italic"),
      strip.text      = element_text(size = 13, face = "bold"),
      panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      plot.margin     = margin(t = 8, r = 140, b = 8, l = 8),
      legend.position = "none"
    ) +
    coord_cartesian(clip = "off")

  dir_strip_df <- data.frame(
    gene      = factor(rev(gene_order), levels = rev(gene_order)),
    x         = 1L,
    direction = dir_lookup[rev(gene_order)]
  )
  p_strip <- ggplot(dir_strip_df, aes(x = x, y = gene, fill = direction)) +
    geom_tile(width = 0.7, height = 0.85, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = group_colors, name = "Highest\nExpr in") +
    scale_y_discrete(limits = rev(gene_order)) +
    scale_x_continuous(limits = c(0.35, 1.65), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(legend.position = "right",
          plot.margin = margin(t = 26, r = 5, b = 4, l = 4))

  p_dir <- patchwork::wrap_plots(
    p + theme(plot.margin = margin(t = 8, r = 5, b = 8, l = 8)),
    p_strip,
    widths = c(n_subtypes * n_groups, 1)
  )

  ph <- max(6, n_genes * 0.17 + 3)
  pw <- n_subtypes * 0.55 * n_groups + 7

  ggsave(file.path(output_dir, sprintf("GO_functional_dotplot_%s.pdf", scope_name)),
         p, width = pw, height = ph, units = "in",
         device = "pdf", useDingbats = FALSE, limitsize = FALSE)
  ggsave(file.path(output_dir, sprintf("GO_functional_dotplot_dir_%s.pdf", scope_name)),
         p_dir, width = pw + 1, height = ph, units = "in",
         device = "pdf", useDingbats = FALSE, limitsize = FALSE)
  write.csv(gene_to_term %>% filter(gene %in% gene_order),
            file.path(output_dir, sprintf("GO_sections_%s.csv", scope_name)),
            row.names = FALSE)
  invisible(p_dir)
}

# Invocation
for (scope_name in names(cell_subsets)) {
  scope_labels <- cell_subsets[[scope_name]]
  so_scope <- if (is.null(scope_labels)) so else {
    subset(so, cells = rownames(so@meta.data)[so@meta.data[[LABEL_COL]] %in% scope_labels])
  }
  so_scope@meta.data[[GROUP_COL]] <- factor(so_scope@meta.data[[GROUP_COL]], levels = GROUPS)

  scope_mk <- all_markers[grep(paste0("^", scope_name, "_"), names(all_markers))]

  make_go_functional_dotplot(
    so_obj       = so_scope,
    scope_markers = scope_mk,
    scope_name   = scope_name,
    label_col    = LABEL_COL,
    group_col    = GROUP_COL,
    groups       = GROUPS,
    group_colors = GROUP_COLORS,
    label_order  = LABEL_ORDER,
    label_colors = LABEL_COLORS,
    output_dir   = OUTPUT_DIR
  )
}
```

---

## Output File Naming

Per comparison (main loop):
```
{output_dir}/{scope}_{comp_label}_DE_full.csv
{output_dir}/{scope}_{comp_label}_volcano.pdf
{output_dir}/{scope}_{comp_label}_heatmap_overall.pdf        (ComplexHeatmap)
{output_dir}/{scope}_{comp_label}_heatmap_functional.pdf     (ComplexHeatmap; skipped if no functional_gene_sets)
{output_dir}/{scope}_{comp_label}_anatomical_dotplot.pdf
{output_dir}/{scope}_{comp_label}_functional_dotplot.pdf     (skipped if no functional_gene_sets)
{output_dir}/{scope}_{comp_label}_functional_dotplot_dir.pdf (skipped if no functional_gene_sets)
{output_dir}/{scope}_{comp_label}_pathway_barplot.pdf
```

Per scope (GO functional loop):
```
{output_dir}/GO_functional_dotplot_{scope}.pdf
{output_dir}/GO_functional_dotplot_dir_{scope}.pdf
{output_dir}/GO_sections_{scope}.csv
```

---

## Figure Sizing Reference

| Plot | Width | Height |
|---|---|---|
| Volcano | 9 | 7 |
| Anatomical dotplot (N groups) | `n_subtypes * 0.55 * N + 7` | `max(6, n_genes * 0.25 + 3)` |
| Functional dotplot (N groups, with direction) | above + 1 | same |
| GO functional dotplot | `n_subtypes * 0.55 * N + 7` | `max(6, n_genes * 0.17 + 3)` |
| Overall heatmap | 12 | 10 |
| Pathway barplot | 14 | `max(4, n_bars * 0.28 + 2)` |

Always add `limitsize = FALSE` to all `ggsave` calls. Plots taller than 50 in error without it.

---

## Project-Specific Values (Stage for Phase 4 examples/)

`examples/bonemarrow_3site_anatomical.md` must define:

- `GROUPS = c("Vertebrae", "Iliac Crest", "Femoral Head")` (the 3 anatomical sites)
- `GROUP_COL = "site"`, `LABEL_COL = "unified_label"`
- `STROMA_ORDER`: canonical ordered vector of bone marrow stromal cell type names.
  This was referenced in v1 anatomical_de_analysis.md (lines 68 and 107) but **was never
  defined** in any v1 library file. Must be curated from BoneMarrowStroma project records
  before running `make_anatomical_dotplot`. Suggested entries: Adipo-MSC, APOD+ MSC,
  Osteoblast — exact order pending confirmation.
- `LABEL_COLORS`: named color vector keyed by `STROMA_ORDER` entries (undefined in v1;
  must be created for the examples/ file)
- `functional_gene_sets`: named list of bone marrow stromal gene sets (MSC fate,
  osteogenic, adipogenic, fibro-inflammatory programs; biology-specific)
- `GROUP_COLORS`: named vector for the 3 anatomical site values (for direction strip)
- `comparisons`: 3 pairwise comparisons: GroupA_vs_GroupB, GroupA_vs_GroupC, GroupB_vs_GroupC
- Validated run: BoneMarrowStroma, 2026-03-16
- Validated script filenames: 1.5_AnatomicalComparison.R, 1.6_AnatomicalDE_CellTypeResolved.R,
  1.7_SiteDE_GOFunctionalDotplot.R
