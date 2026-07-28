---
requires_context:
  palettes:
    - subtype_colors   # optional — named vector for cell subtype labels in UMAP/violin/split-violin
    - group_colors     # optional — named vector for group/tissue colors in violin and bar plots
  metadata_columns:
    required:
      - subtype_col    # cell subtype label column (for violin, split violin, TF dot plot)
    optional:
      - group_col      # grouping column (tissue, condition, etc.) for group-stratified plots
  brief_keys:
    required:
      - output_dir
      - downstream_analyses.bulk_concordance.mode
      - downstream_analyses.bulk_concordance.bulk_csv
      - downstream_analyses.bulk_concordance.experiment_label
      - downstream_analyses.bulk_concordance.subtype_col
    optional:
      - downstream_analyses.bulk_concordance.group_col
      - downstream_analyses.bulk_concordance.exclude_subtypes
      - downstream_analyses.bulk_concordance.target_tf
references:
  - "@primitives/differential_expression.md"
  - "@primitives/visualization.md"
  - "@primitives/seurat_v5_rules.md"
---

# Module: Bulk × scRNA-seq Concordance Analysis

Cross-references a bulk RNA-seq perturbation or overexpression experiment with in-vivo
single-cell data to ask: does the perturbation signature map onto a specific cell
identity, tissue condition, or transcription factor program?

Two analytical modes address this question at different resolutions:

- **Mode 1 — Signature score concordance** (`mode: signature_score`): assigns every cell
  a scalar score (AddModuleScore up-genes minus down-genes), then visualizes score
  distributions across subtypes, tissues, and optionally a TF program (Part 3). Best for
  exploratory analysis, global pattern detection, and integration with trajectory data.

- **Mode 2 — Parallel log2FC heatmap** (`mode: parallel_lfc`): independently computes
  bulk DESeq2 LFC and scRNA-seq FindMarkers LFC for the same curated gene set, then
  places them side by side in a single pheatmap. Best for publication-ready figures
  showing specific gene-level evidence of concordance.

---

## Brief Schema

```yaml
downstream_analyses:
  bulk_concordance:
    mode: project_specific           # REQUIRED: signature_score | parallel_lfc
    bulk_csv: project_specific       # REQUIRED: path to bulk DE results CSV
                                     # Required columns: gene (or symbol), log2FoldChange, padj
    experiment_label: project_specific  # REQUIRED: short ASCII label for axis/title use
                                     # e.g. "NKX2-3 OE" — no Unicode, no em-dash
    subtype_col: project_specific    # REQUIRED: metadata column for cell subtype labels
    group_col: project_specific      # optional: metadata column for grouping (tissue, condition)
                                     # set to null to skip group-stratified plots
    exclude_subtypes: []             # optional: character vector of subtype values to exclude
                                     # e.g. low-quality clusters or off-target cell types
    target_tf: null                  # Mode 1 only: TF of interest for Part 3 TF analysis
                                     # set to null to skip Part 3
```

---

## Choose This When

```
Want per-cell scores and a global landscape of concordance?
  → Mode 1 (signature_score)
    Operates on the live Seurat object (AddModuleScore at runtime)
    Good for: trajectory integration, exploring which clusters are most concordant
    Outputs: UMAP, violin, split violin, mean±SE bar, per-cell score CSV
    Optional: Part 3 TF analysis when target_tf is set in the brief

Want a compact publication figure for a specific gene program?
  → Mode 2 (parallel_lfc)
    Operates on cached DE tables — Seurat object not needed at viz time
    Good for: main figures demonstrating gene-level concordance for reviewers
    Output: single pheatmap (PDF + PNG) with bulk LFC + scRNA-seq LFC columns
    Requires: curated biology_gene_sets (caller-provided — see Phase 4 examples/)

Running both? Run Mode 1 first (exploration), then Mode 2 (publication figure).
```

---

## Mode Comparison

| | Mode 1 — Signature score | Mode 2 — Parallel LFC |
|---|---|---|
| **Unit of analysis** | Cell (per-cell score) | Gene (per-gene LFC) |
| **Bulk DE method** | Gene list → up/down sets | Full DESeq2 re-run |
| **scRNA-seq method** | AddModuleScore (no new DE) | FindMarkers (independent DE) |
| **Seurat object at viz time?** | Yes | No (cached tables) |
| **Main output** | UMAP, violins, bar plots, score CSV | pheatmap (single figure) |
| **Best for** | Exploratory, trajectory integration | Publication-ready gene-level evidence |
| **Part 3 TF analysis** | Yes (conditional on target_tf) | No |

**Shared between modes:** bulk_csv input, experiment_label, subtype_col, group_col,
exclude_subtypes, output directory structure, BULK_LFC_CUT threshold.

---

## Critical Constraints

| ❌ Don't | ✅ Do instead | Why |
|---|---|---|
| Include all bulk genes regardless of effect size | Apply `BULK_LFC_CUT` on top of padj filtering | Small-effect bulk genes dilute module scores and weaken concordance figures |
| Use inconsistent LFC thresholds across parts | Define `BULK_LFC_CUT` once in the config block and use it everywhere | Inconsistent thresholds make concordant-gene counts incomparable across plots |
| Run `AddModuleScore` on genes absent from the object | `intersect(bulk_genes, rownames(obj))` before scoring | Silently incorrect scores if the function cannot find genes |
| Use Unicode or emoji in PDF titles | ASCII only: `-` not `—`, `>=` not `≥` | `mbcsToSbcs` / `grid.Call.graphics` rendering errors on this R/font stack |
| Combine bulk and scRNA-seq LFCs on a shared color scale (Mode 2) | Show raw log2FC from each independently | scRNA-seq LFCs are smaller in magnitude; shared z-scoring collapses real differences |
| Run on excluded low-quality subtypes | Filter by `EXCLUDE_SUBTYPES` before any scoring | Abnormal expression profiles distort module scores and marker tests |
| Fill NA in scRNA-seq LFC matrix with 0 (Mode 2) | Keep NA; set `na_col = "grey88"` in pheatmap | 0 = "no change"; grey = "not tested" — these are different |
| Use SCT assay for FindMarkers in Mode 2 | RNA assay after `NormalizeData` | SCT multi-model throws `unequal library sizes` on some Seurat v5 builds |
| Negate LFC silently | Add a `negate_invivo` parameter to scatter functions | If direction is flipped, the intent must be traceable in the code |
| Use base R | `@primitives/r_environment.md` execution rules | sp package broken in base R |

---

## R Packages Required

```r
# Both modes
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

# Mode 1 (additional)
library(rstatix)
library(ggpubr)
library(ComplexHeatmap)   # Part 3 TF LFC heatmap only
library(circlize)          # Part 3 TF LFC heatmap colorRamp2

# Mode 2 (additional)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(tibble)
```

---

## Shared Configuration Block

```r
# ── CONFIG (shared across both modes) ────────────────────────────────────────
MODE             <- "project_specific"   # REPLACE: "signature_score" or "parallel_lfc"
BULK_DE_PATH     <- "project_specific"   # REPLACE: path to bulk DE results CSV
EXPERIMENT_LABEL <- "project_specific"   # REPLACE: short ASCII label (no Unicode)
                                          # Used in axis titles and score column names
SUBTYPE_COL      <- "project_specific"   # REPLACE: metadata column for cell subtype labels
GROUP_COL        <- "project_specific"   # REPLACE: metadata column for grouping
                                          # set to NULL to skip group-stratified plots
RDS_IN           <- "project_specific"   # REPLACE: path to Seurat object
OUTPUT_DIR       <- file.path("output", "bulk_concordance")

EXCLUDE_SUBTYPES <- c()          # REPLACE: character vector of subtypes to exclude
                                  # leave as c() to keep all subtypes

# Thresholds
PADJ_CUT     <- 0.05             # in-vivo DE significance cutoff
LFC_CUT      <- 0.5              # in-vivo |log2FC| minimum
BULK_LFC_CUT <- 1.0              # bulk |log2FC| minimum (on top of padj filtering)
                                  # raise to 1.5 for noisy bulk experiments

# Color palettes — provide from brief context_overrides.palettes or define here
subtype_colors <- c("project_specific" = "#project_specific")  # REPLACE: one entry per subtype
group_colors   <- c("project_specific" = "#project_specific")  # REPLACE: one entry per group

set.seed(42)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Load and categorize bulk DE ───────────────────────────────────────────────
bulk_raw <- read.csv(BULK_DE_PATH, stringsAsFactors = FALSE)
colnames(bulk_raw)[colnames(bulk_raw) == "symbol"] <- "gene"   # normalise if needed

# Bulk DE category labels use EXPERIMENT_LABEL for generality — no hardcoded condition names
BULK_UP_LABEL   <- paste0("Up (", EXPERIMENT_LABEL, ")")
BULK_DOWN_LABEL <- paste0("Down (", EXPERIMENT_LABEL, ")")

bulk_colors <- c(
  setNames("#D73027", BULK_UP_LABEL),
  setNames("#4575B4", BULK_DOWN_LABEL),
  "NS" = "#CCCCCC"
)

bulk_df <- bulk_raw %>%
  filter(!is.na(padj), padj < PADJ_CUT) %>%
  mutate(bulk_cat = case_when(
    log2FoldChange >  BULK_LFC_CUT ~ BULK_UP_LABEL,
    log2FoldChange < -BULK_LFC_CUT ~ BULK_DOWN_LABEL,
    TRUE                           ~ "NS"
  ))

bulk_gene_cat <- setNames(bulk_df$bulk_cat, bulk_df$gene)
bulk_up       <- bulk_df$gene[bulk_df$bulk_cat == BULK_UP_LABEL]
bulk_down     <- bulk_df$gene[bulk_df$bulk_cat == BULK_DOWN_LABEL]
bulk_sig      <- c(bulk_up, bulk_down)

get_bulk_status <- function(genes) {
  status <- ifelse(genes %in% names(bulk_gene_cat), bulk_gene_cat[genes], "NS")
  factor(status, levels = c(BULK_UP_LABEL, BULK_DOWN_LABEL, "NS"))
}
```

---

## Mode 1: Signature Score Concordance

### Step 1: Load Object and Apply Exclusions

```r
obj <- readRDS(RDS_IN)

# Apply subtype exclusions if specified
if (length(EXCLUDE_SUBTYPES) > 0 && SUBTYPE_COL %in% colnames(obj@meta.data)) {
  keep_cells <- colnames(obj)[!obj@meta.data[[SUBTYPE_COL]] %in% EXCLUDE_SUBTYPES]
  obj <- obj[, keep_cells]
  message(sprintf("After exclusion: %d cells remaining", ncol(obj)))
}

# See @primitives/seurat_v5_rules.md Rule 1 — JoinLayers required before FindMarkers
obj <- JoinLayers(obj)
```

### Step 2: Per-Cell Concordance Score (AddModuleScore)

```r
bulk_up_use   <- intersect(bulk_up,   rownames(obj))
bulk_down_use <- intersect(bulk_down, rownames(obj))

message(sprintf("Genes in object — up: %d/%d, down: %d/%d",
  length(bulk_up_use), length(bulk_up),
  length(bulk_down_use), length(bulk_down)))

# Derive column names from EXPERIMENT_LABEL — make.names() handles spaces/special chars
SCORE_UP_BASE <- make.names(paste0(tolower(gsub(" ", "_", EXPERIMENT_LABEL)), "_up_score"))
SCORE_DN_BASE <- make.names(paste0(tolower(gsub(" ", "_", EXPERIMENT_LABEL)), "_down_score"))
SCORE_COL     <- make.names(paste0(tolower(gsub(" ", "_", EXPERIMENT_LABEL)), "_concordance"))

# AddModuleScore appends "1" to the name argument — retrieve via paste0(name, "1")
obj <- AddModuleScore(obj, features = list(bulk_up_use),
                      name = SCORE_UP_BASE, nbin = 24, seed = 42)
obj <- AddModuleScore(obj, features = list(bulk_down_use),
                      name = SCORE_DN_BASE, nbin = 24, seed = 42)

obj@meta.data[[SCORE_COL]] <- obj@meta.data[[paste0(SCORE_UP_BASE, "1")]] -
                               obj@meta.data[[paste0(SCORE_DN_BASE, "1")]]

# Sanity check — a well-powered signature should have SD > 0.03
message(sprintf("Score — mean: %.3f, SD: %.3f, range: [%.3f, %.3f]",
  mean(obj@meta.data[[SCORE_COL]]), sd(obj@meta.data[[SCORE_COL]]),
  min(obj@meta.data[[SCORE_COL]]),  max(obj@meta.data[[SCORE_COL]])))
```

### Step 3: UMAP Colored by Concordance Score

```r
# Detect UMAP reduction name — see @primitives/seurat_v5_rules.md Rule 5
umap_key <- names(obj@reductions)[grepl("umap", names(obj@reductions), ignore.case = TRUE)][1]
umap_df  <- as.data.frame(Embeddings(obj, umap_key))
colnames(umap_df) <- c("umap_1", "umap_2")
umap_df$score <- obj@meta.data[[SCORE_COL]]

score_lim <- quantile(umap_df$score, c(0.02, 0.98))   # clip outliers for colour scale

p_umap <- ggplot(umap_df, aes(x = umap_1, y = umap_2, colour = score)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_colour_gradientn(
    colours = c("#4575B4", "#FFFFBF", "#D73027"),
    limits  = score_lim,
    oob     = scales::squish,
    name    = paste0(EXPERIMENT_LABEL, "\nconcordance\nscore")
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank())

ggsave(file.path(OUTPUT_DIR, "concordance_umap.pdf"),
       p_umap, width = 6, height = 5, device = "pdf", useDingbats = FALSE)
```

### Step 4: Violin by Subtype and Group

```r
meta_score <- obj@meta.data[, c(SUBTYPE_COL, GROUP_COL), drop = FALSE]
meta_score$score <- obj@meta.data[[SCORE_COL]]

make_concordance_violin <- function(df, x_col, fill_col, colors, output_path,
                                    plot_width = 7, plot_height = 5.5) {
  df_work <- df
  colnames(df_work)[colnames(df_work) == x_col] <- "group_var"

  stat_res <- tryCatch({
    df_work %>%
      rstatix::wilcox_test(score ~ group_var, p.adjust.method = "BH") %>%
      rstatix::add_significance() %>%
      filter(p.adj < 0.05) %>%
      rstatix::add_xy_position(x = "group_var")
  }, error = function(e) NULL)

  p <- ggplot(df_work, aes(x = group_var, y = score, fill = !!sym(fill_col))) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85, colour = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 colour = "grey30", linewidth = 0.5) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_x_discrete(guide = guide_axis(angle = 35)) +
    labs(y = paste0(EXPERIMENT_LABEL, " concordance score"), x = NULL) +
    theme_classic(base_size = 14)

  if (!is.null(stat_res) && nrow(stat_res) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(stat_res, label = "p.adj.signif",
                                         tip.length = 0.01, step.increase = 0.06, size = 3.5)
  }

  ggsave(output_path, p, width = plot_width, height = plot_height,
         device = "pdf", useDingbats = FALSE)
}

make_concordance_violin(
  meta_score, SUBTYPE_COL, SUBTYPE_COL, subtype_colors,
  file.path(OUTPUT_DIR, "concordance_violin_subtype.pdf"),
  plot_width = 7, plot_height = 5.5
)

if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  make_concordance_violin(
    meta_score, GROUP_COL, GROUP_COL, group_colors,
    file.path(OUTPUT_DIR, "concordance_violin_group.pdf"),
    plot_width = 6.5, plot_height = 5.5
  )
}
```

### Step 5: Split Violin (Subtype × Group)

```r
if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  meta_cross <- meta_score

  # Exclude subtype × group combinations with < 20 cells (avoids empty violin shapes)
  cell_counts <- meta_cross %>% count(!!sym(SUBTYPE_COL), !!sym(GROUP_COL))
  keep_pairs  <- cell_counts %>% filter(n >= 20)
  meta_cross  <- meta_cross %>% semi_join(keep_pairs, by = c(SUBTYPE_COL, GROUP_COL))

  n_subtypes <- length(unique(meta_cross[[SUBTYPE_COL]]))
  split_w    <- max(9, n_subtypes * 1.6 + 2)

  p_split <- ggplot(meta_cross, aes(x = !!sym(SUBTYPE_COL), y = score,
                                     fill = !!sym(GROUP_COL))) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8,
                position = position_dodge(0.85), colour = NA) +
    geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white",
                 colour = "grey40", linewidth = 0.4,
                 position = position_dodge(0.85)) +
    scale_fill_manual(values = group_colors, name = "Group") +
    scale_x_discrete(guide = guide_axis(angle = 35)) +
    labs(y = paste0(EXPERIMENT_LABEL, " concordance score"), x = NULL) +
    theme_classic(base_size = 13)

  ggsave(file.path(OUTPUT_DIR, "concordance_split_violin.pdf"),
         p_split, width = split_w, height = 4.75, device = "pdf", useDingbats = FALSE)
}
```

### Step 6: Mean ± SE Bar Plot by Group

```r
if (!is.null(GROUP_COL) && GROUP_COL %in% colnames(obj@meta.data)) {
  summ_df <- meta_score %>%
    group_by(!!sym(GROUP_COL)) %>%
    summarise(mean_score = mean(score), se = sd(score) / sqrt(n()),
              n = n(), .groups = "drop") %>%
    arrange(desc(mean_score))
  summ_df[[GROUP_COL]] <- factor(summ_df[[GROUP_COL]], levels = summ_df[[GROUP_COL]])

  min_y <- min(summ_df$mean_score - summ_df$se) * 1.3

  p_bar <- ggplot(summ_df, aes(x = !!sym(GROUP_COL), y = mean_score,
                                fill = !!sym(GROUP_COL))) +
    geom_col(width = 0.65, alpha = 0.9) +
    geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                  width = 0.2, colour = "grey30", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_text(aes(label = format(n, big.mark = ",")),
              y = min_y, size = 3, colour = "grey40") +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_x_discrete(guide = guide_axis(angle = 35)) +
    labs(y = paste0(EXPERIMENT_LABEL, " concordance\n(mean +/- SE)"), x = NULL) +
    theme_classic(base_size = 14)

  n_groups <- nrow(summ_df)
  ggsave(file.path(OUTPUT_DIR, "concordance_barplot.pdf"),
         p_bar, width = max(5, n_groups * 0.9 + 2), height = 5,
         device = "pdf", useDingbats = FALSE)
}
```

### Step 7: Save Per-Cell Scores to CSV

```r
export_cols <- intersect(c(SUBTYPE_COL, GROUP_COL, SCORE_COL), colnames(obj@meta.data))

write.csv(
  cbind(data.frame(cell_id = colnames(obj)), obj@meta.data[, export_cols, drop = FALSE]),
  file.path(OUTPUT_DIR, "concordance_scores.csv"),
  row.names = FALSE
)
```

---

### Part 1 (Mode 1): Bulk Status Annotation Overlay

Annotates existing DE outputs (from @primitives/differential_expression.md volcano/dot plots)
with a color strip indicating which genes are up, down, or NS in the bulk experiment.

```r
# Build a narrow bulk-status tile strip for attaching to dot plots
# gene_order: character vector of gene names in y-axis display order (top to bottom)
# top_margin_pt: pixel height of the main plot's top axis labels (~35pt for 45-degree labels)
make_bulk_strip <- function(gene_order, top_margin_pt = 0) {
  strip_df <- data.frame(
    gene     = factor(gene_order, levels = rev(gene_order)),
    bulk_cat = get_bulk_status(gene_order),
    x        = 1L
  )
  ggplot(strip_df, aes(x = x, y = gene, fill = bulk_cat)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = bulk_colors, name = "Bulk status") +
    scale_x_continuous(breaks = 1, labels = "Bulk") +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x    = element_text(size = 8, angle = 45, hjust = 0, vjust = 0),
      axis.text.y    = element_blank(),
      axis.title     = element_blank(),
      panel.grid     = element_blank(),
      legend.text    = element_text(size = 8),
      legend.title   = element_text(size = 9),
      plot.margin    = margin(t = top_margin_pt, r = 2, b = 4, l = 2)
    )
}

# Usage: append to an existing dot plot where gene_order is the y-axis gene vector
# p_combined <- p_main + make_bulk_strip(gene_order, top_margin_pt = 35) +
#   plot_layout(widths = c(1, 0.08), guides = "collect")
```

Four-way colored volcano (in-vivo × bulk):

```r
four_way_colors <- c(
  "Concordant up"   = "#B2182B",
  "Concordant down" = "#2166AC",
  "In-vivo only"    = "#4DAC26",
  "NS"              = "#CCCCCC"
)

# de_df: data frame from @primitives/differential_expression.md with gene, avg_log2FC, p_val_adj
annotate_four_way <- function(de_df) {
  de_df %>%
    mutate(
      invivo_sig   = p_val_adj < PADJ_CUT & abs(avg_log2FC) > LFC_CUT,
      bulk_up_gene = gene %in% bulk_up,
      bulk_dn_gene = gene %in% bulk_down,
      four_way = factor(case_when(
        invivo_sig & avg_log2FC >  0 & bulk_up_gene ~ "Concordant up",
        invivo_sig & avg_log2FC <= 0 & bulk_dn_gene ~ "Concordant down",
        invivo_sig                                   ~ "In-vivo only",
        TRUE                                         ~ "NS"
      ), levels = names(four_way_colors))
    )
}

# Bulk LFC scatter: in-vivo LFC (y) vs bulk LFC (x)
# negate_invivo = TRUE if ident1/ident2 ordering conflicts with biological direction
make_bulk_scatter <- function(de_full, negate_invivo = FALSE, output_path) {
  scatter_df <- de_full %>%
    left_join(bulk_df %>% select(gene, bulk_lfc = log2FoldChange), by = "gene") %>%
    filter(!is.na(bulk_lfc)) %>%
    mutate(
      lfc_plot   = if (negate_invivo) -avg_log2FC else avg_log2FC,
      concordant = (lfc_plot > 0 & bulk_lfc > 0) | (lfc_plot < 0 & bulk_lfc < 0),
      sig_both   = p_val_adj < PADJ_CUT & abs(lfc_plot) > LFC_CUT & concordant
    )

  p <- ggplot(scatter_df, aes(x = bulk_lfc, y = lfc_plot)) +
    geom_point(aes(colour = sig_both), alpha = 0.6, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    scale_colour_manual(values = c("TRUE" = "#F9A825", "FALSE" = "#CCCCCC"),
                        labels = c("TRUE" = "Concordant", "FALSE" = "Other"), name = NULL) +
    labs(x = paste0("Bulk log2FC (", EXPERIMENT_LABEL, ")"), y = "In-vivo log2FC") +
    theme_classic(base_size = 14)

  ggsave(output_path, p, width = 6.5, height = 5.5, device = "pdf", useDingbats = FALSE)
}
```

---

### Part 3 (Mode 1): TF-Focused Analysis

**Conditional on `TARGET_TF` being non-null.**

This section requires a caller-provided TF gene list appropriate to your cell type biology.
For EC-biology TF lists (Lambert et al. 2018, Cell), see examples/. For immune, stromal,
or other contexts, provide a curated list for your biology in the project CLAUDE.md.

```r
TARGET_TF <- NULL   # REPLACE from brief: downstream_analyses.bulk_concordance.target_tf
                     # Set to NULL to skip Part 3 entirely

if (!is.null(TARGET_TF)) {

# ── TF list — caller-provided; biology-specific ───────────────────────────────
# EC biology: ETS family, KLF/SP, SOX, FOX, nuclear receptors, AP-1, STAT, NF-kB, etc.
# Immune: IRF, STAT, NF-kB, ETS
# Fibroblast: AP-1, KLF, ZEB/SNAI
# See examples/ for validated biology-specific TF lists
human_tfs <- c("project_specific")   # REPLACE: curated TF list for your cell type

# ── Filter to expressed TFs (>= 1% of cells in at least one subtype) ─────────
expr_mat   <- GetAssayData(obj, layer = "data")
tfs_in_obj <- intersect(human_tfs, rownames(expr_mat))

tfs_use <- tfs_in_obj[sapply(tfs_in_obj, function(g) {
  any(tapply(as.numeric(expr_mat[g, ]) > 0, obj@meta.data[[SUBTYPE_COL]], mean) >= 0.01)
})]

bulk_tf_up   <- tfs_use[tfs_use %in% bulk_up]
bulk_tf_down <- tfs_use[tfs_use %in% bulk_down]

# ── TF dot plot: ordered by bulk status (up first, then down, then NS) ────────
gene_order <- c(
  sort(bulk_tf_up[bulk_tf_up   %in% tfs_use]),
  sort(bulk_tf_down[bulk_tf_down %in% tfs_use]),
  sort(tfs_use[!tfs_use %in% c(bulk_tf_up, bulk_tf_down)])
)
# Attach make_bulk_strip(gene_order) from Part 1 on the right
n_tfs      <- length(tfs_use)
n_subtypes <- length(unique(obj@meta.data[[SUBTYPE_COL]]))
dot_h      <- n_tfs * 0.22 + 2.5
dot_w      <- n_subtypes * 0.75 + 4

# ── One-vs-all TF markers per subtype ────────────────────────────────────────
# See @primitives/seurat_v5_rules.md Rule 1 — JoinLayers before FindMarkers
Idents(obj) <- SUBTYPE_COL

marker_list <- lapply(unique(obj@meta.data[[SUBTYPE_COL]]), function(ct) {
  tryCatch({
    m <- FindMarkers(obj, ident.1 = ct,
                     features        = tfs_use,
                     min.pct         = 0.05,
                     logfc.threshold = 0.1,
                     test.use        = "wilcox",
                     verbose         = FALSE)
    m$gene    <- rownames(m)
    m$subtype <- ct
    m
  }, error = function(e) NULL)
})
markers_all <- bind_rows(marker_list)

markers_sig <- markers_all %>%
  filter(p_val_adj < PADJ_CUT, avg_log2FC > 0) %>%
  mutate(bulk_cat = factor(
    ifelse(gene %in% bulk_tf_up,   BULK_UP_LABEL,
      ifelse(gene %in% bulk_tf_down, BULK_DOWN_LABEL, "NS")),
    levels = c(BULK_UP_LABEL, BULK_DOWN_LABEL, "NS")
  ))

write.csv(markers_sig,
          file.path(OUTPUT_DIR, paste0("tf_markers_", SUBTYPE_COL, ".csv")),
          row.names = FALSE)

# ── TF LFC heatmap: in-vivo comparisons + bulk column (ComplexHeatmap) ───────
# Rows = TFs significant in >= 1 in-vivo comparison; columns = in-vivo LFCs + bulk LFC
# Assemble mat from your cached DE CSVs (produced by @primitives/differential_expression.md):
#   invivo_de_files <- list(
#     "Comparison A" = file.path(OUTPUT_DIR, "comp_a_de.csv"),
#     ...
#   )
#   mat[, invivo_cols] — in-vivo LFC matrix (rows = TFs)
#   mat[, "Bulk LFC"]  — bulk LFC column (same rows)
#
# Two Heatmap() objects concatenated with + so bulk column has its own colour scale:
#   ht <- Heatmap(
#     mat[, invivo_cols], name = "In-vivo\nlog2FC",
#     col = colorRamp2(c(-lfc_lim, 0, lfc_lim), c("#4575B4", "white", "#D73027")),
#     column_split = scope_factor, ...
#   ) + Heatmap(
#     mat[, "Bulk LFC", drop = FALSE], name = "Bulk\nlog2FC",
#     col = colorRamp2(c(-bulk_lim, 0, bulk_lim), c("#4575B4", "white", "#D73027")),
#     width = unit(0.8, "cm"), ...
#   )
#
# ComplexHeatmap exception — pdf()/dev.off() required (see CONVENTIONS.md §4 exception #1)
#   pdf(file.path(OUTPUT_DIR, "tf_heatmap.pdf"),
#       width = 9, height = max(5, n_tfs * 0.22 + 3))
#   draw(ht)
#   dev.off()

}  # end if (!is.null(TARGET_TF))
```

---

## Mode 2: Parallel log2FC Heatmap

### Architecture: Two-Stage with Caching

```
Stage 1 — Compute (runs once, then cached):
  cache_bulk_DE.Rds   <- DESeq2 res_df + sample_df
  cache_sc_DE.Rds     <- named list of FindMarkers data.frames

Stage 2 — Visualise (always runs from cache):
  gene selection from caller-provided biology_gene_sets
  LFC matrix assembly (bulk col + scRNA-seq cols)
  pheatmap -> PDF + PNG (pheatmap filename= writes directly; no pdf()/dev.off() needed)
```

To force recomputation: `rm cache_bulk_DE.Rds` or `rm cache_sc_DE.Rds` (caches are independent).

---

### Part 1 (Mode 2): Bulk DE (DESeq2)

```r
library(DESeq2)

BULK_CACHE <- file.path(OUTPUT_DIR, "cache_bulk_DE.Rds")

if (file.exists(BULK_CACHE)) {
  message("Loading cached bulk DE...")
  b         <- readRDS(BULK_CACHE)
  res_df    <- b$res_df
  sample_df <- b$sample_df
} else {
  # ── Provide from bulk experiment design ───────────────────────────────────────
  counts_mat <- "project_specific"    # REPLACE: load bulk count matrix (genes x samples)
  sample_df  <- "project_specific"    # REPLACE: load sample metadata with condition + batch cols

  condition_samples <- c("project_specific")   # REPLACE: perturbation sample names
  control_samples   <- c("project_specific")   # REPLACE: matched control sample names
  # Sample selection: include only batches that have BOTH condition AND control samples

  keep_samples <- c(condition_samples, control_samples)
  counts_mat   <- counts_mat[, keep_samples]
  sample_df    <- sample_df[sample_df$sample_name %in% keep_samples, ]

  sample_df$condition <- factor(
    ifelse(sample_df$sample_name %in% condition_samples, "Perturbation", "Control"),
    levels = c("Control", "Perturbation")
  )
  sample_df$batch <- as.factor(sample_df$batch_col)   # REPLACE batch_col with actual column name

  # Pre-filter: keep genes with >= 10 counts in >= 3 samples
  counts_mat <- counts_mat[rowSums(counts_mat >= 10) >= 3, ]

  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = sample_df,
    design    = ~ batch + condition
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "Perturbation", "Control"), alpha = 0.05)

  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj, desc(log2FoldChange))

  write.csv(res_df,
            file.path(OUTPUT_DIR,
                      paste0("DESeq2_", make.names(EXPERIMENT_LABEL), "_vs_Control.csv")),
            row.names = FALSE)
  saveRDS(list(res_df = res_df, sample_df = sample_df), BULK_CACHE)
  message(sprintf("Bulk DE complete: %d genes, cached.", nrow(res_df)))
}
```

---

### Part 2 (Mode 2): scRNA-seq DE (FindMarkers)

**Always use RNA assay — see Critical Constraints.**

```r
SC_CACHE <- file.path(OUTPUT_DIR, "cache_sc_DE.Rds")

if (file.exists(SC_CACHE)) {
  sc_de <- readRDS(SC_CACHE)
} else {
  # Assign sc_group column: target subtypes vs a reference population
  # Positive LFC = enriched in target vs reference
  obj$sc_group <- NA_character_

  # Reference group: cells NOT in the target cell population (use broad label or Idents)
  ref_cells <- colnames(obj)[!grepl("project_specific", as.character(Idents(obj)))]
  obj$sc_group[ref_cells] <- "Reference"

  # Target subtypes: assign from the high-resolution metadata column
  # Repeat for each target subtype:
  # obj$sc_group[obj@meta.data[[SUBTYPE_COL]] == "SubtypeA"] <- "SubtypeA"
  # obj$sc_group[obj@meta.data[[SUBTYPE_COL]] == "SubtypeB"] <- "SubtypeB"

  obj_sub <- subset(obj, cells = colnames(obj)[!is.na(obj$sc_group)])
  DefaultAssay(obj_sub) <- "RNA"
  obj_sub <- JoinLayers(obj_sub)   # see @primitives/seurat_v5_rules.md Rule 1
  obj_sub <- NormalizeData(obj_sub, verbose = FALSE)
  Idents(obj_sub) <- "sc_group"

  target_groups <- setdiff(unique(obj_sub$sc_group), "Reference")
  cat("Group sizes:\n"); print(sort(table(obj_sub$sc_group), decreasing = TRUE))

  sc_de <- lapply(target_groups, function(grp) {
    cat(sprintf("  FindMarkers: %s vs Reference...\n", grp))
    df <- FindMarkers(obj_sub,
                      ident.1             = grp,
                      ident.2             = "Reference",
                      assay               = "RNA",
                      test.use            = "wilcox",
                      min.pct             = 0.05,
                      logfc.threshold     = 0,
                      max.cells.per.ident = 8000,   # downsample large reference for speed
                      verbose             = FALSE)
    df %>% tibble::rownames_to_column("gene") %>% mutate(group = grp)
  })
  names(sc_de) <- target_groups
  saveRDS(sc_de, SC_CACHE)
}
```

---

### Part 3 (Mode 2): Gene Selection from Biology Sets

```r
# biology_gene_sets: caller-provided named list of gene vectors by functional category
# Project-specific — define in project CLAUDE.md or brief; see examples/ for biology-specific sets
# Structure: list("Category Name" = c("GENE1", "GENE2", ...), ...)
biology_gene_sets <- list("project_specific" = c("project_specific"))  # REPLACE

MIN_BULK_LFC <- 0.5     # minimum bulk effect size (raise to 1.0 for publication-tight figures)
MAX_PER_CAT  <- 6       # max genes per category (4-6 for main; up to 8 for supplementary)

sig_bulk <- res_df %>% filter(padj < 0.05, log2FoldChange > MIN_BULK_LFC)

selected_list <- lapply(names(biology_gene_sets), function(cat) {
  hits <- sig_bulk %>%
    filter(gene %in% biology_gene_sets[[cat]]) %>%
    arrange(desc(log2FoldChange)) %>%
    slice_head(n = MAX_PER_CAT)
  if (nrow(hits) == 0) return(NULL)
  hits %>% mutate(category = cat)
})

sel_df <- bind_rows(selected_list) %>%
  mutate(category = factor(category, levels = names(biology_gene_sets))) %>%
  arrange(category, desc(log2FoldChange)) %>%
  filter(!duplicated(gene))   # gene in multiple categories: keep first (category order controls)

heatmap_genes <- sel_df$gene
```

---

### Part 4 (Mode 2): Heatmap Assembly

```r
# Column label for bulk — ASCII only (no Unicode)
BULK_COL_LABEL <- paste0(EXPERIMENT_LABEL, "\nvs Control")

lfc_mat <- data.frame(
  gene             = heatmap_genes,
  !!BULK_COL_LABEL := res_df$log2FoldChange[match(heatmap_genes, res_df$gene)],
  check.names = FALSE
)

# scRNA-seq LFC columns — one per target subtype; column headers describe the comparison
sc_col_labels <- setNames(paste0(names(sc_de), "\nvs Reference"), names(sc_de))

for (grp in names(sc_de)) {
  cn  <- sc_col_labels[grp]
  sub <- sc_de[[grp]] %>%
    filter(gene %in% heatmap_genes) %>%
    select(gene, avg_log2FC) %>%
    rename(!!cn := avg_log2FC)
  lfc_mat <- lfc_mat %>% left_join(sub, by = "gene")
}

# Keep NAs — never replace with 0 (see Critical Constraints)
lfc_mat <- lfc_mat %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()
lfc_mat <- lfc_mat[heatmap_genes, ]   # enforce category-sorted row order

# Row annotation: functional category
row_ann <- data.frame(
  Category  = as.character(sel_df$category[match(rownames(lfc_mat), sel_df$gene)]),
  row.names = rownames(lfc_mat))

# Column annotation: data source
col_ann <- data.frame(
  Data      = c("Bulk RNAseq", rep("scRNA-seq", ncol(lfc_mat) - 1)),
  row.names = colnames(lfc_mat),
  check.names = FALSE)

# Category colors — set based on biology_gene_sets category names
cat_colors <- setNames(
  RColorBrewer::brewer.pal(max(3, length(biology_gene_sets)), "Set1")[
    seq_len(length(biology_gene_sets))],
  names(biology_gene_sets)
)
data_colors <- c("Bulk RNAseq" = "#2c2c54", "scRNA-seq" = "#218c74")
ann_colors  <- list(Category = cat_colors, Data = data_colors)

LFC_MAX     <- 3
lfc_clamped <- pmax(pmin(lfc_mat, LFC_MAX), -LFC_MAX)

cat_rle  <- rle(row_ann$Category)
gap_rows <- cumsum(cat_rle$lengths)[-length(cat_rle$lengths)]

# pheatmap writes directly via filename= parameter — no pdf()/dev.off() needed
make_lfc_heatmap <- function(filename, w, h, res = NULL) {
  args <- list(
    mat               = lfc_clamped,
    color             = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    breaks            = seq(-LFC_MAX, LFC_MAX, length.out = 101),
    annotation_row    = row_ann,
    annotation_col    = col_ann,
    annotation_colors = ann_colors,
    cluster_cols      = FALSE,      # column order is meaningful; never cluster
    cluster_rows      = FALSE,      # category grouping is explicit; never cluster
    gaps_row          = gap_rows,
    na_col            = "grey88",   # not tested = grey; 0 = white (different meanings)
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 12,
    fontsize_col      = 13,
    fontsize          = 12,
    border_color      = "grey85",
    cellwidth         = 55,
    cellheight        = 18,
    legend_breaks     = c(-3, -2, -1, 0, 1, 2, 3),
    legend_labels     = c("-3", "-2", "-1", "0", "1", "2", "3"),
    main              = paste0(EXPERIMENT_LABEL,
                               " Program - Concordance\n(log2 fold change)"),
    filename          = filename,
    width             = w,
    height            = h
  )
  if (!is.null(res)) args$res <- res
  do.call(pheatmap::pheatmap, args)
}

# Width: cellwidth=55 per column + ~3.5 inches for row labels and annotations
# Height: cellheight=18 per row + ~2.5 inches for column labels; for >40 genes drop fontsize_row to 10
hm_w <- ncol(lfc_clamped) * 1.6 + 3.5
hm_h <- nrow(lfc_clamped) * 0.6 + 3.0

make_lfc_heatmap(file.path(OUTPUT_DIR, "LFC_Concordance_Heatmap.pdf"), w = hm_w, h = hm_h)
make_lfc_heatmap(file.path(OUTPUT_DIR, "LFC_Concordance_Heatmap.png"), w = hm_w, h = hm_h,
                 res = 200)
```

---

## Output File Summary

### Mode 1 (signature_score) — `output/bulk_concordance/`

| File | Description |
|---|---|
| `concordance_umap.pdf` | UMAP colored by per-cell concordance score (clipped 2–98th percentile) |
| `concordance_violin_subtype.pdf` | Violin + boxplot by subtype; Wilcoxon BH p-values |
| `concordance_violin_group.pdf` | Violin + boxplot by group (if GROUP_COL set) |
| `concordance_split_violin.pdf` | Split violin: subtype × group (if both set; pairs <20 cells excluded) |
| `concordance_barplot.pdf` | Mean ± SE bar plot by group, sorted descending |
| `concordance_scores.csv` | Per-cell scores; load for downstream use without reloading Seurat object |
| `tf_markers_{subtype_col}.csv` | One-vs-all TF markers per subtype (Part 3 only; if TARGET_TF set) |
| `tf_heatmap.pdf` | TF LFC heatmap: in-vivo comparisons + bulk column (Part 3 only; ComplexHeatmap) |

### Mode 2 (parallel_lfc) — `output/bulk_concordance/`

| File | Description |
|---|---|
| `cache_bulk_DE.Rds` | DESeq2 results + sample_df (do not delete; expensive to recompute) |
| `cache_sc_DE.Rds` | Named list of FindMarkers data.frames (delete to recompute sc DE) |
| `DESeq2_{label}_vs_Control.csv` | Full DESeq2 results table |
| `LFC_Concordance_Heatmap.pdf` | Main figure: curated gene rows, bulk + scRNA-seq LFC columns |
| `LFC_Concordance_Heatmap.png` | Same figure at 200 dpi |

---

## Key Pitfalls

### Mode 1
- **Score column naming:** `AddModuleScore` appends "1" to the `name` argument. The module
  derives score column names from `EXPERIMENT_LABEL` via `make.names()`. If `EXPERIMENT_LABEL`
  contains spaces, `make.names()` converts them to dots — verify that the retrieved column
  `paste0(SCORE_UP_BASE, "1")` exists in `obj@meta.data` before proceeding.
- **Split violin pairs < 20 cells:** excluded to avoid empty violin shapes. Check the
  `cell_counts` table if expected groups are missing from the plot.

### Mode 2
- **Grey vs white in heatmap:** `na_col = "grey88"` means not tested; white = 0 = no change.
  Do not replace NA with 0.
- **Negative scRNA-seq LFC for a bulk-upregulated gene:** may reflect high expression in the
  reference group (not absent in target). These rows are data — keep them and note in figure legend.
- **Row order:** `heatmap_genes <- sel_df$gene` after `arrange(category, ...)` — `rle()` gaps
  only produce correct category boundaries if rows are already sorted by category.
- **Garbled PDF title:** `EXPERIMENT_LABEL` must be ASCII-only — no em-dash, Unicode
  subscripts, or `≥`.

