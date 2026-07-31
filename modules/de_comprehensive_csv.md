---
requires_context:
  palettes: []
  metadata_columns:
    required:
      - group_col    # comparison group column (tissue, condition, treatment, etc.)
    optional:
      - label_col    # cell subtype label column (for top_subtype annotation within top group)
  brief_keys:
    required:
      - output_dir
      - downstream_analyses.de_comprehensive_csv.group_col
    optional:
      - downstream_analyses.de_comprehensive_csv.label_col
      - downstream_analyses.de_comprehensive_csv.group_order
      - downstream_analyses.de_comprehensive_csv.subtype_labels
references:
  - "@primitives/differential_expression.md"
  - "@primitives/seurat_v5_rules.md"
---

# Module: Comprehensive Multi-Group DE Statistics Table

**is_confound() vs is_ambient() — why this module uses is_confound():**
This module uses `is_confound()` (broader filter than `is_ambient()`) because comprehensive
multi-group tables benefit from excluding sex-linked genes (XIST, Y-chromosome markers),
HLA class II (HLA-DR/DP/DQ), histones (HIST*), and unannotated identifiers (ENSG, ENSMUSG,
AC/AL/AP/LINC) in addition to the ambient RNA and mitochondrial/ribosomal genes covered by
`is_ambient()`. Both functions are defined in `@primitives/differential_expression.md`.
Use `is_ambient()` for standard volcano/DE labeling; use `is_confound()` here.

Generates a comprehensive CSV covering all significantly variable genes across multiple groups
simultaneously (tissue types, conditions, treatments, etc.), annotated with:

- top-group enrichment (group with highest mean expression)
- log2FC of top group vs all others (cell-count-weighted)
- top cell subtype within the top group
- Kruskal-Wallis H statistic, p-value, BH-adjusted p-value, epsilon-squared effect size
- mean expression per group (natural scale via expm1)

This is the natural complement to pairwise DE from `@primitives/differential_expression.md`:
pairwise DE gives comparison-specific results; this gives a single ranked table across all
groups simultaneously. Designed as a reference table the user can sort, filter, and share.

---

## When to Use

- You have ≥ 3 groups and want a single ranked gene list with group annotations
- User asks for an "Excel sheet", "comprehensive table", or "annotated gene list"
- As a supplement to pairwise DE: covers all genes, not just top hits from specific comparisons
- You need per-group mean expression in natural scale for sharing with collaborators

---

## Brief Schema

```yaml
downstream_analyses:
  de_comprehensive_csv:
    enabled: true
    group_col: project_specific        # REQUIRED: metadata column holding group labels
    label_col: project_specific        # optional: cell subtype column for top_subtype annotation
                                       # if null, top_subtype column will contain NA
    group_order: project_specific      # optional: character vector setting group display order
                                       # if null, uses sort(unique(group values))
    subtype_labels: project_specific   # optional: character vector of subtype label values
                                       # required when label_col is set
```

---

## Statistical Approach: Kruskal-Wallis

Non-parametric one-way test across all groups at the cell level.

**Why KW instead of ANOVA:** scRNA-seq expression is zero-inflated and non-normal.
KW tests rank distributions rather than means, which is robust to this.

**Effect size:** epsilon-squared = `(H - k + 1) / (n - k)` where H = KW statistic,
k = number of groups, n = total cells. Ranges 0–1.

**Performance:** A naive `lapply(genes, kruskal.test)` loop on 12k+ genes with 37k+
cells takes 30+ minutes. Always use the vectorized implementation in Step 2.

---

## Critical Constraints

| ❌ Don't | ✅ Do instead | Why |
|---|---|---|
| `lapply(genes, kruskal.test)` | Vectorized `rowRanks` + matrix math | Per-gene loop takes 30+ min on 12k genes × 37k cells |
| Skip the KW cache | Always save and load `kw_results_cache.Rds` | Recomputing from scratch is wasteful when only filters or output change |
| Report mean in log-normalized space | `expm1(mean_mat)` for CSV columns | Natural scale is more interpretable for users |
| Use unweighted mean for "other" groups | Weight by cell count per group | Unequal group sizes bias logFC toward small groups |
| Append per-group logFC columns for every group | Emit only `logFC_vs_others_top_group` | Compact table; additional logFC columns are rarely consulted |
| Use `is_ambient()` | Use `is_confound()` from @primitives/differential_expression.md | Comprehensive tables need the broader confound filter — see module header |
| Use base R | `@primitives/r_environment.md` execution rules | sp package broken in base R |

---

## R Packages Required

```r
library(Seurat)
library(dplyr)
library(matrixStats)
```

---

## Configuration Block

```r
# ── CONFIG ───────────────────────────────────────────────────────────────────
GROUP_COL    <- "project_specific"   # REPLACE: metadata column holding group labels
                                      # (brief: downstream_analyses.de_comprehensive_csv.group_col)

LABEL_COL    <- "project_specific"   # REPLACE: cell subtype label column for top_subtype annotation
                                      # set to NULL to omit top_subtype from output
                                      # (brief: downstream_analyses.de_comprehensive_csv.label_col)

GROUP_ORDER  <- NULL   # REPLACE: optional ordered character vector of group levels
                        # e.g. c("Control", "Treatment_A", "Treatment_B")
                        # if NULL, uses sort(unique(meta[[GROUP_COL]]))

SUBTYPE_LABELS <- NULL  # REPLACE: character vector of subtype label values (levels of LABEL_COL)
                         # required if LABEL_COL is set; ignored if LABEL_COL is NULL

RDS_IN       <- "project_specific"   # REPLACE: path to Seurat subset object
OUTPUT_DIR   <- file.path("output", "de_comprehensive")

PADJ_CUT     <- 0.05   # BH-adjusted KW p-value cutoff for significance
EPS          <- 1e-4   # pseudocount for log2FC computation (avoids log2(0))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
```

---

## Step 1: Load Object and Prepare Data

```r
obj <- readRDS(RDS_IN)
# See @primitives/seurat_v5_rules.md Rule 1 — JoinLayers before GetAssayData
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"

meta     <- obj@meta.data
data_mat <- GetAssayData(obj, layer = "data")   # log-normalized expression matrix

group_order <- if (!is.null(GROUP_ORDER)) GROUP_ORDER else sort(unique(meta[[GROUP_COL]]))
k           <- length(group_order)

# Retain only genes expressed in at least 1 cell
genes_keep <- rownames(data_mat)[rowSums(data_mat > 0) > 0]
message(sprintf("Genes to test: %d", length(genes_keep)))
```

---

## Step 2: Vectorized Kruskal-Wallis (Cached)

```r
KW_CACHE <- file.path(OUTPUT_DIR, "kw_results_cache.Rds")

if (file.exists(KW_CACHE)) {
  message("Loading cached KW results...")
  kw_df <- readRDS(KW_CACHE)
} else {
  message("Computing vectorized KW (5-30 min depending on object size)...")

  # Materialise dense matrix once — avoids repeated sparse row extraction
  expr_dense <- as.matrix(data_mat[genes_keep, ])   # genes × cells

  # Rank each gene across all cells (C-level, handles ties)
  message("Computing ranks...")
  rank_mat <- rowRanks(expr_dense, ties.method = "average")   # genes × cells

  group_vec <- meta[[GROUP_COL]]
  group_int <- as.integer(factor(group_vec, levels = group_order))
  n  <- ncol(expr_dense)
  nj <- tabulate(group_int, nbins = k)

  # Group rank sums: genes × k matrix
  Rj_mat <- do.call(cbind, lapply(seq_len(k), function(j)
    rowSums(rank_mat[, group_int == j, drop = FALSE])
  ))

  # Kruskal-Wallis H statistic and p-values
  H_vec <- 12 / (n * (n + 1)) * rowSums(sweep(Rj_mat^2, 2, nj, "/")) - 3 * (n + 1)
  p_vec <- pchisq(H_vec, df = k - 1, lower.tail = FALSE)

  kw_df <- data.frame(gene = genes_keep, H = H_vec, p_val = p_vec)
  kw_df$p_adj  <- p.adjust(kw_df$p_val, method = "BH")
  kw_df$eps_sq <- pmax((kw_df$H - k + 1) / (n - k), 0)
  kw_df        <- kw_df %>% arrange(desc(H))

  saveRDS(kw_df, KW_CACHE)
  message(sprintf("KW complete: %d genes tested, cached to %s", nrow(kw_df), KW_CACHE))
}
# Delete the cache only if groups change or you need a fresh run on new data.
```

---

## Step 3: Filter Significant Genes (is_confound from primitives)

```r
# is_confound() is defined in @primitives/differential_expression.md
# It applies a superset of is_ambient() filters — see module header for the distinction.
# The raw KW cache (kw_results_cache.Rds) remains unfiltered so you can reuse it.

sig_genes <- kw_df %>%
  filter(p_adj < PADJ_CUT, !is_confound(gene)) %>%
  arrange(p_adj) %>%
  pull(gene)

message(sprintf("Significant genes after is_confound() filter: %d", length(sig_genes)))
```

---

## Step 4: Mean Expression per Group

```r
mean_mat <- sapply(group_order, function(g) {
  cells <- rownames(meta)[meta[[GROUP_COL]] == g]
  rowMeans(as.matrix(data_mat[sig_genes, cells, drop = FALSE]))
})
rownames(mean_mat) <- sig_genes   # genes × groups
```

---

## Step 5: LogFC Each Group vs All Others (Cell-Count-Weighted)

Weighted mean of other groups by cell count prevents unequal group sizes from
biasing the fold change toward small groups.

```r
n_cells_total <- ncol(obj)
n_cells_g     <- sapply(group_order, function(g) sum(meta[[GROUP_COL]] == g))

logfc_mat <- sapply(seq_along(group_order), function(i) {
  n_others    <- n_cells_total - n_cells_g[i]
  mean_others <- (rowSums(sweep(mean_mat, 2, n_cells_g, "*")) -
                    mean_mat[, i] * n_cells_g[i]) / n_others
  log2((mean_mat[, i] + EPS) / (mean_others + EPS))
})
colnames(logfc_mat) <- group_order
rownames(logfc_mat) <- sig_genes
```

---

## Step 6: Top Group Annotation + Top Subtype Within Top Group

```r
ranked_groups <- t(apply(mean_mat, 1, function(x)
  group_order[order(x, decreasing = TRUE)]
))

top_group_1 <- ranked_groups[, 1]
top_group_2 <- ranked_groups[, 2]

logfc_top1 <- logfc_mat[cbind(seq_along(sig_genes), match(top_group_1, group_order))]
logfc_top2 <- logfc_mat[cbind(seq_along(sig_genes), match(top_group_2, group_order))]

# Populate top_group2 only when second group is also genuinely enriched (logFC > 0.25)
two_tops <- logfc_top2 > 0.25

# ── Top subtype within top_group (optional — requires LABEL_COL and SUBTYPE_LABELS) ─
if (!is.null(LABEL_COL) && !is.null(SUBTYPE_LABELS) &&
    LABEL_COL %in% colnames(meta)) {

  mean_mat_subtype <- lapply(group_order, function(g) {
    sapply(SUBTYPE_LABELS, function(st) {
      cells <- rownames(meta)[meta[[GROUP_COL]] == g & meta[[LABEL_COL]] == st]
      if (length(cells) == 0) return(rep(0, length(sig_genes)))
      rowMeans(as.matrix(data_mat[sig_genes, cells, drop = FALSE]))
    })
  })
  names(mean_mat_subtype) <- group_order

  top_subtype <- sapply(seq_along(sig_genes), function(i) {
    g         <- top_group_1[i]
    sub_means <- mean_mat_subtype[[g]][i, ]
    SUBTYPE_LABELS[which.max(sub_means)]
  })

} else {
  top_subtype <- rep(NA_character_, length(sig_genes))
}
```

---

## Step 7: Assemble and Save CSV

```r
mean_exp_cols <- as.data.frame(expm1(mean_mat))   # natural scale for interpretability
colnames(mean_exp_cols) <- paste0("mean_exp_", gsub(" ", "_", group_order))

out_df <- data.frame(
  gene                      = sig_genes,
  top_group                 = top_group_1,
  logFC_vs_others_top_group = round(logfc_top1, 3),
  top_subtype               = top_subtype,
  top_group2                = ifelse(two_tops, top_group_2,    NA_character_),
  logFC_top2                = ifelse(two_tops, round(logfc_top2, 3), NA_real_),
  stringsAsFactors = FALSE
) %>%
  left_join(kw_df %>% select(gene, KW_H = H, p_val, p_adj, eps_sq), by = "gene") %>%
  arrange(p_adj) %>%
  cbind(mean_exp_cols[match(.$gene, sig_genes), , drop = FALSE])

# Do NOT append per-group logFC columns for every group — only top_group logFC is included.
# This keeps the table compact and avoids columns users never consult.

write.csv(out_df,
          file.path(OUTPUT_DIR, "comprehensive_DE_stats.csv"),
          row.names = FALSE)

message(sprintf("Output: %s (%d genes)",
  file.path(OUTPUT_DIR, "comprehensive_DE_stats.csv"), nrow(out_df)))
```

---

## Output Column Reference

| Column | Description |
|---|---|
| `gene` | Gene symbol |
| `top_group` | Group with highest mean expression |
| `logFC_vs_others_top_group` | log2FC of top_group vs all other groups (cell-count-weighted) |
| `top_subtype` | Highest-expressing cell subtype within top_group; NA if LABEL_COL not set |
| `top_group2` | Second group if also genuinely enriched (logFC > 0.25); else NA |
| `logFC_top2` | log2FC of top_group2 vs all others; else NA |
| `KW_H` | Kruskal-Wallis H statistic |
| `p_val` | Raw KW p-value |
| `p_adj` | BH-adjusted p-value |
| `eps_sq` | Epsilon-squared effect size (0–1) |
| `mean_exp_{group}` | Mean expression per group (natural/expm1 scale) |

Table is sorted by `p_adj` ascending (most significant first).

