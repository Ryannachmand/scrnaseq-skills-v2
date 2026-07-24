---
# Pseudobulk edgeR — Canonical Recipe — v2
# Authored in v2 Phase 4+ (new method)
# Written for: ThymusEC re-analysis (redmondlab4, 2026-07-23)
# FLAGGED FOR LIBRARY REVIEW — new primitive, not yet in Phase 1-3 library
#
# Handles: edgeR glmQLFit on a pre-built DGEList (e.g. from pseudobulk aggregation).
# Does NOT build the DGEList — assumes caller provides it (e.g. loaded from RDS).
# Design note: intended for designs with replicated pseudobulk samples (n >= 2 per group).
# For n=1 designs, the resulting DE is DESCRIPTIVE ONLY — see design_warning below.
requires_context:
  palettes:
    - group_colors    # optional — named vector for volcano/heatmap coloring
  metadata_columns:
    required:
      - sample_col    # column in DGEList$samples identifying pseudobulk samples
      - group_col     # column in DGEList$samples holding the comparison group variable
    optional:
      - subtype_col   # column in DGEList$samples for EC subtype / cell type (for per-subtype loops)
  brief_keys:
    required:
      - output_dir
    optional:
      - downstream_analyses.pseudobulk_edger.contrasts   # list of contrast specs
      - downstream_analyses.pseudobulk_edger.fdr_cut
      - downstream_analyses.pseudobulk_edger.lfc_cut
      - downstream_analyses.pseudobulk_edger.design_formula   # e.g. "~ 0 + group_col + covariate"
---

# Pseudobulk edgeR — Canonical Recipe

Run edgeR glmQLFit on a pre-built DGEList object of pseudobulk aggregated counts.
Handles: filtering low-count genes, TMM normalization, dispersion estimation,
model fitting, contrasts, and FDR-controlled output.

**Design-replicated tier (inference-valid):** requires n >= 2 biological samples per
comparison group. Produces p-values and FDR with the usual glmQLFit interpretation.

**Design-unreplicated tier (descriptive-only):** for n=1 sort designs (e.g. single
EC-sort library per age × timepoint), biological replication is absent. Results are
logFC estimates and tagwise dispersions only — no valid p-values. All outputs must be
labeled `[DESCRIPTIVE — no replication; p-values unreliable]`.

---

## Critical Constraints

| Do NOT | Do Instead | Why |
|---|---|---|
| Interpret p-values from n=1 designs | Label as descriptive; use logFC magnitude only | edgeR uses dispersion estimates from pooled data when n=1 — p-values are anti-conservative |
| Run calcNormFactors after subsetting | Always call on full DGEList before subsetting | Norm factors are computed on all samples jointly for comparability |
| Use LRT instead of QLFit | Always `glmQLFit` + `glmQLFTest` | QL approach is more conservative and appropriate for bulk-like pseudobulk data |
| Hardcode column names in model matrix | Always use `model.matrix(~ 0 + .)` with dynamically constructed formula | Brittle to column reordering in samples metadata |

---

## Configuration Block

```r
# ── Thresholds ────────────────────────────────────────────────────────────────
FDR_CUT     <- 0.05
LFC_CUT     <- 0.5        # |log2FC| minimum
MIN_CPM     <- 1          # minimum CPM for gene filtering
MIN_SAMPLES <- 2          # gene must pass MIN_CPM in at least this many samples
set.seed(42)

# ── DGEList loading ───────────────────────────────────────────────────────────
DGELIST_PATH  <- "project_specific"   # REPLACE: path to pre-built DGEList RDS
OUTPUT_DIR    <- "output/pseudobulk_edgeR"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Design ────────────────────────────────────────────────────────────────────
SAMPLE_COL  <- "project_specific"   # REPLACE: column in dgl$samples (e.g. "library")
GROUP_COL   <- "project_specific"   # REPLACE: column in dgl$samples (e.g. "irradiated_day")
SUBTYPE_COL <- "project_specific"   # REPLACE: EC subtype column (e.g. "cell_type_subset"), or NULL for whole
DESIGN_FORMULA <- NULL              # NULL = "~ 0 + GROUP_COL"; set formula string to override

# ── Contrasts ─────────────────────────────────────────────────────────────────
# Each contrast: list(label, numerator, denominator)
# numerator / denominator = values of GROUP_COL to compare
contrasts_list <- list(
  list(label = "project_specific", numerator = "project_specific", denominator = "project_specific")
  # Add additional entries for additional comparisons
)

# ── Design replication status — MUST set before running ──────────────────────
# Set TRUE when n >= 2 biological replicates per group exist (safe to report p-values)
# Set FALSE for n=1-per-group designs (descriptive logFC only; p-values unreliable)
REPLICATED_DESIGN <- TRUE   # REPLACE: FALSE for thymus EC sort radiation axis
```

---

## Step 1: Load and Filter

```r
library(edgeR)
library(ggplot2)
library(dplyr)

dgl <- readRDS(DGELIST_PATH)

# Confirm sample metadata columns exist
stopifnot(SAMPLE_COL  %in% colnames(dgl$samples))
stopifnot(GROUP_COL   %in% colnames(dgl$samples))
cat("Samples:", nrow(dgl$samples), "\n")
cat("Genes (pre-filter):", nrow(dgl$counts), "\n")
cat("Group distribution:\n")
print(table(dgl$samples[[GROUP_COL]]))

# Filter low-count genes (MIN_CPM in at least MIN_SAMPLES samples)
keep <- filterByExpr(dgl, group = dgl$samples[[GROUP_COL]],
                     min.count = MIN_CPM, min.total.count = 5)
dgl_filt <- dgl[keep, , keep.lib.sizes = FALSE]
cat("Genes post-filter:", nrow(dgl_filt$counts), "\n")

# TMM normalization (always on filtered object)
dgl_filt <- calcNormFactors(dgl_filt, method = "TMM")
```

---

## Step 2: Dispersion and Model Fit

```r
run_pseudobulk_edgeR <- function(dgl_obj, group_col, design_formula = NULL) {
  # dgl_obj:        filtered DGEList (post-filterByExpr + calcNormFactors)
  # group_col:      column in dgl_obj$samples for group factor
  # design_formula: optional string formula; default = "~ 0 + group_col"
  #
  # Returns: list(dgl_fit, design_mat, contrast_names)

  grp_factor <- factor(dgl_obj$samples[[group_col]])
  grp_factor <- droplevels(grp_factor)

  if (is.null(design_formula)) {
    design_mat <- model.matrix(~ 0 + grp_factor)
    colnames(design_mat) <- gsub("grp_factor", "", colnames(design_mat))
  } else {
    design_mat <- model.matrix(as.formula(design_formula), data = dgl_obj$samples)
  }

  cat("Design matrix columns:", paste(colnames(design_mat), collapse = ", "), "\n")
  cat("Rank:", Matrix::rankMatrix(design_mat), "\n")

  dgl_obj <- estimateDisp(dgl_obj, design_mat, robust = TRUE)
  dgl_fit <- glmQLFit(dgl_obj, design_mat, robust = TRUE)
  gc()

  list(dgl_fit = dgl_fit, design_mat = design_mat)
}

fit_obj <- run_pseudobulk_edgeR(dgl_filt, GROUP_COL, DESIGN_FORMULA)
```

---

## Step 3: Contrast Testing

```r
get_edger_contrast <- function(fit_obj, numerator, denominator, label,
                                output_dir, fdr_cut = FDR_CUT, lfc_cut = LFC_CUT,
                                replicated = TRUE) {
  # fit_obj:      list from run_pseudobulk_edgeR (dgl_fit, design_mat)
  # numerator:    value in GROUP_COL that is "up" (positive logFC)
  # denominator:  value in GROUP_COL that is "down" (negative logFC)
  # label:        filename-safe label for outputs
  # replicated:   FALSE = descriptive design; suppress p-value interpretation
  #
  # Returns: data frame with gene, logFC, logCPM, F (or LR), PValue, FDR, direction

  dgl_fit <- fit_obj$dgl_fit
  design_mat <- fit_obj$design_mat

  # Build contrast vector
  col_num <- which(colnames(design_mat) == numerator)
  col_den <- which(colnames(design_mat) == denominator)
  if (length(col_num) == 0 || length(col_den) == 0) {
    stop(sprintf("Contrast columns not found: '%s' vs '%s'. Available: %s",
                 numerator, denominator, paste(colnames(design_mat), collapse = ", ")))
  }

  contrast_vec <- numeric(ncol(design_mat))
  contrast_vec[col_num] <-  1
  contrast_vec[col_den] <- -1

  qlf <- glmQLFTest(dgl_fit, contrast = contrast_vec)
  results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  results$gene <- rownames(results)
  results$direction <- ifelse(results$logFC > 0, numerator, denominator)
  rownames(results) <- NULL

  # Label unreplicated results
  design_label <- if (!replicated) "[DESCRIPTIVE - no replication; p-values unreliable]" else ""

  out_prefix <- file.path(output_dir, label)
  write.csv(results, paste0(out_prefix, "_edgeR_full.csv"), row.names = FALSE)

  n_sig <- sum(results$FDR < fdr_cut & abs(results$logFC) > lfc_cut)
  cat(sprintf("%s %s vs %s: %d genes FDR<%g |logFC|>%g  %s\n",
              label, numerator, denominator, n_sig, fdr_cut, lfc_cut, design_label))

  results
}
```

---

## Step 4: MDS Plot (QC)

```r
make_edger_mds_plot <- function(dgl_obj, group_col, sample_col,
                                 group_colors = NULL, output_file) {
  # dgl_obj:      DGEList post-normalization
  # group_col:    column for color coding
  # sample_col:   column for point labels
  # group_colors: optional named vector; if NULL uses Okabe-Ito

  mds <- plotMDS(dgl_obj, plot = FALSE)
  mds_df <- data.frame(
    sample = rownames(mds$distance.matrix.squared),
    Dim1   = mds$x,
    Dim2   = mds$y
  )
  mds_df[[group_col]] <- dgl_obj$samples[[group_col]][match(mds_df$sample, rownames(dgl_obj$samples))]
  mds_df[[sample_col]] <- dgl_obj$samples[[sample_col]][match(mds_df$sample, rownames(dgl_obj$samples))]

  PALETTE_OKABE_ITO <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#000000")
  if (is.null(group_colors)) {
    grps <- unique(mds_df[[group_col]])
    group_colors <- setNames(PALETTE_OKABE_ITO[seq_len(length(grps))], grps)
  }

  p <- ggplot(mds_df, aes_string(x = "Dim1", y = "Dim2", color = group_col, label = sample_col)) +
    geom_point(size = 3.5, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.4) +
    scale_color_manual(values = group_colors) +
    labs(title = "Pseudobulk MDS — leading log-FC", x = "Dim 1", y = "Dim 2") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(output_file, plot = p, width = 8, height = 6,
         units = "in", device = "pdf", useDingbats = FALSE)
  invisible(p)
}
```

---

## Step 5: Volcano Plot (edgeR results)

```r
make_edger_volcano <- function(results, comp_label, subset_name,
                                numerator, denominator,
                                fdr_cut = FDR_CUT, lfc_cut = LFC_CUT,
                                replicated = TRUE,
                                n_label = 20, output_file) {
  # results: data frame from get_edger_contrast()
  # replicated: if FALSE, adds [DESCRIPTIVE] tag to title and removes FDR axis interpretation

  df <- results %>%
    mutate(
      neg_log_fdr = -log10(FDR + 1e-300),
      sig  = FDR < fdr_cut & abs(logFC) > lfc_cut,
      dir  = case_when(
        sig & logFC >  0 ~ "up_num",
        sig & logFC <= 0 ~ "up_den",
        TRUE             ~ "ns"
      ),
      label = ""
    )

  top_num <- df %>% filter(dir == "up_num") %>% arrange(FDR) %>% head(n_label) %>% pull(gene)
  top_den <- df %>% filter(dir == "up_den") %>% arrange(FDR) %>% head(n_label) %>% pull(gene)
  df$label <- ifelse(df$gene %in% c(top_num, top_den), df$gene, "")

  colors <- c(up_num = "#B2182B", up_den = "#2166AC", ns = "grey78")
  title_tag <- if (!replicated) " [DESCRIPTIVE ONLY]" else ""

  p <- ggplot(df, aes(x = logFC, y = neg_log_fdr, color = dir, size = dir)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey50") +
    ggrepel::geom_text_repel(aes(label = label), size = 3.5, fontface = "italic",
                              max.overlaps = 30, box.padding = 0.4, segment.size = 0.3) +
    scale_color_manual(values = colors, guide = "none") +
    scale_size_manual(values = c(up_num = 1.3, up_den = 1.3, ns = 0.7), guide = "none") +
    labs(
      title    = sprintf("%s | %s%s", subset_name, gsub("_", " ", comp_label), title_tag),
      subtitle = sprintf("Up in %s: %d  |  Up in %s: %d\n(FDR < %g, |logFC| > %g)",
                         numerator, sum(df$dir == "up_num"),
                         denominator, sum(df$dir == "up_den"), fdr_cut, lfc_cut),
      x = "log2 Fold Change", y = expression(-log[10](FDR))
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title    = element_text(hjust = 0.5, size = 17, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey35"),
          panel.grid.major = element_blank())

  ggsave(output_file, plot = p, width = 9, height = 8,
         units = "in", device = "pdf", useDingbats = FALSE)
  invisible(p)
}
```

---

## Per-Subtype Loop Pattern

```r
# When SUBTYPE_COL is defined: loop over EC subtypes, subsetting dgl_filt each time.
# When NULL: run on all pseudobulk samples together.

subtypes_to_run <- if (!is.null(SUBTYPE_COL) && SUBTYPE_COL %in% colnames(dgl_filt$samples)) {
  unique(dgl_filt$samples[[SUBTYPE_COL]])
} else {
  "AllEC"
}

for (subtype in subtypes_to_run) {
  cat("\n=== Subtype:", subtype, "===\n")

  dgl_sub <- if (subtype == "AllEC") {
    dgl_filt
  } else {
    keep_samp <- dgl_filt$samples[[SUBTYPE_COL]] == subtype
    if (sum(keep_samp) < 2) {
      message(sprintf("SKIPPED %s — fewer than 2 samples", subtype)); next
    }
    dgl_filt[, keep_samp, keep.lib.sizes = FALSE]
  }

  # Re-normalize and re-fit within subtype
  dgl_sub <- calcNormFactors(dgl_sub, method = "TMM")
  fit_sub  <- run_pseudobulk_edgeR(dgl_sub, GROUP_COL, DESIGN_FORMULA)

  out_sub <- file.path(OUTPUT_DIR, gsub("[: ]", "_", subtype))
  dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)

  # MDS plot
  make_edger_mds_plot(dgl_sub, group_col = GROUP_COL, sample_col = SAMPLE_COL,
                       output_file = file.path(out_sub, "mds_plot.pdf"))

  # Run each contrast
  for (ct in contrasts_list) {
    results <- get_edger_contrast(fit_sub, ct$numerator, ct$denominator, ct$label,
                                   output_dir = out_sub,
                                   fdr_cut = FDR_CUT, lfc_cut = LFC_CUT,
                                   replicated = REPLICATED_DESIGN)
    if (is.null(results)) next

    make_edger_volcano(results, ct$label, subtype,
                        ct$numerator, ct$denominator,
                        fdr_cut = FDR_CUT, lfc_cut = LFC_CUT,
                        replicated = REPLICATED_DESIGN,
                        output_file = file.path(out_sub, paste0(ct$label, "_volcano.pdf")))
  }
}
```

---

## Output File Naming

```
{output_dir}/{subtype}/{label}_edgeR_full.csv
{output_dir}/{subtype}/{label}_volcano.pdf
{output_dir}/{subtype}/mds_plot.pdf
```

For descriptive (unreplicated) results, the CSV filename gains suffix `_DESCRIPTIVE`:
```
{output_dir}/{subtype}/{label}_edgeR_full_DESCRIPTIVE.csv
```

Callers should add `_DESCRIPTIVE` to the label argument when `REPLICATED_DESIGN = FALSE`.

---

## Design Warning (ALWAYS READ before using on n=1 data)

When `REPLICATED_DESIGN = FALSE`, edgeR estimates dispersion by borrowing across genes
(estimateDisp with robust=TRUE), but there is only one observation per group. The
resulting test statistics and p-values are not valid for inference. Report ONLY:
- logFC magnitude and direction (biologically meaningful)
- logCPM (expression level evidence)
- MDS plot (shows overall sample structure)

Never report FDR < threshold from an unreplicated design as evidence of significant
differential expression. Use the pseudobulk results for exploration and hypothesis
generation only. The aging-axis comparisons (with multiple libraries per age group)
ARE replicated and CAN support inference.

---

## Brief Configuration

```yaml
downstream_analyses:
  pseudobulk_edger:
    enabled: true
    dgelist_path: project_specific       # path to pre-built DGEList RDS
    sample_col: project_specific         # column in DGEList$samples identifying samples
    group_col: project_specific          # grouping column for contrasts
    subtype_col: project_specific        # optional: EC subtype column (null = whole-object)
    replicated_design: true              # FALSE for n=1-per-group designs (descriptive only)
    design_formula: null                 # null = "~ 0 + group_col"
    fdr_cut: 0.05
    lfc_cut: 0.5
    contrasts:
      - label: project_specific
        numerator: project_specific
        denominator: project_specific
```
