# Doublet Removal Primitive

Uses scDblFinder for per-sample doublet detection on raw counts.
Applies in LargeDataset pipeline only. NOT for IntegratePublicData
(public atlas data is assumed pre-QC'd and would be processed data, not raw counts).

## Required R packages

Suspend startup messages and load:
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BiocParallel)

All available in the r-env conda environment.

## Function

Function: run_doublet_removal(obj, sample_col = "sample_id",
                               score_threshold = 0.7,
                               output_dir = "output/doublets/")

Arguments:
  obj             Seurat object with raw counts in RNA assay, merged across all
                  in-house samples. Must contain sample_col in metadata.
  sample_col      Metadata column identifying per-capture samples. Default "sample_id".
  score_threshold Score cutoff for confident-doublet removal (default 0.7).
  output_dir      Where to write the doublet annotation CSV and decision log entries.

Behavior:
  1. Extract raw counts via GetAssayData(obj, assay="RNA", layer="counts").
  2. Build a SingleCellExperiment and run scDblFinder with samples=obj[[sample_col]][,1].
     Use BPPARAM=SerialParam() for reproducibility on shared infrastructure.
  3. Write scDblFinder.score (numeric) and scDblFinder.class (singlet|doublet) into
     obj metadata. Always write both, regardless of removal decision.
  4. Apply SOFT REMOVAL POLICY: remove only cells where
       scDblFinder.class == "doublet" AND scDblFinder.score > score_threshold
     This is "Option B" -- belt-and-suspenders. The class label alone is not
     sufficient because its calibration threshold (~0.5) is the classifier's
     decision boundary, not a high-confidence call.
  5. Log per-sample doublet counts and removal counts to
     <output_dir>/decision_log.txt (append).
  6. Save the full per-cell annotation CSV (sample_id, score, class,
     removed_yes_or_no) to <output_dir>/doublet_annotations.csv.
  7. If any sample has >30% of cells flagged class=="doublet" (regardless of
     score), log a WARNING -- this suggests sample quality issues. Do not block.
  8. Return the filtered Seurat object.

## Reference implementation skeleton

```r
run_doublet_removal <- function(obj,
                                sample_col = "sample_id",
                                score_threshold = 0.7,
                                output_dir = "output/doublets/") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(output_dir, "decision_log.txt")

  counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
  samples <- obj[[sample_col, drop = TRUE]]

  sce <- scDblFinder(
    SingleCellExperiment(list(counts = counts_mat)),
    samples = samples,
    BPPARAM = SerialParam(),
    verbose = TRUE
  )

  obj$scDblFinder.score <- sce$scDblFinder.score
  obj$scDblFinder.class <- sce$scDblFinder.class

  # Soft removal policy: class == doublet AND score > threshold
  keep <- !(obj$scDblFinder.class == "doublet" &
            obj$scDblFinder.score > score_threshold)

  # Per-sample summary
  summary_df <- data.frame(
    sample = unique(samples),
    total = as.numeric(table(samples)),
    class_doublet = as.numeric(table(samples, obj$scDblFinder.class == "doublet")[, "TRUE"]),
    removed = as.numeric(table(samples, !keep)[, "TRUE"]),
    stringsAsFactors = FALSE
  )
  summary_df$doublet_class_pct <- 100 * summary_df$class_doublet / summary_df$total
  summary_df$removed_pct <- 100 * summary_df$removed / summary_df$total

  cat(sprintf("[%s] scDblFinder per-sample summary:\n", Sys.time()),
      file = log_file, append = TRUE)
  cat(paste(capture.output(print(summary_df)), collapse = "\n"),
      "\n", file = log_file, append = TRUE)

  # High-doublet-rate warning
  high_rate <- summary_df$sample[summary_df$doublet_class_pct > 30]
  if (length(high_rate) > 0) {
    cat(sprintf(
      "[%s] WARNING: samples with >30%% doublet class rate: %s\n",
      Sys.time(), paste(high_rate, collapse = ", ")),
      file = log_file, append = TRUE)
  }

  # Per-cell CSV
  write.csv(data.frame(
    cell_barcode = colnames(obj),
    sample = samples,
    scDblFinder.score = obj$scDblFinder.score,
    scDblFinder.class = obj$scDblFinder.class,
    removed = !keep
  ), file.path(output_dir, "doublet_annotations.csv"), row.names = FALSE)

  obj <- obj[, keep]
  cat(sprintf("[%s] Doublet removal complete. Kept %d / %d cells.\n",
              Sys.time(), sum(keep), length(keep)),
      file = log_file, append = TRUE)

  return(obj)
}
```

## Where to call this

**LargeDataset pipeline:** Insert after sample merge and basic QC (Stage 1 or 2),
before HVG selection, PCA, integration, or clustering. The merged object must have
the per-cell sample_id metadata column populated.

**Cell-type subclustering:** Call on the cell-type subset immediately after
subsetting and before NormalizeData. The subset retains raw counts and must have
the per-cell sample_col populated. This second pass catches doublets that are only
detectable within a single cell-type context (e.g. EC+EC from distinct functional
states that were masked in whole-object space).

## DO NOT call this in IntegratePublicData

The atlas data side of an integration is already-QC'd published data, not raw
counts. scDblFinder on the atlas side would (a) operate on normalized data and
give wrong results, and (b) second-guess the atlas's own QC pipeline. The
in-house side of an IntegratePublicData job that originates from a LargeDataset
handoff has already been doublet-filtered upstream, so it does not need a
second pass.
