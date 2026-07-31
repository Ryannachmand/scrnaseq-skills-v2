---
# Module: Label Harmonization
# Pipeline: IntegratePublicData
# Migrated from ~/claude-skills/pipelines/IntegratePublicData/methods/label_harmonization.md
# broad_label default parameterized; result_distribution example values documented in examples/
requires_context:
  palettes: []
  metadata_columns:
    required: []
    optional:
      - unified_label_col  # destination column for harmonized labels (default: "unified_label")
  brief_keys:
    required:
      - label_col   # dict: dataset_name → source metadata column holding original labels
    optional:
      - label_transfer_reference    # list of dataset names to use as transfer reference
      - label_transfer_score_high   # confidence threshold: high-confidence label (default: 0.75)
      - label_transfer_score_low    # confidence threshold: broad-label floor (default: 0.40)
      - label_transfer_broad_label  # broad fallback label for low-confidence cells — project_specific
      - known_label_mappings        # dict of pre-approved label mappings from brief
references:
  - "@primitives/seurat_v5_rules.md"
---

# Module: Label Harmonization

Harmonizes heterogeneous cell type labels from multiple datasets into a single
`unified_label` column (or a column named by `brief$unified_label_col`).

Two parts:
- **3A** — YAML-based mapping for datasets that already have cell type labels
- **3B** — Seurat label transfer for datasets with no existing labels

Set at top of script:
```r
UNIFIED_LABEL_COL <- "unified_label"  # brief: metadata.unified_label_col
OUTPUT_DIR        <- "output/label_harmonization"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
```

---

## Part 3A — YAML Mapping

### Generating the YAML

For each dataset, extract unique values from the designated label column:
```r
for (name in names(object_list)) {
  so   <- object_list[[name]]
  col  <- brief$label_col[[name]]   # source label column per dataset
  lbls <- unique(na.omit(so@meta.data[[col]]))
  cat("\n", name, "—", col, ":\n")
  print(sort(lbls))
}
```

Write `output/label_harmonization.yaml` pre-filled with known mappings from
`brief$known_label_mappings`. Mark entries not in the brief as `REVIEW`.

### Applying approved mappings

```r
# After user approves the YAML:
yaml_data <- yaml::read_yaml(file.path(OUTPUT_DIR, "label_harmonization.yaml"))

for (name in names(object_list)) {
  so      <- object_list[[name]]
  mapping <- yaml_data[[name]]$labels
  col     <- yaml_data[[name]]$source_column

  if (!is.null(col) && col != "") {
    orig  <- as.character(so@meta.data[[col]])
    mapped <- mapping[orig]
    # Unmapped values → Unknown (see Common Pitfalls below)
    mapped[is.na(mapped)] <- "Unknown"
    so@meta.data[[UNIFIED_LABEL_COL]] <- unname(mapped)
  }
  object_list[[name]] <- so
}
```

---

## Part 3B — Seurat Label Transfer

Use when public datasets have no usable existing labels.
Requires a trusted reference object specified via `brief$label_transfer_reference`.

### Building the reference

```r
ref_objects <- object_list[brief$label_transfer_reference]
ref_merged  <- merge(ref_objects[[1]], y = ref_objects[-1])
ref_merged  <- JoinLayers(ref_merged, assay = "RNA")  # seurat_v5_rules.md Rule 1
ref_merged  <- NormalizeData(ref_merged) %>%
               FindVariableFeatures() %>%
               ScaleData(features = VariableFeatures(.)) %>%
               RunPCA()
```

### Running transfer

```r
# broad_label: biology-appropriate fallback label for cells with moderate confidence.
# Must be defined in brief$label_transfer_broad_label.
# Example: "Mesenchymal" for stromal datasets, "Epithelial" for epithelial datasets.
# No default — caller must specify; leaving this blank will produce "Unknown" for all
# moderate-confidence cells, losing potentially informative low-confidence assignments.
broad_label <- brief$label_transfer_broad_label
if (is.null(broad_label)) {
  stop("brief$label_transfer_broad_label is required. ",
       "Set to a biology-appropriate broad label for moderate-confidence cells ",
       "(e.g., 'Mesenchymal', 'Stromal', 'Epithelial') in the project brief.")
}

for (name in unlabeled_datasets) {
  query <- object_list[[name]]

  # Skip if too few cells — TransferData will fail with k.weight error
  if (ncol(query) < 50) {
    query[[UNIFIED_LABEL_COL]] <- "Unknown"
    message("Skipping transfer for ", name, " — too few cells: ", ncol(query))
    object_list[[name]] <- query
    next
  }

  query <- NormalizeData(query) %>% FindVariableFeatures() %>%
           ScaleData(features = VariableFeatures(.)) %>% RunPCA()

  anchors <- FindTransferAnchors(
    reference           = ref_merged,
    query               = query,
    dims                = 1:30,
    reference.reduction = "pca"
  )
  preds <- TransferData(anchorset = anchors,
                        refdata   = ref_merged[[UNIFIED_LABEL_COL]],
                        dims      = 1:30)
  query <- AddMetaData(query, metadata = preds)

  # Confidence thresholds from brief (or lab defaults)
  score_high <- brief$label_transfer_score_high %||% 0.75
  score_low  <- brief$label_transfer_score_low  %||% 0.40

  query[[UNIFIED_LABEL_COL]] <- dplyr::case_when(
    query$prediction.score.max >= score_high ~ query$predicted.id,
    query$prediction.score.max >= score_low  ~ broad_label,
    TRUE                                     ~ "Unknown"
  )

  object_list[[name]] <- query
  cat(name, "label transfer complete:\n")
  print(table(query[[UNIFIED_LABEL_COL]]))
}
```

### Document transfer results in YAML

After running transfer, add a `transferred` section to the YAML for traceability:
```yaml
<dataset_name>:
  source_column: transferred
  notes: "Seurat label transfer from <reference_dataset_1> + <reference_dataset_2>"
  score_high: 0.75
  score_low: 0.40
  result_distribution:
    <label_1>: <count>
    <label_2>: <count>
    Unknown: <count>
```
Replace angle-bracket placeholders with actual dataset names and observed counts.

---

## Controlled Vocabulary Rules

The unified label vocabulary is always defined in the brief. General principles:

1. **Specific before broad** — use specific labels when confidence is high
2. **Unknown, not blank** — never leave `unified_label` as NA after Stage 3
3. **Contamination** — non-target cells that passed the filter; distinct from Unknown
4. **RNAlo / QC pool survivors** — cells with low RNA content that survived targeted QC;
   label communicates identity without discarding them

### Special label situations

**Sort-strategy-based labeling:**
When cell type labels are unavailable but sort strategy is known, use the sort strategy
as the label (e.g., "CD14 Isolated Stroma", "PE Isolated Stroma").
These are informative for provenance even without transcriptomic annotation.

**Passaged or cultured cells:**
When all cells in a dataset come from a single experimental condition with no subclustering,
label all as a single informative label (e.g., "Passaged Stroma", "Cultured <CellType>").
Subclusters can be resolved after integration.

**Derived display columns:**
For heatmaps and cohort plots, create a derived `heatmap_label` column rather than
modifying `unified_label`. Example: prepend "Cultured " and append "-like" suffix.
Never modify `unified_label` for display purposes — create a separate column.

*Note: worked examples of these harmonization patterns live in the
examples/ directory; see any `*_label_harmonization.md` file for a
project-specific application.*

---

## Common Pitfalls

### Named vector mapping fails silently

When using `mapping[orig]` to remap labels, NAs are produced for any `orig`
value not in `mapping`. Always check and replace with "Unknown":
```r
mapped <- mapping[orig]
mapped[is.na(mapped)] <- "Unknown"
```

### Seurat v5 metadata assignment

Use direct slot assignment when assigning a vector of the same length as `ncol(so)`.
See @primitives/seurat_v5_rules.md Rule 3 for the named-vector keyed-by-cluster case:
```r
so@meta.data[[UNIFIED_LABEL_COL]] <- mapped_vector   # correct
# so$unified_label <- mapped_vector                  # risky in v5 with named vectors
```
