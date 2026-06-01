---
# Large Dataset Pipeline — scRNAseq from assembled multi-sample files
# Pipeline: large_dataset
# Version: v2
# Authored: Phase 3
requires_context:
  always_required:
    - "@primitives/r_environment.md"
    - "@primitives/seurat_v5_rules.md"
    - "@primitives/aesthetics.md"
    - "@primitives/visualization.md"
    - "@primitives/harmony_integration.md"
    - "@context/lab_context.md"
    - "@context/color_palettes.md"
  when_downstream_de:
    - "@primitives/differential_expression.md"
references:
  stage_1_load:
    - "@modules/load_formats.md"
  stage_8_modules:
    - "@modules/feature_umap_plot.md"
    - "@modules/celltype_subclustering.md"
    - "@modules/metabolic_profile.md"
    - "@modules/trajectory_monocle3.md"
    - "@modules/pyscenic_regulons.md"
    - "@modules/bulk_concordance.md"
    - "@modules/de_comprehensive_csv.md"
    - "@modules/cellchat.md"
    - "@modules/multi_group_de_analysis.md"
---

# Large Dataset Pipeline

Multi-sample scRNAseq analysis pipeline for datasets assembled from .h5, .rds, or MEX
files. Covers assembly through annotation, subset analysis, and optional downstream
modules dispatched from the brief. Expect 9 stages and 5 checkpoints from raw data to
complete analysis.

---

## Universal Rules — apply at every stage and in every continuation session

- Before writing any plot code (UMAP, violin, dotplot, heatmap, volcano, bar,
  chord, or any other visualization), re-read @primitives/aesthetics.md and
  @primitives/visualization.md. Do this even in continuation sessions and even
  when regenerating a plot already produced earlier in this job.
- Before writing any R script, re-read @primitives/r_environment.md and
  @primitives/seurat_v5_rules.md.
- These re-reads are mandatory at every stage. Do not rely on recalled content
  from earlier in the session or from a prior continuation session.

---

## Stage 1 — Assemble

**Depends on:** Nothing (pipeline start)
**Brief keys consumed:** `inputs.type`, `inputs.paths`, `inputs.metadata_csv`
**Reference:** @modules/load_formats.md

Load each file listed in `brief.inputs.paths` using the format-appropriate recipe from
@modules/load_formats.md: `Read10X_h5()` for `h5`, `Read10X()` for `mex`, direct load
for `rds`. If any path yields a list (multi-modal h5), apply the `is.list(counts)` check
from load_formats before creating the Seurat object.

Parse metadata from filenames or from `inputs.metadata_csv`. If `metadata_csv` is
provided, merge the CSV by sample identifier after object creation. If absent, infer
sample IDs from file path components per the project naming convention.

Merge all per-sample objects into one Seurat object with `merge()`. Do NOT call
`JoinLayers()` here — layers remain split by sample until Stage 3 (see
@primitives/seurat_v5_rules.md Rule 1).

Compute `percent.mt` from MT gene patterns (`^MT-` for human, `^mt-` for mouse per
`brief.organism`). Save the unfiltered merged object to `output/qc/raw.rds`.

---

## Stage 2 — QC

**Depends on:** Stage 1 complete
**Brief keys consumed:** `pipeline_params.qc_thresholds`
**Checkpoint:** CP1 — User reviews violin plots and approves filter thresholds

Load `output/qc/raw.rds`. Generate three violin plots — `nFeature_RNA`, `nCount_RNA`,
`percent.mt` — split by `pipeline_params.batch_correction_var`. Save to `output/qc/`.

Compute proposed thresholds using the percentiles in `pipeline_params.qc_thresholds`:
- `nFeature_RNA_min` at `n_feature_min_percentile` (default 2nd percentile)
- `nFeature_RNA_max` at `n_feature_max_percentile` (default 98th percentile)
- `percent_mt_max` from the fixed value (default 20)

Write proposed thresholds plus per-sample summary statistics (n_cells, median nFeature,
median nCount) to `output/qc/thresholds_proposed.yaml`.

**CP1:** Present violin plots and `thresholds_proposed.yaml`. The user approves or
modifies the proposed thresholds. After approval, apply the filters and save the filtered
object to `output/qc/filtered.rds`. Do not proceed to Stage 3 until CP1 is approved.

---

## Stage 3 — Process whole object

**Depends on:** Stage 2 complete, CP1 approved, `output/qc/filtered.rds` exists
**Brief keys consumed:** `pipeline_params.n_variable_features`, `pipeline_params.n_pcs`,
  `pipeline_params.clustering_resolution`, `pipeline_params.batch_correction_var`
**Reference:** @primitives/harmony_integration.md
**Checkpoint:** CP2 — User reviews elbow plot and UMAPs

Load `output/qc/filtered.rds`. Execute the canonical Seurat+Harmony processing
sequence from @primitives/harmony_integration.md using `pipeline_params.*` values:

1. `NormalizeData()`
2. `FindVariableFeatures(nfeatures = n_variable_features)` (default 4000)
3. `JoinLayers()` — required before ScaleData on multi-sample objects (Rule 1)
4. `ScaleData(features = VariableFeatures(so))` — variable features only (Rule 2)
5. `RunPCA(npcs = n_pcs)` (default 30)
6. `RunHarmony(group.by.vars = batch_correction_var)` → reduction `harmony`
7. `FindNeighbors(reduction = "harmony", dims = 1:n_pcs)`
8. `FindClusters(resolution = clustering_resolution)` (default 0.8)
9. `RunUMAP(reduction = "harmony", dims = 1:n_pcs)` → stored as `harmony.umap`

Save processed object to `output/processed/processed.rds`.

**CP2:** Generate elbow plot (PCA variance explained) and two diagnostic UMAPs:
(A) colored by Seurat cluster, (B) colored by `batch_correction_var`. Present for
user review. Confirm that Harmony removed batch structure while preserving biological
clusters. The user may request re-processing with adjusted `n_pcs` or
`clustering_resolution` before proceeding.

---

## Stage 4 — Annotation

**Depends on:** Stage 3 complete, CP2 approved, `output/processed/processed.rds` exists
**Brief keys consumed:** `annotation.top_n_markers`
**References:** @primitives/differential_expression.md (for `is_ambient()`)
**Checkpoint:** CP3 — User provides cluster-to-cell-type label mapping

Load `output/processed/processed.rds`. Run `FindAllMarkers()` across all clusters
(`min.pct = 0.25`, `logfc.threshold = 0.25`, `test.use = "wilcox"`). Apply
`is_ambient()` from @primitives/differential_expression.md to exclude ambient
contamination genes from the output marker table.

Write the top `annotation.top_n_markers` (default 25) markers per cluster to
`output/annotation/markers.yaml`, organized by cluster number with gene name,
log2FC, adjusted p-value, and pct.1/pct.2.

**CP3:** Present `markers.yaml` to the user. The user supplies a cluster-to-label
mapping (YAML dict: cluster_id → cell_type_name). Apply the mapping to Seurat
metadata under the `label_col` column (default column name: `cell_type`). Save
the annotated object to `output/annotation/annotated.rds`.

---

## Stage 5 — Whole-object visualization

**Depends on:** Stage 4 complete, CP3 labels applied, `output/annotation/annotated.rds` exists
**Brief keys consumed:** `context_overrides.palettes.group_colors`,
  `context_overrides.metadata_columns.group_col`,
  `context_overrides.metadata_columns.sample_col`
**Reference:** @primitives/visualization.md
**Checkpoint:** CP4 — User reviews plot aesthetics and approves style

Load `output/annotation/annotated.rds`. Generate whole-object visualizations using
functions from @primitives/visualization.md:

- `make_labeled_umap()` — UMAP colored by cell type label with LabelClusters overlay
- `make_stacked_violin()` — canonical marker genes per cell type
- `make_proportion_plot()` — cell type proportions by sample and by group

Use color vectors from @context/color_palettes.md for unlabeled groups; override with
`context_overrides.palettes.group_colors` for project-specific group coloring. The
`group_col` from `context_overrides.metadata_columns` (if set) is used for group
stratification in the proportion plot.

Save all plots to `output/visualization/`.

**CP4:** Present all whole-object plots. The user confirms colors, label placement,
and plot dimensions are acceptable. The user may supply updated `group_colors` or
adjusted axis orders before sign-off.

---

## Stage 6 — Subset and re-process

**Depends on:** Stage 4 complete (CP3 labels applied)
**Brief keys consumed:** `subset.target_celltypes`, `subset.subset_name`, `subset.pipeline_params`
**Reference:** @primitives/harmony_integration.md
**Checkpoint:** CP5 — User provides subtype label mapping

If `brief.subset.target_celltypes` is null or absent, skip Stages 6 and 7 (proceed
directly to Stage 8 or Stage 9 if Stage 8 is also empty).

Load `output/annotation/annotated.rds`. Subset cells whose `label_col` value is in
`subset.target_celltypes`. Use factor-safe subsetting:
`so[, so@meta.data[[label_col]] %in% target_celltypes]`.

Re-run the full Seurat+Harmony processing pipeline on the subset using parameters from
`subset.pipeline_params` (if provided) or default subset values from
@primitives/harmony_integration.md (n_variable_features: 2850, n_pcs: 40,
clustering_resolution: 0.39). See @primitives/harmony_integration.md for the complete
canonical sequence. Save to `output/subset/subset_processed.rds`.

Run `AddModuleScore()` for marker gene sets to guide initial subtype annotation. Run
`FindAllMarkers()` on subset clusters and apply `is_ambient()`. Write top markers per
cluster to `output/subset/markers.yaml`.

**CP5:** Present subset UMAP (colored by cluster), `AddModuleScore` feature UMAPs, and
`markers.yaml`. The user provides a cluster-to-subtype label mapping and optionally a
list of clusters to exclude. Apply labels under `subtype_col` (default: `subtype`).
Save labeled subset to `output/subset/subset_annotated.rds`.

---

## Stage 7 — Core subset analysis

**Depends on:** Stage 6 complete, CP5 labels applied, `output/subset/subset_annotated.rds` exists
**Brief keys consumed:** `context_overrides.palettes.subtype_colors`,
  `context_overrides.metadata_columns.subtype_col`,
  `context_overrides.metadata_columns.group_col`
**References:** @primitives/visualization.md, @primitives/differential_expression.md

Load `output/subset/subset_annotated.rds`. Generate core analytical plots using
functions from @primitives/visualization.md:

- `make_canonical_dotplot()` — curated marker genes organized by functional section,
  colored by average expression (PALETTE_EXPRESSION), sized by percent expressing,
  shape 21; gene sections are caller-provided from context or established during CP5
- `make_tf_diamond_plot()` — TF fold changes as diamond points (shape 23), produced
  only when a TF gene list is available from context or from a dispatched
  `multi_group_de` module; requires `group_col` to be set
- `make_stacked_violin()` — marker genes per subtype, dynamically sized
- `make_diff_abundance_heatmap()` — differential abundance by group, produced only
  when `group_col` is set in `context_overrides.metadata_columns`

Save all outputs to `output/subset_analysis/`. Stage 7 produces the core manuscript
figures for the subset. No additional checkpoint follows Stage 7.

---

## Stage 8 — Optional downstream analyses

**Depends on:** Stage 7 complete (or Stage 6 if Stage 7 produced no subset object)
**Brief keys consumed:** `downstream_analyses.*`
**References:** All modules listed in `references.stage_8_modules`

If `brief.downstream_analyses` is absent or empty (`{}`), skip Stage 8 and proceed to
Stage 9.

### 8.1 Brief key → module file mapping

The brief key in `downstream_analyses` maps to the module file as follows:

| Brief key | Module file |
|---|---|
| `multi_group_de` | @modules/multi_group_de_analysis.md |
| `feature_umap_plot` | @modules/feature_umap_plot.md |
| `subclustering` | @modules/celltype_subclustering.md |
| `metabolic_profile` | @modules/metabolic_profile.md |
| `trajectory` | @modules/trajectory_monocle3.md |
| `pyscenic` | @modules/pyscenic_regulons.md |
| `bulk_concordance` | @modules/bulk_concordance.md |
| `de_comprehensive_csv` | @modules/de_comprehensive_csv.md |
| `cellchat` | @modules/cellchat.md |
| `deg` | @primitives/differential_expression.md [1] |

[1] `deg` dispatches to a **primitive** (not a module). It is invoked via an explicit
`optional_analyses.deg` configuration block in the brief. Load
`@primitives/differential_expression.md` before executing any comparison in this key.

### 8.2 Three-axis DE pattern

When both a condition column (`group_col`) and a cluster/subtype column (`label_col`) exist
on the object, the default DE pattern covers three axes. All three should run by default
unless the brief explicitly suppresses one:

- **Axis A (cluster markers):** Produced by `@modules/celltype_subclustering.md` during
  annotation, or by Stage 4 for top-level clustering. FindAllMarkers, condition-agnostic,
  one-vs-rest per cluster. NOT produced by the `deg` primitive.
- **Axis B (condition global):** Via `deg.axes.condition_global: true`. Runs the primary
  condition comparison across all cells (or scopes) using `@primitives/differential_expression.md`.
- **Axis C (condition per cluster):** Via `deg.axes.condition_per_cluster: true` with
  `label_col` set. Iterates over every unique value of `label_col` and runs the condition
  comparison within each subset. Gated by 100 cells per side per cluster; underpowered
  clusters are skipped and logged. Output: `DE_full_{comp_label}_within_{label_value}.csv`.

### 8.3 Gene-set dotplot enforcement

Any gene-set dotplot in any module or job script must call `make_functional_dotplot()`
from `@primitives/differential_expression.md`. This function filters by
`p_val_adj < PADJ_CUT` AND `|avg_log2FC| > LFC_CUT` before plotting.

Inline reinventions that skip the DE filter (e.g., a `make_geneset_dotplot()` that plots
all genes in the set regardless of significance) are **forbidden**. If a job script
requires gene-set-on-expression visualization that genuinely should not be DE filtered,
the brief must explicitly request it and the deployment agent must approve -- this is not
a default-available capability.

If `make_functional_dotplot()` returns NULL (zero genes pass the filter), log the decision
to `output/decision_log.txt` and do not render an empty plot.

### 8.4 Validation

Before executing any module:
- Confirm each key in `downstream_analyses` appears in the mapping table above
- Confirm the module's required brief keys are present in the config block
- If validation fails for a module, log the error to `output/<module_name>/error.log`
  and skip that module (continue with others)

### 8.5 Input object

By default, Stage 8 modules operate on `output/subset/subset_annotated.rds` (the
annotated subset). If no subset was produced (Stage 6 was skipped), modules operate
on `output/annotation/annotated.rds`. A module may specify `input_object` in its
config block to override the default input path.

### 8.6 Execution order

Execute in this dependency-aware order:

1. `multi_group_de` — first, because `bulk_concordance` with annotation_overlay may
   read DE outputs from this module
2. `bulk_concordance` — second, if `mode: signature_score` with annotation_overlay;
   otherwise it can run in the independent group below
3. All remaining independent modules (no inter-module dependencies):
   `de_comprehensive_csv`, `feature_umap_plot`, `subclustering`, `metabolic_profile`,
   `trajectory`, `pyscenic`, `cellchat`, `deg`

### 8.7 Error handling

If a module raises an exception, save the full error trace to
`output/<module_name>/error.log`. Continue with the next module. Report all module
outcomes in Stage 9.

---

## Stage 9 — Pipeline complete

**Depends on:** Stage 8 complete (or the last completed stage if 8 was skipped)

Print a summary of all outputs generated: each output file, its location, and the stage
that produced it. Report any Stage 8 modules that logged errors, with pointers to their
error logs. Offer to document new analysis patterns discovered during this project in
the library modules/ directory.
