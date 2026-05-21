---
# Integrate Public Data Pipeline — Multi-source heterogeneous dataset integration
# Pipeline: integrate_public_data
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
  when_stage_1_contains_geo:
    - "@primitives/geo_download.md"
    - "@primitives/file_downloads.md"
  when_downstream_de:
    - "@primitives/differential_expression.md"
references:
  stage_1_load:
    - "@modules/load_formats.md"
  stage_3_harmonize:
    - "@modules/label_harmonization.md"
  stage_7_cohort:
    - "@modules/cohort_plots.md"
  stage_8_modules:
    - "@modules/multi_group_de_analysis.md"
    - "@modules/atlas_co_umap.md"
    - "@modules/cross_dataset_dotplot.md"
    - "@modules/bulk_concordance.md"
    - "@modules/de_comprehensive_csv.md"
    - "@modules/cellchat.md"
    - "@modules/feature_umap_plot.md"
    - "@modules/celltype_subclustering.md"
    - "@modules/metabolic_profile.md"
    - "@modules/trajectory_monocle3.md"
    - "@modules/pyscenic_regulons.md"
---

# Integrate Public Data Pipeline

Integration of heterogeneous scRNAseq datasets from multiple sources — in-house RDS,
public GEO downloads, or other formats — into a single harmonized Seurat object. Covers
loading through Harmony integration, optional cohort characterization plots, and
Stage 8 post-integration downstream analyses. Expect 8 stages and up to 6 checkpoints.

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

### Standard metadata columns produced by this pipeline

After Stage 3, all cells in the merged object carry these columns:

| Column | Content |
|---|---|
| `dataset` | Source dataset name (= `name` from `inputs.datasets` entry) |
| `unified_label` | Harmonized cell type label (produced by Stage 3) |
| `sample_id` | Sample identifier (from source metadata or assigned in Stage 1) |
| `condition` | Experimental condition if available from source; NA otherwise |
| `organ` | Anatomical origin if available from source; NA otherwise |
| `data_type` | Sequencing modality (e.g., 10x_v3, 10x_v2) if available; NA otherwise |

Downstream modules reference these columns by default. Override with
`context_overrides.metadata_columns` if a project uses different column names.

---

## Stage 1 — Load

**Depends on:** Nothing (pipeline start)
**Brief keys consumed:** `inputs.datasets`
**References:** @modules/load_formats.md; @primitives/geo_download.md (for `type: "geo"`)
**Checkpoint:** CP1 — User reviews loaded object summaries and confirms metadata resolved

For each entry in `brief.inputs.datasets`, dispatch to the format-appropriate recipe
in @modules/load_formats.md based on `type`:
- `rds`: direct load; verify layer structure; apply @primitives/seurat_v5_rules.md Rule 5
  for dynamic UMAP reduction name detection
- `h5`: `Read10X_h5()` → `CreateSeuratObject()`; check `is.list(counts)` for
  multi-modal h5 files per load_formats guidance
- `mex`: `Read10X()` from the directory
- `geo`: follow @primitives/geo_download.md acquisition priority rules (prefer
  processed count matrix over raw reads; check normalization status before assuming
  raw counts; download h5 or MEX format)

Assign the standard `dataset` column to all cells in each loaded object (value =
`name` from the brief entry). Assign `organ`, `data_type`, and `condition` if these
are available in the source metadata; leave as NA otherwise.

Save per-dataset Seurat objects to `output/loaded/<dataset_name>.rds`. Write a
summary table to `output/loaded/summary.yaml`: dataset name, n_cells, n_genes,
normalization_status (raw_counts | log_normalized | unknown), metadata_columns_found.
Flag any pre-normalized datasets for handling in Stage 6.

**CP1:** Present the summary table. The user confirms all datasets loaded correctly,
cell counts match expectations, and the metadata columns resolved as intended.

---

## Stage 2 — Filter target cells

**Depends on:** Stage 1 complete, CP1 approved
**Brief keys consumed:** `filtering.target_markers`, `filtering.threshold_method`,
  `filtering.reference_object`
**Checkpoint:** CP2 — User reviews cell counts before/after filtering

If `brief.filtering` is null, skip Stage 2 (proceed directly to Stage 3).

**Module-score filtering** (`threshold_method: "module_score"`): For each dataset
listed under `filtering.target_markers.positive`, compute `AddModuleScore()` using
the provided positive marker genes. Auto-suggest a filter threshold at the top quartile
of the module score distribution. Retain cells above the threshold.

For any dataset with `filtering.target_markers.negative` entries, filter cells that
express any negative marker gene above background (mean + 2 SD of that gene's
expression in the dataset).

**Label-transfer filtering** (`threshold_method: "label_transfer"`): If
`filtering.reference_object` is provided (path to a well-annotated reference RDS),
run Seurat label transfer to identify the target cell population. Retain cells with
prediction score ≥ 0.5 for the target label.

Save filtered per-dataset objects to `output/filtered/<dataset_name>.rds`.

**CP2:** Present cell counts before/after filtering per dataset. The user confirms
the filter retained the expected target population size and no major cell type was
inadvertently removed.

---

## Stage 3 — Label harmonization

**Depends on:** Stage 2 complete (or Stage 1 if Stage 2 skipped)
**Brief keys consumed:** `label_harmonization.*`
**Reference:** @modules/label_harmonization.md
**Checkpoint:** CP3 — User reviews harmonization YAML and resolves REVIEW-flagged entries

Dispatch to @modules/label_harmonization.md for both sub-stages:

**Stage 3A — YAML-based mapping:** For each labeled dataset, look up the source
label column from `label_harmonization.label_col` (a dict: dataset_name → column).
Map raw labels to `unified_label` using `label_harmonization.known_mappings`
(source_label → unified_label). Labels not present in `known_mappings` are
flagged as `REVIEW` in the output YAML. Assign `label_harmonization.broad_label`
to cells with ambiguous or low-confidence assignments.

**Stage 3B — Label transfer:** For unlabeled datasets listed in
`label_harmonization.reference_datasets`, run Seurat label transfer using a reference
built from the already-labeled datasets. Apply `score_high` (default 0.75) and
`score_low` (default 0.40) thresholds. Cells below `score_low` receive `broad_label`.

Write the full harmonization mapping — including REVIEW-flagged entries — to
`output/harmonization/harmonization.yaml`.

**CP3:** Present `harmonization.yaml`. The user resolves all REVIEW entries by
providing the correct `unified_label` for each flagged source label. After approval,
apply the finalized `unified_label` column to all dataset objects.

**Invariant:** `unified_label` must be non-null for all cells after CP3. Any remaining
NA values indicate a gap in the mapping that must be resolved before Stage 4.

---

## Stage 4 — Downsampling

**Depends on:** Stage 3 complete, CP3 approved
**Brief keys consumed:** `downsampling.ceiling`, `downsampling.strategy`
**Checkpoint:** CP4 — User reviews cell counts after downsampling

If `brief.downsampling.ceiling` is null, skip Stage 4 (proceed to Stage 5).

Perform stratified downsampling to at most `downsampling.ceiling` cells per stratum:
- `strategy: "by_sample"`: downsample each (unified_label × sample_id) combination
  independently to at most `ceiling` cells
- `strategy: "by_dataset"`: downsample each (unified_label × dataset) combination
  independently to at most `ceiling` cells

Save downsampled per-dataset objects. Write a summary to
`output/downsampled/summary.yaml`: n_cells before/after per dataset and per
unified_label category.

**CP4:** Present the downsampling summary. The user confirms cell counts are
acceptable. If any cell type is severely underrepresented after downsampling, the
user may adjust `ceiling` before proceeding.

---

## Stage 5 — Merge

**Depends on:** Stage 4 complete (or Stage 3 if downsampling was skipped)
**Reference:** @primitives/seurat_v5_rules.md (Rule 1: JoinLayers after merge)
**Checkpoint:** CP5 — User reviews merged object summary

Find the common gene set across all datasets (intersection of rownames). Subset each
object to the common genes before merging to ensure matrix dimensions are compatible.

Merge all filtered/downsampled objects into one Seurat object with `merge()`.
Call `JoinLayers()` immediately after merge per @primitives/seurat_v5_rules.md Rule 1,
preparing the object for `ScaleData()` in Stage 6.

Propagate all standard metadata columns: `dataset`, `unified_label`, `sample_id`,
`condition`, `organ`, `data_type`.

Save merged object to `output/merged/merged.rds`. Write summary: total n_cells,
n_genes in common gene set, n_cells per dataset, n_cells per unified_label.

**CP5:** Present the merged object summary. The user confirms the merge is correct
(expected cell counts, metadata intact). Flag if the common gene set is substantially
smaller than expected — this often indicates protocol or species mismatches between
datasets.

---

## Stage 6 — Integration

**Depends on:** Stage 5 complete, CP5 approved, `output/merged/merged.rds` exists
**Brief keys consumed:** `integration.batch_correction_var`, `integration.n_variable_features`,
  `integration.n_pcs`, `integration.clustering_resolution`
**Reference:** @primitives/harmony_integration.md
**Checkpoint:** CP6 — User reviews integration UMAP

Load `output/merged/merged.rds`. Execute the canonical Seurat+Harmony processing
sequence from @primitives/harmony_integration.md using `integration.*` parameters.
Note: for IntegratePublicData, `JoinLayers()` was already called in Stage 5, so the
sequence begins at `NormalizeData()`:

1. `NormalizeData()`
2. `FindVariableFeatures(nfeatures = integration.n_variable_features)` (default 2250)
3. `ScaleData(features = VariableFeatures(so))` (variable features only, Rule 2)
4. `RunPCA(npcs = integration.n_pcs)` (default 30)
5. `RunHarmony(group.by.vars = integration.batch_correction_var)` → reduction `harmony`
6. `FindNeighbors(reduction = "harmony", dims = 1:n_pcs)`
7. `FindClusters(resolution = integration.clustering_resolution)` (default 0.3)
8. `RunUMAP(reduction = "harmony", dims = 1:n_pcs)` → stored as `harmony.umap`

Note: if any source dataset was pre-normalized in Stage 1, handle normalization
compatibility before this step (do not re-normalize pre-normalized data; merge
on the raw counts layer).

Save integrated object to `output/integrated/integrated.rds`.

**CP6:** Generate three diagnostic UMAPs: (A) colored by Seurat cluster, (B) colored
by `dataset`, (C) colored by `unified_label`. Present for user review. Confirm that
Harmony removed batch structure between datasets while preserving biological variation
(cell type clusters should be mixed across datasets). The user may request re-integration
with adjusted parameters before proceeding.

---

## Stage 7 — Cohort plots

**Depends on:** Stage 6 complete, CP6 approved
**Brief keys consumed:** `cohort_plots`
**Reference:** @modules/cohort_plots.md

If `brief.cohort_plots` is false or absent, skip Stage 7.

Dispatch to @modules/cohort_plots.md for all three plot types:

- **Plot A — Correlation heatmap:** HVG z-scored expression averages per cell type
  and dataset, displayed as a Pearson correlation heatmap (ComplexHeatmap). Uses the
  `pdf()` + `ComplexHeatmap::draw()` exception (documented in CONVENTIONS.md §4).
- **Plot B — Chord diagram:** Dataset × cell type composition using circlize
  `chordDiagramFromMatrix()`. Uses the `pdf()` + `dev.off()` exception for circlize.
- **Plot C — Proportion heatmap:** Cell type percentages per sample with annotation
  strips for `organ`, `condition`, and `data_type` (strips shown only for columns
  where the value is non-NA).

The module requires `dataset_colors` and `ct_colors` from `context_overrides.palettes`
(or defaults from @context/color_palettes.md), plus `dataset_col`, `label_col`, and
`sample_col` from `context_overrides.metadata_columns`.

Save all outputs to `output/cohort_plots/`. No checkpoint follows Stage 7; cohort
plots are informational and do not gate Stage 8.

---

## Stage 8 — Post-integration downstream analyses

**Depends on:** Stage 7 complete (or Stage 6 if cohort_plots was skipped)
**Brief keys consumed:** `downstream_analyses.*`, `atlas` (if set)
**References:** All modules listed in `references.stage_8_modules`

If `brief.downstream_analyses` is absent or empty (`{}`), skip Stage 8 and print
the pipeline completion summary.

### 8.1 Known-Atlas Convention

If `brief.atlas` is set to a recognized atlas name string, look up the atlas in
@context/known_atlases.md before dispatching any atlas-aware modules. This allows the
brief to use a short name (e.g., `atlas: "Tabula Sapiens"`) in place of specifying
full atlas parameters in every module block.

The registry in @context/known_atlases.md provides per-atlas values for:
- `source_col` — the metadata column identifying dataset source (standard: `dataset`)
- `atlas_label` — short atlas identifier for plot titles (e.g., `"TS"`)
- `atlas_group_col` — atlas cell type column (e.g., `cell_ontology_class` for TS)
- `recommended_downsample_n` — validated DOWNSAMPLE_N for this atlas

These values are injected into the `atlas_co_umap` and `cross_dataset_dotplot`
config blocks as defaults before execution. Explicit values in the brief's module
config blocks take precedence over registry defaults.

The following `atlas_co_umap` and `cross_dataset_dotplot` parameters remain
project-specific even when `brief.atlas` is set; they must be specified in the
module config block:
- `atlas_co_umap`: `inhouse_label`, `subtype_col`, `highlight_inhouse`,
  `highlight_atlas`, `coords_csv`, and optionally `subtype_colors`
- `cross_dataset_dotplot`: `sources_order`, `marker_genes`, `col_group_col`,
  `inhouse_label` (as the source_col value for in-house data), and optionally
  `sec1_force`

### 8.2 Brief key → module file mapping

| Brief key | Module file |
|---|---|
| `multi_group_de` | @modules/multi_group_de_analysis.md |
| `atlas_co_umap` | @modules/atlas_co_umap.md |
| `cross_dataset_dotplot` | @modules/cross_dataset_dotplot.md |
| `bulk_concordance` | @modules/bulk_concordance.md |
| `de_comprehensive_csv` | @modules/de_comprehensive_csv.md |
| `cellchat` | @modules/cellchat.md |
| `feature_umap_plot` | @modules/feature_umap_plot.md |
| `subclustering` | @modules/celltype_subclustering.md |
| `metabolic_profile` | @modules/metabolic_profile.md |
| `trajectory` | @modules/trajectory_monocle3.md |
| `pyscenic` | @modules/pyscenic_regulons.md |

### 8.3 Validation

Before executing any module:
- Confirm each key in `downstream_analyses` appears in the mapping table above
- Confirm the module's required brief keys are present in the config block
- For atlas-aware modules (`atlas_co_umap`, `cross_dataset_dotplot`): if `brief.atlas`
  is set, confirm the name resolves in @context/known_atlases.md; otherwise require
  explicit `atlas_group_col`, `source_col`, and `atlas_label` in the config block
- If validation fails for a module, log the error to `output/<module_name>/error.log`
  and skip that module (continue with others)

### 8.4 Input object

All Stage 8 modules operate on `output/integrated/integrated.rds` by default. A module
may specify `input_object` in its config block to override the default. The `subclustering`
module produces its own subset object at
`output/celltype_subclustering/<subset_name>_annotated.rds`, which can be used as
`input_object` for any subsequent module targeting that subset.

### 8.5 Execution order

Execute modules in this dependency-aware order:

1. `multi_group_de` — first; `bulk_concordance` with annotation_overlay may read
   DE outputs from this module
2. `bulk_concordance` — second if `mode: signature_score` and requires DE annotation;
   otherwise runs in the independent group below
3. All remaining independent modules (no inter-module dependencies):
   `atlas_co_umap`, `cross_dataset_dotplot`, `de_comprehensive_csv`, `cellchat`,
   `feature_umap_plot`, `subclustering`, `metabolic_profile`, `trajectory`, `pyscenic`

### 8.6 Error handling and completion

If a module raises an exception, save the full error trace to
`output/<module_name>/error.log`. Do not halt the pipeline. Continue with the next
module in the execution order.

Stage 8 runs silently to completion with no intermediate checkpoints. The cleanup
agent's biological interpretation at job completion is the review gate for all Stage 8
outputs. If a project requires reviewing DE results before dispatching downstream
pathway analysis, structure the initial brief with only `multi_group_de` in
`downstream_analyses`, review the output, then resubmit a continuation brief with
`bulk_concordance` and other DE-dependent modules.

Print the pipeline completion summary after Stage 8: all output files, their locations,
and any module error logs.
