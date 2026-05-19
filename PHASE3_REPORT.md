---
# Phase 3 Report — Pipeline Orchestration Layer
# Written: 2026-05-19
# Author: Phase 3 migration agent
---

# Phase 3 Report

Covers all five Phase 3 deliverables: two pipeline.md files, two brief_template.txt
files, and the Known-Atlas Convention registry. Includes brief schema validation,
self-check results, cross-file findings, and open items for Phase 4/6.

---

## 1. Files Created

| File | Lines | Notes |
|---|---|---|
| `pipelines/large_dataset/pipeline.md` | 295 | 9-stage orchestration, 5 CPs, Stage 8 dispatch table |
| `pipelines/large_dataset/brief_template.txt` | 164 | v2 schema; all HumanFat defaults removed |
| `pipelines/integrate_public_data/pipeline.md` | 372 | 8-stage orchestration, 6 CPs, new Stage 8 with atlas convention |
| `pipelines/integrate_public_data/brief_template.txt` | 225 | New file; first formal schema for this pipeline |
| `context/known_atlases.md` | 177 | Tabula Sapiens entry; Phase 6 placeholder entries |
| **Total** | **1,233** | |

---

## 2. LargeDataset Pipeline

### 2.1 Stage-by-stage summary

| Stage | What it does | Modules dispatched | Checkpoint |
|---|---|---|---|
| 1 — Assemble | Load .h5/rds/mex, parse metadata, merge, compute percent.mt, save raw RDS | @modules/load_formats.md | — |
| 2 — QC | Violin plots, propose 2nd/98th percentile thresholds, write thresholds_proposed.yaml | — | **CP1** user approves/modifies thresholds |
| 3 — Process | NormalizeData → FindVariableFeatures → JoinLayers → ScaleData → PCA → Harmony → Cluster → UMAP | @primitives/harmony_integration.md | **CP2** elbow + UMAPs |
| 4 — Annotation | FindAllMarkers + is_ambient() filter, write markers.yaml | @primitives/differential_expression.md | **CP3** user provides cluster→label mapping |
| 5 — Visualization | make_labeled_umap, make_stacked_violin, make_proportion_plot | @primitives/visualization.md | **CP4** style approval |
| 6 — Subset | Subset target_celltypes, re-process (Harmony), AddModuleScore, FindAllMarkers | @primitives/harmony_integration.md | **CP5** user provides subtype mapping |
| 7 — Core subset | make_canonical_dotplot, make_tf_diamond_plot, make_stacked_violin, make_diff_abundance_heatmap | @primitives/visualization.md, @primitives/differential_expression.md | — |
| 8 — Downstream | Dispatch to downstream_analyses modules in dependency order | Per-module @modules/ files | — |
| 9 — Complete | Print output summary, report errors | — | — |

### 2.2 Brief schema fields by stage

| Stage | Brief key | Required | Default |
|---|---|---|---|
| Identity | `project_name`, `organism`, `tissue`, `output_dir` | Yes | — |
| Stage 1 | `inputs.type`, `inputs.paths`, `inputs.metadata_csv` | Yes (paths required) | metadata_csv: null |
| Stage 2 | `pipeline_params.qc_thresholds.*` | No | 2nd/98th pct, percent_mt: 20 |
| Stage 3 | `pipeline_params.n_variable_features`, `n_pcs`, `clustering_resolution`, `batch_correction_var` | Yes | 4000, 30, 0.8, sample_id |
| Stage 4 | `annotation.top_n_markers` | No | 25 |
| Stage 6 | `subset.target_celltypes`, `subset.subset_name`, `subset.pipeline_params` | No (skipped if null) | subset_name: "subset" |
| Stage 7 | `context_overrides.palettes.subtype_colors`, `metadata_columns.subtype_col`, `metadata_columns.group_col` | No | lab defaults |
| Stage 8 | `downstream_analyses.*` | No | {} |

### 2.3 Changes from v1

| v1 content | v2 change | Reason |
|---|---|---|
| `batch_correction_var: "source_file"` shown as example | Replaced with `"sample_id"` (lab default per lab_context.md) | source_file is HumanFat-specific |
| Adipose tissue defaults throughout | No tissue defaults | Tissue is always project_specific |
| EC-specific palette names (ec_colors, tissue_colors, adipose_type_colors) in stage descriptions | Replaced with generic context_overrides.palettes.group_colors | HumanFat-specific |
| Inline R code in stage descriptions | All code removed; @references point to primitives/modules | Orchestration vs. implementation separation |
| Stage 8 "optional extras" listed as ad-hoc tasks | Stage 8 dispatch table with exact brief keys and module mapping | Explicit routing |
| No module dispatch table | Table mapping brief key → module file | Required for deployment agent routing |
| JoinLayers timing not specified | Explicit: after merge, before ScaleData (Rule 1) | Seurat v5 compliance |

---

## 3. IntegratePublicData Pipeline

### 3.1 Stage-by-stage summary

| Stage | What it does | Modules dispatched | Checkpoint |
|---|---|---|---|
| 1 — Load | Load per-dataset objects (rds/h5/mex/geo), assign standard metadata columns | @modules/load_formats.md, @primitives/geo_download.md | **CP1** load summary review |
| 2 — Filter | AddModuleScore or label-transfer filtering; skip if filtering: null | — | **CP2** before/after cell counts |
| 3 — Harmonize | YAML mapping (3A) + Seurat label transfer (3B) → unified_label | @modules/label_harmonization.md | **CP3** user resolves REVIEW entries |
| 4 — Downsample | Stratified downsample by (unified_label × sample_id or dataset); skip if ceiling: null | — | **CP4** post-downsample counts |
| 5 — Merge | Common gene intersection, merge(), JoinLayers | @primitives/seurat_v5_rules.md | **CP5** merge summary |
| 6 — Integrate | NormalizeData → FindVariableFeatures → ScaleData → PCA → Harmony → Cluster → UMAP | @primitives/harmony_integration.md | **CP6** integration UMAPs |
| 7 — Cohort plots | Correlation heatmap (A), chord diagram (B), proportion heatmap (C); skip if cohort_plots: false | @modules/cohort_plots.md | — |
| 8 — Downstream | Known-Atlas Convention expansion, then dispatch to downstream_analyses modules | Per-module @modules/ files | — |

### 3.2 Stage 8 specification

**Trigger:** `brief.downstream_analyses` is present and non-empty.

**Known-Atlas Convention:**
1. If `brief.atlas` is set, look up the name in @context/known_atlases.md
2. On match, inject registry values into `atlas_co_umap` and `cross_dataset_dotplot`
   config blocks as defaults: `source_col`, `atlas_label`, `atlas_group_col`,
   `n_per_subtype` / `recommended_downsample_n`
3. Explicit values in the module config block take precedence over registry defaults

**Validation:** For each module key:
- Confirm it appears in the brief key → module file mapping table
- Confirm required brief keys are present
- For atlas-aware modules: if `brief.atlas` is unrecognized, require explicit params
- On failure: write `output/<module_name>/error.log`, skip module, continue

**Execution order:**
1. `multi_group_de` — first (other modules may annotate against DE results)
2. `bulk_concordance` — second if annotation_overlay requires DE outputs
3. All others in parallel (no inter-module dependencies)

**Error handling:** No hard stop on module failure. Each failed module writes its
error trace to `output/<module_name>/error.log`. Subsequent modules continue.

**No Stage 8 checkpoint:** Stage 8 runs silently to completion. Rationale: the cleanup
agent's biological interpretation at job completion is the review gate. Multi-pass
brief strategy: run `multi_group_de` first, review DE output, then resubmit with
dependent modules.

### 3.3 Brief schema fields by stage

| Stage | Brief key | Required | Default |
|---|---|---|---|
| Identity | `project_name`, `organism`, `tissue`, `output_dir` | Yes | — |
| Stage 1 | `inputs.datasets[*].{name, type, path}` | Yes | — |
| Stage 2 | `filtering.*` | No | null (skip stage) |
| Stage 3 | `label_harmonization.{label_col, broad_label, known_mappings}` | Yes | score_high: 0.75, score_low: 0.40 |
| Stage 4 | `downsampling.{ceiling, strategy}` | No | ceiling: null (skip) |
| Stage 6 | `integration.*` | Yes | n_variable_features: 2250, n_pcs: 30, resolution: 0.3, batch_var: sample_id |
| Stage 7 | `cohort_plots` | No | false |
| Stage 8 | `atlas`, `downstream_analyses.*` | No | atlas: null, analyses: {} |

### 3.4 Changes from v1

| v1 state | v2 change |
|---|---|
| Pipeline ends at CP6 with optional Stage 7 cohort plots | New Stage 8 dispatches all downstream analyses |
| No brief_template.txt | New file created (225 lines) |
| No formal brief schema | Full YAML schema with defaults documented |
| No atlas shorthand | `atlas: "Tabula Sapiens"` → Known-Atlas Convention |
| No dispatch routing | Explicit brief key → module file mapping table |
| Stage 1 loads with no metadata column conventions | Standard columns (dataset, unified_label, sample_id, condition, organ, data_type) codified in pipeline header |
| No pre-normalized data handling | Explicit note in Stage 6 to handle pre-normalized datasets before NormalizeData |
| atlas_co_umap and cross_dataset_dotplot as orphan files | Wired into Stage 8 dispatch with Known-Atlas Convention |

---

## 4. Known-Atlas Convention

### 4.1 Registry contents

**Registered atlas:** Tabula Sapiens

| Parameter | Value |
|---|---|
| `display_names` | "Tabula Sapiens", "tabula sapiens", "TS", "Tabula_Sapiens", "TabulaSapiens" |
| `citation` | Tabula Sapiens Consortium, Science 2022 |
| `source_col` | "dataset" (standard IntegratePublicData column) |
| `atlas_label` | "TS" |
| `atlas_group_col` | "cell_ontology_class" |
| `recommended_downsample_n` | 1500 (validated TabulaSapiensComparison project) |
| `validated_date` | 2026-04-15 |

**Placeholder entries (not validated):** Human Cell Atlas (HCA), Mouse Cell Atlas (MCA).

### 4.2 Cross-references in modules

| Module | Fields consumed from registry |
|---|---|
| @modules/atlas_co_umap.md | `source_col`, `atlas_label`, `atlas_group_col`, `n_per_subtype` |
| @modules/cross_dataset_dotplot.md | `source_col`, `atlas_group_col` |

### 4.3 What the deployment agent needs in Phase 6

1. Parse `brief.atlas` string against all `display_names` lists (case-insensitive match)
2. On match: read the registry entry, inject default values into `atlas_co_umap` and
   `cross_dataset_dotplot` config blocks
3. Log which values came from registry vs. explicit brief for traceability
4. On no match: emit a warning and require explicit atlas parameters in module config blocks
5. Registry should remain static YAML data; expansion logic belongs in agent code

---

## 5. Brief Schema Validation

For each module dispatchable from Stage 8, confirmation that its required brief keys
are satisfiable from the pipeline brief template + `downstream_analyses` block:

| Module (file) | Brief key | Required fields | All satisfiable? | Notes |
|---|---|---|---|---|
| multi_group_de_analysis.md | `multi_group_de` | group_col, label_col, groups, comparisons | ✓ | `functional_gene_sets` labeled REQUIRED but module code gates functional plots conditionally; effectively optional for base DE |
| atlas_co_umap.md | `atlas_co_umap` | source_col, inhouse_label, atlas_label, subtype_col, atlas_group_col, highlight_inhouse, highlight_atlas, coords_csv | ✓ | 4 of 8 fields auto-populated by Known-Atlas Convention when brief.atlas is set; 4 remain project-specific |
| cross_dataset_dotplot.md | `cross_dataset_dotplot` | source_col, subtype_col, sources_order, marker_genes, col_group_col | ✓ | source_col and atlas_group_col auto-populated by convention; sources_order, marker_genes, col_group_col always project-specific |
| bulk_concordance.md | `bulk_concordance` | mode, bulk_csv, experiment_label, subtype_col | ✓ | mode values: signature_score \| parallel_lfc |
| de_comprehensive_csv.md | `de_comprehensive_csv` | group_col | ✓ | label_col optional |
| cellchat.md | `cellchat` | label_col | ✓ | group_col optional; null skips Script 6 |
| feature_umap_plot.md | `feature_umap_plot` | output_dir (top-level), genes (top-level) | ⚠️ | **Mismatch:** module uses `brief$genes` and `brief$output_dir` at the top level, not under `downstream_analyses.feature_umap_plot`. Stage 8 dispatch must map `downstream_analyses.feature_umap_plot.genes` → `brief$genes`. Documented as P3-1. |
| celltype_subclustering.md | `subclustering` | label_col (from metadata), batch_correction_var | ✓ | label_col comes from context_overrides.metadata_columns, not the subclustering block; Stage 8 must pass it from context |
| metabolic_profile.md | `metabolic_profile` | gene_sets, label_col, group_col | ✓ | gene_sets is biology-specific and must be in the module config block |
| trajectory_monocle3.md | `trajectory` | label_col, start_cluster | ✓ | end_cluster is informational; Monocle3 infers endpoints |
| pyscenic_regulons.md | `pyscenic` | label_col, scenic_python_path, database_dir, tf_list | ✓ | All machine-specific paths; no registry defaults |

**One mismatch identified (P3-1):** `feature_umap_plot.md` uses `brief$genes` at the
top level. All other modules use `brief$downstream_analyses.<module_key>.*` for their
parameters. The Stage 8 dispatch implementation must handle this special case by lifting
`downstream_analyses.feature_umap_plot.genes` to the top-level `genes` key before
executing the module's script template.

---

## 6. Self-Check Grep Results

Terms checked across all five Phase 3 files:

| Term | Hits in pipelines/ or context/known_atlases.md | Classification |
|---|---|---|
| HumanFat | 0 | Clean |
| ec_colors, tissue_colors, adipose_type_colors | 0 | Clean |
| mylabel, STROMA_ORDER | 0 | Clean |
| Vertebrae, Iliac Crest, Femoral Head | 0 | Clean |
| Myelolipoma | 0 | Clean |
| CapEC, VenEC, AEC, Fat EC, Kidney EC | 0 | Clean |
| PPARG | 0 | Clean |
| source_file | 0 | Clean |
| NKXSpleen | 0 | Clean |
| ec_subtype_colors | 0 | Clean |
| tissue_type | 0 | Clean |
| patient_id | 0 | Clean |
| EC-GC, EC-AEA | 0 | Clean |
| LiposuctionFat, BreastFat | 0 | Clean |
| KidneyNew | 0 | Clean |
| BoneMarrowStroma | 0 | Clean |
| TabulaSapiensComparison | 1 (known_atlases.md:85) | **Expected** — notes block documenting validation context; atlas registry explicitly names source projects |
| SEC1_FORCE, EC_SUBTYPES | 0 | Clean |
| "Source Sans 3", "Playfair Display" | 1 each (known_atlases.md:94) | **Expected** — font-restriction warning in Tabula Sapiens notes block; explicitly listed in task brief as acceptable hits in known_atlases.md |

All three hits are in `context/known_atlases.md` in the `notes` documentation block for
the Tabula Sapiens entry. All are expected per task brief guidance ("Hits in
known_atlases.md mentioning Tabula Sapiens are EXPECTED and correct"). No hits in
pipeline.md or brief_template.txt files.

---

## 7. Cross-File Findings

### 7.1 Brief key names differ from module file names

Four modules use abbreviated brief keys that do not match their file names:

| Module file | Brief key | Confirmed from |
|---|---|---|
| multi_group_de_analysis.md | `multi_group_de` | module frontmatter line 17 |
| celltype_subclustering.md | `subclustering` | module frontmatter line 18 |
| trajectory_monocle3.md | `trajectory` | module frontmatter line 13 |
| pyscenic_regulons.md | `pyscenic` | module frontmatter line 13 |

Both pipeline.md dispatch tables and both brief_template.txt files use the correct
abbreviated keys matching the module frontmatter. The task brief's proposed template
would have used longer names (`multi_group_de_analysis`, etc.); the actual module
schemas take precedence.

### 7.2 bulk_concordance mode values corrected

The task brief proposed `mode: 1` and `mode: 2`. The actual module schema uses
`mode: "signature_score"` and `mode: "parallel_lfc"`. Both brief templates use
the string values.

### 7.3 feature_umap_plot brief key structure is non-standard

All other modules store their parameters under `downstream_analyses.<key>.*`.
`feature_umap_plot.md` uses `brief$genes` and `brief$output_dir` at the top level.
This is a Phase 2A design decision that predates the Stage 8 dispatch layer. Both
brief templates document the module config block under `downstream_analyses.feature_umap_plot`
for consistency, with a note that Stage 8 dispatch must lift `genes` to the top level.

### 7.4 Cohort plots module requires dataset_colors

`cohort_plots.md` requires `dataset_colors` and `ct_colors` from context. The
`integrate_public_data/brief_template.txt` includes `context_overrides.palettes.dataset_colors`
to address this. The `large_dataset` pipeline does not dispatch cohort_plots (it is
an IntegratePublicData-specific module), so no corresponding key is needed there.

### 7.5 JoinLayers timing is different between the two pipelines

- **large_dataset:** `JoinLayers()` is called in Stage 3 (before ScaleData), immediately
  after loading the Stage 2 filtered RDS. This is correct per Rule 1 for multi-sample
  objects loaded in Stage 1.
- **integrate_public_data:** `JoinLayers()` is called in Stage 5 (Merge), immediately
  after `merge()`. Stage 6 (Integration) starts at `NormalizeData()` and skips
  `JoinLayers()` (already joined). Both pipelines correctly follow Rule 1.

---

## 8. Open Items for Phase 4 and Phase 6

### P3-1. feature_umap_plot top-level brief key (Phase 6)

When the deployment agent implements Stage 8 dispatch, it must handle the special case
where `feature_umap_plot.md` reads `brief$genes` and `brief$output_dir` at the top
level. The dispatch layer should copy `downstream_analyses.feature_umap_plot.genes`
to the top-level `genes` key before running the module's script template. The module
itself needs no modification.

### P3-2. SKILL.md routing table needs updating (Phase 4)

`SKILL.md` still references v1 module names and paths from the original routing table
(e.g., `modules/interactome_cellchat.md` (Phase 2)). The routing table in SKILL.md
should be updated to reference the v2 module files and the new pipeline brief templates.
This is a library maintenance task for Phase 4.

### P3-3. Atlas registry — validated_date for Tabula Sapiens (Phase 4)

The `validated_date: "2026-04-15"` in `known_atlases.md` is derived from the Phase
report notes. If the actual TabulaSapiensComparison run used a different date, update
the registry entry. Also: the `typical_cell_count` of 483152 is an approximation —
confirm against the actual TS download used in that project.

### P3-4. HCA and MCA placeholder entries (Phase 6)

`Human Cell Atlas` and `Mouse Cell Atlas` appear as commented-out placeholder entries
in `known_atlases.md`. These should be validated and populated when a project using
one of these atlases runs atlas_co_umap or cross_dataset_dotplot for the first time.

### P3-5. cross_dataset_dotplot inhouse_label not in Known-Atlas Convention (Phase 6)

The `inhouse_label` field (source_col value for in-house data) cannot be auto-populated
from the atlas registry because it depends on the project's dataset naming convention.
The deployment agent must always prompt for this value when `brief.atlas` is set and
`cross_dataset_dotplot` is in `downstream_analyses`. Document this in the deployment
agent's known-atlas expansion logic.

### P3-6. SKILL.md pipeline entries still say "Phase 2 (not yet written)" (Phase 4)

SKILL.md's Available Pipelines table and the pipelines/ directory entry both indicate
Phase 2 as the target. Now that Phase 3 has created the pipeline files, SKILL.md should
be updated to reflect that pipelines/ is complete and point to the new files. This
is a library maintenance update for Phase 4 alongside P3-2.

### P3-7. multi_group_de functional_gene_sets labeled REQUIRED but conditional (Phase 4)

In `multi_group_de_analysis.md` Brief Schema, `functional_gene_sets` is labeled
REQUIRED but the module code conditionally executes functional plots only when it is
non-null. The module schema should be updated to label this as optional (with a note
that functional dotplots are skipped if null). This is a Phase 4 documentation cleanup
for the module itself — no change needed in Phase 3 files.

---

## 9. Anything Unexpected

### U1 — IntegratePublicData pipeline had no brief_template.txt at all in v1

The entire IntegratePublicData pipeline relied on the shared CONVENTIONS.md brief
schema plus informal project CLAUDE.md files. Creating a formal schema revealed that
the pipeline has substantially more configuration surface than the LargeDataset pipeline
— 8 stages vs 9, but Stage 3 (label harmonization) alone requires a complex multi-key
schema (label_col dict, reference_datasets list, score thresholds, known_mappings dict).
The new brief_template.txt is 225 lines vs 164 for LargeDataset, reflecting this
complexity.

### U2 — atlas_co_umap requires coords_csv (not documented in task brief)

The task brief's proposed brief schema for `atlas_co_umap` did not include `coords_csv`
(path to save/load UMAP coordinate CSV). Reading the actual module frontmatter revealed
`coords_csv` is REQUIRED. This was added to the integrate_public_data brief_template.txt.
The large_dataset pipeline does not dispatch atlas_co_umap (it is available only in
IntegratePublicData Stage 8), so no change needed there.

### U3 — cross_dataset_dotplot requires sources_order and marker_genes

The task brief's proposed brief schema did not include these fields. Reading the module
schema revealed both are REQUIRED. These cannot be auto-populated from the Known-Atlas
Convention — they are biology-specific inputs. Both brief_template.txt files now
document these as required fields in the commented `cross_dataset_dotplot` block.

### U4 — Two pipelines share the same Stage 8 module set except atlas-aware modules

LargeDataset Stage 8 offers: `multi_group_de`, `feature_umap_plot`, `subclustering`,
`metabolic_profile`, `trajectory`, `pyscenic`, `bulk_concordance`, `de_comprehensive_csv`,
`cellchat`. IntegratePublicData Stage 8 adds: `atlas_co_umap`, `cross_dataset_dotplot`.
The atlas-aware modules (`atlas_co_umap`, `cross_dataset_dotplot`) are not in
LargeDataset because the pipeline does not produce a multi-dataset `dataset` column,
making atlas comparison semantically undefined. This was a deliberate scoping decision.
