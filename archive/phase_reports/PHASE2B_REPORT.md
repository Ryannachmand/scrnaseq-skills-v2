---
# Phase 2 Group B Report — celltype_subclustering, metabolic_profile, trajectory_monocle3,
#   pyscenic_regulons, bulk_concordance, de_comprehensive_csv
# Written: 2026-05-15 (retroactive — covers all six modules)
# Author: Phase 2B resume agent (four modules authored by Phase 2B original agent;
#         two modules authored by Phase 2B resume agent)
---

# Phase 2 Group B Report

Human review document. Covers all six Group B modules: per-module parameterization tables,
context dependency summaries, primitives referenced, deprecated pattern fixes, and Phase 4
staging notes. Also covers punch-list resolution, self-check grep results, pdf()/dev.off()
exceptions, cross-file findings, and open items.

---

## 1. Files Created

| File | v2 Lines | v1 Source(s) | v1 Lines | Delta |
|---|---|---|---|---|
| `modules/celltype_subclustering.md` | 566 | `LargeDataset/methods/celltype_subclustering.md` | 521 | +45 |
| `modules/metabolic_profile.md` | 454 | `LargeDataset/methods/metabolic_profile.md` | 257 | +197 |
| `modules/trajectory_monocle3.md` | 292 | `LargeDataset/methods/trajectory_monocle3.md` | 183 | +109 |
| `modules/pyscenic_regulons.md` | 517 | `LargeDataset/methods/pyscenic_regulon_analysis.md` | 428 | +89 |
| `modules/bulk_concordance.md` | 945 | `LargeDataset/methods/bulk_concordance.md` (759) + `bulk_lfc_concordance_heatmap.md` (583) | 1342 | -397 |
| `modules/de_comprehensive_csv.md` | 359 | `LargeDataset/methods/de_comprehensive_csv.md` | 253 | +106 |
| **Total** | **3,133** | | **2,984** | **+149** |

**Note on bulk_concordance.md:** Two v1 files (1342 combined lines) were merged into a
single module (945 lines). The reduction (-397 lines) reflects elimination of duplicate
framing/documentation that existed in both v1 files and consolidation of shared config patterns.

---

## 2. Per-Module Details

---

### 2.1 modules/celltype_subclustering.md

#### Parameterization Table

| v1 hardcoded value | v2 parameter |
|---|---|
| `HARMONY_BY <- "source_file"` (HumanFat batch column) | `BATCH_CORRECTION_VAR` from brief `pipeline_params.batch_correction_var` |
| `"tissue_type"` metadata column | `GROUP_COL` config variable (from brief `metadata.group_col`) |
| `"source_file"` as SAMPLE_COL | `SAMPLE_COL` config variable (defaults to `BATCH_CORRECTION_VAR`) |
| `adipose_type_colors` (6 adipose depot palette) | Removed; `group_cols` from brief `context_overrides.palettes.group_colors` |
| `c("Macrophage", "Monocyte")` example label | `c("project_specific")` blank template |
| `mylabel` implicit subtype column | `subtype_label` (new column name assigned during Phase 2) |
| `pdf()` + `dev.off()` for Phase 1 diagnostic plots | All plots use `ggsave()` — P2-7 fix |
| Adipose depot tissue color definitions inline | Removed; sentinel `subtype_cols` / `group_cols` |

#### Context Dependency Declaration Summary

- **Palettes (optional):** `group_colors`, `subtype_colors`
- **Metadata columns (required):** `label_col` (for subset from whole object)
- **Metadata columns (optional):** `group_col`, `sample_col`, `batch_col`
- **Brief keys (required):** `output_dir`, `pipeline_params.batch_correction_var`
- **Brief keys (optional):** `downstream_analyses.subclustering.*` (target_celltypes, resolution, n_pcs, keep_clusters, exclude_groups)

#### Primitives Referenced

- `@primitives/harmony_integration.md` — re-clustering recipe (referenced, not inlined)
- `@primitives/visualization.md` — color scale and sizing conventions
- `@primitives/seurat_v5_rules.md` — Rule 1 (JoinLayers before FindAllMarkers), Rule 3 (metadata assignment)
- `@context/color_palettes.md` — Paul Tol muted palette for cluster exploration

#### Deprecated v1 Patterns Fixed

- **P2-7 Fix:** v1 used `pdf()` + `dev.off()` for Phase 1 diagnostic plots (lines 108–110).
  v2 uses `ggsave(..., device = "pdf", useDingbats = FALSE)` for all ggplot outputs.
  The fix is documented inline with a "P2-7 fix:" comment.

#### Project-Specific Values Staged for Phase 4

- `examples/humanfat_mac_subclustering.md`:
  - `CELLTYPE_LABELS` for Macrophage/Monocyte subset
  - `adipose_type_colors` palette
  - `source_file` as `BATCH_CORRECTION_VAR` and `SAMPLE_COL`
  - `tissue_type` as `GROUP_COL`
  - Validated 2026-03-24 cluster → label mapping

---

### 2.2 modules/metabolic_profile.md

#### Parameterization Table

| v1 hardcoded value | v2 parameter |
|---|---|
| OXPHOS gene set with invalid range syntax `"NDUFA1"-"NDUFA10"` | Not carried forward — caller-provided `gene_sets` argument |
| FA_Synthesis gene set with invalid range syntax `"ELOVL1"-"ELOVL7"` | Not carried forward |
| PPARG_Targets gene set (biology-specific) | Not carried forward — Phase 4 examples/ |
| `tissue_type != "Myelolipoma"` exclusion filter | Not carried forward — caller responsibility |
| `ec_subtype_colors` referenced in output code | Not carried forward — `subtype_colors` from brief |
| Script names referencing "HumanFat" | Generic output filenames (`metabolic_scored.Rds`, etc.) |
| `SUBSAMPLE_CEILING = 2500` hardcoded | Config variable with documented default |
| 3-stop blue-yellow-red diverging palette | `PALETTE_DIVERGING_6` from `@context/color_palettes.md` (6-stop) |

#### Context Dependency Declaration Summary

- **Palettes (optional):** `subtype_colors`, `group_colors`
- **Metadata columns (required):** `label_col`
- **Metadata columns (optional):** `subtype_col`, `group_col`
- **Brief keys (required):** `output_dir`, `downstream_analyses.metabolic_profile.gene_sets`
- **Brief keys (optional):** `metabolic_profile.label_col`, `.subtype_col`, `.group_col`, `.subsample_ceiling`

#### Primitives Referenced

- `@primitives/aucell_scoring.md` — `run_aucell()` and `add_auc_to_seurat()` (called, not redefined)
- `@primitives/seurat_v5_rules.md` — Rule 1 (GetAssayData layer), Rule 5 (UMAP detection)
- `@context/color_palettes.md` — `PALETTE_DIVERGING_6` for AUCell UMAP coloring

#### Deprecated v1 Patterns Fixed

- **Invalid gene set range syntax:** v1 used `"NDUFA1"-"NDUFA10"` (not valid R). These gene sets
  are not in v2 at all — gene sets are caller-provided.
- **3-stop diverging palette drift:** v1 used a unique blue-yellow-red 3-stop. v2 uses the
  lab canonical `PALETTE_DIVERGING_6` per `duplication_report.md §3.1`.
- **HumanFat assumptions in gene sets:** PPARG_Targets and all explicit gene vectors removed.

#### pdf()/dev.off() Usage

Two ComplexHeatmap calls — both exception #1 from CONVENTIONS.md §4:
1. `make_gene_heatmap()`: `pdf(...); draw(ht); dev.off()`
2. AUC z-score heatmap: `pdf(...); draw(ht); dev.off()`

#### Project-Specific Values Staged for Phase 4

- `examples/humanfat_metabolic_profile.md`:
  - Five metabolic pathway gene sets (Glycolysis, Beta_Oxidation, TCA_Cycle, OXPHOS, FA_Synthesis)
    with explicit gene vectors replacing the invalid range syntax from v1
  - PPARG_Targets gene set
  - Myelolipoma exclusion rationale and filter
  - `ec_subtype_colors` from HumanFat validated_examples.yaml

---

### 2.3 modules/trajectory_monocle3.md

#### Parameterization Table

| v1 hardcoded value | v2 parameter |
|---|---|
| `LABEL_COL = "ec_subtype"` | `LABEL_COL <- "project_specific"` from brief `downstream_analyses.trajectory.label_col` |
| `EXCLUDE_PATTERN = "^AEC"` | `EXCLUDE_PATTERN <- NULL` (optional; from brief `.trajectory.exclude_pattern`) |
| `START_CLUSTER = "CapEC2"` | `START_CLUSTER <- "project_specific"` from brief `.trajectory.start_cluster` |
| `SUBSET_NAME = "ec_traj"` | `SUBSET_NAME <- "celltype_traj"` (generic prefix) |
| Machine-specific UMAP reduction name assumption | `umap_key <- names(obj@reductions)[grepl(...)]` — Rule 5 detection |
| HumanFat-specific validated note in comments | Removed from module body; staged to Phase 4 |

#### Context Dependency Declaration Summary

- **Palettes (optional):** `subtype_colors`
- **Metadata columns (required):** `label_col`
- **Metadata columns (optional):** `group_col`
- **Brief keys (required):** `output_dir`, `downstream_analyses.trajectory.label_col`, `.trajectory.start_cluster`
- **Brief keys (optional):** `.trajectory.exclude_pattern`, `.trajectory.end_cluster`, `.trajectory.supplementary_plots`, `.trajectory.group_col`

#### Primitives Referenced

- `@primitives/visualization.md` — color scale references
- `@primitives/seurat_v5_rules.md` — Rule 1 (JoinLayers note for GetAssayData), Rule 5 (UMAP detection)
- `@context/color_palettes.md` — `PALETTE_PSEUDOTIME` (magma, NOT inferno — validated convention)

#### Deprecated v1 Patterns Fixed

- **SeuratWrappers dependency removed:** v1 assumed `SeuratWrappers::as.cell_data_set()`.
  v2 explicitly bypasses this (unavailable in r-env) and uses `new_cell_data_set()` manually.
- **`preprocess_cds()` removed:** causes `as_cholmod_sparse` ABI error. v2 skips it.

#### Project-Specific Values Staged for Phase 4

- `examples/humanfat_trajectory.md`:
  - `LABEL_COL = "ec_subtype"`, `EXCLUDE_PATTERN = "^AEC"`, `START_CLUSTER = "CapEC2"`
  - HumanFat EC subtype colors
  - Validated on 53,487 EC cells, 13 subtypes (HumanFat_Yang run 2026-03-05)

---

### 2.4 modules/pyscenic_regulons.md

#### Parameterization Table

| v1 hardcoded value | v2 parameter |
|---|---|
| `/media/david/Mayo/ryan/scRNAseq/Sinusoid/GeneLists/` (P2-8 machine path) | `database_dir` brief field |
| `/home/ryannachman/scenicplus/scenicplus/resources/allTFs_hg38.txt` (P2-8) | `tf_list` brief field |
| `/home/ryannachman/anaconda3/envs/scenicenv/bin/python` (P2-8) | `scenic_python_path` brief field |
| `n_workers = 16` hardcoded | `n_workers` brief field with default 16 |
| `LABEL_COL = "ec_subtype"` | `LABEL_COL <- "project_specific"` from brief |
| NKXSpleen subtype order in Python script | `SUBTYPE_ORDER = sorted(...)` + `project_specific` REPLACE comment |

#### Context Dependency Declaration Summary

- **Palettes (optional):** `subtype_colors`
- **Metadata columns (required):** `label_col`
- **Metadata columns (optional):** `sample_col`
- **Brief keys (required):** `output_dir`, `downstream_analyses.pyscenic.label_col`, `.pyscenic.scenic_python_path`, `.pyscenic.database_dir`, `.pyscenic.tf_list`
- **Brief keys (optional):** `.pyscenic.n_workers`, `.pyscenic.min_pct_cells`, `.pyscenic.min_mean_expr`

#### Primitives Referenced

- `@primitives/seurat_v5_rules.md` — Rule 1 (JoinLayers before GetAssayData), Rule 5 (UMAP name detection)
- `@primitives/r_environment.md` — execution environment rules; `nohup` usage note

#### Deprecated v1 Patterns Fixed

- **P2-8 Fix:** All three machine-specific absolute paths replaced with `project_specific` brief fields.
  The fix is documented in a dedicated "P2-8 Fix" section with a table comparing v1 and v2.
- **`pyscenic aucell` CLI bug documented:** v1 used the CLI. v2 explicitly marks this as broken
  (numpy `np.float` removal) and mandates Python API via `scenic_01b_aucell.py`.

#### Project-Specific Values Staged for Phase 4

- `examples/nkxspleen_pyscenic.md`: NKXSpleen project paths, spleen EC subtype ordering,
  validated 317 regulons from KidneyNew run
- `examples/humanfat_pyscenic.md`: HumanFat EC subtype configuration, 6 subtypes, 320 regulons

---

### 2.5 modules/bulk_concordance.md ← Merge of two v1 files

#### Source Files

| v1 file | Lines | Validation context |
|---|---|---|
| `bulk_concordance.md` | 759 | HumanFat PPARG OE, 2026-03-12 |
| `bulk_lfc_concordance_heatmap.md` | 583 | NKXSpleen, 2026-03-13 |

#### Parameterization Table

| v1 hardcoded value | v2 parameter |
|---|---|
| `BULK_CSV` → HumanFat-specific path | `BULK_DE_PATH <- "project_specific"` (brief: `bulk_csv`) |
| `ec_colors` (6-subtype EC palette) | `subtype_colors <- c("project_specific" = "#project_specific")` |
| `mylabel` subtype column | `SUBTYPE_COL <- "project_specific"` (brief: `subtype_col`) |
| `tissue_type` group column | `GROUP_COL <- "project_specific"` (brief: `group_col`) |
| `tissue_colors` (adipose tissue palette) | `group_colors <- c("project_specific" = "#project_specific")` |
| `RibHighEC` exclusion | `EXCLUDE_SUBTYPES <- c()` (brief: `exclude_subtypes`) |
| `PPARG` as TF of interest | `TARGET_TF <- NULL` (brief: `target_tf`; NULL skips Part 3) |
| `pparg_concordance` as score column | Derived: `make.names(paste0(EXPERIMENT_LABEL, "_concordance"))` |
| `"PPARG1_ROSI vs Control_ROSI HUVECs"` in comments | `EXPERIMENT_LABEL <- "project_specific"` (brief: `experiment_label`) |
| `"Up (PPARG)"` / `"Down (PPARG)"` category labels | `paste0("Up (", EXPERIMENT_LABEL, ")")` — generic, driven by brief |
| `output3/` directory | `OUTPUT_DIR <- file.path("output", "bulk_concordance")` |
| Condition/control sample names (Mode 2) | `condition_samples` / `control_samples` as `project_specific` sentinels |
| NKXSpleen biology_gene_sets (chemokines, etc.) | `biology_gene_sets <- list("project_specific" = c(...))` |
| EC-biology TF list (~180 TFs, Lambert 2018) | `human_tfs <- c("project_specific")` with examples/ reference |
| `output/` directory (Mode 2 v1) | Same `OUTPUT_DIR` from brief |

#### Two-Mode Comparison

| | Mode 1 — Signature score | Mode 2 — Parallel LFC |
|---|---|---|
| **Unit of analysis** | Cell (per-cell score) | Gene (per-gene LFC) |
| **Bulk DE method** | Gene list → up/down sets | Full DESeq2 re-run |
| **scRNA-seq method** | AddModuleScore (no new DE) | FindMarkers (independent DE) |
| **Seurat object at viz time?** | Yes | No (cached tables only) |
| **Main outputs** | UMAP, violins, bar plots, score CSV | pheatmap PDF + PNG |
| **Part 3 TF analysis** | Yes (conditional on target_tf) | No |
| **pdf()/dev.off() needed?** | Yes (Part 3 ComplexHeatmap only) | No (pheatmap filename= parameter) |

**Shared between modes:** `BULK_DE_PATH`, `EXPERIMENT_LABEL`, `SUBTYPE_COL`, `GROUP_COL`,
`EXCLUDE_SUBTYPES`, `OUTPUT_DIR`, `BULK_LFC_CUT`, `PADJ_CUT`, `LFC_CUT`, and the
`get_bulk_status()` helper.

#### Context Dependency Declaration Summary

- **Palettes (optional):** `subtype_colors`, `group_colors`
- **Metadata columns (required):** `subtype_col`
- **Metadata columns (optional):** `group_col`
- **Brief keys (required):** `output_dir`, `bulk_concordance.mode`, `.bulk_csv`, `.experiment_label`, `.subtype_col`
- **Brief keys (optional):** `.group_col`, `.exclude_subtypes`, `.target_tf`

#### Primitives Referenced

- `@primitives/differential_expression.md` — `is_ambient()` in annotation overlay; `make_volcano()` referenced as source for `de_df`
- `@primitives/visualization.md` — plot conventions
- `@primitives/seurat_v5_rules.md` — Rule 1 (JoinLayers), Rule 5 (UMAP detection)

#### Cross-Mode Design Decisions

1. **pheatmap filename= in Mode 2:** The v1 bulk_lfc_concordance_heatmap.md used pheatmap's
   built-in `filename=` parameter to write directly to file. This does NOT require pdf()/dev.off()
   (pheatmap manages its own device). v2 preserves this compliant pattern and explicitly notes
   "no pdf()/dev.off() needed" in the architecture comment.

2. **Part 3 TF analysis conditional:** v1 unconditionally ran TF analysis tailored to PPARG.
   v2 gates all of Part 3 on `!is.null(TARGET_TF)`. Projects with non-TF-OE bulk experiments
   (drug treatments, disease vs control) skip Part 3 entirely without modification.

3. **Score column naming:** v1 used `pparg_concordance` as a hardcoded column name. v2 derives
   it dynamically from `EXPERIMENT_LABEL` via `make.names()`. The sanity check (SD > 0.03) is
   preserved from v1.

4. **bulk_lfc_concordance_heatmap.md Mode 2 biology_gene_sets:** The NKXSpleen-specific gene
   categories (chemokines, adhesion, angiocrine, ECM, sinusoidal markers, TF program) are not
   in the v2 module. The `biology_gene_sets` list is `project_specific` sentinel with structure
   documented; the NKXSpleen instantiation goes in examples/nkxspleen_bulk_lfc.md.

#### pdf()/dev.off() Usage

- Mode 1 Part 3 TF heatmap: ComplexHeatmap — exception #1. The code is provided as a
  commented skeleton (since the calling context — which DE CSVs to load — is project-specific).
- Mode 2 pheatmap: No pdf()/dev.off(). pheatmap's `filename=` parameter handles device management.

#### Project-Specific Values Staged for Phase 4

- `examples/pparg_bulk_concordance.md`:
  - Mode 1 context: PPARG1_ROSI vs Control_ROSI HUVECs
  - HumanFat tissue color palette (6 adipose depots) for group_colors
  - EC subtype colors (6 subtypes) for subtype_colors
  - `RibHighEC` exclusion in EXCLUDE_SUBTYPES with rationale
  - EC-biology TF list (~180 TFs, Lambert 2018 basis) for human_tfs
  - Validated figure dimensions from 2026-03-12 run

- `examples/nkxspleen_bulk_lfc.md`:
  - Mode 2 context: 80,140 EC cells, 10 subtypes, 12 bulk samples, 6 Control / 6 NKX2-3 OE
  - biology_gene_sets for spleen EC biology (chemokines, adhesion, angiocrine, ECM, sinusoidal, TF)
  - 2-batch DESeq2 sample selection
  - Validated 2026-03-13 figure

---

### 2.6 modules/de_comprehensive_csv.md

#### Parameterization Table

| v1 hardcoded value | v2 parameter |
|---|---|
| `is_confound()` phantom (never defined anywhere) | Wired to `@primitives/differential_expression.md` — P2-3 fix |
| `SUBTYPE_LABELS` as hard dependency | `SUBTYPE_LABELS <- NULL` config variable; optional when `LABEL_COL` is set |
| `SUBTYPE_COL` as hard dependency | `LABEL_COL <- "project_specific"` config variable; set to NULL to omit `top_subtype` |
| `GROUP_COL` hardcoded implicit | `GROUP_COL <- "project_specific"` with REPLACE sentinel |
| `ec` as object name | `obj` (generic) |
| `meta[[group_col]]` (already parameterized in v1 Step 2) | Preserved; `GROUP_COL` now a config variable |
| `kw_results_cache.Rds` filename | Preserved in `OUTPUT_DIR` |
| `comprehensive_DE_stats.csv` fixed output name | `file.path(OUTPUT_DIR, "comprehensive_DE_stats.csv")` |
| No `group_order` parameter | `GROUP_ORDER <- NULL` config variable (optional ordered vector) |

#### P2-3 Punch-List Resolution

**Status: RESOLVED.**

v1 called `is_confound()` at Step 2 (`filter(!is_confound(gene))`). This function was a phantom —
called but never defined anywhere in the v1 library (per per_file_inventory.md §de_comprehensive_csv).
The inventory confirmed is_confound() was meant to be a superset of is_ambient() that also excludes
sex-linked genes, HLA class II, histones, and unannotated identifiers.

v2 resolution:
1. `is_confound()` is defined in `@primitives/differential_expression.md` (Phase 1 F3 fix).
2. v2 de_comprehensive_csv.md references `@primitives/differential_expression.md` in YAML frontmatter.
3. The module header documents the `is_confound()` vs `is_ambient()` distinction.
4. The Critical Constraints table explicitly says "Use `is_confound()` from @primitives/..."
5. A comment at the call site says "is_confound() is defined in @primitives/differential_expression.md".

No phantom function remains. The module is fully wired to the real primitive.

#### Context Dependency Declaration Summary

- **Palettes:** none required
- **Metadata columns (required):** `group_col`
- **Metadata columns (optional):** `label_col`
- **Brief keys (required):** `output_dir`, `downstream_analyses.de_comprehensive_csv.group_col`
- **Brief keys (optional):** `.de_comprehensive_csv.label_col`, `.group_order`, `.subtype_labels`

#### Primitives Referenced

- `@primitives/differential_expression.md` — `is_confound()` function (primary reason for reference)
- `@primitives/seurat_v5_rules.md` — Rule 1 (JoinLayers before GetAssayData)

#### Deprecated v1 Patterns Fixed

- **Phantom is_confound() (P2-3):** Wired to primitive. See P2-3 section above.
- **`ec` object name:** Renamed to generic `obj`.
- **`SUBTYPE_LABELS` / `SUBTYPE_COL` as hard dependencies:** Made optional; `top_subtype` column
  is NA when `LABEL_COL` is not set.

#### Analytical Contributions Preserved

- Vectorized KW implementation (rowRanks + matrix math) — preserved exactly from v1
- Two-stage caching pattern (kw_results_cache.Rds) — preserved
- Cell-count-weighted logFC vs others — preserved
- top_group2 populated only when logFC > 0.25 — preserved
- expm1() conversion for natural-scale mean expression — preserved

#### Project-Specific Values Staged for Phase 4

- `examples/humanfat_de_comprehensive.md`:
  - GROUP_COL = tissue group column, 6 tissue groups
  - LABEL_COL = EC subtype column
  - SUBTYPE_LABELS = 6 EC subtype names
  - GROUP_ORDER for tissue group display ordering
  - Validated output: ~37,692 EC cells; ~253 significant genes after is_confound(); 2026-04-23

---

## 3. Punch-List Resolution Across Group B

| Item | Status | Module(s) that resolve it |
|---|---|---|
| **P2-3** — `is_confound()` phantom in de_comprehensive_csv.md | ✅ RESOLVED | `de_comprehensive_csv.md` — wired to `@primitives/differential_expression.md` |
| **P2-7** — `celltype_subclustering.md` uses pdf()/dev.off() for ggplot outputs | ✅ RESOLVED | `celltype_subclustering.md` — all plots use `ggsave()` |
| **P2-8** — pySCENIC machine-specific paths | ✅ RESOLVED | `pyscenic_regulons.md` — all three paths replaced with `project_specific` brief fields |

---

## 4. Self-Check Grep Results — ALL SIX Modules

**Leakage terms checked:** HumanFat, ec_colors, tissue_colors, adipose_type_colors, mylabel,
STROMA_ORDER, Vertebrae, Iliac Crest, Femoral Head, Myelolipoma, CapEC, PPARG, TabulaSapiens,
WCM, source_file, NKXSpleen, ec_subtype_colors, tissue_type, patient_id, RibHighEC, BULK_CSV,
output3

**Result: ZERO matches in executable code.**

All 14 matches found are in Phase 4 staging sections at the bottom of each module file
("Project-Specific Values (Stage for Phase 4 examples/)"). Classification:

| Match | File | Location | Classification |
|---|---|---|---|
| `adipose_type_colors`, `source_file`, `tissue_type` | `celltype_subclustering.md` line 565–566 | Phase 4 staging section | **Acceptable** |
| `PPARG_Targets` | `metabolic_profile.md` line 65 | Brief schema comment showing example `gene_sets` key | **Acceptable** |
| `PPARG_Targets`, `Myelolipoma`, `ec_subtype_colors`, `HumanFat` | `metabolic_profile.md` lines 451–453 | Phase 4 staging section | **Acceptable** |
| `ec_subtype` (LABEL_COL example), `CapEC2`, `HumanFat` | `trajectory_monocle3.md` lines 290–292 | Phase 4 staging section | **Acceptable** |
| `NKXSpleen`, `HumanFat` | `pyscenic_regulons.md` lines 512–516 | Phase 4 staging section | **Acceptable** |
| `PPARG1_ROSI...`, `HumanFat`, `RibHighEC`, `NKXSpleen` | `bulk_concordance.md` lines 931–940 | Phase 4 staging section | **Acceptable** |
| `humanfat_de_comprehensive.md` | `de_comprehensive_csv.md` line 353 | Phase 4 staging section | **Acceptable** |

**CLEAN — zero leakage in executable code across all six modules.**

---

## 5. pdf()/dev.off() Exceptions Across All Six Modules

| Module | Exception # | Library | Context |
|---|---|---|---|
| `metabolic_profile.md` | #1 — ComplexHeatmap | `ComplexHeatmap::draw()` | AUC z-score heatmap + per-gene expression heatmap (2 calls) |
| `bulk_concordance.md` | #1 — ComplexHeatmap | `ComplexHeatmap::draw()` | Part 3 TF LFC heatmap (commented skeleton; conditional on TARGET_TF) |

All other modules: zero pdf()/dev.off() calls. All usages fall under CONVENTIONS.md §4 exception #1.

**Exception #2 (pheatmap + grid.text) and Exception #3 (circlize):** Neither appears in Group B.
Both were added to CONVENTIONS.md in Phase 2B pre-work (commit 92306fe) based on Group A findings.

**Mode 2 pheatmap in bulk_concordance.md:** Uses pheatmap's built-in `filename=` parameter for
direct file output. This does NOT use pdf()/dev.off() — pheatmap manages its own device when
`filename` is specified. This is compliant and requires no exception.

---

## 6. Cross-File Findings

### 6.1 UMAP Detection Rule 5 Used Consistently

All four modules that read Seurat embeddings use the dynamic UMAP detection pattern:
```r
umap_key <- names(obj@reductions)[grepl("umap", names(obj@reductions), ignore.case = TRUE)][1]
```
This is referenced to `@primitives/seurat_v5_rules.md` Rule 5 in every module. Consistent.

### 6.2 JoinLayers Rule 1 Annotated Consistently

Every module that calls `GetAssayData()` or `FindMarkers()` includes a `JoinLayers(obj)` call
annotated with a cross-reference to `@primitives/seurat_v5_rules.md` Rule 1. Consistent.

### 6.3 bulk_concordance.md — Mode Choice Mirrors the Decision Made in v1

The v1 self-documented the distinction between the two files at lines 1–35 of
`bulk_lfc_concordance_heatmap.md`. v2 preserves this decision framework as the "Choose This When"
section. The two-mode design accurately reflects the genuine analytical question difference:
cell-level (Mode 1) vs gene-level (Mode 2).

### 6.4 de_comprehensive_csv.md — Vectorized KW Not a Named Function

The v1 per_file_inventory.md noted: "the vectorized KW implementation is a code block, not a named
function — would benefit from being a named function." In v2 the implementation remains a code block
(consistent with the existing module structure where analysis steps are blocks, not functions).
Naming it `run_vectorized_kw()` would require it to be in a primitive or returned/usable, which
is not the pattern for module-level analysis steps. Left as-is by design.

### 6.5 is_confound vs is_ambient Decision Documented at the Call Site

The is_confound() usage in de_comprehensive_csv.md is the only place in Group B where this
function is called. The distinction is documented:
- In the module header (why is_confound is appropriate for comprehensive tables)
- In the Critical Constraints table
- As an inline comment at the call site

This is the only Group B module that uses is_confound(). All other Group B modules that do DE
filtering (e.g. bulk_concordance Part 1 annotation overlay) use is_ambient() from
`@primitives/differential_expression.md`, which is appropriate for the standard volcano context.

### 6.6 No New Primitives Created in Group B

No new primitives were authored in Group B. All six modules reference existing primitives.
The one new analytical pattern introduced (bulk concordance score via AddModuleScore subtraction)
is a caller-side helper within the module, not a primitive — consistent with the rule that
"caller-side helpers that do not belong in primitives" stay in modules.

---

## 7. Open Items for Group C and Group D

### From Group B findings

- **TODO-BC-1 (Part 3 TF heatmap skeleton):** The TF LFC heatmap in bulk_concordance.md Part 3
  is provided as commented pseudocode because the matrix assembly step requires project-specific
  DE CSV paths that cannot be generalized. The Phase 4 examples/ for HumanFat should provide the
  full working implementation as a reference pattern.

- **TODO-BC-2 (Mode 2 sc_group assignment pattern):** The `sc_group` column assignment in Mode 2
  Part 2 includes a pattern for splitting target vs reference cells, but the target subtype
  assignment lines are left as comments with `project_specific`. The Phase 4 NKXSpleen example
  should demonstrate this concretely.

- **TODO-BC-3 (Mode 1 + Mode 2 cross-validation):** The two modes have not been validated
  to produce coherent results on the same dataset (e.g. confirming that cells with high Mode 1
  scores in a given subtype correspond to the genes showing high Mode 2 concordance for that
  subtype). This cross-validation is a Phase 3/4 item.

### Carried forward from Phase 1 and Phase 2A

- **P2-1 (STROMA_ORDER undefined):** Still open. Not addressed in Group B.
- **P2-2 (IntegratePublicData DE injection):** Still open.
- **P2-4 (make_topgene_dotplot name collision):** Still open.
- **P2-5 (IntegratePublicData brief_template):** Still open.
- **P2-6 (CellChat projectData() bug):** Still open.
- **P2-9 (RenameFeatures P2-9 fix in load_formats):** Resolved in Group A.
- **P2-10 (doublet detection gap):** Still open.
- **G1 (anatomical_de_analysis.md orphan):** Still open.
- **G2 (atlas_co_umap):** Still open.
- **G3 (cross_atlas_dotplot):** Still open.

---

## 8. Anything Unexpected

### U1 — bulk_concordance.md merge produces a significant size reduction

The two v1 files (759 + 583 = 1342 lines) merged into 945 lines (-397). The reduction came
from: (a) eliminating duplicate "distinction from the other method" framing in both files,
(b) consolidating the shared config block (thresholds, bulk data loading, get_bulk_status()),
and (c) the mode comparison table replacing two separate overview sections.

### U2 — pheatmap filename= compliance clarifies a potential false alarm

The v1 `bulk_lfc_concordance_heatmap.md` used `pheatmap::pheatmap(filename = ...)`. This is
NOT a pdf()/dev.off() usage — pheatmap handles its own device when filename is specified.
The check confirms this is compliant without exception. Documented explicitly in both the
module and in this report to prevent future confusion.

### U3 — Mode 1 Part 3 TF list is the largest project-specific staging item

The EC-biology TF list in v1 bulk_concordance.md Part 3 contains ~180 genes organized by TF
family. This is the single largest body of biology-specific content in Group B. It was not
carried into v2 executable code (per leakage rules) and must be in examples/. The examples/
entry for this list is particularly important because the list organization (by family, with
emphasis on EC-relevant families) is itself analytical guidance, not just data.

### U4 — de_comprehensive_csv.md top_subtype is now genuinely optional

The v1 implementation required `SUBTYPE_LABELS` and `SUBTYPE_COL` as hard dependencies (would
fail silently or with an error if not set). v2 makes top_subtype fully optional — the column
is present in the output but contains NA when LABEL_COL is not set. This matches the module's
stated use case (≥ 3 groups; not every project has a named subtype column).
