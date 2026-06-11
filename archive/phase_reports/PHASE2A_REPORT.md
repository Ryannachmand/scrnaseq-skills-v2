---
# Phase 2 Group A Report — load_formats, cohort_plots, label_harmonization, feature_umap_plot
# Written: 2026-05-14
# Author: Phase 2 Group A migration agent
---

# Phase 2 Group A Report

Human review document. Covers what was created, parameterization decisions, deprecated
pattern fixes, project-specific values staged for Phase 4 examples/, uncertainties left
as TODOs, self-check results, and open items for subsequent groups.

---

## 1. Files Created

| File | v2 Lines | v1 Lines | Delta |
|---|---|---|---|
| `modules/load_formats.md` | 381 | 341 | +40 |
| `modules/cohort_plots.md` | 323 | 235 | +88 |
| `modules/label_harmonization.md` | 227 | 181 | +46 |
| `modules/feature_umap_plot.md` | 187 | 131 | +56 |
| **Total** | **1,118** | **888** | **+230** |

Delta is positive in all cases due to added frontmatter, configuration blocks,
parameterization scaffolding, and expanded documentation of exceptions.

---

## 2. Per-Module Details

---

### 2.1 modules/load_formats.md

#### Parameterization Table

| v1 hardcoded value | v2 parameter / change |
|---|---|
| `"_raw_gene_bc_matrices_h5\\.h5$"` Leimkuler file pattern | Replaced with generic comment: "adjust pattern to match your dataset's actual filenames" |
| `GSM4423510_MNC_op/2.2.filtered_feature_bc_matrix/` Wang path | Replaced with generic nested-directory example |
| `GSM4731560_CK80_gli1_barcodes.tsv.gz` Leimkuler file names | Replaced with `paste0(accession, "_barcodes.tsv.gz")` pattern + REPLACE comment |
| `GSE190965_barcodes.tsv.gz` Li spliced/unspliced file names | Replaced with `paste0(accession, "_barcodes.tsv.gz")` pattern + REPLACE comment |
| `project = "Li"` in CreateSeuratObject | Replaced with `project = dataset_name` |
| `isolatedstroma` object reference in comment | Removed (BoneMarrowStroma-specific) |
| `FemoralHead stores UMAP as UMAP_dim30` | Removed; general guidance references seurat_v5_rules.md Rule 5 |
| `CLEC4G + LYVE1 > 100 FPKM confirms liver sinusoidal ECs` | Preserved as example but annotated: "substitute appropriate markers for your tissue" |
| `rownames(so@assays$RNA@counts) <- sym` (deprecated slot assignment) | Replaced with `RenameFeatures(so, new.names = sym, assay = "RNA")` with fallback documented (see TODO below) |
| "Li dataset finding: 9 distinct sample prefixes" note | Removed (dataset-specific note) |
| "Last validated: TFEC Expression Atlas run" header text | Removed (project name in header) |

#### Context Dependency Declaration Summary

- **Palettes:** none required
- **Metadata columns (required):** none (columns are created during loading)
- **Metadata columns (optional):** `sample_col`
- **Brief keys (required):** `inputs.paths`, `inputs.type`
- **Brief keys (optional):** `project_name`

#### Primitives Referenced

- `@primitives/seurat_v5_rules.md` — Rule 1 (JoinLayers), Rule 5 (dynamic UMAP name)
- `@primitives/r_environment.md` — execution environment rules

#### Deprecated v1 Patterns Fixed

- **P2-9 Fix:** `rownames(so@assays$RNA@counts) <- sym` replaced with `RenameFeatures()` call.
  The v1 direct slot assignment is deprecated in Seurat v5 because the assay object structure
  changed from S4 slots to a layered architecture.

#### Project-Specific Values Staged for Phase 4

- `examples/load_formats_leimkuler.md` — Leimkuler file pattern `_raw_gene_bc_matrices_h5\.h5$`, GSE156644 context
- `examples/load_formats_wang.md` — Wang nested directory structure, GSM4423510 accession
- `examples/load_formats_li.md` — Li spliced/unspliced format, GSE190965, 9-sample splitting result
- `examples/load_formats_bonemarrowstroma.md` — FemoralHead UMAP_dim30 reduction name, isolatedstroma context

#### TODOs

**TODO-LF-1 (RenameFeatures API):** The v2 module calls `RenameFeatures(so, new.names = sym, assay = "RNA")`.
This is the documented Seurat v5 approach per the PHASE1_REPORT.md P2-9 note. However, the exact
signature has not been verified against the specific Seurat version in r-env. The fallback using
`GetAssayData()` / `SetAssayData()` per layer is documented inline in the module. Validate during
the next actual IntegratePublicData run with an Ensembl-ID object.

---

### 2.2 modules/cohort_plots.md

#### Parameterization Table

| v1 hardcoded value | v2 parameter / change |
|---|---|
| `so$dataset` column reference | `so[[DATASET_COL]]` — config variable, default "dataset" |
| `so$unified_label` column reference | `so[[LABEL_COL]]` — config variable, default "unified_label" |
| `so$sample_id` column reference | `so[[SAMPLE_COL]]` — config variable, default "sample_id" |
| `so$condition`, `so$organ`, `so$data_type` references | `CONDITION_COL`, `ORGAN_COL`, `DATA_TYPE_COL` — config variables |
| `ann_col` data.frame with hardcoded column names `Dataset`, `Organ`, `Condition`, `Data.Type` | `ann_track_cols` character vector — caller specifies which metadata columns to include |
| `dataset_colors` and `ct_colors` (assumed to exist in scope) | Explicit `project_specific` sentinels with REPLACE comments in config block |
| `dataset_gap_positions` using `col_order$dataset` | `col_order[[DATASET_COL]]` (parameterized) |
| BoneMarrowStroma validated figure size `13.25 × 7 inches` | Replaced with dynamic formula: `fig_w = max(n_samples * 0.3 + 4, 8)`, `fig_h = max(n_celltypes * 0.35 + 2, 5)` + REVIEW comment |
| `brief$label_col[[name]]` pattern (3A) | Already parameterized in v1; preserved as-is |
| `n_hvg` from brief | `N_HVG` config variable, default 2250 |

#### Context Dependency Declaration Summary

- **Palettes (required):** `dataset_colors`, `ct_colors`
- **Metadata columns (required):** `dataset_col`, `label_col`, `sample_col`
- **Metadata columns (optional):** `condition_col`, `organ_col`, `data_type_col`
- **Brief keys (required):** `project_name`, `output_dir`
- **Brief keys (optional):** `n_hvg`

#### Primitives Referenced

- `@primitives/visualization.md` — color scale, sizing conventions
- `@primitives/aesthetics.md` — typography rules, output conventions
- `@context/color_palettes.md` — diverging blue→white→red palette values
- `@primitives/seurat_v5_rules.md` — Rule 6 (AverageExpression dgCMatrix conversion)

#### Deprecated v1 Patterns Fixed

- **Correlation heatmap pdf()/dev.off():** The v1 uses `pdf()` + `dev.off()` for the correlation
  heatmap to support `grid.draw()` title overlay. In v2:
  - Simple output (no styled title): uses `ggsave(plot = ph$gtable, ...)` — compliant.
  - Styled title + subtitle overlay: retains `pdf()/dev.off()` as a documented exception.
    This is analogous to the ComplexHeatmap exception in PHASE1_REPORT.md §4.
    The technical reason is identical: pheatmap renders via grid and requires an active device
    for `grid.text()` title overlay. ggsave cannot intercept mid-render grid calls.

- **Chord diagram pdf()/dev.off():** circlize uses base R graphics and has no ggsave path.
  Documented as a third intentional exception (after ComplexHeatmap and pheatmap+title).

- **JoinLayers inline note:** The v1 (line 111) has an inline JoinLayers reminder in the
  checklist. In v2 this is replaced with a cross-reference to @primitives/seurat_v5_rules.md Rule 1.

#### Project-Specific Values Staged for Phase 4

- `examples/bone_marrow_stroma_cohort_plots.md` — validated figure size 13.25 × 7 in,
  specific annotation strip configuration (Dataset/Organ/Condition/Data.Type)

#### TODOs

**TODO-CP-1 (ann_colors structure):** The `ann_colors` variable is set to `project_specific`
sentinel. In practice this is a named list where each element is a named color vector matching
a column in `ann_col`. The Phase 4 BoneMarrowStroma example should document the specific
structure used in that run so future projects have a concrete reference.

**TODO-CP-2 (chord diagram output convention):** The chord diagram uses `pdf()/dev.off()` (third
exception). Consider whether this should be noted in CONVENTIONS.md §4 alongside the two existing
exceptions (ComplexHeatmap, pheatmap+grid-title). Left for human review.

---

### 2.3 modules/label_harmonization.md

#### Parameterization Table

| v1 hardcoded value | v2 parameter / change |
|---|---|
| `so@meta.data$unified_label <- unname(mapped)` hardcoded column | `so@meta.data[[UNIFIED_LABEL_COL]] <- unname(mapped)` — config variable |
| `ref_merged$unified_label` in TransferData | `ref_merged[[UNIFIED_LABEL_COL]]` — parameterized |
| `query$unified_label` assignment | `query[[UNIFIED_LABEL_COL]]` — parameterized |
| `"Mesenchymal"` as hardcoded broad_label in case_when | `brief$label_transfer_broad_label` — required brief field; stop() if not set |
| "Vertebrae + FemoralHead reference" in YAML example notes | Replaced with `<reference_dataset_1> + <reference_dataset_2>` placeholders |
| BoneMarrowStroma result_distribution values (Adipo-MSC: 1823, Osteo-MSC: 912) | Replaced with `<label_1>: <count>` placeholders |
| `score_high: 0.75`, `score_low: 0.40` in code | Preserved as defaults via `%||%` operator; configurable from brief |

#### Context Dependency Declaration Summary

- **Palettes:** none
- **Metadata columns (optional):** `unified_label_col`
- **Brief keys (required):** `label_col` (dict: dataset_name → source column)
- **Brief keys (optional):** `label_transfer_reference`, `label_transfer_score_high`, `label_transfer_score_low`, `label_transfer_broad_label`, `known_label_mappings`

#### Primitives Referenced

- `@primitives/seurat_v5_rules.md` — Rule 1 (JoinLayers in ref_merged construction), Rule 3 (metadata assignment safety)

#### Deprecated v1 Patterns Fixed

None — the v1 file had no deprecated Seurat v5 slot assignment issues or pdf()/dev.off() patterns.
The `dplyr::case_when()` namespace-qualified call was already correct in v1 and is preserved.

#### Project-Specific Values Staged for Phase 4

- `examples/bone_marrow_stroma_label_harmonization.md`:
  - `broad_label = "Mesenchymal"` for low-confidence stromal cells
  - `label_transfer_reference = ["Vertebrae", "FemoralHead"]` datasets used as reference
  - Result distribution: Adipo-MSC (1823), Osteo-MSC (912), Unknown (75) from Li transfer
  - Sort-strategy labeling patterns: "CD14 Isolated Stroma", "PE Isolated Stroma"
  - Passaged cells labeling: "Passaged Stroma"
  - Derived display column: "Cultured Adipo-like MSC" from "Adipo-MSC"

#### TODOs

**TODO-LH-1 (broad_label nullability design):** The v2 module requires `brief$label_transfer_broad_label`
and throws a stop() if it is missing. An alternative design would use "Unknown" as a universal
fallback, which is safer but less informative. The current design prioritizes correctness (forcing
explicit biology-appropriate labeling) over convenience. Human review: is the stop() too strict?
If many projects will use this without a meaningful broad_label, consider changing to a warning.

---

### 2.4 modules/feature_umap_plot.md

#### Parameterization Table

| v1 hardcoded value | v2 parameter / change |
|---|---|
| `so <- readRDS("output/your_object.Rds")` | `RDS_PATH <- "project_specific"` config variable with REPLACE comment |
| `genes <- c("GENE1", "GENE2", "GENE3")` | `GENES <- brief$genes %||% c("GENE1", "GENE2", "GENE3")` — pulls from brief |
| `ggsave("output/feature_umap.pdf", ...)` hardcoded output path | `file.path(OUTPUT_DIR, "feature_umap.pdf")` — OUTPUT_DIR from config |
| No subset filtering option | Added `SUBSET_COL` / `SUBSET_VAL` config variables (both NULL by default) |
| No mention of subset compatibility | Added header note: "compatible with any Seurat object from LargeDataset or IntegratePublicData pipelines, and any subset object" |

#### Context Dependency Declaration Summary

- **Palettes:** none (uses fixed `c("lightgrey", "blue")` per Seurat FeaturePlot default)
- **Metadata columns (optional):** `label_col` (for subsetting)
- **Brief keys (required):** `output_dir`
- **Brief keys (optional):** `genes`

#### Primitives Referenced

- `@primitives/visualization.md` — color scale references, general plot structure
- `@primitives/aesthetics.md` — UMAP Feature Plots section (lightgrey→blue scale, raster=FALSE rule)
- `@primitives/seurat_v5_rules.md` — Rule 5 (dynamic UMAP reduction name)

#### Deprecated v1 Patterns Fixed

- None — v1 already used `ggsave()` correctly with `device = "pdf", useDingbats = FALSE`.
- Dynamic UMAP detection code is kept inline (critical to the recipe) with a cross-reference
  to @primitives/seurat_v5_rules.md Rule 5 noting the canonical location.
  This is consistent with the duplication_report.md §1.3 recommendation: "Feature_umap has extra
  context; keep both but add cross-ref."

#### Project-Specific Values Staged for Phase 4

- None. The v1 file had no project-specific values beyond the generic template placeholders.
  `examples/coculture_feature_umap.md` (Phase 4) can document the specific CoCulture gene list
  and the run context if needed, but nothing in v1 required migration.

#### TODOs

**TODO-FU-1 (AutoPointSize internal API):** The recipe calls `Seurat:::AutoPointSize(...)` using
the triple-colon internal accessor. This is the only documented way to match Seurat's exact
point size computation (validated against native FeaturePlot in the CoCulture run). However,
internal Seurat functions can change between versions without warning. The verification block
at the bottom of the module catches any regression. If a future Seurat update breaks this,
compare pt.size to a native FeaturePlot and refit the formula.

---

## 3. Cross-File Findings

### Duplication handled

**JoinLayers (duplication_report.md §1.2):** All three IntegratePublicData modules that
previously contained inline JoinLayers reminders now reference @primitives/seurat_v5_rules.md
Rule 1 instead of repeating the code. This eliminates the duplication noted across:
- load_formats.md (v1 line 159)
- cohort_plots.md (v1 line 111)
- label_harmonization.md (v1 lines referenced)

**Dynamic UMAP detection (duplication_report.md §1.3):** The inline code block in
feature_umap_plot.md is kept (it provides critical context for the recipe) and cross-referenced
to seurat_v5_rules.md Rule 5. This follows the inventory recommendation.

**AverageExpression dgCMatrix note:** cohort_plots.md references seurat_v5_rules.md Rule 6 for
the `as.matrix()` wrapper requirement — consistent with the duplication_report.md §2.4 pattern.

### Annotation metadata column standard

cohort_plots.md, label_harmonization.md, and load_formats.md all reference the same set of
IntegratePublicData standard columns (`dataset`, `unified_label`, `sample_id`, `condition`,
`organ`, `data_type`). In v2 these are all parameterized via config variables with defaults
matching the pipeline standard. This is consistent across all three modules.

### pdf()/dev.off() exceptions — new occurrences

Phase 1 documented one exception: ComplexHeatmap.
Phase 2 Group A adds two more:

| Exception | File | Reason |
|---|---|---|
| ComplexHeatmap | differential_expression.md | No ggsave path for draw() |
| pheatmap + grid title overlay | cohort_plots.md Part A | grid.text() requires active device |
| circlize chord diagram | cohort_plots.md Part B | Base R graphics, no ggsave path |

**Recommendation for human review:** Consider adding these two exceptions to CONVENTIONS.md §4.
Currently CONVENTIONS.md §4 says "Never use `pdf()` / `dev.off()` in generated scripts" with
ComplexHeatmap as "the only exception" — but two new exceptions are now documented.

---

## 4. Self-Check Grep Results

Grep run on all four module files in `/home/ryannachman/claude-skills-v2/modules/` for:
`HumanFat, ec_colors, tissue_colors, adipose_type_colors, mylabel, STROMA_ORDER, Vertebrae,
Iliac Crest, Femoral Head, Myelolipoma, CapEC, PPARG, TabulaSapiens, WCM, source_file`

**Result: ZERO matches in executable code.**

Secondary grep for additional sentinel values:
`ec_subtype, tissue_type, patient_id, "Li", Adipo-MSC, Osteo-MSC, BoneMarrowStroma, Mesenchymal`

**Result:** Four matches found, all in documentation comments:
1. `label_harmonization.md` line 5 — migration header comment (not code)
2. `label_harmonization.md` line 108 — comment showing example values for the `broad_label` parameter (not code)
3. `label_harmonization.md` line 115 — stop() error message example text (not code)
4. `label_harmonization.md` line 204 — Phase 4 cross-reference note (not code)

All four are in comment/documentation context. No project-specific values appear in
executable code blocks in any of the four modules.

Tertiary grep for GEO accession patterns:
`GSM, GSE, CK80, gli1, MNC_op, UMAP_dim30, isolatedstroma`

**Result:** Two matches — `GSM` and `GSE` appear in the "GSM vs GSE Accessions" section of
`load_formats.md` as generic GEO terminology (not specific accession numbers). CLEAN.

---

## 5. Open Items for Group B/C/D

### For whoever handles Group B (IntegratePublicData pipeline + remaining methods)

- **P2-5:** `IntegratePublicData/brief_template.md` is still missing. The CONVENTIONS.md brief
  schema covers both pipelines, but a pipeline-specific template with the IntegratePublicData
  standard metadata columns (`dataset`, `unified_label`, `condition`, `organ`, `data_type`) would
  be valuable. These columns are referenced in load_formats, cohort_plots, and label_harmonization
  but never formally declared as pipeline standards.

- **P2-2:** The IntegratePublicData pipeline still has no direct injection of DE functions
  (`run_findmarkers`, `make_volcano`, etc.) for the `anatomical_de_analysis.md` orphan. When
  that file is migrated (Group B or C), it must explicitly reference
  `@primitives/differential_expression.md`.

- **IntegratePublicData standard columns:** Group A's three modules all default to
  `dataset`, `unified_label`, `sample_id`, `condition`, `organ`, `data_type` as metadata column
  names. These should be codified as the IntegratePublicData pipeline standard in either
  `context/lab_context.md` (add a pipeline_metadata_defaults section) or a new
  `context/integratePublicData_defaults.md`. Currently each module declares them independently.

### For whoever handles the pipeline.md migration

- **cohort_plots.md routing:** SKILL.md currently shows `modules/cohort_plots.md` under the
  `downstream_analyses.cohort_plots` key. This is correct. The module is ready for injection.

- **pdf()/dev.off() exception list:** CONVENTIONS.md §4 needs updating to list three exceptions
  (ComplexHeatmap, pheatmap+grid-title, circlize). This is a documentation-only change, not a code change.

### For Phase 4 (examples/ authoring)

Based on Phase 2 Group A migration, the following project-specific instantiation files are needed:

| File | Content |
|---|---|
| `examples/load_formats_bonemarrowstroma.md` | FemoralHead UMAP_dim30, isolatedstroma context, GSM file patterns for Wang/Li/Leimkuler |
| `examples/bone_marrow_stroma_cohort_plots.md` | Validated figure size 13.25 × 7, ann_colors structure |
| `examples/bone_marrow_stroma_label_harmonization.md` | Mesenchymal broad_label, Vertebrae+FemoralHead reference, result distributions |
| `examples/coculture_feature_umap.md` | Gene list from CoCulture run (optional) |
| `examples/ec_functional_gene_sets.md` | EC biology gene sets from v1 differential_expression.md (referenced in Phase 1 as F4 item) |

---

## 6. Anything Unexpected

### U1 — pheatmap and circlize require two more pdf()/dev.off() exceptions

CONVENTIONS.md currently says ComplexHeatmap is "the only exception" to the ggsave rule.
But cohort_plots.md uses pheatmap + grid for the correlation heatmap title overlay, and
circlize for the chord diagram. Both have legitimate technical reasons (no ggsave path).
Human review should decide whether to:
(a) Update CONVENTIONS.md to list three exceptions, or
(b) Restructure the correlation heatmap to avoid the grid overlay (simpler title in pheatmap's main= parameter, accepting the limitation of unstyled text)

### U2 — label_harmonization.md `%||%` operator dependency

The module uses `%||%` (null-coalescing operator from rlang) in several places:
```r
score_high <- brief$label_transfer_score_high %||% 0.75
GENES <- brief$genes %||% c("GENE1", ...)
```
This operator is available when rlang (a Seurat dependency) is loaded, but it is not
explicitly documented as available in r_environment.md. If a future script runs without
loading Seurat, `%||%` will fail. A safer alternative is the base R pattern
`if (!is.null(x)) x else default`. Flagged here for awareness; not fixed because rlang
is always loaded in these scripts via Seurat.

### U3 — feature_umap_plot.md aesthetics discrepancy (pre-existing)

The v2 module correctly uses `c("lightgrey", "blue")` per @primitives/aesthetics.md UMAP
Feature Plots section. However, primitives/visualization.md uses
`c("#F5F5F5","#FFF9C4","#FFB300","#E53935")` for gene expression in dot plots. These are
different plot types with intentionally different color scales — this is correct behavior.
The potential confusion: a developer reading visualization.md might apply the warm palette
to feature UMAPs. The UMAP feature plots section of aesthetics.md correctly documents the
lightgrey→blue scale. No fix needed; noted for awareness.
