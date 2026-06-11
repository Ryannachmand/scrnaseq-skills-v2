# Phase 2 Group C Report

**Modules authored:** `multi_group_de_analysis.md`, `atlas_co_umap.md`, `cross_dataset_dotplot.md`
**Commit:** Phase 2 Group C (commit bb19b01)
**Date:** 2026-05-15

---

## 1. multi_group_de_analysis.md

**File path:** `modules/multi_group_de_analysis.md`
**v2 line count:** 687
**Source v1 file:** `~/claude-skills/pipelines/IntegratePublicData/methods/anatomical_de_analysis.md`
**v1 line count:** 347

### Parameterization Table

| v1 hardcoded value | v2 parameter name |
|---|---|
| `c("Vertebrae", "Iliac Crest", "Femoral Head")` factor levels | `GROUPS` (character vector, N ≥ 2) |
| `"site"` as the group column | `GROUP_COL` |
| `"unified_label"` as the cell type column | `LABEL_COL` |
| `STROMA_ORDER` for x-axis cell type ordering | `LABEL_ORDER` (optional; NULL = alphabetical) |
| `label_colors` named vector (coupled to STROMA_ORDER) | `LABEL_COLORS` (optional; NULL = grey30) |
| `"Femoral Head"` for section label pinning | `groups[length(groups)]` (auto-derived: last group = rightmost panel) |
| `ncol = 3` in `facet_wrap` | `ncol = length(groups)` |
| `is_ambient()` (patterns-only, local definition) | `is_ambient()` from `@primitives/differential_expression.md` (full dual-filter) |
| EC-biology `functional_gene_sets` | `functional_gene_sets` caller-provided (project_specific sentinel) |
| `run_findmarkers`, `make_volcano`, etc. (undefined cross-pipeline) | Referenced from `@primitives/differential_expression.md` |
| `N_GO_TERMS = 8`, `N_GENES_PER_SECT = 8` | `N_GO_TERMS`, `N_GENES_PER_SECT` config constants |

### Context Dependency Declarations

```yaml
requires_context:
  palettes: [group_colors, label_colors (optional)]
  metadata_columns:
    required: [label_col, group_col]
    optional: [label_order]
  brief_keys:
    required: [output_dir, groups, group_col, label_col, comparisons, functional_gene_sets]
    optional: [label_order, label_colors, n_each, padj_cut, lfc_cut]
```

### Primitives Referenced

- `@primitives/differential_expression.md` — `run_findmarkers`, `make_volcano`, `make_overall_heatmap`, `make_functional_heatmap`, `make_functional_dotplot`, `make_pathway_barplot`, `is_ambient`
- `@primitives/visualization.md` — dotplot and direction strip conventions
- `@primitives/aesthetics.md` — typography and color philosophy
- `@primitives/seurat_v5_rules.md` — Rule 1 (JoinLayers), Rule 5 (UMAP detection)

### Risk Handling

**G1 Risk #1 — cross-pipeline injection:** RESOLVED. Module declares `references: ["@primitives/differential_expression.md"]` in frontmatter. The deployment agent injects this primitive for all jobs using this module, making `run_findmarkers`, `make_volcano`, etc. available. This resolves the original VertebraeMyelolipoma failure where those functions were called but never injected.

**G1 Risk #2 — N-group plot layout:** SCOPED. Module explicitly documents: "Comparisons are always pairwise — each comparison is a `list(label, ident1, ident2)` specifying two groups." The N-group structure is for the Seurat object (N groups present); DE and visualizations are pairwise within that object. `make_anatomical_dotplot` and `make_functional_dotplot` (from primitive) both use `facet_wrap(ncol = length(groups))` so all N panels render from the same pairwise comparison. N-way pairwise across all group combinations is NOT attempted.

**G1 Risk #3 — label_colors coupling:** RESOLVED. `LABEL_COLORS` is optional. When `NULL`, all x-axis labels render in grey30. The `make_anatomical_dotplot` function handles the `NULL` case with `rep("grey30", length(ct_present))`.

**G1 Risk #4 — AMBIENT_EXPLICIT absence:** RESOLVED. Module references `is_ambient()` from `@primitives/differential_expression.md` which contains the full dual-filter (`AMBIENT_PATTERNS` + `AMBIENT_EXPLICIT` with MALAT1, platelet markers, etc.). The local v1 patterns-only definition is not carried forward.

### TODOs in Module

- `LABEL_ORDER` (STROMA_ORDER): The v1 file references `STROMA_ORDER` at lines 68 and 107 but it is never defined in any v1 library file. The module supports any `label_order` vector; the BoneMarrowStroma example (Phase 4) must curate and provide this. Documented with `# TODO:` comment at the config variable.

### N-Group Scope Statement

Confirmed: the module supports **N-group objects with pairwise comparisons**. The `GROUPS` vector defines all N group levels for factor ordering and panel display. Each entry in `comparisons` is a pairwise `list(label, ident1, ident2)`. N-way simultaneous multi-group testing is explicitly out of scope — use `@modules/de_comprehensive_csv.md` for that.

### Function Name Collision Resolution

**make_topgene_dotplot:** The v1 anatomical_de_analysis.md defined `make_topgene_dotplot()` with a different signature (hardcoded 3-panel layout, no `groups` argument). The v2 primitive's `make_topgene_dotplot` is already a parameterized N-group version. This module defines **`make_anatomical_dotplot()`** — a new function that adds `groups`, `label_order`, and `label_colors` capabilities. There is no collision because the names are different. The primitive's `make_topgene_dotplot` remains available for use in other modules.

**make_functional_dotplot:** The v1 anatomical version hardcoded `"Femoral Head"` as the pin-to panel. The v2 primitive's `make_functional_dotplot` already handles N groups dynamically, pinning section labels to `tail(levels(factor(dot_df$group_var)), 1)`. Since `GROUP_COL` factor levels are set to `GROUPS` before the call, the rightmost panel is always `GROUPS[length(GROUPS)]`. No rename needed — the primitive version is used directly.

### Project-Specific Values Staged for Phase 4

`examples/bonemarrow_3site_anatomical.md`:
- `GROUPS = c("Vertebrae", "Iliac Crest", "Femoral Head")`
- `STROMA_ORDER` (requires curation from project records — never defined in v1)
- `LABEL_COLORS` keyed by `STROMA_ORDER` (undefined in v1; must be created)
- `functional_gene_sets` for bone marrow stromal biology
- `GROUP_COLORS` for 3 anatomical site groups (direction strip)
- 3 pairwise `comparisons` entries
- Validated run: BoneMarrowStroma, 2026-03-16

---

## 2. atlas_co_umap.md

**File path:** `modules/atlas_co_umap.md`
**v2 line count:** 367
**Source v1 file:** `~/claude-skills/pipelines/IntegratePublicData/methods/co_umap_embedding.md`
**v1 line count:** 166

### Parameterization Table

| v1 hardcoded value | v2 parameter name |
|---|---|
| `"WCM"` as in-house dataset label | `INHOUSE_LABEL` |
| `"TS"` as atlas dataset label | `ATLAS_LABEL` |
| `"dataset"` as source column | `SOURCE_COL` |
| `"CapEC"` as highlighted in-house subtype | `HIGHLIGHT_INHOUSE` |
| `"Fat EC"` as highlighted atlas group | `HIGHLIGHT_ATLAS` |
| `CB_CAPEC = "#E69F00"` (orange) | `COLOR_INHOUSE` (default `"#E69F00"`, Okabe-Ito orange) |
| `CB_FAT_EC = "#0072B2"` (blue) | `COLOR_ATLAS` (default `"#0072B2"`, Okabe-Ito blue) |
| `ec_colors` (HumanFat palette) | `subtype_colors` (optional, passed as `SUBTYPE_COL` color vector) |
| `N_PER_SUBTYPE = 1500` (constant) | `DOWNSAMPLE_N` (config constant, tunable with documented rationale) |
| `ec_subtype` / `organ` column names | `SUBTYPE_COL`, `ATLAS_GROUP_COL` |
| `output2/` directory | `OUTPUT_DIR` |
| `ts_harmony_umap_presample_coords.csv` filename | `COORDS_CSV` |

### Context Dependency Declarations

```yaml
requires_context:
  palettes: [subtype_colors (optional)]
  metadata_columns:
    required: [source_col, subtype_col, atlas_group_col]
    optional: []
  brief_keys:
    required: [output_dir, source_col, inhouse_label, atlas_label, subtype_col,
               atlas_group_col, highlight_inhouse, highlight_atlas, coords_csv]
    optional: [n_per_subtype, color_inhouse, color_atlas, subtype_colors]
```

### Primitives Referenced

- `@primitives/harmony_integration.md` — preprocessing recipe for the co-embedded object
- `@primitives/visualization.md` — UMAP recipes and panel construction
- `@primitives/aesthetics.md` — typography and color rules
- `@primitives/seurat_v5_rules.md` — Rule 1 (JoinLayers), Rule 5 (UMAP reduction detection)

### Risk Handling

**G2 Risk #1 — n_per_subtype tuning:** DOCUMENTED. Module includes a "Tuning the downsampling rate" subsection with the formula `target_n = floor(0.25 * n_atlas_cells / n_subtypes)` for recalculating with different size ratios. The Phase 4 examples/ file documents the specific WCM/TS ratio analysis (in-house ~38k, atlas ~100k → DOWNSAMPLE_N = 1500 gives atlas = 91% of UMAP input).

**G2 Risk #2 — UMAP geometry sensitivity:** HANDLED. A runtime warning fires if any subtype yields < 500 cells after downsampling: `warning(sprintf("DOWNSAMPLE_N = %d produces < 500 cells for subtype(s): %s ...", ...))`.

**G2 Risk #3 — Custom fonts:** RESOLVED. Module defaults to system fonts. The v1's `base_family = "Source Sans 3"` is documented in the Risk section and Phase 4 staging as a known issue (font unavailable in r-env). No font arguments appear in the module's plot code.

**G2 Risk #4 — cairo_pdf:** DOCUMENTED. Module explicitly uses `device = cairo_pdf` for patchwork output and documents: (a) why cairo is required (only first panel renders with base `pdf()`), (b) that `useDingbats = FALSE` must be omitted with `cairo_pdf`, (c) what to check if cairo is unavailable.

### TODOs in Module

None. All four risks have been addressed with documented behavior or explicit warnings.

### Project-Specific Values Staged for Phase 4

`examples/tabula_sapiens_co_umap.md`:
- `INHOUSE_LABEL = "WCM"`, `ATLAS_LABEL = "TS"`, `SOURCE_COL = "dataset"`
- `HIGHLIGHT_INHOUSE = "CapEC"`, `HIGHLIGHT_ATLAS = "Fat EC"`
- `DOWNSAMPLE_N = 1500` with WCM/TS ratio analysis (in-house ~38k, atlas ~100k; 1500/subtype × 6 subtypes = ~9k in-house vs ~100k atlas → atlas 91%)
- `COORDS_CSV = "output2/ts_harmony_umap_presample_coords.csv"`
- Validated script filenames: `06f_ts_umap_presample.R`, `06g_ts_umap_replot.R`
- Validated figure dimensions: `width = 21, height = 5.5`
- Known font issue: v1 used `"Source Sans 3"` (unavailable in r-env)
- Validated run: TabulaSapiensComparison, 2026-04-15

---

## 3. cross_dataset_dotplot.md

**File path:** `modules/cross_dataset_dotplot.md`
**v2 line count:** 443
**Source v1 file:** `~/claude-skills/pipelines/IntegratePublicData/methods/cross_atlas_dotplot.md`
**v1 line count:** 232

### Parameterization Table

| v1 hardcoded value | v2 parameter name |
|---|---|
| `"WCM"` / `"TS"` dataset labels | `SOURCES_ORDER` character vector; individual sources referenced via `SOURCE_COL` |
| `EC_SUBTYPES` (WCM subtype vector) | `subtype_col` — subtypes derived from metadata, not hardcoded |
| `"WCM Fat CapEC"` column label convention | Constructed via `paste0(SOURCE_COL_value, " ", SUBTYPE_COL_value)` — configurable |
| `SEC1_FORCE` (hardcoded TF gene vector) | `FORCE_GENES_SEC1` (optional, NULL default; requires per-project justification) |
| Background colors (warm yellow WCM, light red TS Fat) | `GROUP_COLORS` named vector (optional) |
| Custom fonts `"Source Sans 3"`, `"Playfair Display"` | Removed; system fonts as default |
| `z_cap = ±2` | `Z_CAP` config constant |
| `n_variable = 50` (centroid correlation) | `N_VARIABLE` config constant |
| `output2/` directory | `OUTPUT_DIR` |

### Context Dependency Declarations

```yaml
requires_context:
  palettes: [group_colors (optional)]
  metadata_columns:
    required: [source_col, subtype_col, col_group_col]
    optional: []
  brief_keys:
    required: [output_dir, source_col, subtype_col, sources_order, marker_genes, col_group_col]
    optional: [n_variable, z_cap, group_colors, sec1_force]
```

### Primitives Referenced

- `@primitives/visualization.md` — dotplot recipes
- `@primitives/aesthetics.md` — typography, color philosophy
- `@context/color_palettes.md` — diverging scale for z-scored expression
- `@primitives/seurat_v5_rules.md` — Rule 1 (JoinLayers), Rule 6 (AverageExpression dgCMatrix)

### Risk Handling

**G3 Risk #1 — depth correction comparability assumption:** DOCUMENTED PROMINENTLY. The module header contains a dedicated "Depth Correction Assumption" section (before any code) warning that within-dataset z-scoring is only valid when both datasets use comparable normalization and gene panels. Includes guidance for detecting failures (genes consistently at z_cap in one direction).

**G3 Risk #2 — SEC1_FORCE as project hack:** RESOLVED. The v1 `SEC1_FORCE` constant becomes `FORCE_GENES_SEC1` in the v2 config block (NULL default). The module header documents it as "a project-specific analysis decision that must be justified per project — it is not a general-purpose feature." The Phase 4 staging section explicitly requires per-gene justification in examples/.

**G3 Risk #3 — Custom fonts:** RESOLVED. `"Source Sans 3"` and `"Playfair Display"` are documented as unavailable in r-env (Risk section + Phase 4 staging). No `family` arguments appear in the module's plot code. System fonts used throughout.

**G3 Risk #4 — Column ordering fragility:** DOCUMENTED AND HANDLED. The module includes `colnames(cor_mat) <- gsub("\\.", " ", colnames(cor_mat))` for CSV round-trip fix. The Critical Constraints table documents the whitespace-to-dot issue. The column ordering code uses `sub(" .*", "", col_labels)` to extract source prefix — this is robust to special characters in subtype names as long as the source prefix doesn't contain spaces.

### TODOs in Module

None. All four risks have been addressed.

### Project-Specific Values Staged for Phase 4

`examples/tabula_sapiens_dotplot.md`:
- Dataset source labels and `SOURCE_COL = "dataset"`
- Column label convention: in-house as `"{inhouse_label} {subtype}"`, atlas organs as `"{organ_tissue} EC"` (toTitleCase)
- `marker_genes`: 3-section named list with HumanFat-specific content (Section 1 shared, Section 2 in-house-specific, Section 3 atlas-specific)
- `FORCE_GENES_SEC1`: specific TF gene vector with rationale (TFs present in in-house but below atlas threshold due to depth)
- `GROUP_COLORS`: warm yellow (`"#FFF8E1"`) for in-house, light red (`"#FFEBEE"`) for primary atlas column
- Validated figure dimensions from 2026-03-24 run
- Known font issue: v1 used `"Source Sans 3"` (x-axis) and `"Playfair Display"` (y-axis gene labels) — both unavailable

---

## Known-Atlas Convention

### Requirement

The user has indicated Tabula Sapiens will be the most frequently-used atlas and asked
that the platform recognize "Tabula Sapiens" in a user brief and configure the relevant
modules appropriately. This is a **Phase 3 concern** (deployment agent + context registry
+ brief schema) and is NOT implemented in any Group C module.

### Proposed Design

A `known_atlases:` block in `context/validated_examples.yaml` (or a new
`context/known_atlases.md`) that documents Tabula Sapiens as a named atlas with its
canonical parameters:

```yaml
known_atlases:
  tabula_sapiens:
    display_name: "Tabula Sapiens"
    source_col: dataset
    atlas_label: "TS"
    atlas_group_col: null          # TODO: confirm column name from TS download
    typical_cell_count: ~100000
    recommended_downsample_n: 1500  # validated for in-house ~38k, 6 subtypes
    common_organ_labels:            # values appearing in atlas_group_col
      - "Fat EC"
      - "Kidney EC"
      - "Liver EC"
      # etc. — expand from TabulaSapiensComparison project records
    notes: |
      Official name: "Tabula Sapiens" (Tabula Sapiens Consortium, Science 2022).
      Download format: h5ad. Cell type column confirmed as 'cell_ontology_class'.
      Validated 2026-04-15 (TabulaSapiensComparison project).
```

When the deployment agent reads a brief containing `atlas: "Tabula Sapiens"` (or a
recognized alias), it looks up this entry and injects:
- `ATLAS_LABEL <- "TS"` into `atlas_co_umap.md` CONFIG
- `atlas_group_col` into `atlas_co_umap.md` CONFIG
- `SOURCES_ORDER` second element into `cross_dataset_dotplot.md` CONFIG
- `DOWNSAMPLE_N` recommendation into `atlas_co_umap.md` CONFIG

The modules remain atlas-agnostic — the deployment agent performs the injection, not the
module code itself.

### Cross-References in Modules

Both `atlas_co_umap.md` and `cross_dataset_dotplot.md` include a **Known-Atlas Convention**
section at the top of the module (after the YAML frontmatter) that reads:

> "When a brief names a recognized atlas (e.g. 'Tabula Sapiens'), the deployment agent
> can auto-populate `source_col`, `atlas_label`, `atlas_group_col`, and `n_per_subtype`
> from the context registry (`context/known_atlases` — proposed, not yet implemented).
> See PHASE2C_REPORT.md for the design. This module is atlas-agnostic; all atlas identity
> flows through parameters."

---

## Punch-List Items Resolved in Group C

| Item | Status |
|---|---|
| **P2-1** (STROMA_ORDER undefined in v1 library) | STAGED — module supports `label_order` argument; `STROMA_ORDER` is a Phase 4 examples/ concern for BoneMarrowStroma. Not resolvable in the module — requires human curation of cell type ordering from project records. TODO documented in module. |
| **G1** (anatomical_de_analysis cross-pipeline injection) | RESOLVED — `@primitives/differential_expression.md` in `references:` frontmatter |
| **G2** (co_umap_embedding generalization) | RESOLVED — `atlas_co_umap.md` authored with 7 parameterized values, 4 risks addressed |
| **G3** (cross_atlas_dotplot generalization) | RESOLVED — `cross_dataset_dotplot.md` authored with 8 parameterized values, 4 risks addressed |

---

## Self-Check Grep Results

All three modules passed the self-check grep. Classification of all remaining hits:

| Term | Location | Classification |
|---|---|---|
| `STROMA_ORDER` | multi_group_de_analysis.md:151 | R comment in CONFIG block (`# TODO:`) — not executable code |
| `STROMA_ORDER` | multi_group_de_analysis.md:674, 679 | Phase 4 staging section — expected |
| `Vertebrae`, `Iliac Crest`, `Femoral Head` | multi_group_de_analysis.md:672 | Phase 4 staging section — expected |
| `CapEC`, `Fat EC` | atlas_co_umap.md:357 | Phase 4 staging section — expected |
| `TabulaSapiensComparison` | atlas_co_umap.md:355, 365; cross_dataset_dotplot.md:425, 439 | Phase 4 staging section — expected |
| `WCM`, `TS` | atlas_co_umap.md:354 | Phase 4 staging section — expected |
| `Tabula Sapiens` | atlas_co_umap.md:42; cross_dataset_dotplot.md:38 | Known-Atlas Convention prose header — required by task |
| `SEC1_FORCE` | cross_dataset_dotplot.md:433 | Phase 4 staging section — expected |
| `N_PER_SUBTYPE` | atlas_co_umap.md:358 | Phase 4 staging section — expected |
| `Source Sans 3` | atlas_co_umap.md:80, 366; cross_dataset_dotplot.md:440 | Risk documentation prose and Phase 4 staging — expected |
| `Playfair Display` | cross_dataset_dotplot.md:440, 443 | Phase 4 staging — expected |

**Zero hits in executable R code** for any prohibited term.

---

## pdf()/dev.off() Exceptions Used

| Module | Usage | CONVENTIONS.md §4 exception |
|---|---|---|
| multi_group_de_analysis.md | `pdf(...); ComplexHeatmap::draw(ht); dev.off()` (overall and functional heatmaps) | Exception #1 — ComplexHeatmap::draw() |
| cross_dataset_dotplot.md | `pdf(...); ComplexHeatmap::draw(ht); dev.off()` (centroid correlation heatmap) | Exception #1 — ComplexHeatmap::draw() |

No non-excepted `pdf()/dev.off()` patterns used. `cairo_pdf` is passed as a device to
`ggsave()` (not a `pdf()/dev.off()` pattern) — this is permitted.

---

## Cross-File Findings

1. **`make_functional_dotplot` from primitive is N-group compatible.** The primitive's
   implementation already uses `facet_wrap(~ group_var, ncol = n_groups)` and pins section
   labels to `tail(levels(factor(dot_df$group_var)), 1)`. Setting group factor levels to
   `GROUPS` on the object before calling propagates the canonical ordering. No wrapper function
   was needed.

2. **`cairo_pdf` + `useDingbats` incompatibility** documented in both atlas_co_umap.md and
   cross_dataset_dotplot.md. The v1 cross_atlas_dotplot.md already documented this pitfall;
   v2 modules carry it forward.

3. **Column name whitespace→dot conversion** is a cross-cutting issue affecting any module
   that round-trips through CSV. The fix (`gsub("\\.", " ", ...)`) is documented in
   cross_dataset_dotplot.md and should be noted in any future module that writes/reads
   correlation matrices.

4. **`rlang::sym()` used for NSE in `dplyr` / `ggplot2` code.** Both atlas_co_umap.md and
   cross_dataset_dotplot.md use `!!rlang::sym(COL_NAME)` for column-variable indirection.
   This requires `rlang` to be available (it is — it's a tidyverse dependency). Primitive
   modules that use the `rename to fixed column name before pivoting` pattern avoid this; the
   choice between the two is context-dependent.

---

## Open Items for Group D and Phase 3

1. **STROMA_ORDER curation** — must be done before running `make_anatomical_dotplot` in the
   BoneMarrowStroma project. Human review of BoneMarrowStroma project records required.
   Suggested examples/ file: `examples/bonemarrow_3site_anatomical.md`.

2. **Known-Atlas Convention implementation** — the proposed `context/known_atlases` design
   described above is a Phase 3 deployment-agent concern. It requires:
   - A new `context/known_atlases.md` or YAML block in `validated_examples.yaml`
   - Deployment agent logic to recognize atlas names in briefs
   - Population of canonical parameters for Tabula Sapiens from TabulaSapiensComparison records

3. **`TabulaSapiensComparison` project TODO fields** — `validated_examples.yaml` has several
   `null` / TODO fields for this project (data paths, subtype_col, last_validated date).
   These should be filled before Phase 4 examples/ files are written.

4. **Multi-source depth correction** (cross_dataset_dotplot.md Part B) — the current
   implementation handles 2-dataset cases cleanly (split into inhouse/atlas, z-score each,
   combine). The comment documents the pattern for >2 sources (group_by source + gene, scale
   within each source). A future version could formalize this into a `run_depth_correction()`
   helper in a primitive.

5. **`make_anatomical_dotplot` vs `make_topgene_dotplot` convergence** — the two functions
   have overlapping capabilities. `make_anatomical_dotplot` adds `groups` + `label_order` +
   `label_colors` to the primitive's dotplot pattern. Future consideration: fold `label_order`
   and `label_colors` into the primitive's `make_topgene_dotplot` so only one function exists.
   Left to a future refactor — don't merge without confirming that all callers of the primitive
   version are compatible.

---

## Unexpected Findings

1. The v1 `cross_atlas_dotplot.md` uses `make_tf_diamond_plot()` referenced from
   `@primitives/visualization.md`. This function IS defined in the v2 primitives (visualization.md
   defines `make_tf_diamond_plot()`). The Group C module's TF diamond section (Part C) does not
   call this function — instead it patches the dotplot with `shape = 23`. Future consideration:
   refactor Part C to use `make_tf_diamond_plot()` if that primitive version is equivalent.

2. The v1 `cross_atlas_dotplot.md` output directory was switched from `output/` to `output2/`
   on 2026-03-24. The v2 module uses `OUTPUT_DIR` parameter (defaults to `output/cross_dataset_dotplot`).
   Projects that used `output2/` (like TabulaSapiensComparison) will set `OUTPUT_DIR` to
   `"output2"` in their examples/ file.

3. `make_go_functional_dotplot` in `multi_group_de_analysis.md` is a newly authored function
   not present in any primitive. If it proves useful across projects, it is a candidate for
   promotion to `@primitives/differential_expression.md` in a future refactor. For now it lives
   in the module.
