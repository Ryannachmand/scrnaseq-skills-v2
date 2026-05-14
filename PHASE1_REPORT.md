---
# Phase 1 Report — Foundation: Primitives and Context
# Written: 2026-05-14
# Author: Phase 1 migration agent
---

# Phase 1 Report

Human review document before Phase 2 begins. Covers what was created, parameterization decisions,
fix verification, self-check results, inventory discrepancies, and open items for Phase 2.

---

## 1. Files Created

| File | Lines | Provenance |
|---|---|---|
| `SKILL.md` | 96 | Authored in v2 Phase 1 (analogous to v1 SKILL.md but restructured) |
| `CONVENTIONS.md` | 210 | Authored in v2 Phase 1 (no v1 equivalent) |
| `PROJECT_CLAUDE_TEMPLATE.md` | 29 | Migrated from v1 PROJECT_CLAUDE_TEMPLATE.md; updated @reference |
| `primitives/r_environment.md` | 61 | Migrated from v1 shared/r_environment.md — F1 fix applied |
| `primitives/seurat_v5_rules.md` | 102 | Migrated from v1 shared/seurat_v5_rules.md |
| `primitives/file_downloads.md` | 31 | Migrated from v1 shared/file_downloads.md |
| `primitives/geo_download.md` | 34 | Migrated from v1 shared/geo_download.md |
| `primitives/aesthetics.md` | 165 | Migrated from v1 shared/aesthetics.md |
| `primitives/visualization.md` | 351 | Migrated from v1 LargeDataset/methods/aesthetics.md — F5 fix applied |
| `primitives/differential_expression.md` | 743 | Migrated from v1 LargeDataset/methods/differential_expression.md — F2/F3/F4 fixes applied |
| `primitives/harmony_integration.md` | 87 | Authored in v2 Phase 1 (G8 generalization candidate — no v1 file) |
| `primitives/aucell_scoring.md` | 110 | Migrated from v1 metabolic_profile.md lines 62–86 — F6 fix applied |
| `context/lab_context.md` | 139 | Migrated from v1 lab_context.md — project-specific content stripped |
| `context/color_palettes.md` | 144 | Authored in v2 Phase 1 (no v1 equivalent) |
| `context/validated_examples.yaml` | 298 | Migrated from v1 lab_context/validated_examples.yaml — schema enriched |
| **Total** | **2,600** | |

---

## 2. Parameterization Decisions

### primitives/r_environment.md

| v1 hardcoded | v2 parameter / change |
|---|---|
| `conda run -n r-env Rscript` | `conda run --no-capture-output -n r-env Rscript` (F1 fix) |
| Incomplete package list | Expanded: added AUCell, ComplexHeatmap, circlize, CellChat, clusterProfiler, etc. |

### primitives/visualization.md (F5)

| v1 hardcoded | v2 parameter |
|---|---|
| `patient_id` column (proportion plot) | `group_col` argument |
| `cell_fraction` column | `subtype_col` argument |
| `source_file` column | `sample_col` argument |
| `ec_colors` | `group_colors` argument (named vector, caller-provided) |
| `tissue_colors` | `group_colors` argument |
| `adipose_type_colors` | `group_colors` argument |
| `"EC Functional Dot Plot"` function | `make_canonical_dotplot()` — no EC assumptions |
| Stacked violin EC-subset size comment `7.85` | Removed (HumanFat-specific) |
| Age/BMI/sex metadata track hardcoded | `meta_tracks` argument (optional named list) |
| `annotate("text")` for section labels | `geom_text(data=section_label_df)` pinned to one group |

### primitives/differential_expression.md (F2/F3/F4)

| v1 hardcoded | v2 parameter / change |
|---|---|
| `is_ambient()` regex-only version (anatomical_de_analysis.md) | Unified to full dual-filter (F2) |
| `is_confound()` phantom — never defined | Authored as real function (F3) |
| `functional_gene_sets` EC biology list hardcoded | Placeholder comment; caller must provide (F4) |
| `ec_sub$mylabel` hardcoded in heatmap annotation | `label_col` argument |
| `ec_colors` in heatmap column annotation | `subtype_colors` argument |
| `tissue_colors` in heatmap column annotation | `group_colors` argument |
| `expr$mylabel` in dot plot | `expr$subtype_var` (fixed column name to avoid !!sym failures) |
| `ec_colors[ec_order]` for x-axis colors | `subtype_colors[levels(factor(...))]` |
| `comp$col` references assume tissue_type column | `group_col` argument; `comp$col` preserved in comparisons list |

### primitives/harmony_integration.md (new)

No v1 parameterization needed — this is a new file. v1 duplicated this recipe
in two pipeline files (LargeDataset Stage 3, IntegratePublicData Stage 6) with
slight parameter drift. v2 provides one canonical reference.

### primitives/aucell_scoring.md (F6)

| v1 hardcoded | v2 parameter / change |
|---|---|
| OXPHOS gene set with invalid range syntax `"NDUFA1"-"NDUFA10"` | Not carried forward |
| FA_Synthesis gene set with invalid range syntax `"ELOVL1"-"ELOVL7"` | Not carried forward |
| PPARG target signature | Not carried forward |
| `ec_subtype_colors` referenced in output code | Not carried forward (primitive is pure compute) |
| Gene sets defined in the same file | `gene_sets` argument (caller-provided) |
| `tissue_type != "Myelolipoma"` filter | Not carried forward |
| Script names referencing HumanFat | Not carried forward |

### context/lab_context.md

| v1 content removed | Reason |
|---|---|
| Per-project notes (fat-scrnaseq-continued, KidneyNew, EyePublicData) | Moved to validated_examples.yaml |
| EC subtype color palettes (ec_subtype_colors_adipose, ec_subtype_colors_kidney) | Moved to validated_examples.yaml |
| adipose_type_colors, organ_tissue_colors, celltype_colors | Moved to validated_examples.yaml |
| Commented-out bone marrow stroma entry | Removed (D3 deletion candidate) |
| Commented-out unified_label_vocabulary and known_aliases | Removed (D4 deletion candidate) |

### context/validated_examples.yaml

| v1 schema field | v2 schema field |
|---|---|
| `job_id` | `project_name` (key) |
| `description_summary` | Incorporated into `notes` |
| `pipeline_variant` | Incorporated into `inputs.notes` |
| `geo_ids` | `inputs.paths` (typed) |
| `non_default_params` | `context_defaults` (richer structure) |
| `key_gotchas` | Consolidated into `notes` |
| `notes` | `notes` (expanded) |
| (missing) | `status`, `validated_outputs`, `last_validated`, `context_defaults.palettes` |

---

## 3. F1–F6 Fix Verification

### F1 — conda run --no-capture-output
**Applied:** `primitives/r_environment.md` line 13.
The command reads: `conda run --no-capture-output -n r-env Rscript <script.R> > <logfile> 2>&1`
A CRITICAL warning is included explaining what happens without the flag (hanging).
The v1 command (`conda run -n r-env Rscript`) is not present anywhere in v2 primitives.

### F2 — is_ambient() full dual-filter
**Applied:** `primitives/differential_expression.md`.
The function checks BOTH `AMBIENT_PATTERNS` (14 regex patterns) AND `AMBIENT_EXPLICIT`
(16 explicit gene names: ALB, FGA, FGB, FGG, FGL1, PF4, PPBP, GP9, GP1BA, GP1BB, ITGA2B,
MALAT1, NEAT1, TPSAB1, TPSB2, CPA3). The patterns-only version from anatomical_de_analysis.md
is explicitly NOT carried forward. A header comment marks this as F2 fix.

### F3 — is_confound() authored
**Applied:** `primitives/differential_expression.md`.
The function is newly authored with:
- All AMBIENT_PATTERNS (14 regex patterns, superset of is_ambient)
- CONFOUND_PATTERNS adds: `^HLA-D[RPQ]`, `^HIST[0-9]`, `^ENSG[0-9]+`, `^ENSMUSG[0-9]+`,
  `^AC[0-9]+\.[0-9]+`, `^AL[0-9]+\.[0-9]+`, `^AP[0-9]+\.[0-9]+`, `^LINC[0-9]+`
- CONFOUND_EXPLICIT adds: XIST, RPS4Y1, DDX3Y, EIF1AY, KDM5D, NLGN4Y, RPS4Y2, USP9Y, UTY, ZFY
- File header documents: is_confound is broader (comprehensive DE sweeps); is_ambient is narrower
  (standard volcano/DE labeling). This answers coverage_gaps.md §2.1.

### F4 — functional_gene_sets not hardcoded
**Applied:** `primitives/differential_expression.md`.
The v1 EC-specific gene sets (Adhesion & Immune Trafficking, Signaling & Angiocrine,
Extracellular Matrix, Metabolic & Specialized, ~80 genes) are NOT in this primitive.
The file contains a placeholder comment block explaining that functional_gene_sets must
be caller-provided, with a YAML brief schema showing where to define it. All downstream
functions (make_functional_heatmap, make_functional_dotplot, make_volcano) accept
functional_gene_sets as an explicit argument and document it as REQUIRED/optional accordingly.

### F5 — visualization.md sanitized of HumanFat assumptions
**Applied:** `primitives/visualization.md`.
Verified clean:
- Proportion plot: `patient_id`, `cell_fraction`, `source_file` → `group_col`, `subtype_col`, `sample_col`
- `ec_colors` → `group_colors` argument
- `tissue_colors` → `group_colors` argument
- `adipose_type_colors` → `group_colors` argument
- "EC Functional Dot Plot" → `make_canonical_dotplot()` with shape=21, no EC assumptions
- EC size comment in stacked violin (7.85) removed
- Stacked violin `make_stacked_violin()` takes generic `label_col`, `gene_col`, `value_col` args
- `annotate("text")` pattern replaced with `geom_text` + pinned data frame (Critical Constraints)

### F6 — aucell_scoring.md no hardcoded gene sets / no invalid R syntax
**Applied:** `primitives/aucell_scoring.md`.
The invalid range syntax (`"NDUFA1"-"NDUFA10"`, `"ELOVL1"-"ELOVL7"`, `"PLIN1"-"PLIN5"`)
is NOT in the primitive — these gene sets are not carried into v2 primitives at all.
The primitive provides `run_aucell(seurat_obj, gene_sets, assay)` where `gene_sets`
is a caller-provided named list. The PPARG target signature is also not in the primitive.
A note documents that gene sets belong in examples/ (Phase 4).

---

## 4. Self-Check Results

Each primitives file was checked for: project-specific column names, hardcoded color vectors,
hardcoded gene sets, use of `pdf()/dev.off()` instead of `ggsave()`.

| File | Project-specific columns | Hardcoded colors | Hardcoded gene sets | pdf()/dev.off() |
|---|---|---|---|---|
| r_environment.md | None | None | None | N/A |
| seurat_v5_rules.md | None | None | None | N/A |
| file_downloads.md | None | None | None | N/A |
| geo_download.md | None | None | None | N/A |
| aesthetics.md | None | References color_palettes.md | None | Rule stated |
| visualization.md | None (all args) | PALETTE_EXPRESSION inline constants only | None | All ggsave |
| harmony_integration.md | None (all args) | None | None | N/A |
| aucell_scoring.md | None (all args) | None | None | N/A |
| differential_expression.md | None (all args) | COLOR_UP_IDENT1/2 constants | None (placeholder) | Note: ComplexHeatmap uses pdf()/draw() — documented |

**One intentional exception:** `make_overall_heatmap()` and `make_functional_heatmap()` use
`pdf()` + `ComplexHeatmap::draw()` for output. This is unavoidable — ComplexHeatmap's `draw()`
cannot be piped into `ggsave()`. The comment in the primitive documents this as intentional.
This exception applies only to ComplexHeatmap outputs; all ggplot functions use `ggsave()`.

**Color constants in visualization.md:** The canonical palette color values
(e.g., `c("#F5F5F5", "#FFF9C4", "#FFB300", "#E53935")`) appear in the ggplot scale calls.
These are palette constants (not project-specific colors) and are acceptable in the primitive.
They correspond to the entries in `context/color_palettes.md`.

---

## 5. Discrepancies Between Inventory and v1 Source Files

### D1 — AMBIENT_EXPLICIT list count
The inventory (function_index.md Table 3.1 and SUMMARY.md) mentions "13 specific genes"
in AMBIENT_EXPLICIT. Direct reading of v1 differential_expression.md lines 94–99 shows
**16 genes**: ALB, FGA, FGB, FGG, FGL1, PF4, PPBP, GP9, GP1BA, GP1BB, ITGA2B,
MALAT1, NEAT1, TPSAB1, TPSB2, CPA3. The task brief specifies 16. v2 uses 16 (correct).

### D2 — Package list incompleteness
v1 r_environment.md listed only: Seurat, harmony, ggplot2, dplyr, tidyr, patchwork, scales,
monocle3, SingleCellExperiment, cowplot, ggrepel. The inventory confirmed ~12 additional
packages are used across methods files but unlisted. v2 r_environment.md adds the missing
packages based on methods file usage.

### D3 — metabolic_profile.md gene sets
v1 metabolic_profile.md shows **three** gene set families with invalid range syntax:
OXPHOS, FA_Synthesis, AND PPARG (PLIN1-PLIN5 range syntax also invalid). The inventory
mentioned only OXPHOS and FA_Synthesis. All three were excluded from v2 primitive.

### D4 — validated_examples.yaml v1 projects
The task brief specified migrating: HumanFat, KidneyNew, BoneMarrowStroma, CoCulture,
eye-test-01, TabulaSapiensComparison. The v1 validated_examples.yaml contained:
eye-test-01 (stub), fat-scrnaseq-continued, KidneyNew, EyePublicData.
BoneMarrowStroma, CoCulture, and TabulaSapiensComparison were NOT in v1 validated_examples.yaml
— they are documented in per_file_inventory.md and SUMMARY.md. v2 creates entries for all
six projects (eye-test-01 excluded per D2 deletion candidate); BoneMarrowStroma/CoCulture/
TabulaSapiensComparison entries have null values where information was unavailable.

### D5 — is_confound() sex-linked gene list
The task brief listed 10 sex-linked genes for is_confound(). The coverage_gaps.md §2.1
mentions a shorter list (XIST, UTY, DDX3Y, EIF1AY, KDM5D, USP9Y). v2 implements
all 10 from the task brief (adds NLGN4Y, RPS4Y1, RPS4Y2, ZFY to the coverage_gaps list).

---

## 6. Open Items for Phase 2

### P2-1. STROMA_ORDER is undefined
`STROMA_ORDER` is referenced in v1 `anatomical_de_analysis.md` lines 68, 107 but never defined
in any library file. Phase 2 must establish the canonical BoneMarrowStroma cell type ordering
before the anatomical DE module can work. Likely lived in a project-specific CLAUDE.md.
**Action needed:** Human provides the correct cell type ordering for BoneMarrowStroma.

### P2-2. IntegratePublicData pipeline modules need injection of DE functions
The v1 orphan files (`anatomical_de_analysis.md`, `co_umap_embedding.md`, `cross_atlas_dotplot.md`)
are NOT injected by the v1 IntegratePublicData pipeline.md. Phase 2 modules must:
- Create `modules/multi_group_de_analysis.md` (G1 generalization candidate)
- Wire the module into a conditional IntegratePublicData pipeline stage
- Ensure `run_findmarkers()`, `make_volcano()` etc. from the DE primitive are accessible

### P2-3. de_comprehensive_csv.md now has is_confound() defined
v2 adds `is_confound()` in `primitives/differential_expression.md`. When Phase 2 creates
the de_comprehensive_csv module, it should replace the phantom `is_confound()` reference
with the real function from the primitive. No further authoring needed — just wire it in.

### P2-4. make_topgene_dotplot name collision
v1 has two `make_topgene_dotplot` implementations with incompatible signatures
(function_index.md Table 3.2). v2 primitives use `make_topgene_dotplot` for the general
2-group version. Phase 2 module for anatomical/multi-site DE should use `make_anatomical_dotplot`
(as recommended by function_index.md) to avoid the name collision.

### P2-5. IntegratePublicData brief_template missing
v1 has no `IntegratePublicData/brief_template.txt` (coverage_gaps.md §1.1). Phase 2 must
create this. Until then, the CONVENTIONS.md brief schema covers both pipelines.

### P2-6. CellChat projectData() bug still in v1
v1 `interactome_cellchat.md` line 71 contains `projectData(cellchat, PPI.human)` which
was removed in CellChat v2 (D1 deletion candidate). When Phase 2 creates the CellChat module,
this call must be deleted from the inference pipeline skeleton.

### P2-7. celltype_subclustering.md uses pdf()/dev.off()
v1 `celltype_subclustering.md` lines 108–110 uses `pdf() + dev.off()` for Phase 1 diagnostic
plots instead of `ggsave()`. Phase 2 must fix this when creating the subclustering module.

### P2-8. pySCENIC machine-specific paths
v1 `pyscenic_regulon_analysis.md` contains three hardcoded absolute paths on a specific machine:
`/media/david/Mayo/ryan/scRNAseq/Sinusoid/GeneLists/`, `/home/ryannachman/scenicplus/...`,
`/home/ryannachman/anaconda3/envs/scenicenv/bin/python`. Phase 2 pySCENIC module must
replace these with `project_specific` sentinels or configurable brief fields.

### P2-9. Load formats RDS pattern uses deprecated Seurat v5 call
v1 `load_formats.md` line 185–186 uses direct slot assignment for Ensembl ID rename:
`rownames(so@assays$RNA@counts) <- sym`. This is deprecated in Seurat v5 and should
use `RenameFeatures()`. Phase 2 must fix when creating the load_formats module.

### P2-10. Doublet detection gap
Neither pipeline has a doublet detection step (coverage_gaps.md §5.4). Phase 2 should
add a DoubletFinder or scDblFinder option to the LargeDataset QC stage.

---

## 7. Anything Unexpected

### U1 — EyePublicData not in task brief's project list
The task brief listed: HumanFat, KidneyNew, BoneMarrowStroma, CoCulture, eye-test-01,
TabulaSapiensComparison. EyePublicData exists in v1 validated_examples.yaml but was not
in the task brief's list. EyePublicData was included in v2 validated_examples.yaml because
it is the most complete public-data integration project in the library with well-documented
gotchas. Human review: should EyePublicData be kept or archived?

### U2 — NKXSpleen project referenced but not in registry
v1 `bulk_lfc_concordance_heatmap.md` was "validated NKXSpleen run" and `pyscenic_regulon_analysis.md`
references NKXSpleen project structure. NKXSpleen does not appear in v1 validated_examples.yaml
or any registry. It is referenced only in methods files as a validation context.
Phase 2/4: consider adding NKXSpleen as an examples/ entry or validated_examples.yaml entry
if it represents a distinct completed project.

### U3 — make_functional_dotplot section label pinning
The v1 implementation pins section labels to `comp$ident2` (the right panel in a 2-panel facet).
v2 implements this as pinning to `tail(levels(factor(dot_df$group_var)), 1)` (last factor level),
which is equivalent. However, the exact yintercept calculation for section label midpoints
uses `as.numeric(factor(gene, levels = rev(gene_order)))` which may differ slightly from
v1's implementation depending on how genes are ordered after the sig_df filter.
This should be validated against a real run in Phase 3/4.

### U4 — color_palettes.md diverging palette recommendation
The `PALETTE_DIVERGING_6` values (`c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#D6604D","#B2182B")`)
match v1 `shared/aesthetics.md`. However, `metabolic_profile.md` in v1 uses a different
3-stop palette (`c("#4575B4","#FFFFBF","#D73027")`) for AUCell UMAPs. This is genuine drift
(duplication_report.md §3.1). v2 consolidates to the 6-stop palette as canonical. The
metabolic profile module in Phase 2 should use `PALETTE_DIVERGING_6` rather than the
blue-yellow-red variant.
