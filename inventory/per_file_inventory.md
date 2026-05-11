# Per-File Inventory — ~/claude-skills/
*Generated: 2026-05-11 by inventory agent*

---

## Ordering: shared/ → LargeDataset/ → IntegratePublicData/ → lab_context/

---

## shared/aesthetics.md
**Lines:** 180  
**Status:** SHARED_INFRASTRUCTURE

**Summary:** Defines generalizable visual *principles* for all pipelines — typography tables, color philosophy, continuous/categorical palettes, layout/spacing rules, figure sizing formulas for four plot types, and UMAP feature plot standards. Contains no executable R code whatsoever; every entry is a rule or table.

**Functions defined:** None

**Functions called from outside this file:** None

**Project-specific hard dependencies:** None. All content is abstract principles. Color values listed (`#F5F5F5` → `#E53935`, `#2166AC` → `#B2182B`, etc.) are principles, not project bindings.

**Generalizable content:** Everything. This is the most generalizable file in the library.

**Project-specific content:** None.

**Content that duplicates other files:**
- Typography tables duplicate/conflict with specific sizes embedded inline in `LargeDataset/methods/aesthetics.md` (e.g., shared says axis text 13–14pt; LargeDataset/aesthetics hardcodes 9.85pt for EC violin y-axis). These are instances of the same principle at different specificity levels.
- Color palette values (gene expression `#F5F5F5`→`#E53935`, TF diverging blue/red, differential abundance `steelblue`/`firebrick`) are repeated verbatim in `LargeDataset/methods/aesthetics.md`, `differential_expression.md`, `cross_atlas_dotplot.md`, and `cohort_plots.md`.
- Figure sizing principles are repeated and specialized in every methods file.

**Inconsistencies with other files:**
- States "For plot-specific implementations, see each pipeline's `methods/aesthetics.md`" (line 4) but only `LargeDataset/` has a `methods/aesthetics.md`; `IntegratePublicData/` has none.
- UMAP feature plot section (lines 122–146) overlaps with but is less specific than `feature_umap_plot.md`; neither references the other.
- Proportion plot sizing formula (`max(10/1.65, n_samples*0.38/1.65+3.5)`) matches `LargeDataset/methods/aesthetics.md` but the formula origin is not documented.

---

## shared/seurat_v5_rules.md
**Lines:** 95  
**Status:** SHARED_INFRASTRUCTURE

**Summary:** Six permanently-validated Seurat v5 bug fixes applied in every generated script: JoinLayers after merge, ScaleData on variable features only, metadata assignment via `@meta.data$col <-` not `$col <-`, external CSV metadata merge pattern, dynamic UMAP reduction name detection, and matrix conversion for log2OR.

**Functions defined:** None (code snippets, not named functions)

**Functions called from outside this file:** None

**Project-specific hard dependencies:** None.

**Generalizable content:** Everything. All six rules are universal Seurat v5 behavior corrections.

**Project-specific content:** None.

**Content that duplicates other files:**
- JoinLayers rule repeated inline in `load_formats.md` (line 159), `IntegratePublicData/pipeline.md` Stage 5 (line 264), `cohort_plots.md` (line 111), `pyscenic_regulon_analysis.md` (line 84).
- Dynamic UMAP name detection code duplicated in `feature_umap_plot.md` (lines 52–54).
- `@meta.data$col <-` rule referenced (not repeated verbatim) in `label_harmonization.md` (lines 177–181) and `celltype_subclustering.md` (line 288).

**Inconsistencies:** The validated date is 2026-04-10 ("TFEC Expression Atlas run") but the file header rule 5 (dynamic UMAP detection) references `Seurat:::AutoPointSize` usage pattern documented more fully in `feature_umap_plot.md`. Minor: Rule 6 says "never `tibble::column_to_rownames()`" but does not note that `tibble` IS available in r-env (only listed packages don't include it in `r_environment.md`). Unclear if tibble is actually absent.

---

## shared/r_environment.md
**Lines:** 37  
**Status:** SHARED_INFRASTRUCTURE

**Summary:** Mandates use of `r-env` conda environment for all R execution. Documents why (sp package ABI incompatibility), lists available packages, and specifies script structure conventions (self-contained, gc() calls, versioned RDS saves, ggsave with useDingbats=FALSE).

**Functions defined:** None

**Functions called from outside this file:** None

**Project-specific hard dependencies:** Machine-specific: `/home/ryannachman/` implied by conda environment location. The r-env name is machine-specific.

**Generalizable content:** Script structure conventions (lines 28–37) are general. Package list is general.

**Project-specific content:** The conda command itself is machine-specific.

**Content that duplicates other files:**
- The "always use r-env" rule is re-stated in almost every methods file header and in the Critical Constraints tables of `differential_expression.md` (line 36), `bulk_concordance.md`, `celltype_subclustering.md`, `metabolic_profile.md`, and `pyscenic_regulon_analysis.md`.
- `ggsave(useDingbats=FALSE)` rule duplicated in `shared/aesthetics.md` File Output section and multiple methods files.

**Inconsistencies:**
- States "SeuratWrappers NOT installed" but does not mention AUCell, circlize, ComplexHeatmap, CellChat, clusterProfiler, org.Hs.eg.db, matrixStats, data.table, readxl — all used in methods files. Package list is incomplete.
- Does not mention the `scenicenv` Python environment used by `pyscenic_regulon_analysis.md`.
- The CLAUDE.md for this project specifies `conda run --no-capture-output -n r-env Rscript` (with `--no-capture-output`) but this file shows `conda run -n r-env Rscript` (without it). The CLAUDE.md flag is the validated form.

---

## shared/file_downloads.md
**Lines:** 24  
**Status:** SHARED_INFRASTRUCTURE

**Summary:** Describes how to convert Box.com shared-link URLs from `/s/<token>` (preview) to `/shared/static/<token>` (direct download) for use with `curl -L`.

**Functions defined:** None

**Functions called from outside this file:** None

**Project-specific hard dependencies:** None (URL pattern is general to Box.com).

**Generalizable content:** Everything.

**Project-specific content:** None.

**Content that duplicates other files:** None.

**Inconsistencies:** None.

---

## shared/geo_download.md
**Lines:** 18  
**Status:** SHARED_INFRASTRUCTURE

**Summary:** Defines priority order for GEO data acquisition (processed matrix > per-sample files > supplementary > FASTQ never), efficiency rules (inspect before downloading, prefer merged matrices, flag >20GB), and prohibits FASTQ downloads without approval.

**Functions defined:** None

**Functions called from outside this file:** None

**Project-specific hard dependencies:** None.

**Generalizable content:** Everything.

**Project-specific content:** None.

**Content that duplicates other files:**
- GEO download discussion in `load_formats.md` section "GEO Bulk Supplementary Files" (lines 238–310) provides far more detail on specific GEO file formats than this file but doesn't repeat its priority-order rules.

**Inconsistencies:** None.

---

## pipelines/LargeDataset/pipeline.md
**Lines:** 222  
**Status:** SHARED_INFRASTRUCTURE

**Summary:** Nine-stage checkpoint-driven pipeline definition for multi-sample scRNAseq from .h5 files. Covers assemble → QC → process → annotate → visualize → subset → core analysis → optional extras → complete. Defines five checkpoints and a brief parameters reference table.

**Functions defined:** None (pipeline logic, not R functions)

**Functions called:** None directly; references method files via Stage 8 `optional_analyses` keys.

**Project-specific hard dependencies:**
- `source_file` appears in Brief Parameters Reference as the default `batch_correction_var` — this is from the HumanFat project (fat-scrnaseq-continued).
- `metadata_display_cols` default includes `tissue_type`, `additional_notes` — HumanFat-specific column names.
- Brief template defaults (`n_variable_features: 4000`, `n_pcs: 30`, resolution 0.5) come from validated HumanFat parameters.

**Generalizable content:** The 9-stage structure, checkpoint pattern, and optional_analyses key mapping are fully general. All processing steps are standard Seurat pipeline.

**Project-specific content:** Default column names (`source_file`, `tissue_type`, `additional_notes`). These override lab_context.md defaults implicitly.

**Content that duplicates other files:**
- Stage 3 Seurat processing sequence duplicated in `IntegratePublicData/pipeline.md` Stage 6.
- UMAP recipe in Stage 5 description partially duplicates `LargeDataset/methods/aesthetics.md`.
- Stage 8 interactome entry says "(🔲 placeholder)" — this is a stale placeholder marker indicating the CellChat integration was pending. `interactome_cellchat.md` now exists but the pipeline.md was not updated.

**Inconsistencies:**
- Stage 8 `interactome` entry: "🔲 placeholder" contradicts SKILL.md which shows `interactome_cellchat.md` as "✅ validated KidneyNew run".
- Brief Parameters says `batch_correction_var` default is `source_file` but `lab_context.md` sets `batch_correction_var: sample_id`. These conflict.

---

## pipelines/LargeDataset/brief_template.txt
**Lines:** 166  
**Status:** SHARED_INFRASTRUCTURE

**Summary:** YAML template for LargeDataset project briefs with every parameter field annotated as REQUIRED, CONFIRM, or OPTIONAL. Covers project identity, data source, metadata, processing parameters, QC thresholds, annotation, subset configuration, core analysis, and optional_analyses block.

**Functions defined:** None

**Project-specific hard dependencies:**
- Default `batch_correction_var: source_file` is HumanFat-specific (see conflict with lab_context.md).
- `metadata_columns` schema (`patient_id`, `sex`, `age`, `bmi`, `diagnosis`, `additional_notes`) reflects HumanFat metadata structure.
- `metadata_display_cols` default list includes `tissue_type`, `additional_notes`, `age`, `bmi`, `sex` — HumanFat metadata columns.
- `deg` block `group_var: tissue_type` — HumanFat-specific field name.

**Generalizable content:** All structural elements. The field annotations (REQUIRED/CONFIRM/OPTIONAL) are highly useful.

**Project-specific content:** Column name defaults, metadata schema. These should be cleared/generalized for v2.

**Content that duplicates other files:** The optional_analyses block partially mirrors `differential_expression.md` brief configuration section (lines 778–807).

**Inconsistencies:** `batch_correction_var: source_file` conflicts with `lab_context.md`'s `batch_correction_var: sample_id`.

---

## pipelines/LargeDataset/methods/aesthetics.md
**Lines:** 207  
**Status:** ACTIVELY_USED (referenced by pipeline.md Stages 5 and 7)

**Summary:** Plot-specific R code for six visualization types: UMAP (DimPlot/LabelClusters), stacked violin, cell-type proportion plot, EC functional dot plot, EC TF diamond plot, and differential abundance heatmap plus trajectory plots. Contains complete ggplot code with validated aesthetics from HumanFat_Yang run.

**Functions defined:** None (code recipes, not named functions)

**Functions called from outside this file:**
- `LabelClusters()` — Seurat function, called from pipeline.md Stage 5 description
- The complete dot plot and TF plot recipes here are what `differential_expression.md` functions `make_topgene_dotplot` and `make_functional_dotplot` implement — partial conceptual overlap

**Project-specific hard dependencies:**
- Cell Type Proportion Plot (lines 69–96): references `patient_id`, `cell_fraction`, `source_file` as metadata column names — HumanFat-specific.
- The proportion plot metadata track colors (age gradient `#FFF5EB→#FD8D3C→#7F2704`, BMI blue `#EFF3FF→#6BAED6→#084594`, sex `#E07A9F`/`#45B7D1`) match HumanFat metadata fields.
- Stacked violin y-axis comment "For EC subset (more subtypes), use size = 7.85 instead" — HumanFat EC-specific.
- Section label positioning using `annotate("text", ...)` in EC dot plot (line 122) is correct for non-faceted dot plots but inconsistent with `differential_expression.md` which forbids `annotate("text")` in faceted contexts.

**Generalizable content:**
- UMAP recipe (lines 1–23) — fully general.
- Stacked violin recipe (lines 27–66) — mostly general; HumanFat EC comment should be removed.
- Differential abundance heatmap (lines 170–190) — fully general.
- Trajectory plot recipes (lines 192–207) — fully general.
- Dot plot and TF plot structural patterns (shape, color scale, size range, stroke) — general; only the section label annotations using `annotate()` are non-faceted-specific.

**Project-specific content:** Proportion plot (lines 69–96) metadata column references. Stacked violin EC size note.

**Content that duplicates other files:**
- Differential abundance heatmap code (lines 170–190) is semantically the same operation described in `differential_expression.md` Step 3 header; they target different specific heatmaps but use the same `scale_fill_gradient2(low="steelblue",high="firebrick")`.
- UMAP aesthetics partially duplicated in `co_umap_embedding.md` panel design.
- Dot plot recipe here is the "EC Functional Dot Plot" — a related but simpler version appears in `anatomical_de_analysis.md` as `make_topgene_dotplot()` (different function signature and 3-panel vs 1-panel design).
- Trajectory plot recipe at lines 192–207 is verbatim identical to `trajectory_monocle3.md` lines 111–168 (different scope: this is a 2-plot summary, trajectory_monocle3.md has the full implementation).
- Color palettes repeated from `shared/aesthetics.md`.

**Inconsistencies:**
- EC dot plot section labels use `annotate("text", ...)` (line 122) but `differential_expression.md` lines 364–370 explicitly forbid `annotate("text")` in faceted contexts — these are different contexts (this file is non-faceted) but the inconsistency creates risk of misapplication.
- No file references this for IntegratePublicData, leaving that pipeline without any equivalent aesthetic code file.

---

## pipelines/LargeDataset/methods/differential_expression.md
**Lines:** 807  
**Status:** ACTIVELY_USED (referenced by pipeline.md Stage 8 `deg` key)

**Summary:** Complete DE visualization suite for between-group comparisons within a subset object. Defines configuration block, ambient filter with both PATTERNS and EXPLICIT lists, functional gene sets (EC-specific biology), and seven named functions: run_findmarkers, make_volcano, make_overall_heatmap, make_functional_heatmap, make_topgene_dotplot, make_functional_dotplot (with direction strip), make_pathway_barplot. Includes main loop structure, figure sizing reference, and brief configuration.

**Functions defined:**
- `is_ambient(genes)` — line 100. Checks genes against `AMBIENT_PATTERNS` regex list AND `AMBIENT_EXPLICIT` named list. Returns logical vector.
- `run_findmarkers(ec_obj, comp)` — line 187. Runs `FindMarkers()` with wilcox test, logfc.threshold=0 (returns all genes), min.pct=0.05; skips if <10 cells per group.
- `make_volcano(markers, comp, subset_name, n_label = N_LABEL_VOLCANO)` — line 221. Full ggplot volcano with ggrepel labels, significance coloring, ambient filter, biological gene prioritization.
- `make_overall_heatmap(ec_obj, markers, comp, subset_name, n_cells = N_CELLS_HM, n_genes = N_TOP_GENES_HM)` — line 285. ComplexHeatmap cell-level z-score heatmap with column split by group.
- `make_functional_heatmap(...)` — line 333 (described in Step 4). Uses `functional_gene_sets` + `row_split` for section gaps.
- `make_topgene_dotplot(...)` — line 415 (Step 5 description). % expressing + avg expression per EC subtype, faceted by tissue group.
- `make_functional_dotplot(..., show_direction = FALSE)` — line 483. Organizes dot plot into functional sections with optional patchwork direction strip.
- `make_pathway_barplot(markers, comp, subset_name, universe_genes, n_pathways = 15)` — line 594. GO:BP ORA with `enrichGO()` + `simplify()`, diverging bar plot.

**Functions called from outside this file:**
- All 7 functions above called from `anatomical_de_analysis.md` main loop (lines 316–333) — **UNRESOLVED**: anatomical_de_analysis.md does NOT include or source differential_expression.md.
- `is_ambient()` referenced in `de_comprehensive_csv.md` as `is_confound()` — different name, different file, possibly referring to the same concept (see Phantom Functions).

**Project-specific hard dependencies:**
- `functional_gene_sets` (lines 119–168): EC biology specific (Adhesion & Immune Trafficking, Signaling & Angiocrine, Extracellular Matrix, Metabolic & Specialized). Gene lists are curated for EC biology (PPARG, FABP4, KDR, NOTCH1, COL4A1, etc.).
- `ec_sub$mylabel` hardcoded in `make_overall_heatmap` column annotation (line 313) — HumanFat/KidneyNew EC label column name.
- `ec_colors` and `tissue_colors` used in heatmap annotations — project-specific color dictionaries not defined in this file.
- Configuration thresholds (PADJ_CUT=0.05, LFC_CUT=0.5, N_CELLS_HM=300) are reasonable defaults but embedded as constants.
- Volcano direction colors `#B2182B`/`#2166AC` match general diverging palette but are hardcoded.

**Generalizable content:**
- `run_findmarkers()`, `make_volcano()`, `make_overall_heatmap()`, `make_pathway_barplot()` are fully general (just need configuration block and ambient filter).
- Critical Constraints table (lines 29–37) is highly generalizable architectural guidance.
- Configuration block pattern (comparisons + ec_subsets as named lists) is general.
- `is_ambient()` with both PATTERNS and EXPLICIT is general (explicit list is human-generic contamination).

**Project-specific content:**
- `functional_gene_sets` — EC biology specific; must be replaced for non-EC contexts.
- `mylabel` column reference in heatmap annotation.
- Section label annotations assume 2-panel faceting (ident1/ident2).

**Content that duplicates other files:**
- `is_ambient()` definition duplicated with divergence in `anatomical_de_analysis.md` (lines 340–347).
- `make_functional_dotplot()` direction strip pattern is partially duplicated in `anatomical_de_analysis.md` direction strip (lines 186–209) with minor differences (3-site vs 2-site).
- `make_topgene_dotplot()` described in Step 5 but a different named implementation `make_topgene_dotplot()` exists in `anatomical_de_analysis.md` with different behavior.
- GO:BP ORA pattern repeated conceptually in `anatomical_de_analysis.md` `1.7_SiteDE_GOFunctionalDotplot.R` description.
- Ambient filter AMBIENT_PATTERNS list appears in both files.

**Inconsistencies:**
- `make_overall_heatmap` uses `ec_sub$mylabel` (line 313) — hardcoded column name, not parameterized.
- `make_topgene_dotplot` uses `ec_sub$mylabel` for x-axis (line 350) — same issue.
- `ec_colors` and `tissue_colors` referenced but not defined in this file (must come from the CLAUDE.md context block).

---

## pipelines/LargeDataset/methods/de_comprehensive_csv.md
**Lines:** 253  
**Status:** ACTIVELY_USED (referenced by SKILL.md as validated)

**Summary:** Vectorized Kruskal-Wallis multi-group statistics table using `matrixStats::rowRanks`. Produces a comprehensive CSV with KW statistics, per-group mean expression, top-group annotation, and log2FC vs others. Includes a cache pattern to avoid recomputation.

**Functions defined:**
- `run_vectorized_kw()` — not named explicitly; the vectorized KW implementation in Step 1 (lines 63–92) is a code block, not a named function. **Would benefit from being a named function.**
- `get_gene_fpkm()` — NOT in this file. (This is in `load_formats.md`.)

**Functions called from outside this file:**
- `is_confound(gene)` — called at line 116 in Step 2. **NEVER DEFINED ANYWHERE IN THE LIBRARY.** This is a phantom function. The intent is likely `is_ambient()` (defined in `differential_expression.md`) plus additional filters (sex-linked genes, HLA class II, histone genes, unannotated ENSG IDs — mentioned in the "Exclusion filters" section at line 215). **Critical gap.**

**Project-specific hard dependencies:**
- `SUBTYPE_LABELS` and `SUBTYPE_COL` at lines 168–172 — must be set per project (e.g., `c("Capillary","Venous","Arterial")` for HumanFat EC).
- `GROUP_COL` at line 171 — project-specific metadata column.
- Output filename `comprehensive_DE_stats.csv` is fixed (not parameterized).

**Generalizable content:** The entire statistical approach (vectorized KW, top-group annotation, log2FC weighting) is project-agnostic. The cache pattern is highly reusable.

**Project-specific content:** Subtype label column name references, group column name. Otherwise general.

**Content that duplicates other files:**
- Conceptual overlap with `is_ambient()` through the phantom `is_confound()` reference.

**Inconsistencies:**
- References `is_confound()` which is never defined. The "Exclusion filters" section (line 215) says "Apply `is_confound()` (defined in `differential_expression.md`)" — but `differential_expression.md` defines `is_ambient()`, not `is_confound()`. These may be the same function or `is_confound()` may be an extended version that was never written.
- Note at line 216 says "Exclusion filters" include "sex-linked genes, HLA class II, histone genes, unannotated ENSG IDs" — none of these are in `is_ambient()`. So `is_confound()` would need to be a superset — but it doesn't exist.

---

## pipelines/LargeDataset/methods/bulk_concordance.md
**Lines:** 759  
**Status:** ACTIVELY_USED (referenced by SKILL.md)

**Summary:** Three-part analysis cross-referencing in-vivo scRNA-seq with a bulk RNA overexpression experiment: (1) annotation overlay re-labeling volcano/dot/heatmap outputs with bulk concordance status, (2) per-cell AUCell module scoring showing which cells resemble the bulk signature, (3) TF-focused analysis with TF dot plots, TF markers, TF volcanos, and a TF LFC heatmap. Validated on HumanFat PPARG OE experiment.

**Functions defined:**
- Multiple unnamed helper blocks for: bulk-status annotation (`mutate(bulk_status = ...)`), module score logic, split violin recipe, TF marker extraction, TF LFC heatmap construction. None are named functions — all inline code.

**Functions called from outside this file:**
- `is_ambient()` — called in annotation overlay section; resolves to `differential_expression.md`.
- `make_volcano()`, `make_topgene_dotplot()`, `make_functional_dotplot()` — conceptually expected but not explicitly called; file re-implements annotation overlay versions of these rather than calling the originals.

**Project-specific hard dependencies:**
- `BULK_CSV` points to `DE_significant_Group4_vs_Group2.csv` (PPARG OE vs HUVEC control) — HumanFat specific.
- `ec_colors`, `mylabel` column, `tissue_type` column, `tissue_colors` — HumanFat/EC project specific.
- `RibHighEC` exclusion pattern — HumanFat specific cluster name.
- `PPARG` as the TF of interest in Part 3 — project biology specific.
- Bulk experiment described as "PPARG1_ROSI vs Control_ROSI HUVECs" — very specific experimental context.

**Generalizable content:**
- The three-part conceptual framework (annotation overlay + module score + TF focus) is generalizable to any bulk perturbation cross-reference.
- Configuration block pattern is good.
- Per-cell module score UMAP/violin recipe is general.
- Split violin recipe (comparing two groups with ggplot violin) is general.

**Project-specific content:** Substantial. The PPARG biology, HumanFat tissue types, specific bulk experiment columns, `RibHighEC` filter, and TF selection criteria are all project-specific. Generalizing this file requires parameterizing the bulk input, the cell subtype column, the tissue column, and removing any hardcoded PPARG references.

**Content that duplicates other files:**
- Module score approach (AddModuleScore → UMAP overlay → violin) shares conceptual overlap with `IntegratePublicData/pipeline.md` Stage 2 cell type pre-filtering.
- TF dot plot recipe partially overlaps with EC TF Diamond Plot in `LargeDataset/methods/aesthetics.md`.
- The distinction from `bulk_lfc_concordance_heatmap.md` is clearly documented within the file (lines 1–35) — good self-documentation.

**Inconsistencies:**
- `make_volcano()` function signature used in annotation overlay context may differ from `differential_expression.md` definition (same name, potentially different caller context).

---

## pipelines/LargeDataset/methods/bulk_lfc_concordance_heatmap.md
**Lines:** 583  
**Status:** ACTIVELY_USED (referenced by SKILL.md; validated NKXSpleen run)

**Summary:** Parallel log2FC heatmap comparing DESeq2 bulk LFC and scRNA-seq FindMarkers LFC side by side for the same gene set. Uses a two-stage caching architecture. Validated on NKXSpleen project (80k EC cells, 12 bulk samples). Clear distinction from `bulk_concordance.md` documented in the file.

**Functions defined:**
- `get_gene_fpkm()` — This function appears in `load_formats.md` (line 264), not here. This file has no named function definitions; uses inline code blocks with caching architecture.

**Functions called from outside this file:**
- `JoinLayers()`, `FindMarkers()` — standard Seurat.
- No custom library functions called.

**Project-specific hard dependencies:**
- NKXSpleen biology: target-organ EC subtypes (10 subtypes), 20 reference-organ EC types, bulk RNA-seq design (6 Control / 6 TF overexpression).
- Two-batch structure in DESeq2 formula (`~batch + condition`).
- Curated gene categories (chemokines, adhesion molecules, angiocrine factors, ECM, cell-type markers, TF program) — project biology specific.
- `na_col = "grey88"` rule for missing scRNA-seq genes — good general pattern.

**Generalizable content:**
- Two-stage caching architecture (compute once, iterate visualization) is an excellent general pattern.
- `pheatmap` NA handling rule (`na_col = "grey88"`) is universally applicable.
- The distinction between bulk LFC and scRNA-seq LFC (independent scales, never z-scored together) is an important methodological note for any similar analysis.

**Project-specific content:** Gene category curation, subtype labels, batch variable, specific experiment design. The figure is fundamentally publication-specific.

**Content that duplicates other files:**
- `bulk_concordance.md` distinction documented internally — good.
- Two-stage caching pattern similar to `de_comprehensive_csv.md` KW cache pattern.

**Inconsistencies:** None identified.

---

## pipelines/LargeDataset/methods/celltype_subclustering.md
**Lines:** 521  
**Status:** ACTIVELY_USED (referenced by SKILL.md; validated HumanFat Interactions run)

**Summary:** Two-phase subclustering methodology: Phase 1 produces exploration outputs (UMAP, feature plots, FindAllMarkers) for user to decide labels; Phase 2 applies user-provided labels and generates endpoint figures (labeled UMAP, proportions, curated dotplot). Documents Seurat v5 factor/unname gotchas.

**Functions defined:** None (code recipes, not named functions)

**Functions called from outside this file:** None

**Project-specific hard dependencies:**
- `source_file` column used throughout (HumanFat sample column name) — lines 55, 196, 204, 205, 209, 391–401.
- `tissue_type` column used throughout — lines 105, 157, 163, 165, 168, etc.
- `adipose_type_colors` referenced in diagnostic plot (line 105) — HumanFat color dictionary.
- Tissue colors defined inline as an adipose depot palette (lines 332–341: Subcutaneous Fat, Visceral Fat, Lipoma, Liposuction Fat, Breast Fat, Orbital Fat).
- Paul Tol muted colorblind palette hardcoded for cluster exploration (lines 190–191).
- Phase 1 HARMONY_BY hardcoded to `"source_file"` (line 55).
- CELLTYPE_LABELS example `c("Macrophage", "Monocyte")` — from HumanFat Interactions run.

**Generalizable content:**
- The two-phase exploration→endpoint framework is fully general.
- Seurat v5 gotchas table (lines 24–31) is general and valuable.
- Critical Constraints table has good general fixes (factor/unname issue).
- Cluster coloring with Paul Tol colorblind palette is general.
- FindAllMarkers → top100 CSV handoff pattern is general.
- Iteration patterns section (lines 509–521) is general.

**Project-specific content:** `source_file`, `tissue_type`, `adipose_type_colors`, `tissue_cols` with adipose depot values, HARMONY_BY defaults.

**Content that duplicates other files:**
- UMAP + proportion plot aesthetic code overlaps with `LargeDataset/methods/aesthetics.md`.
- `JoinLayers()` call (line 222) repeats `seurat_v5_rules.md` Rule 1.

**Inconsistencies:**
- Phase 1 code hardcodes `HARMONY_BY <- "source_file"` (line 55) but the general configuration variable should come from the brief.
- Phase 1 uses `pdf()` + `dev.off()` (lines 107–109, 140–142) for some plots instead of the mandated `ggsave(..., device="pdf", useDingbats=FALSE)` — inconsistent with `shared/r_environment.md` and `shared/aesthetics.md` conventions.

---

## pipelines/LargeDataset/methods/differential_expression.md *(see above)*

---

## pipelines/LargeDataset/methods/feature_umap_plot.md
**Lines:** 131  
**Status:** ACTIVELY_USED (referenced by SKILL.md; validated CoCulture run)

**Summary:** Reproduces Seurat FeaturePlot behavior exactly using ggplot. Documents seven specific gotchas discovered during CoCulture run (AutoPointSize formula, oob=squish, gradientn vs gradient, q95 limit compression, etc.). Includes verification code to compare against native FeaturePlot.

**Functions defined:** None (one complete inline template)

**Functions called from outside this file:**
- `Seurat:::AutoPointSize()` — internal Seurat function, called with `data.frame(x=seq_len(ncol(so)))` pattern documented here.

**Project-specific hard dependencies:** None beyond the general Seurat object structure.

**Generalizable content:** Everything. The gotcha table is the most valuable content (7 specific debugging lessons).

**Project-specific content:** None.

**Content that duplicates other files:**
- Dynamic UMAP reduction name detection (lines 52–54) is a subset of `seurat_v5_rules.md` Rule 5 — same code.
- Gene expression color scale `c("lightgrey","blue")` referenced in `shared/aesthetics.md` UMAP Feature Plots section (line 124) with less detail.

**Inconsistencies:**
- `shared/aesthetics.md` (line 116) says "UMAP point size: `pt.size = 0.3` (cell type / categorical); `0.8` (feature/expression plots)" but this file correctly shows the AutoPointSize approach produces pt.size ≈ 0.082 for ~19k cells — contradicting the shared principles file. The shared file's `0.8` advice is wrong for large datasets.

---

## pipelines/LargeDataset/methods/interactome_cellchat.md
**Lines:** 1302  
**Status:** ACTIVELY_USED (referenced by pipeline.md Stage 8; validated KidneyNew + HumanFat)

**Summary:** CellChat v2 complete implementation in a 6-script architecture (inference, basic plots, stacked bubbles, bar plots, circos pathway diagrams, tissue comparison). Documents all known v2 bugs. Covers pathway categorization, individual pathway circle/bubble/chord/heatmap/bar plots, and a multi-tissue comparison framework.

**Functions defined:**
- `get_filtered_net(cellchat, paths, min_cells = 10)` — helper described in Scripts 3–4 section; provides filtered interaction table for chord diagrams (bypassing pairLR.use bug). Definition is provided inline in the plotting scripts section.
- `plot_heatmap(cellchat, paths, ...)` — custom ggplot replacement for `netVisual_heatmap()` described in the Critical Constraints and implemented inline. Not a standalone named function; pattern shown as code block.

**Functions called from outside this file:**
- All CellChat v2 functions (`createCellChat`, `computeCommunProb`, `netVisual_circle`, `chordDiagram`, etc.) — external library.
- `plan("sequential")` — futures package, critical fix documented.

**Project-specific hard dependencies:**
- KidneyNew example categories hardcoded (lines 131–149): Chemokine_Adhesion, Angiocrine, ECM, Metabolic categories with specific pathway lists for kidney ECs.
- HumanFat chord color palettes tied to adipose-specific tissue/cell type colors.
- Script architecture references NKXSpleen and KidneyNew project directory structures.
- CellChat color schemes tied to project-specific color palettes.
- Multi-tissue comparison framework (Script 6) compares LiposuctionFat vs BreastFat — HumanFat specific.

**Generalizable content:**
- Critical Constraints table (6 bugs with fixes) — fully general, extremely valuable.
- 6-script architecture pattern — general.
- Pathway categorization strategy (lines 120–130) — general framework.
- `get_filtered_net()` helper — general.
- All CellChat API usage patterns — general for any CellChat v2 project.

**Project-specific content:**
- KidneyNew pathway category definitions (lines 131–180).
- Multi-tissue comparison logic assumes 2 specific tissue types.
- Color palette references.

**Content that duplicates other files:**
- `plan("sequential")` rule mentioned in `lab_context.md` KidneyNew notes and here — consistent.
- Critical v2 bug list overlaps with `lab_context.md` KidneyNew notes (lines 203–208).

**Inconsistencies:**
- Line 72: `cellchat <- projectData(cellchat, PPI.human)` is included in the inference pipeline skeleton, but Critical Constraints says to REMOVE this line ("Function removed in CellChat v2"). This is a direct contradiction within the same file — inference pipeline code still contains the invalid call that the constraints table says to remove.

---

## pipelines/LargeDataset/methods/metabolic_profile.md
**Lines:** 257  
**Status:** ACTIVELY_USED (referenced by SKILL.md; validated HumanFat MetabolicProfile run)

**Summary:** AUCell-based metabolic pathway scoring for 5 core pathways (Glycolysis, Beta_Oxidation, TCA_Cycle, OXPHOS, FA_Synthesis) plus PPARG target signature. Includes gene sets, AUCell helper functions, full output structure for HumanFat (6 scripts, 3 subdirectories), and all plot type recipes.

**Functions defined:**
- `run_aucell(seurat_obj, gene_sets, assay = "RNA")` — line 62. Builds AUCell rankings from counts and returns AUC scores as data frame.
- `add_auc_to_seurat(seurat_obj, auc_df)` — line 70. Adds AUC scores as metadata columns to a Seurat object.

**Functions called from outside this file:** None

**Project-specific hard dependencies:**
- Gene sets: FA_Synthesis and OXPHOS gene sets use invalid R pseudo-notation (`"NDUFA1"-"NDUFA10"`, `"ELOVL1"-"ELOVL7"`) — these are not valid R range expressions; actual gene vectors would need explicit enumeration.
- PPARG gene set (`pparg_gene_set`) — biology-specific to PPARG OE project.
- `tissue_type != "Myelolipoma"` exclusion filter in Key Pitfalls section 5 — HumanFat specific (Myelolipoma is an adipose tissue type in this project).
- Script filenames reference "HumanFat" project structure.
- `ec_subtype_colors` referenced in Faceted Tissue Dot Plot section — project-specific.
- Output directory `MetabolicProfile/` and `output/metabolic/` structure — project-specific.

**Generalizable content:**
- `run_aucell()` and `add_auc_to_seurat()` helper functions — fully general.
- Evidence dot plot recipe (pathway × group, size = fraction genes detected, fill = mean expression) — general.
- AUCell violin + boxplot recipe — general.
- Pathway schematic pattern — general.
- Memory optimization (stratified sampling before AUCell for large objects) — general.

**Project-specific content:** Gene sets (especially PPARG), Myelolipoma filter, specific script names, output subdirectory names. Core pathway gene sets (Glycolysis, Beta_Oxidation, TCA_Cycle, OXPHOS) are biology-general but the PPARG signature is project-specific.

**Content that duplicates other files:**
- PPARG target gene list partially overlaps with `functional_gene_sets` "Metabolic & Specialized" section in `differential_expression.md` (PPARG, CD36, FABP4, FABP5, SCD, LPL, DGAT1, DGAT2 appear in both).
- AUCell concept similar to IntegratePublicData Stage 2 AddModuleScore approach — different functions but same scoring philosophy.

**Inconsistencies:**
- OXPHOS and FA_Synthesis gene sets use invalid R range syntax (`"NDUFA1"-"NDUFA10"`) — these will fail if copied directly. Must be treated as pseudo-notation requiring explicit enumeration.
- `scale_colour_gradientn(colours = c("#4575B4","#FFFFBF","#D73027"))` in AUCell UMAP (line 151) is blue-yellow-red but the shared aesthetics file specifies `c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#D6604D","#B2182B")` for diverging data — slightly different palettes for similar use cases.

---

## pipelines/LargeDataset/methods/pyscenic_regulon_analysis.md
**Lines:** 428  
**Status:** ACTIVELY_USED (referenced by SKILL.md; validated KidneyNew + HumanFat)

**Summary:** pySCENIC TF regulon inference using 4 scripts across two conda environments (r-env for export, scenicenv for Python pipeline). Documents critical pySCENIC 0.12.1 bugs (NumPy float removal, load_motifs DataFrame type error). Covers GRN → ctx → AUCell pipeline with resource requirements, database file paths, gene filtering strategy, and analysis/visualization scripts.

**Functions defined:** None (code recipes in R and Python across 4 scripts)

**Functions called from outside this file:** None

**Project-specific hard dependencies:**
- **Machine-specific absolute paths**: `/media/david/Mayo/ryan/scRNAseq/Sinusoid/GeneLists/` (database files) and `/home/ryannachman/scenicplus/scenicplus/resources/allTFs_hg38.txt` (TF list) and `/home/ryannachman/anaconda3/envs/scenicenv/bin/python` (Python binary). All three are non-portable.
- Python environment `scenicenv` — machine-specific installation.
- `ec_subtype` metadata column — project-specific label column name.
- `UMAP_1`/`UMAP_2` coordinate export — works but may need adjustment for objects using non-standard UMAP names.
- `nohup` approach references Claude Code session ending behavior — machine-specific operational context.

**Generalizable content:**
- Bug fix documentation for pySCENIC 0.12.1 — fully general for all users of this version.
- Gene filtering thresholds for tractable runtime (≥5% cells, mean ≥0.05) — general.
- Script architecture (export → prep → GRN → AUCell → plots) — general.
- AUCell Python API usage pattern — general.

**Project-specific content:** All absolute paths, Python environment name. The analysis scripts themselves are general except for the label column name.

**Content that duplicates other files:**
- `JoinLayers()` call (line 84) repeats `seurat_v5_rules.md` Rule 1.
- "Always run with r-env" note repeats `r_environment.md`.

**Inconsistencies:**
- Machine-specific paths are hardcoded throughout and will break on any other machine. These should be documented as placeholders, not literal paths.

---

## pipelines/LargeDataset/methods/trajectory_monocle3.md
**Lines:** 183  
**Status:** ACTIVELY_USED (referenced by pipeline.md Stage 8; validated HumanFat_Yang run)

**Summary:** Monocle3 pseudotime trajectory analysis bypassing SeuratWrappers and preprocess_cds (both unavailable/broken in r-env). Injects Seurat embeddings directly into CDS object. Produces pseudotime UMAP, cell-type UMAP, density plot, and violin plot outputs.

**Functions defined:** None (inline code)

**Functions called from outside this file:** None

**Project-specific hard dependencies:**
- `ec_traj` object name and `ec_subtype` column name (lines 35, 73, 101, 125, 128) — KidneyNew/HumanFat EC specific.
- Brief configuration example uses `exclude_pattern: "^AEC"`, `start_cluster: CapEC2`, `end_cluster: VenEC3` — HumanFat EC subtype names.
- `brief$trajectory$exclude_pattern` and `brief$trajectory$start_cluster` — assumes a parsed brief object exists.

**Generalizable content:**
- The SeuratWrappers bypass pattern (Steps 2–4) — fully general for any r-env user.
- CDS construction and Seurat embedding injection — general.
- `use_partition=FALSE` rule — general.
- All output plot recipes — general (pseudotime palette, violin/boxplot/median-point pattern).

**Project-specific content:** Variable names, specific subtype references in configuration example. Easy to parameterize.

**Content that duplicates other files:**
- Pseudotime plot recipe (lines 111–131) is repeated verbatim in `LargeDataset/methods/aesthetics.md` lines 192–207. The trajectory file has the full implementation including supplementary plots; aesthetics.md has a 2-plot summary.

**Inconsistencies:** None significant.

---

## pipelines/IntegratePublicData/pipeline.md
**Lines:** 416  
**Status:** SHARED_INFRASTRUCTURE

**Summary:** Seven-stage adaptive integration pipeline for heterogeneous multi-source datasets. Defines 6 checkpoints, adaptive behavior rules, `stratified_downsample()` function, iterative QC patterns, and brief parameters reference. Validated on BoneMarrowStroma run.

**Functions defined:**
- `stratified_downsample(so, ceiling, label_col = "unified_label")` — line 220. Stratified random sampling preserving cell type proportions. Pure R implementation.

**Functions called from outside this file:**
- `JoinLayers()` — line 243, 264; resolves to `seurat_v5_rules.md` Rule 1.
- `stratified_downsample()` — self-referential in Stage 4 example.

**Project-specific hard dependencies:**
- Default metadata column values in the checkpoint 1 summary table example (line 107): `Vertebrae`, `IliacCrest`, `FemoralHead`, `Li`, `Wang` — BoneMarrowStroma dataset names.
- Iterative QC "targeted pool QC" example files `QC_cells.txt`, `Scrub_cells.txt` — generic names but the pattern is from a specific run.
- Stage 6 diagnostic UMAPs list six fixed metadata columns: `unified_label`, `dataset`, `condition`, `data_type`, `organ`, `sample_id` — these are pipeline standards, but `data_type` and `organ` may not always be present.

**Generalizable content:** The full pipeline framework is general. `stratified_downsample()` is general.

**Project-specific content:** Example data in checkpoint summary table, iterative QC file examples.

**Content that duplicates other files:**
- Stage 6 Seurat processing sequence (NormalizeData → FindVariableFeatures → ScaleData → RunPCA → RunHarmony → etc.) duplicates LargeDataset Stage 3.
- JoinLayers calls duplicate `seurat_v5_rules.md`.

**Inconsistencies:**
- "Always Read Before Starting" (lines 31–38) lists `methods/cohort_plots.md` as conditional on `cohort_plots: true` — but does NOT list `methods/anatomical_de_analysis.md`, `co_umap_embedding.md`, or `cross_atlas_dotplot.md`. These are injected via SKILL.md "read all files in methods/" but are not in the pipeline's own read list. This is the root cause of the VertebraeMyelolipoma failure.
- No Stage 8 for post-integration DE analysis — pipeline ends at CP6 with "user specifies next steps."
- No `brief_template.txt` file — IntegratePublicData has only inline Brief Parameters Reference.

---

## pipelines/IntegratePublicData/methods/load_formats.md
**Lines:** 341  
**Status:** ACTIVELY_USED (referenced by pipeline.md Stage 1)

**Summary:** Format decision tree + loading recipes for all supported formats: 10x HDF5, standard MEX directory, loose MEX files, spliced/unspliced (scVelo), in-house RDS (standard, Ensembl IDs, non-standard reductions, CDS/Monocle3, HTO-demultiplexed). Also covers GEO bulk supplementary formats (Cufflinks tmap, Excel FPKM matrices). Post-loading checklist.

**Functions defined:**
- `get_gene_fpkm(tbl, gene)` — line 264. Extracts gene-level FPKM from Cufflinks tmap table, preferring exact-match (`class_code == "="`) transcripts.

**Functions called from outside this file:** None

**Project-specific hard dependencies:**
- Leimkuler-specific file pattern (lines 46–50): `"_raw_gene_bc_matrices_h5\\.h5$"` regex and GSM4731560 prefixes.
- Wang dataset directory structure (lines 62–70): `GSM4423510_MNC_op/2.2.filtered_feature_bc_matrix/`.
- Li dataset barcodes (lines 82–89) referencing specific GSM accession numbers.
- `isolatedstroma` object reference (line 165 comment) — BoneMarrowStroma specific.
- FemoralHead UMAP reduction name `"UMAP_dim30"` (lines 193–198) — BoneMarrowStroma specific.
- TFEC Expression Atlas run (header line 2) — most recent validation.
- FPKM validation ranges (lines 275–285) are human poly-A general, but the tissue marker example (CLEC4G + LYVE1 for liver sinusoidal ECs, line 285) is TFEC Atlas specific.

**Generalizable content:**
- Format decision tree — fully general.
- HDF5 and standard MEX loading — general.
- Loose MEX files pattern — general.
- Spliced/unspliced handling with transpose check — general.
- JoinLayers after multi-layer detection — general (Rule 1 from seurat_v5_rules.md).
- Ensembl ID conversion via biomaRt — general.
- HTO demultiplexing detection — general.
- Post-loading checklist — general.

**Project-specific content:** Specific GSM accession file patterns, FemoralHead reduction name, isolatedstroma comments.

**Content that duplicates other files:**
- JoinLayers instruction (line 159) duplicates `seurat_v5_rules.md` Rule 1.
- Dynamic UMAP reduction name note (lines 193–198) duplicates `seurat_v5_rules.md` Rule 5.
- CDS to Seurat conversion (lines 202–208) uses `SingleCellExperiment` — same library documented in `shared/r_environment.md`.

**Inconsistencies:**
- The "RDS Objects — In-House" example code at lines 183–197 has `biomaRt` direct slot assignment `rownames(so@assays$RNA@counts) <- sym` (line 185–186) — this is a deprecated Seurat v5 pattern. Should use `RenameFeatures()` or equivalent.

---

## pipelines/IntegratePublicData/methods/label_harmonization.md
**Lines:** 181  
**Status:** ACTIVELY_USED (referenced by pipeline.md Stage 3)

**Summary:** Two-part label harmonization: 3A (YAML-based mapping for labeled datasets), 3B (Seurat label transfer for unlabeled datasets). Documents confidence threshold approach, controlled vocabulary rules, and special label situations from BoneMarrowStroma (sort-strategy labeling, passaged cells, heatmap_label derived columns).

**Functions defined:** None (code recipes)

**Functions called from outside this file:**
- `yaml::read_yaml()` — external package.
- `FindTransferAnchors()`, `TransferData()`, `AddMetaData()` — Seurat functions.
- `dplyr::case_when()` — explicitly namespaced (good practice given Seurat masking).

**Project-specific hard dependencies:**
- Broad label "Mesenchymal" at line 110 — BoneMarrowStroma specific default for low-confidence label transfer.
- BoneMarrowStroma result_distribution example (lines 129–133): Adipo-MSC, Osteo-MSC counts.
- Special label situations section describes BoneMarrowStroma-specific patterns: CD14/PE sort strategy labels, "Passaged Stroma" label, "Adipo-like MSC" suffix convention.

**Generalizable content:**
- YAML mapping pattern (3A) — fully general.
- Label transfer with confidence thresholds (3B) — general.
- Controlled vocabulary rules (lines 138–146) — general.
- Common pitfalls (named vector NA handling, v5 metadata assignment) — general.

**Project-specific content:** Default broad label in case_when, BoneMarrowStroma-specific special situations. The special situations section is more example-documentation than general rule.

**Content that duplicates other files:**
- `@meta.data` assignment note (lines 177–181) cross-references `seurat_v5_rules.md` Rule 3.
- `dplyr::select()` namespace note in `anatomical_de_analysis.md` Critical Constraints (line 47) — similar explicit namespace awareness, different function.

**Inconsistencies:** None significant.

---

## pipelines/IntegratePublicData/methods/cohort_plots.md
**Lines:** 235  
**Status:** ACTIVELY_USED (referenced by pipeline.md Stage 7)

**Summary:** Three cohort plot types with complete R code: (A) Pearson correlation heatmap using scaled within-dataset averages (z-score per dataset before cross-dataset average — corrects batch without limma), (B) circlize chord diagram of dataset × cell type composition, (C) pheatmap proportion heatmap of cell types per sample with metadata annotation strips.

**Functions defined:** None (code recipes)

**Functions called from outside this file:**
- `AverageExpression()` — Seurat; critical note about `as.matrix()` wrapper for Seurat v5 dgCMatrix output.
- `pheatmap::pheatmap()`, `circlize::chordDiagramFromMatrix()`, `grid::grid.draw()`.

**Project-specific hard dependencies:**
- BoneMarrowStroma validated figure size: "13.25 × 7 inches" for proportion heatmap (line 224) — project-specific.
- Chord diagram `gap.after` pattern (lines 155–156) uses list structure that implies the number of datasets/cell types.
- Annotation strip fields (`Dataset`, `Organ`, `Condition`, `Data.Type`) — match IntegratePublicData pipeline metadata columns.

**Generalizable content:** Everything except the hardcoded figure size. The `adjustcolor` matrix fix, `pivot_wider` vs for-loop bug avoidance, `pheatmap` factor level mismatch fix, `as.matrix()` after `AverageExpression()` — all fully general.

**Project-specific content:** Only the one hardcoded figure size for BoneMarrowStroma.

**Content that duplicates other files:**
- `as.matrix()` after `AverageExpression()` rule duplicated with `seurat_v5_rules.md` Rule 6 conceptually (same theme: dgCMatrix coercion needed).
- Diverging blue-white-red correlation palette matches `cross_atlas_dotplot.md` Heatmap color scheme and `differential_expression.md` z-score heatmap color.

**Inconsistencies:**
- Chord diagram sector labels use `circos.track()` and the file comments "NOT `circos.trackPlotRegion`" — the latter function may or may not exist in the circlize version in r-env; this is undocumented.

---

## pipelines/IntegratePublicData/methods/anatomical_de_analysis.md
**Lines:** 347  
**Status:** ORPHAN (in methods/ but NOT referenced in pipeline.md "Always Read Before Starting" or any pipeline stage)

**Summary:** Anatomical DE analysis framework for BoneMarrowStroma 3-site comparison (Vertebrae, Iliac Crest, Femoral Head). Defines `make_topgene_dotplot()` locally (3-panel version), describes `make_functional_dotplot()` with 3-site direction strip, and `1.7_SiteDE_GOFunctionalDotplot.R` (data-driven GO sections). Calls `run_findmarkers()`, `make_volcano()`, `make_overall_heatmap()`, `make_functional_heatmap()`, `make_pathway_barplot()` without defining them.

**Functions defined:**
- `make_topgene_dotplot(obj, markers, comp, subset_name, n_each = 12)` — line 77. 3-panel (Vertebrae | Iliac Crest | Femoral Head) dot plot. Different signature and behavior from `differential_expression.md`'s version.
- `make_functional_dotplot(...)` — line 142 (described, partial implementation). 3-panel version with "Femoral Head" pinning vs 2-panel with `comp$ident2` pinning in `differential_expression.md`.
- `is_ambient(genes)` — line 344. Uses only `AMBIENT_PATTERNS`; lacks `AMBIENT_EXPLICIT` list. Different implementation from `differential_expression.md`.

**Functions called from outside this file (as orphan, these would be called within script context):**
- `run_findmarkers(so_scope, comp)` — line 323 in main loop. **UNRESOLVED: not defined in this file or any IntegratePublicData methods file.**
- `make_volcano(markers, comp, scope_name)` — line 325. **UNRESOLVED.**
- `make_overall_heatmap(so_scope, markers, comp, scope_name)` — line 326. **UNRESOLVED.**
- `make_functional_heatmap(so_scope, markers, comp, scope_name)` — line 327. **UNRESOLVED.**
- `make_pathway_barplot(markers, comp, scope_name, universe_genes)` — line 331. **UNRESOLVED.**

**Project-specific hard dependencies:**
- Factor levels `c("Vertebrae","Iliac Crest","Femoral Head")` hardcoded in: `make_topgene_dotplot()` (line 104), Critical Constraints table (line 44), `make_functional_dotplot()` section label pinning (line 168), main loop comparisons (lines 312–314).
- `STROMA_ORDER` referenced but not defined in this file (line 68, 107).
- Section label pinned to `"Femoral Head"` (line 168) — always the rightmost of the 3 sites.
- `unified_label` as the cell type column (line 92, 311) — IntegratePublicData standard.
- Site column called `"site"` in comparisons (lines 312–314) — BoneMarrowStroma specific.

**Generalizable content:**
- The conceptual approach (cell-type-resolved dot plots to disentangle composition from expression) is excellent and general.
- Data-driven GO section algorithm (1.7 section, lines 217–258) is fully generalizable.
- Critical Constraints table (lines 40–51) has general fixes.

**Project-specific content:** All site name hardcoding, STROMA_ORDER, 3-panel assumption, "Femoral Head" rightmost panel pinning. The fundamental design limitation: this file was written for exactly 3 specific sites.

**Content that duplicates other files:**
- `is_ambient()` duplicate (diverged version) — see inconsistencies.
- `make_topgene_dotplot()` duplicate definition vs `differential_expression.md`.
- `make_functional_dotplot()` partial duplicate vs `differential_expression.md`.
- Same output file naming convention as `differential_expression.md` (lines 280–299).

**Inconsistencies (critical):**
1. **Phantom functions**: calls `run_findmarkers()`, `make_volcano()`, `make_overall_heatmap()`, `make_functional_heatmap()`, `make_pathway_barplot()` without defining them. These are defined in `LargeDataset/methods/differential_expression.md` but that file is never injected for IntegratePublicData jobs.
2. **`is_ambient()` diverged**: this file's version (line 344) uses only `AMBIENT_PATTERNS` and implements `Reduce(`|`, lapply(AMBIENT_PATTERNS, function(p) grepl(p, genes)))` — returns a plain logical. `differential_expression.md`'s version also checks `genes %in% AMBIENT_EXPLICIT` — a richer filter. The simpler version will miss explicit ambient genes (ALB, MALAT1, platelet markers, etc.).
3. **`make_topgene_dotplot()` name collision**: defines a function with the same name as the one conceptually described in `differential_expression.md` but with different signature (`n_each=12` vs no explicit argument), different behavior (3-panel facet vs single-panel), and different output aesthetics.

---

## pipelines/IntegratePublicData/methods/co_umap_embedding.md
**Lines:** 166  
**Status:** ORPHAN (in methods/ but NOT referenced in pipeline.md stages)

**Summary:** Pre-UMAP downsampling strategy for co-embedding in-house subtypes with a public atlas (Tabula Sapiens). Documents 4-panel figure design (dataset source, in-house subtypes highlighted, atlas organs highlighted, colorblind-safe comparison panel). Validates WCM vs TS downsampling ratio (1500/subtype).

**Functions defined:** None

**Functions called from outside this file:** None

**Project-specific hard dependencies:**
- `"WCM"` and `"TS"` as dataset names throughout — TabulaSapiensComparison specific.
- `"CapEC"` as the highlighted subtype in Panel D — HumanFat EC subtype name.
- `"Fat EC"` as the highlighted atlas organ — TS organ label.
- Script names `06f_ts_umap_presample.R` and `06g_ts_umap_replot.R` — project-specific.
- Output directory `output2/` — project convention.
- Coordinates CSV named `ts_harmony_umap_presample_coords.csv` — project-specific.
- Color constants `CB_CAPEC = "#E69F00"`, `CB_FAT_EC = "#0072B2"` (Okabe-Ito) — general but applied to project-specific groups.
- `N_PER_SUBTYPE = 1500` recommendation tied to the WCM/TS ratio analysis.
- EC subtype colors `ec_colors` — HumanFat palette.

**Generalizable content:**
- The "downsample before UMAP" principle (lines 19–47) is general and valuable for any atlas comparison. The ratio argument (80-85% atlas cells to dominate embedding geometry) is reusable reasoning.
- Save coordinates as CSV for fast replotting — excellent general pattern.
- 4-panel figure design concept (source / in-house / atlas / comparison) — general template.
- `ggsave(..., device = cairo_pdf)` vs `pdf()` + `print()` fix for patchwork — general.
- Okabe-Ito colorblind-safe palette choice — general.

**Project-specific content:** WCM/TS names, CapEC/Fat EC highlighted groups, all hardcoded labels, output filenames.

**Content that duplicates other files:** None significant.

**Inconsistencies:**
- Uses `theme_minimal()` (implied in Panel A ggplot call) but `shared/aesthetics.md` says "Never use default theme_gray()" and implies `theme_classic()` as default. `theme_minimal()` is permissible but the inconsistency is notable.
- Font `base_family = "Source Sans 3"` noted in Common Pitfalls as causing an error — this font is not available in r-env but appears in `cross_atlas_dotplot.md` as well.

---

## pipelines/IntegratePublicData/methods/cross_atlas_dotplot.md
**Lines:** 232  
**Status:** ORPHAN (in methods/ but NOT referenced in pipeline.md stages)

**Summary:** Cross-atlas comparison methods for WCM vs Tabula Sapiens: (A) Harmony centroid correlation heatmap, (B) depth-corrected 3-section dotplot using within-dataset z-scoring, (C) TF diamond variant. Documents sequencing depth batch effect solution (z-score within each dataset before combining).

**Functions defined:** None (code recipes)

**Functions called from outside this file:** None

**Project-specific hard dependencies:**
- `"WCM Fat CapEC"`, `"WCM Fat VenEC1"`, etc. column label convention (lines 31–32) — HumanFat specific.
- `"Fat EC"`, `"Kidney EC"` TS organ label convention (line 32) — TabulaSapiensComparison.
- `EC_SUBTYPES` vector with WCM subtypes (line 120) — HumanFat.
- `"CapEC"` as the marker of interest in WCM (lines 143–145) — HumanFat EC.
- Script names reference TabulaSapiensComparison project.
- `allMarkersEC.csv` and `fat_enriched_markers.Rds` file names — project-specific.
- `capec_markers_current.Rds`, `fat_enriched_markers.Rds` — project-specific cached outputs.
- Custom fonts: `"Source Sans 3"`, `"Playfair Display"` — these are NOT available in standard r-env (confirmed bug in co_umap_embedding.md).
- `SEC1_FORCE` vector for hardcoded Section 1 TF genes — project-specific.

**Generalizable content:**
- Depth-correction approach (z-score within dataset, cap ±2) — fully general for any cross-dataset comparison.
- Background highlight rectangles for column groups (warm yellow for WCM, light red for TS Fat) — general visual pattern.
- `useDingbats` incompatibility with `cairo_pdf` (line 158) — general fix.
- Centroid correlation in Harmony space (Part A) — general approach for any batch-corrected atlas comparison.
- CSV spaces → dots fix for `cor_mat` column names (lines 36–39) — general.

**Project-specific content:** All WCM/TS naming, marker file references, custom fonts, specific column ordering logic.

**Content that duplicates other files:**
- Diverging blue-white-red palette (`#2166AC` → `#B2182B`) duplicated from `shared/aesthetics.md` and `cohort_plots.md`.
- `geom_hline on discrete y-axis` fix (line 152) noted as `geom_hline` works, `annotate("segment")` fails — similar axis issue documented in `differential_expression.md` for `annotate("text")`.
- `ComplexHeatmap::Heatmap` with `row_order`/`column_order` pattern — general usage present in `differential_expression.md` as well.

**Inconsistencies:**
- Uses `theme_minimal()` with `text = element_text(family = "Source Sans 3")` — font unavailable in r-env (documented as bug in co_umap_embedding.md pitfalls).
- y-axis `element_text(family="Playfair Display")` — also unavailable.

---

## lab_context.md
**Lines:** 210  
**Status:** SHARED_INFRASTRUCTURE

**Summary:** Lab-level context file covering default organism, tissues/cell compartments studied, standard QC thresholds (empty — left for data-driven suggestion), label harmonization vocabulary, pipeline convention defaults, color palettes for four contexts (whole object, adipose types, organ/tissue, EC subtypes for adipose + kidney), and per-project notes for fat-scrnaseq-continued, EyePublicData, and KidneyNew.

**Functions defined:** None

**Functions called from outside this file:** None

**Project-specific hard dependencies:**
- Color palettes: `celltype_colors`, `adipose_type_colors`, `organ_tissue_colors`, `ec_subtype_colors_adipose`, `ec_subtype_colors_kidney` — all project-specific from HumanFat and KidneyNew.
- Per-project Notes section: 3 specific project notes (fat-scrnaseq-continued, EyePublicData, KidneyNew) with hardcoded column names, filter rules, and assay conventions.
- `batch_correction_var: sample_id` — intended lab default but conflicts with `LargeDataset/brief_template.txt` default of `source_file`.
- Comment `#- Bone marrow stroma` is commented out from primary_tissues — suggests BoneMarrowStroma was active but is no longer.
- PPARG biology context in fat project notes.

**Generalizable content:** The structural framework (default_organism, tissues, QC thresholds template, label vocabulary, pipeline convention defaults) is general and should be preserved in v2.

**Project-specific content:** All color palettes, all per-project notes, commented-out bone marrow stroma.

**Content that duplicates other files:**
- EC subtype color palettes duplicated in `lab_context/validated_examples.yaml` (fat-scrnaseq-continued entry), `LargeDataset/methods/aesthetics.md` (proportion plot), `cross_atlas_dotplot.md` Section D aesthetics.
- CellChat v2 bug list in KidneyNew notes partially duplicates `interactome_cellchat.md` Critical Constraints.
- ASSAY SWITCHING note in KidneyNew notes partially duplicates seurat_v5_rules.md context.

**Inconsistencies:**
- `batch_correction_var: sample_id` here conflicts with `LargeDataset/brief_template.txt`'s `batch_correction_var: source_file`. The brief_template.txt overrides lab_context.md in practice.
- `unified_label_vocabulary` section (lines 61–69) is all commented out — no standardized cell type vocabulary is actually defined despite this being the intended location.
- `known_aliases` section (lines 71–75) is all commented out — no universal aliases populated.

---

## lab_context/validated_examples.yaml
**Lines:** 162  
**Status:** EXAMPLE_DOCUMENTATION

**Summary:** Four job records (eye-test-01 as template, fat-scrnaseq-continued, KidneyNew, EyePublicData) documenting successful runs with description summaries, pipeline choices, non-default parameters, key gotchas, and notes. The eye-test-01 entry is unfilled (TODO markers).

**Functions defined:** None

**Functions called from outside this file:** None

**Project-specific hard dependencies:** All entries are by definition project-specific. eye-test-01 is an empty template.

**Generalizable content:** The schema (job_id, description_summary, pipeline, organism, tissue, non_default_params, key_gotchas, notes) is the generalizable part. The filled examples are reference documentation.

**Project-specific content:** All filled entries.

**Content that duplicates other files:** 
- fat-scrnaseq-continued key_gotchas duplicates `lab_context.md` Notes section for that project.
- KidneyNew key_gotchas duplicates `lab_context.md` Notes for KidneyNew.
- EyePublicData notes partially duplicate `lab_context.md` Notes for EyePublicData.

**Inconsistencies:**
- `eye-test-01` has `# TODO` placeholders — this is a stub, not a real validated example. Should either be filled or removed in v2.
- fat-scrnaseq-continued `pipeline_variant` field does not appear in the schema defined at the top of the file — schema drift.

---

## SKILL.md
**Lines:** 101  
**Status:** SHARED_INFRASTRUCTURE

**Summary:** Master router that defines the two available pipelines, how to generate project files (read pipeline + all methods + all shared → generate CLAUDE.md + analysis_brief.txt), how to add new methods, and the file location map. The critical "read ALL files in methods/" instruction is the root cause of irrelevant file injection.

**Functions defined:** None

**Functions called from outside this file:** None

**Project-specific hard dependencies:**
- File location map has specific status comments per file (✅/🔲 markers) referencing specific project names.
- `de_comprehensive_csv.md` indented incorrectly in file tree (4 spaces instead of 8, missing `│   │   └─`) — formatting error.

**Generalizable content:** Everything structural.

**Project-specific content:** Status comments, project-name references in comments.

**Inconsistencies:**
- "Step 3: Read all files in the pipeline's methods/ directory" — this blunt instruction injects all methods files including large, irrelevant ones (1302-line interactome_cellchat.md for non-CellChat jobs; 166-line co_umap_embedding.md for 2-dataset integrations that aren't atlas comparisons).
- `interactome` Stage 8 entry in SKILL.md file map shows `interactome_cellchat.md` as "✅ validated KidneyNew run" but `LargeDataset/pipeline.md` Stage 8 still shows "(🔲 placeholder)" for the interactome entry.
- `anatomical_de_analysis.md` appears in SKILL.md file map under IntegratePublicData/methods but is NOT in pipeline.md "Always Read Before Starting" — creating a silent injection via SKILL.md's "read all files" but with no pipeline stage context.

---

## PROJECT_CLAUDE_TEMPLATE.md
**Lines:** 11  
**Status:** SHARED_INFRASTRUCTURE

**Summary:** Minimal template for new analysis project CLAUDE.md files. Contains only `@~/claude-skills/SKILL.md` and a comment block for project-specific overrides.

**Functions defined:** None  
**Project-specific content:** None.  
**Inconsistencies:** None.
