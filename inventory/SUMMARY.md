# Inventory Summary — ~/claude-skills/
*Generated: 2026-05-11 by inventory agent*

TL;DR for the human reading the full inventory reports.

---

## Overall Library Health Assessment

The library is structurally sound and analytically rich. It contains genuinely validated, working content from at least five real projects (HumanFat/fat-scrnaseq-continued, KidneyNew, BoneMarrowStroma, TabulaSapiensComparison, CoCulture). The larger files (`differential_expression.md`, `interactome_cellchat.md`, `bulk_concordance.md`) contain deep domain knowledge and hard-won bug fixes that would take significant effort to reconstruct.

The primary pathology is **structural drift from bottom-up accretion** without cross-file governance: project-specific column names have propagated into shared templates, functions were duplicated as they crossed project boundaries without unification, and three orphan methods files (anatomical_de_analysis.md, co_umap_embedding.md, cross_atlas_dotplot.md) were written for IntegratePublicData but never wired into its pipeline. The library has excellent raw material for a v2 restructure. The restructure task is primarily editorial (parameterize, deduplicate, wire up) rather than content-generative.

**Estimated content that is broadly reusable as-is:** ~65% (shared/ files, most LargeDataset methods, cohort_plots.md, load_formats.md, label_harmonization.md, the core functions in differential_expression.md)

**Estimated content requiring parameterization before v2:** ~25% (column name defaults in brief_template.txt, HumanFat-specific palette references throughout, metabolic gene set pseudo-notation)

**Estimated content to delete without replacement:** ~10% (eye-test-01 stub, stale comments, the `projectData()` invalid call)

---

## Top 5 Most Important Findings

### 1. Three IntegratePublicData Methods Files Are Orphans with Phantom Dependencies

`anatomical_de_analysis.md`, `co_umap_embedding.md`, and `cross_atlas_dotplot.md` appear in the IntegratePublicData methods/ directory and are read by SKILL.md's blanket "read all files" instruction — but none are in the pipeline.md "Always Read Before Starting" list and none have a corresponding pipeline stage. More critically, `anatomical_de_analysis.md` calls five functions (`run_findmarkers`, `make_volcano`, `make_overall_heatmap`, `make_functional_heatmap`, `make_pathway_barplot`) that are defined in `LargeDataset/methods/differential_expression.md` and never injected into IntegratePublicData jobs. Any script generated from `anatomical_de_analysis.md` in a real IntegratePublicData context will fail with "could not find function" errors. This is the most operationally critical finding.

### 2. `is_confound()` Is Called But Never Defined — `de_comprehensive_csv.md` Is Currently Broken

`de_comprehensive_csv.md` calls `is_confound(gene)` at the core of its gene filtering step (line 116). This function does not exist anywhere in the library. The file's own documentation implies it should be a superset of `is_ambient()` adding sex-linked, HLA class II, histone, and unannotated gene exclusions. Until `is_confound()` is written, every script generated from this methods file will error at the filtering step.

### 3. Two Diverged `is_ambient()` Implementations — Silent Results Difference

The library has two definitions of `is_ambient()` with the same name but different filter completeness. The `differential_expression.md` version (the correct one) checks both regex patterns AND an explicit list of 16 ambient genes including MALAT1, NEAT1, ALB, and platelet markers. The `anatomical_de_analysis.md` version checks patterns only. The comment in `anatomical_de_analysis.md` claims it is "same as LargeDataset pipeline" — this is false. DE analyses using the IntegratePublicData version will include MALAT1 and platelet marker genes in their outputs. This is a silent correctness issue with no error message.

### 4. `brief_template.txt` Has HumanFat-Specific Defaults That Override `lab_context.md`

The `LargeDataset/brief_template.txt` default `batch_correction_var: source_file` contradicts `lab_context.md`'s `batch_correction_var: sample_id`. `source_file` is a HumanFat metadata column name. A new project using the template verbatim will attempt to harmonize on a column that doesn't exist. This also affects `pipeline.md`'s Brief Parameters Reference table. The HumanFat project's specific defaults have leaked into what should be a general starting point.

### 5. `interactome_cellchat.md` Contains a Self-Contradicting Code Block

The inference pipeline code skeleton at line 71 calls `projectData(cellchat, PPI.human)`. The same file's Critical Constraints section explicitly says this function was removed in CellChat v2 and must not be called. Any generated CellChat script will contain this invalid call and fail immediately on a v2 CellChat installation. This is a live bug in the most complex file in the library.

---

## Top 5 Surprises

### 1. `skills_gap_analysis.md` Does Not Exist

The inventory task description instructs reading `~/claude-skills/jobs/skills_gap_analysis.md` as orientation. The file — and the entire `jobs/` directory — does not exist in the library. It was presumably written in a project working directory and never committed to the library. The inventory was conducted without it.

### 2. The Library Has No Doublet Detection

For a production-quality scRNAseq library, the complete absence of doublet detection documentation is surprising. Neither pipeline documents DoubletFinder, scDblFinder, or any other doublet removal approach. QC in both pipelines stops at percent_mt and feature count thresholds.

### 3. `get_filtered_net()` in CellChat Is Fixing a v2 Bug with No Bug Report

`interactome_cellchat.md` defines `get_filtered_net()` to work around a `pairLR.use` bug in CellChat v2 chord diagrams. The bug is documented only within this one methods file with no reference to whether it has been reported to the CellChat maintainers or has a known version fix timeline. Users who upgrade CellChat expecting the bug to be fixed may not know to check this file.

### 4. Three Files Reference `ec_colors` / `tissue_colors` Without Defining Them

`bulk_concordance.md`, `cross_atlas_dotplot.md`, and `metabolic_profile.md` all reference named color vectors (`ec_colors`, `tissue_colors`, `ec_subtype_colors`) that are not defined in those files. The implicit assumption is that these are injected via the project CLAUDE.md from `lab_context.md`. But if the CLAUDE.md generation omits the color block, all three methods files will fail with an R `object not found` error. The library has no mechanism to ensure the color context block is always included.

### 5. The r-env `conda run` Command in the Canonical Skills File Is Wrong

`shared/r_environment.md` shows `conda run -n r-env Rscript` without `--no-capture-output`. The user's system CLAUDE.md specifies `conda run --no-capture-output -n r-env Rscript` as the validated form (the flag prevents the command from hanging after the script completes). The canonical library file has the wrong command, and all methods files that cite it as authority carry the same error.

---

## Recommended Order of Operations for the Restructure

### Phase 1 — Critical Fixes (before v2 is usable for any new project)
These are bugs, not design issues. Fix them first so v2 doesn't inherit broken content.

1. **Fix `projectData()` call** in `interactome_cellchat.md` line 71 (D1).
2. **Fix `conda run` command** in `shared/r_environment.md` to add `--no-capture-output` (D11).
3. **Write `is_confound()`** as a proper named function in `shared/` or a new `primitives/` file (Coverage Gap §2.1).
4. **Fix OXPHOS/FA_Synthesis gene sets** in `metabolic_profile.md` to use explicit vectors (D6).
5. **Fix `brief_template.txt`** — clear `batch_correction_var: source_file` and HumanFat metadata column defaults (D8, D9).

### Phase 2 — Structural Fixes (resolve the orphan and cross-pipeline dependency problems)
6. **Add `differential_expression.md` functions to IntegratePublicData's read scope** — either by injecting the file or by moving the core DE functions to a shared location accessible by both pipelines. This resolves all UNRESOLVED function calls in `anatomical_de_analysis.md`.
7. **Create `IntegratePublicData/brief_template.txt`** (Coverage Gap §1.1).
8. **Wire the three orphan methods files** into `IntegratePublicData/pipeline.md` as conditional optional analyses with clear trigger conditions in the brief.
9. **Unify `is_ambient()` into a single canonical version** in `shared/` using the full dual-filter (patterns + explicit list).

### Phase 3 — Parameterization (make shared content truly reusable)
10. **Parameterize `functional_gene_sets`** out of `differential_expression.md` into a brief-supplied configuration block (G5).
11. **Create a shared `primitives/harmony_integration.md`** recipe (G8).
12. **Parameterize `celltype_subclustering.md`** to remove `source_file`/`tissue_type` hardcoding (G7).
13. **Create `primitives/aucell_scoring.md`** from `run_aucell()` and `add_auc_to_seurat()` (G6).

### Phase 4 — Canonical Color Palette (design once, reference everywhere)
14. **Define a single `context/color_palettes.md`** with the canonical 6-stop diverging palette, gene expression continuous scale, and a clear slot for project-specific palettes. Replace scattered inline palette definitions with references (Duplication Report §3.1, §3.2).

### Phase 5 — Convention and Documentation (clean up what's left)
15. **Rename `make_topgene_dotplot` in `anatomical_de_analysis.md`** to `make_anatomical_dotplot` (function_index.md Table 3.2).
16. **Create `IntegratePublicData/methods/aesthetics.md`** (Coverage Gap §1.2).
17. **Create `examples/` directory** with project-specific instantiation files (BoneMarrow, TabulaSapiens, HumanFat, etc.).
18. **Remove deletion candidates D2–D5 and D10** (cleanup).

---

## Open Questions for Human Review

1. **`STROMA_ORDER` — where is the actual definition?** The variable is referenced in `anatomical_de_analysis.md` but never defined in any library file. It presumably existed in the BoneMarrowStroma project CLAUDE.md. Is there a canonical BoneMarrowStroma cell type ordering that should be brought into the library as part of the BoneMarrowStroma example file?

2. **Is `is_confound()` intentionally broader than `is_ambient()`?** The coverage gaps analysis proposes a full list of additional exclusion categories. But the original intent may have been simpler (just a renamed `is_ambient()`). Human confirmation needed before writing the function.

3. **Should `differential_expression.md` functions be accessible to IntegratePublicData?** The structural decision here is: (a) move the DE core functions to a shared `primitives/` directory, (b) cross-inject the LargeDataset file into IntegratePublicData jobs that need DE, or (c) maintain fully separate DE implementations per pipeline. Each has different maintenance implications.

4. **Are the three atlas comparison methods (`co_umap_embedding.md`, `cross_atlas_dotplot.md`) intended only for TabulaSapiens-scale comparisons, or are they general atlas comparison tools?** The specificity of their naming (TS references everywhere) suggests they were written for one project and the generalization hasn't been validated.

5. **What is the timeline for a planned EyePublicData pipeline?** The lab_context.md primary_tissues includes Eye (RPE/Choroid), and `validated_examples.yaml` has an EyePublicData entry with KidneyNew-style pipeline notes. But there is no Eye-specific pipeline or methods in the library. Is this intentional (using IntegratePublicData for eye work) or is an Eye pipeline planned?

6. **Should `validated_examples.yaml` be the project registry?** Currently project notes live in two places (lab_context.md + validated_examples.yaml). If validated_examples.yaml becomes the single source of per-project truth, it needs a richer schema (dates, status, output directory, key column names used). Should this be structured enough to auto-populate CLAUDE.md defaults for returning-to-old-project scenarios?
