# Duplication Report — ~/claude-skills/
*Generated: 2026-05-11 by inventory agent*

Content duplication beyond function definitions. Three categories: verbatim code blocks, conceptual content described differently in 2+ files, and configuration blocks with drift.

---

## 1. Verbatim or Near-Verbatim Code Block Duplications

### 1.1 Trajectory/Pseudotime Plot Recipe

**Locations:**
- `pipelines/LargeDataset/methods/trajectory_monocle3.md` lines 110–131
- `pipelines/LargeDataset/methods/aesthetics.md` lines 196–207

**What's duplicated:** The `plot_cells(cds, color_cells_by = "pseudotime", cell_size = 0.6) + theme_classic(base_size = 16) + scale_color_viridis_c(option = "magma")` pattern and the supplementary violin sorted-by-median recipe (`stat_summary(fun = median, geom = "point", size = 2.5, color = "#E53935")` + `geom_boxplot(width = 0.07, ...)`).

**Identical or drifted?** Near-identical. `trajectory_monocle3.md` has the full implementation including `ggsave()` calls and supplementary density plots; `aesthetics.md` has a 2-plot summary without `ggsave`. The core plot code is functionally identical. The `# magma — NOT inferno` comment appears in both files.

**Recommendation:** `aesthetics.md` should reference `trajectory_monocle3.md` rather than duplicating the code. The summary in `aesthetics.md` provides no additional information.

---

### 1.2 Seurat v5 JoinLayers Call

**Locations:**
- `shared/seurat_v5_rules.md` — Rule 1 (canonical definition)
- `pipelines/IntegratePublicData/methods/load_formats.md` line 159 (inline note)
- `pipelines/IntegratePublicData/pipeline.md` Stage 5 line 264 (pipeline step)
- `pipelines/IntegratePublicData/methods/cohort_plots.md` line 111 (note in checklist)
- `pipelines/LargeDataset/methods/pyscenic_regulon_analysis.md` line 84 (Script 1 prep)
- `pipelines/LargeDataset/methods/celltype_subclustering.md` line 222 (Phase 1 prep)

**What's duplicated:** The instruction to call `JoinLayers()` after merging multiple Seurat objects, because Seurat v5 keeps layers separate until joined.

**Identical or drifted?** Mostly identical instruction, varying levels of explanation. The canonical version in `seurat_v5_rules.md` explains why (layers are kept separate in v5 memory architecture); the inline versions are reminders without explanation. Not drifted — consistent message.

**Recommendation:** In v2, the canonical rule in `seurat_v5_rules.md` is sufficient. Inline reminders can be replaced with a comment `# see seurat_v5_rules.md Rule 1`.

---

### 1.3 Dynamic UMAP Reduction Name Detection

**Locations:**
- `shared/seurat_v5_rules.md` — Rule 5 (canonical)
- `pipelines/LargeDataset/methods/feature_umap_plot.md` lines 52–54 (verbatim code block)

**What's duplicated:** The pattern for detecting the current UMAP reduction name from `so@reductions`:
```r
umap_key <- names(so@reductions)[grepl("umap", names(so@reductions), ignore.case = TRUE)][1]
```

**Identical or drifted?** Feature_umap_plot.md elaborates with more context about why this is needed (objects with non-standard reduction names from import); the rule in seurat_v5_rules.md is more terse. Not drifted.

---

### 1.4 `r-env` Execution Reminder

**Locations:**
- `shared/r_environment.md` — canonical definition
- `pipelines/LargeDataset/methods/differential_expression.md` Critical Constraints line 36
- `pipelines/LargeDataset/methods/bulk_concordance.md` header
- `pipelines/LargeDataset/methods/celltype_subclustering.md` header
- `pipelines/LargeDataset/methods/metabolic_profile.md` header
- `pipelines/LargeDataset/methods/pyscenic_regulon_analysis.md` header
- `pipelines/IntegratePublicData/methods/anatomical_de_analysis.md` Critical Constraints

**What's duplicated:** "Always run with `conda run -n r-env Rscript`" / "Use r-env conda environment".

**Identical or drifted?** One notable drift: `shared/r_environment.md` shows `conda run -n r-env Rscript` (without `--no-capture-output`), but the user's CLAUDE.md mandates `conda run --no-capture-output -n r-env Rscript`. The inline reminders in methods files do not include `--no-capture-output` either. The canonical rule is therefore stale.

---

### 1.5 `ggsave(useDingbats=FALSE)` Rule

**Locations:**
- `shared/aesthetics.md` File Output section (canonical)
- `shared/r_environment.md` Script Structure conventions
- Multiple methods files (differential_expression.md, trajectory_monocle3.md, co_umap_embedding.md, cross_atlas_dotplot.md)

**What's duplicated:** `ggsave(..., useDingbats = FALSE)` rule for PDF output.

**Identical or drifted?** Consistent. One exception: `celltype_subclustering.md` lines 108–110 uses `pdf() + dev.off()` (base R graphics) instead of `ggsave()`, bypassing this rule entirely — an inconsistency, not a duplication.

---

## 2. Same Concept Described Differently in Multiple Files

### 2.1 Ambient Gene Filter: `AMBIENT_PATTERNS` vs `AMBIENT_EXPLICIT`

**Locations:**
- `pipelines/LargeDataset/methods/differential_expression.md` lines 82–104: Two-part definition — `AMBIENT_PATTERNS` (regex) AND `AMBIENT_EXPLICIT` (named list of 16 specific genes)
- `pipelines/IntegratePublicData/methods/anatomical_de_analysis.md` lines 341–346: Defines only `AMBIENT_PATTERNS` (same regex list as above)
- `pipelines/LargeDataset/methods/de_comprehensive_csv.md` line 116: References `is_confound()` — a phantom that was described as covering additional filters (sex-linked genes, HLA class II, histone genes, ENSG IDs) beyond what either `is_ambient()` version covers

**What diverges:** Three distinct ambient filter specifications exist in the library with increasing strictness: (1) patterns-only, (2) patterns+explicit, (3) patterns+explicit+sex-linked+HLA+histones+ENSG (phantom). None of the three explicitly states it supersedes the others. The comment in `anatomical_de_analysis.md` claims "(same as LargeDataset pipeline)" which is false.

**Impact:** Any script generated from `anatomical_de_analysis.md` will silently pass MALAT1 and platelet markers through the ambient filter.

---

### 2.2 Bulk Concordance vs Bulk LFC Heatmap — Self-Documented Distinction

**Locations:**
- `pipelines/LargeDataset/methods/bulk_concordance.md` Overview section lines 10–24
- `pipelines/LargeDataset/methods/bulk_lfc_concordance_heatmap.md` lines 10–34

**What's duplicated:** Both files describe the distinction between themselves in a comparison table. The distinction is accurate and consistent, making this a positive pattern (self-documenting files). The tables are differently structured but convey the same information.

**Verdict:** This is *intentional* duplication for discoverability. Keep both cross-references in v2. Not a cleanup candidate.

---

### 2.3 CellChat v2 Bug List

**Locations:**
- `pipelines/LargeDataset/methods/interactome_cellchat.md` Critical Constraints section (6 bugs)
- `lab_context.md` KidneyNew project notes (partially overlapping list)

**What diverges:** `interactome_cellchat.md` has the authoritative 6-bug list. `lab_context.md` duplicates ~3 of them in prose notes under KidneyNew. The `lab_context.md` version is less complete and less structured. The `interactome_cellchat.md` version should be treated as canonical.

---

### 2.4 Seurat v5 `@meta.data` Assignment Rule

**Locations:**
- `shared/seurat_v5_rules.md` Rule 3 (canonical)
- `pipelines/IntegratePublicData/methods/label_harmonization.md` lines 177–181 (referenced)
- `pipelines/LargeDataset/methods/celltype_subclustering.md` line 288 (referenced)

**What's duplicated:** The rule that `so$col <-` is unsafe in Seurat v5; must use `so@meta.data$col <-`. These are cross-references, not full duplications.

---

### 2.5 Harmony Integration Processing Sequence

**Locations:**
- `pipelines/LargeDataset/pipeline.md` Stage 3: `NormalizeData → FindVariableFeatures → ScaleData → RunPCA → RunHarmony → FindNeighbors → FindClusters → RunUMAP`
- `pipelines/IntegratePublicData/pipeline.md` Stage 6: Same sequence

**What's duplicated:** The standard Seurat+Harmony processing pipeline. Both describe the identical sequence with nearly identical code.

**Drifted?** Slightly — `LargeDataset` uses `n_variable_features: 4000` as default, `IntegratePublicData` does not specify a default. The `RunHarmony` `group.by.vars` source differs (`batch_correction_var` brief field in LargeDataset, inline `"dataset"` in IntegratePublicData). Otherwise equivalent.

**Recommendation:** This is a strong candidate for a shared `primitives/harmony_integration.md` in v2.

---

## 3. Configuration Blocks Appearing in Multiple Files with Drift

### 3.1 Diverging Blue-White-Red Color Palette

**Locations and values:**
| File | Palette values |
|---|---|
| `shared/aesthetics.md` (canonical) | `c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#D6604D","#B2182B")` (6-stop) |
| `differential_expression.md` volcano | `c(up_ident1 = "#B2182B", up_ident2 = "#2166AC", ns = "grey78")` (endpoints only) |
| `cohort_plots.md` correlation heatmap | `c("#2166AC","white","#B2182B")` (3-stop) |
| `cross_atlas_dotplot.md` correlation heatmap | `c("#2166AC","#B2182B")` (2-stop endpoints in a `colorRampPalette`) |
| `metabolic_profile.md` AUCell UMAP | `c("#4575B4","#FFFFBF","#D73027")` (different blue-yellow-red 3-stop) |
| `LargeDataset/methods/aesthetics.md` diff abundance | `scale_fill_gradient2(low="steelblue",high="firebrick")` (named colors) |

**Assessment:** The endpoints (`#2166AC` blue, `#B2182B` red) are consistent across most uses. The midpoint and intermediate stops vary. `metabolic_profile.md` uses a different color family entirely (`#4575B4`/`#D73027` from RColorBrewer RdYlBu vs the RdBu family used elsewhere). `steelblue`/`firebrick` in aesthetics.md is a low-specificity variant. No file references another as the source. This is genuine drift.

**Recommendation:** Define a single canonical 6-stop diverging palette in `primitives/color_palettes.md` in v2. All uses should derive from it.

---

### 3.2 Gene Expression Continuous Color Scale

**Locations:**
| File | Scale |
|---|---|
| `shared/aesthetics.md` UMAP Feature Plots section | `c("lightgrey", "blue")` (2-stop) |
| `feature_umap_plot.md` | `colorRampPalette(c("#F5F5F5","#E53935"))(100)` (100-stop grey→red) |
| `differential_expression.md` dot plots | `scale_color_gradient(low = "grey92", high = "#B2182B")` |
| `anatomical_de_analysis.md` make_topgene_dotplot | `scale_color_gradient(low = "grey92", high = "#B2182B")` |

**Assessment:** `shared/aesthetics.md` specifies `"lightgrey"→"blue"` for UMAP expression overlays, but the validated implementation in `feature_umap_plot.md` uses grey→red — a fundamentally different color scheme. The dot plot usage (grey→red #B2182B) is consistent between files but inconsistent with the shared principles for UMAPs. The shared file's advice appears outdated or intended for a different plot type than what gets generated.

---

### 3.3 `batch_correction_var` Default

**Locations:**
| File | Default value |
|---|---|
| `lab_context.md` Pipeline Conventions | `sample_id` |
| `pipelines/LargeDataset/brief_template.txt` | `source_file` |
| `pipelines/LargeDataset/pipeline.md` Brief Parameters Reference | `source_file` |

**Assessment:** Three files define a default for the same parameter with two different values. `lab_context.md` represents the intended lab standard; `brief_template.txt` and `pipeline.md` are vestiges of the HumanFat project where `source_file` was the batch variable. The HumanFat-specific default has propagated into the template and now contradicts the lab standard. `source_file` should be removed from the template default and replaced with `sample_id` or a blank prompt.

---

### 3.4 AMBIENT_PATTERNS Regex List

**Locations:**
- `differential_expression.md` lines 82–92 (canonical, with comments)
- `anatomical_de_analysis.md` lines 341–342 (copy, no comments)

**Identical or drifted?** The two regex lists are identical: `"^IGK","^IGL","^IGH","^HBA","^HBB","^HBD","^HBE","^HBG","^HBM","^HBQ","^HBZ","^MT-","^RPS","^RPL"`. No drift between these two copies. But the EXPLICIT list is entirely absent from `anatomical_de_analysis.md` — see Section 2.1.

---

### 3.5 EC Color Palettes

**Locations:**
- `lab_context.md` lines ~135–155: `ec_subtype_colors_adipose` and `ec_subtype_colors_kidney` (canonical definitions as R named vectors)
- `lab_context/validated_examples.yaml` fat-scrnaseq-continued entry: Same palette duplicated in YAML key_gotchas
- `pipelines/LargeDataset/methods/aesthetics.md` proportion plot section: `ec_colors` referenced but values inline in some plot code
- `pipelines/IntegratePublicData/methods/cross_atlas_dotplot.md` Section D: References `ec_colors` without definition

**Assessment:** `lab_context.md` is the canonical source. The duplication in `validated_examples.yaml` serves a documentation purpose (what the run used). Cross-file references that assume `ec_colors` is defined in context are fragile — if the project CLAUDE.md omits the lab_context color block, the scripts will fail with `object 'ec_colors' not found`.

---

### 3.6 `functional_gene_sets` EC Biology Genes

**Locations:**
- `differential_expression.md` lines 119–167: Full 4-section EC gene set (~80 genes)
- `metabolic_profile.md` PPARG target gene set: Overlaps with "Metabolic & Specialized" section in differential_expression.md — genes include PPARG, FABP4, FABP5, SCD, LPL, DGAT1, DGAT2, CD36 which appear in both files

**Assessment:** Not verbatim duplication — `metabolic_profile.md` has a focused PPARG signature while `differential_expression.md` has the broader metabolic category. But the overlap (~8 genes) means a migration agent curation session should decide which file is authoritative for the common genes.

---

## Summary Matrix

| # | Duplication | Files involved | Type | Verdict |
|---|---|---|---|---|
| 1.1 | Trajectory plot recipe | trajectory_monocle3.md, aesthetics.md | Verbatim code | Remove from aesthetics.md, add cross-reference |
| 1.2 | JoinLayers call | seurat_v5_rules.md + 5 others | Instruction repeat | Keep in seurat_v5_rules.md only; inline versions → comment |
| 1.3 | UMAP reduction detection | seurat_v5_rules.md, feature_umap_plot.md | Near-verbatim | Feature_umap has extra context; keep both but add cross-ref |
| 1.4 | r-env reminder | r_environment.md + 6 others | Instruction repeat | Stale in r_environment.md (missing --no-capture-output) |
| 1.5 | ggsave useDingbats | aesthetics.md + 4 others | Instruction repeat | Consistent except celltype_subclustering.md uses pdf() |
| 2.1 | Ambient filter 3 variants | differential_expression.md, anatomical_de_analysis.md, de_comprehensive_csv.md | Diverged concept | **High priority fix** — silent result differences |
| 2.2 | Bulk concordance distinction | bulk_concordance.md, bulk_lfc_concordance_heatmap.md | Intentional | Keep as-is |
| 2.3 | CellChat v2 bug list | interactome_cellchat.md, lab_context.md | Partial duplicate | interactome_cellchat.md is canonical |
| 2.4 | @meta.data rule | seurat_v5_rules.md + 2 others | Cross-reference | Fine |
| 2.5 | Harmony processing sequence | LargeDataset/pipeline.md, IntegratePublicData/pipeline.md | Near-identical | Candidate for shared primitive |
| 3.1 | Diverging color palette | 6 files | Drifted config | **High priority** — define single canonical palette |
| 3.2 | Gene expression color scale | shared/aesthetics.md, feature_umap_plot.md, dot plots | Conflicting advice | shared/aesthetics.md advice is outdated |
| 3.3 | batch_correction_var default | lab_context.md, brief_template.txt, pipeline.md | Drifted config | **High priority** — brief_template.txt must be corrected |
| 3.4 | AMBIENT_PATTERNS list | differential_expression.md, anatomical_de_analysis.md | Near-identical copy | Acceptable; EXPLICIT omission is the real issue |
| 3.5 | EC color palettes | lab_context.md, validated_examples.yaml, aesthetics.md, cross_atlas_dotplot.md | Scattered | lab_context.md is canonical; others are references |
| 3.6 | Metabolic gene sets | differential_expression.md, metabolic_profile.md | Partial overlap | Curate shared PPARG gene list |
