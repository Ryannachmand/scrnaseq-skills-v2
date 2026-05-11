# Function Index — ~/claude-skills/
*Generated: 2026-05-11 by inventory agent*

Cross-cutting map of every named function in the library: definitions, call sites, inconsistencies, orphans, and phantoms.

**Scope note:** "Function" here means a named R function created with `function(...)`. Unnamed inline code blocks are excluded. Conceptual descriptions in prose that merely mention a function name (without a code definition or call) are also excluded from Table 1 and Table 2, but may be flagged in Table 3 notes where they create confusion.

---

## Table 1 — Every Function Definition

| Function name | Signature | Defined in | Line | Description |
|---|---|---|---|---|
| `is_ambient` | `is_ambient(genes)` | `pipelines/LargeDataset/methods/differential_expression.md` | 101 | Returns logical vector TRUE for ambient/contamination genes; checks against AMBIENT_PATTERNS regex list AND AMBIENT_EXPLICIT named list. |
| `is_ambient` | `is_ambient(genes)` | `pipelines/IntegratePublicData/methods/anatomical_de_analysis.md` | 344 | Returns logical vector TRUE for ambient genes; checks AMBIENT_PATTERNS regex list ONLY. Missing AMBIENT_EXPLICIT check — see Table 3. |
| `run_findmarkers` | `run_findmarkers(ec_obj, comp)` | `pipelines/LargeDataset/methods/differential_expression.md` | 187 | Runs `FindMarkers()` with wilcox test, `logfc.threshold=0`, `min.pct=0.05`; skips comparisons with <10 cells per group; returns tibble with `gene` column. |
| `make_volcano` | `make_volcano(markers, comp, subset_name, n_label = N_LABEL_VOLCANO)` | `pipelines/LargeDataset/methods/differential_expression.md` | 223 | Full ggplot volcano with significance coloring (`#B2182B`/`#2166AC`), ggrepel labels; prioritises biological genes (from `functional_gene_sets`) in top-N selection; applies `is_ambient()` filter. |
| `make_overall_heatmap` | `make_overall_heatmap(ec_obj, markers, comp, subset_name, n_cells = N_CELLS_HM, n_genes = N_TOP_GENES_HM)` | `pipelines/LargeDataset/methods/differential_expression.md` | 285 | ComplexHeatmap cell-level z-score heatmap; downsamples to n_cells; column split by group; column annotation uses `ec_sub$mylabel` (hardcoded column name). |
| `make_functional_heatmap` | `make_functional_heatmap(ec_obj, markers, comp, subset_name, ...)` | `pipelines/LargeDataset/methods/differential_expression.md` | 333 | ComplexHeatmap restricted to `functional_gene_sets` genes; uses `row_split` for section gap lines. Full signature not visible from header (described in Step 4 prose). |
| `make_topgene_dotplot` | `make_topgene_dotplot(...)` | `pipelines/LargeDataset/methods/differential_expression.md` | ~415 | Single-panel dot plot of top DE genes per direction; % expressing + avg expression per EC subtype; faceted by tissue group. Exact signature not extracted — see Table 3. |
| `make_topgene_dotplot` | `make_topgene_dotplot(obj, markers, comp, subset_name, n_each = 12)` | `pipelines/IntegratePublicData/methods/anatomical_de_analysis.md` | 77 | **3-panel version** faceted by `c("Vertebrae","Iliac Crest","Femoral Head")`; uses `STROMA_ORDER` for x-axis; section divider between DE directions. Different implementation from `differential_expression.md` version — see Table 3. |
| `make_functional_dotplot` | `make_functional_dotplot(..., show_direction = FALSE)` | `pipelines/LargeDataset/methods/differential_expression.md` | 483 | Dot plot organized into functional gene set sections with dashed section dividers; optional patchwork direction strip when `show_direction = TRUE`; ident2 panel pins section labels. |
| `make_functional_dotplot` | `make_functional_dotplot(...)` | `pipelines/IntegratePublicData/methods/anatomical_de_analysis.md` | 142 | **3-panel version** with section labels pinned to `"Femoral Head"` panel (rightmost of 3 sites); different panel geometry from `differential_expression.md` version — see Table 3. |
| `make_pathway_barplot` | `make_pathway_barplot(markers, comp, subset_name, universe_genes, n_pathways = 15)` | `pipelines/LargeDataset/methods/differential_expression.md` | 594 | GO:BP over-representation analysis with `enrichGO()` + `simplify()`; produces diverging bar plot (up/down arms); uses `clusterProfiler` + `org.Hs.eg.db`. |
| `run_aucell` | `run_aucell(seurat_obj, gene_sets, assay = "RNA")` | `pipelines/LargeDataset/methods/metabolic_profile.md` | 62 | Builds AUCell cell rankings from count matrix; calculates AUC scores; returns data frame (cells × gene_set_names). |
| `add_auc_to_seurat` | `add_auc_to_seurat(seurat_obj, auc_df)` | `pipelines/LargeDataset/methods/metabolic_profile.md` | 70 | Adds AUC score columns from `auc_df` to Seurat object metadata via `AddMetaData()`; handles cell barcode intersection gracefully. |
| `stratified_downsample` | `stratified_downsample(so, ceiling, label_col = "unified_label")` | `pipelines/IntegratePublicData/pipeline.md` | 220 | Stratified random sampling of Seurat object preserving cell type proportions; applies ceiling globally; never upsamples; returns subsetted Seurat object. |
| `get_gene_fpkm` | `get_gene_fpkm(tbl, gene)` | `pipelines/IntegratePublicData/methods/load_formats.md` | 264 | Extracts gene-level FPKM from a Cufflinks tmap table; prefers exact-match transcripts (`class_code == "="`); falls back to all transcripts if none found. |
| `get_filtered_net` | `get_filtered_net(cellchat, paths, min_cells = 10)` | `pipelines/LargeDataset/methods/interactome_cellchat.md` | ~Scripts 3–4 section | Filters CellChat interaction table to a specified pathway list with minimum cell count; bypasses the `pairLR.use` bug in CellChat v2. Defined inline (no explicit line number in header). |

---

## Table 2 — Every Function Call Site

| Function name | Called from | File | Approx. line | Resolves to definition in | Status |
|---|---|---|---|---|---|
| `is_ambient` | `make_volcano()` body | `differential_expression.md` | 242 | `differential_expression.md` line 101 | RESOLVED |
| `is_ambient` | Main analysis loop | `differential_expression.md` | ~650 | `differential_expression.md` line 101 | RESOLVED |
| `is_ambient` | Annotation overlay section | `bulk_concordance.md` | ~120 | `differential_expression.md` line 101 | RESOLVED (requires `differential_expression.md` in context) |
| `is_ambient` | `make_topgene_dotplot()` body | `anatomical_de_analysis.md` | 80 | `anatomical_de_analysis.md` line 344 (local) or `differential_expression.md` line 101 | AMBIGUOUS — uses local definition (simpler; missing EXPLICIT list) |
| `is_confound` | `sig_genes` filter | `de_comprehensive_csv.md` | 116 | *Nowhere* | **PHANTOM — never defined** |
| `run_findmarkers` | Main analysis loop | `anatomical_de_analysis.md` | 323 | `differential_expression.md` line 187 | **UNRESOLVED** — `differential_expression.md` not injected for IntegratePublicData jobs |
| `make_volcano` | Main analysis loop | `anatomical_de_analysis.md` | 325 | `differential_expression.md` line 223 | **UNRESOLVED** — same injection gap |
| `make_overall_heatmap` | Main analysis loop | `anatomical_de_analysis.md` | 326 | `differential_expression.md` line 285 | **UNRESOLVED** |
| `make_functional_heatmap` | Main analysis loop | `anatomical_de_analysis.md` | 327 | `differential_expression.md` line 333 | **UNRESOLVED** |
| `make_pathway_barplot` | Main analysis loop | `anatomical_de_analysis.md` | 331 | `differential_expression.md` line 594 | **UNRESOLVED** |
| `stratified_downsample` | Stage 4 example | `IntegratePublicData/pipeline.md` | ~240 | Same file, line 220 | RESOLVED (self-referential) |
| `get_filtered_net` | Scripts 3–4 chord diagram code | `interactome_cellchat.md` | ~Script 3 section | Same file | RESOLVED (self-contained) |
| `run_aucell` | Main metabolic scoring code | `metabolic_profile.md` | ~100 | Same file, line 62 | RESOLVED (self-contained) |
| `add_auc_to_seurat` | Main metabolic scoring code | `metabolic_profile.md` | ~110 | Same file, line 70 | RESOLVED (self-contained) |
| `JoinLayers` | Post-merge step | `load_formats.md` | 159 | Seurat v5 external library | RESOLVED (external) |
| `JoinLayers` | Stage 5 integration step | `IntegratePublicData/pipeline.md` | 264 | Seurat v5 external library | RESOLVED (external) |
| `JoinLayers` | Phase 1 prep | `celltype_subclustering.md` | 222 | Seurat v5 external library | RESOLVED (external) |
| `JoinLayers` | Script 1 export | `pyscenic_regulon_analysis.md` | 84 | Seurat v5 external library | RESOLVED (external) |
| `AUCell_buildRankings` | `run_aucell()` body | `metabolic_profile.md` | 65 | AUCell external library | RESOLVED (external) |
| `enrichGO` | `make_pathway_barplot()` body | `differential_expression.md` | ~610 | clusterProfiler external library | RESOLVED (external) |
| `chordDiagramFromMatrix` | Chord diagram recipe | `cohort_plots.md` | ~135 | circlize external library | RESOLVED (external) |
| `projectData` | Inference pipeline skeleton | `interactome_cellchat.md` | 71 | CellChat v2 (**REMOVED in CellChat v2**) | **INVALID** — function does not exist; contradicts Critical Constraints in same file |

---

## Table 3 — Mismatched / Inconsistent Definitions

### 3.1 `is_ambient(genes)` — Two Diverged Implementations

**Definition A** — `differential_expression.md` lines 101–104:
```r
is_ambient <- function(genes) {
  pattern_hit <- Reduce(`|`, lapply(AMBIENT_PATTERNS, function(p) grepl(p, genes)))
  genes %in% AMBIENT_EXPLICIT | pattern_hit
}
```
Requires both `AMBIENT_PATTERNS` (regex vector) AND `AMBIENT_EXPLICIT` (named character vector of 13 specific genes: ALB, FGA, FGB, FGG, FGL1, PF4, PPBP, GP9, GP1BA, GP1BB, ITGA2B, MALAT1, NEAT1, TPSAB1, TPSB2, CPA3).

**Definition B** — `anatomical_de_analysis.md` lines 344–346:
```r
is_ambient <- function(genes) {
  Reduce(`|`, lapply(AMBIENT_PATTERNS, function(p) grepl(p, genes)))
}
```
Uses AMBIENT_PATTERNS only. Requires only the regex vector. No explicit gene list.

**Diff:** Definition B is a strict subset of Definition A. Any gene in AMBIENT_EXPLICIT (e.g., `ALB`, `MALAT1`, `PPBP`) will NOT be filtered by Definition B's version if its symbol doesn't match any regex. Concretely: `MALAT1` matches no regex in AMBIENT_PATTERNS, so `is_ambient("MALAT1")` returns FALSE under Definition B but TRUE under Definition A. This means `anatomical_de_analysis.md`'s ambient filter will pass through abundant nuclear lncRNAs (MALAT1, NEAT1), plasma proteins (ALB, FGA, etc.), and platelet markers (PF4, PPBP) that Definition A catches. File header comments in `anatomical_de_analysis.md` line 338 say "(same as LargeDataset pipeline)" — this claim is factually incorrect.

**Risk:** Any DE output from `anatomical_de_analysis.md`-based analyses may include MALAT1, ALB, and platelet genes in its labeled points.

---

### 3.2 `make_topgene_dotplot(...)` — Same Name, Incompatible Implementations

**Definition A** — `differential_expression.md` ~line 415 (Step 5 description):
- Signature: not explicitly shown; uses `ec_sub$mylabel` for x-axis grouping
- Layout: single-panel plot; x-axis = EC subtype (`mylabel`); y-axis = genes; faceted by tissue group (`tissue_type`)
- Top-gene selection: top N per comparison direction
- Section divider: between up_ident1 and up_ident2 genes

**Definition B** — `anatomical_de_analysis.md` line 77:
- Signature: `make_topgene_dotplot(obj, markers, comp, subset_name, n_each = 12)`
- Layout: **3-panel** faceted by `c("Vertebrae","Iliac Crest","Femoral Head")`; x-axis = `unified_label` (cell type); no tissue facet
- Top-gene selection: `n_each = 12` top per direction using `PADJ_CUT` and `LFC_CUT`
- Section divider: `geom_hline(yintercept = divider_y, linetype = "dashed")`
- Axis color: `element_text(color = label_colors[ct_present])` — requires `label_colors` named vector

**Diff:** These share a name but are fundamentally different figures. Definition A is a "which subtypes express this gene across tissue conditions?" plot. Definition B is a "which gene is up in which anatomical site?" plot. They have different data pipelines, different faceting strategies, and different x-axis semantics. Both name their section divider with `geom_hline` but the y-intercept calculation differs. A migration agent cannot unify these without losing the unique contribution of each.

**Recommendation:** Rename B to `make_anatomical_dotplot()` in v2 to eliminate the collision.

---

### 3.3 `make_functional_dotplot(...)` — Same Name, Different Panel Count and Label Logic

**Definition A** — `differential_expression.md` line 483:
- Shows functional gene sets as sections; section labels pinned to rightmost panel (ident2 panel in a 2-comparison facet)
- `show_direction = FALSE` flag suppresses the patchwork direction strip
- Section labels avoid `annotate()` (renders in all panels); uses `geom_text` with ident2 pinning

**Definition B** — `anatomical_de_analysis.md` line 142:
- Same section-based structure but for **3 panels** (Vertebrae | Iliac Crest | Femoral Head)
- Section labels pinned to `"Femoral Head"` (rightmost of the 3) using `factor("Femoral Head", levels=c(...))`
- `coord_cartesian(clip="off")` + `plot.margin = margin(t=8, r=120, b=8, l=8)` for right overflow

**Diff:** The core insight (geom_text + panel pinning to avoid annotate()) is shared. The difference is purely structural: 2-panel (ident1/ident2) vs 3-panel (3 sites). Definition A's `show_direction` flag has no equivalent in B. The label-pinning expression differs (`comp$ident2` vs hardcoded `"Femoral Head"`). These could be unified into a single function with a `pin_to` argument and `n_panels` generalization.

---

## Table 4 — Orphan Functions (defined but never called from another file)

| Function name | Defined in | Reason for orphan status |
|---|---|---|
| `run_aucell` | `metabolic_profile.md` line 62 | Only called within the same file's own code block; no other file calls it. The function is well-designed and portable — a strong primitive candidate. |
| `add_auc_to_seurat` | `metabolic_profile.md` line 70 | Same situation; helper to `run_aucell`. Should travel with it in v2. |
| `get_gene_fpkm` | `load_formats.md` line 264 | Only referenced in inline loading examples within the same file. No pipeline or other methods file calls it. |
| `get_filtered_net` | `interactome_cellchat.md` ~Scripts 3–4 section | Defined inline in a code block within the CellChat methods file; used in the same file's chord diagram scripts. Never exposed to other files. |
| `make_functional_heatmap` | `differential_expression.md` line 333 | Called in the same file's main loop. Never exposed to or called from any other file. |
| `make_overall_heatmap` | `differential_expression.md` line 285 | Same: self-contained within `differential_expression.md`'s own loop. Technically called from `anatomical_de_analysis.md` but that call is UNRESOLVED (not injected). |
| `make_pathway_barplot` | `differential_expression.md` line 594 | Same: called from `anatomical_de_analysis.md` but UNRESOLVED. Otherwise orphaned. |

**Note:** `run_findmarkers`, `make_volcano`, `make_topgene_dotplot`, `make_functional_dotplot` from `differential_expression.md` are all called from `anatomical_de_analysis.md` — but those calls are UNRESOLVED because `differential_expression.md` is never injected into IntegratePublicData jobs. So for practical library purposes they are effectively orphans outside their own file.

---

## Table 5 — Phantom Functions (called but never defined)

| Function name | Called from | Line | Note |
|---|---|---|---|
| `is_confound` | `de_comprehensive_csv.md` | 116 | **Critical gap.** The file's own "Exclusion filters" section (line ~215) says "Apply `is_confound()` (defined in `differential_expression.md`)" — but `differential_expression.md` defines `is_ambient()`, not `is_confound()`. The intent appears to be a superset filter including sex-linked genes, HLA class II, histone genes, and unannotated ENSG IDs — none of which are covered by `is_ambient()`. This function either (a) was planned but never written, (b) was written in a project-level CLAUDE.md and never moved to the library, or (c) was an informal alias for `is_ambient()` with a broader name. Until resolved, any script generated from `de_comprehensive_csv.md` will produce an R error at the `sig_genes` filter step. |
| `projectData` | `interactome_cellchat.md` | 71 | Not a phantom in the phantom-function sense — this is a real CellChat function that was **removed in CellChat v2**. It appears in an inference pipeline code block that contradicts the same file's Critical Constraints section which says to REMOVE this call. Calling it on a v2 install will produce `Error: could not find function "projectData"`. |
