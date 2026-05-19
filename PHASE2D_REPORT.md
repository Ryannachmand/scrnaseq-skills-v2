# Phase 2 Group D Report

**Module authored:** `modules/cellchat.md`
**Commit:** Phase 2 Group D (commit b16d9a0)
**Date:** 2026-05-19

---

## 1. File Statistics

| Metric | Value |
|---|---|
| v2 file path | `modules/cellchat.md` |
| v2 line count | 1,320 |
| Source v1 file | `~/claude-skills/pipelines/LargeDataset/methods/interactome_cellchat.md` |
| v1 line count | 1,302 |

---

## 2. Parameterization Table

| v1 hardcoded value | v2 parameter name |
|---|---|
| `"celltype"` / `"ec_subtype"` / `"mylabel"` (cell type column) | `LABEL_COL` (`project_specific` sentinel) |
| `"tissue_type"` (group/condition column) | `GROUP_COL` (optional; `null` = skip Script 6) |
| `CellChatDB.human` (hardcoded organism) | `ORGANISM` ("human" or "mouse") |
| KidneyNew pathway category lists (PTPRM, PECAM1, CD99, etc.) | `PATHWAY_CATEGORIES` (named list; `project_specific` sentinel) |
| KidneyNew `chord_paths` (subset of paths excluding dominant pathway APP) | `chord_paths` key inside `PATHWAY_CATEGORIES[[cat_id]]` |
| KidneyNew category gradient colors (`"#7B3FA0"`, `"#E53935"`, `"#3D70BE"`) | `grad_high` key inside `PATHWAY_CATEGORIES[[cat_id]]` |
| `ec_senders` / `env_receivers` (KidneyNew EC + environment cell vectors) | `SOURCE_TYPES`, `TARGET_TYPES` (`project_specific` sentinels) |
| `cell_colors` named vector (KidneyNew EC + stroma/epithelial/MoMF) | `CELL_COLORS` (`project_specific` sentinel) |
| `ec_colors` named vector (HumanFat EC subtype palette) | `SOURCE_COLORS` (auto-derived as `CELL_COLORS[SOURCE_TYPES]`) |
| HumanFat `mac_types`, `aspc_types`, `aspc16` (partner cell type groupings) | `TARGET_TYPES` (`project_specific` sentinel, split as needed) |
| HumanFat `tissue_type` comparison (LiposuctionFat vs BreastFat) | `group_col` + `GROUP_COL` (`project_specific` sentinel) |
| HumanFat `cat_colors` (Angiocrine="#FF9800", Metabolic="#4CAF50", etc.) | `CAT_COLORS` (named vector; `project_specific` sentinel) |
| `output/` / `output2/` (project-specific output directories) | `OUTPUT_DIR` (`project_specific` sentinel; default `output/cellchat/`) |
| `"output/CellChat_all_interactions.csv"` (hardcoded filename) | `file.path(OUTPUT_DIR, "CellChat_all_interactions.csv")` |
| Script 6 tissue names in `extract_tissue()` calls (LiposuctionFat, BreastFat) | `group_val` argument passed per-call (`project_specific`) |
| `STACKED_CATS` omitting ECM category (HumanFat decision) | `STACKED_CATS` variable (`project_specific`; documented as a per-project exclusion decision) |
| `N_PER_SUBTYPE = 1500` (not relevant; that is atlas_co_umap) | N/A — not present in CellChat |

---

## 3. Context Dependency Declarations

```yaml
requires_context:
  palettes:
    - cell_colors      # required: named vector mapping all cell types to hex colors
    - source_colors    # optional: subset of cell_colors for source cell types
  metadata_columns:
    required:
      - label_col      # cell type label column in Seurat metadata
    optional:
      - group_col      # condition/tissue column (Script 6 only; null = skip)
  brief_keys:
    required:
      - output_dir
      - downstream_analyses.cellchat.label_col
      - downstream_analyses.cellchat.organism
    optional:
      - downstream_analyses.cellchat.group_col
      - downstream_analyses.cellchat.n_workers
      - downstream_analyses.cellchat.signaling_pathways
      - downstream_analyses.cellchat.source_celltypes
      - downstream_analyses.cellchat.target_celltypes
```

---

## 4. Primitives Referenced

| Primitive | Usage |
|---|---|
| `@primitives/seurat_v5_rules.md` | Rule 1 (JoinLayers before GetAssayData), Rule 5 (dynamic UMAP detection) |
| `@primitives/r_environment.md` | conda run with `--no-capture-output`, script structure |
| `@primitives/visualization.md` | Plot aesthetic conventions |
| `@primitives/aesthetics.md` | Typography and color philosophy |
| `@context/color_palettes.md` | Diverging palette for comparison heatmaps |

---

## 5. Six-Script Structure

All six scripts are present with their full analytical contributions preserved:

| Script | Filename | Analytical contribution |
|---|---|---|
| 1 | `01_cellchat_inference.R` | Creates CellChat object, runs full inference pipeline (identifyOverExpressedGenes → computeCommunProb → filterCommunication → aggregateNet), saves Rds and all-interactions CSV. Writes pathway count checkpoint output. |
| 2 | `02_cellchat_plots.R` | Loads saved object, runs per-category plots (bubble, circle, chord, heatmap) + global plots (all-pathway circle, rankNet, signaling role scatter, outgoing/incoming heatmaps). |
| 3 | `03_cellchat_bubbles.R` | Stacked composite bubble figure; rebuilds from `p$data` to override CellChat's internal scales; reorders LR pairs by dominant source cell type; patchwork-stacked with dynamic height sizing. |
| 4 | `04_cellchat_bars.R` | Pathway category distribution bar plots (count + strength versions) across all source × target partner combinations; fixed y-axis scale for comparability. |
| 5 | `05_cellchat_circos.R` | `make_circos()` function + 2 standard calls (source-to-target, target-to-source); cairo_pdf for Illustrator-editable text; two-panel layout (circos + legend); top-N pathways by interaction count. |
| 6 | `06_cellchat_comparison.R` | Conditional (only when `group_col` set); memory-efficient tissue extraction with caching; strict merge-before-centrality ordering; compareInteractions, netVisual_heatmap, rankNet, scatter; top-differential pathway filtering; same-axis N-condition scatter. |

---

## 6. CellChat v2 Bugs Preserved

All documented CellChat v2 bugs from v1 are carried forward with their workarounds:

| Bug | v1 documentation location | Workaround in v2 module |
|---|---|---|
| `plan("multisession")` crashes chord diagram rendering | Critical Constraints table | `plan("sequential")` at the very top of every script, before any library loads. Documented as Critical Constraint and in every script's CONFIG preamble. |
| `netVisual_aggregate(signaling=NULL)` crashes | Critical Constraints table | Use `netVisual_circle()` on `@net$count` instead. Documented in Critical Constraints and enforced in Script 2 global plot code. |
| `netVisual_heatmap()` annotation dimension bug when filtering pathways | Critical Constraints table | `plot_heatmap()` custom ggplot replacement for all per-category heatmaps. (Merged-object comparison heatmap in Script 6 is safe — no pathway subset applied.) |
| `pairLR.use=` validation precedence bug | Critical Constraints table | `get_filtered_net()` builds a pre-filtered `net` data frame; passed via `net=` parameter (never `pairLR.use=`). Module-level helper function. |
| `top=` parameter does not exist in `netVisual_chord_gene` | Critical Constraints table | `get_filtered_net(top_n=...)` pre-filters to top N LR pairs; result passed via `net=` parameter. |
| `cell.level = TRUE` removed from `mergeCellChat()` | Script 6 compare_pair code comment | Parameter removed; `compare_pair()` never passes `cell.level = TRUE`. Documented in Critical Constraints table. |
| `netAnalysis_computeCentrality()` throws `rowSums` on merged object | Script 6 compare_pair code comment | Run on each individual object inside the `cc_list` loop before calling `mergeCellChat()`. Enforced in `compare_pair()`. |
| `subsetCellChat()` required before merging (non-conformable arrays in netVisual_heatmap) | Script 6 pitfall table | `compare_pair()` finds common cell types and subsets each object before merge. |
| JoinLayers required before GetAssayData (Seurat v5 layer architecture) | Script 6 extract_tissue comment | `JoinLayers(sub_obj)` called in `extract_tissue()` before `GetAssayData()`. References `@primitives/seurat_v5_rules.md` Rule 1. |
| Empty / 3.6 KB PDF when saving ComplexHeatmap without `draw()` | Script 6 saving section comment | `ComplexHeatmap::draw(ht1 + ht2)` inside `pdf()` / `dev.off()`. Exception #1 documented. |

---

## 7. `get_filtered_net()` Placement

`get_filtered_net()` is defined as a **module-level helper** (not promoted to a primitive) per the inventory's classification as "RESOLVED (self-contained)." Reasons:

1. It exists specifically to bypass the `pairLR.use` validation precedence bug in CellChat v2. It has no utility outside of CellChat contexts.
2. The function depends on CellChat-specific object slots (`cellchat@LR$LRsig`, `slot(cellchat, "net")$prob`). Promoting it to a primitive would create a CellChat dependency in primitives/.
3. The inventory function_index.md Table 2 explicitly classifies it as "RESOLVED (self-contained)" — meaning it lives in the same file that uses it.

All other CellChat helper functions (`subset_net`, `plot_heatmap`, `reorder_bubble_by_ec`, `rebuild_bubble`, `assign_category`, `make_bar_data`, `make_bar_plot`, `save_both`, `make_circos`, `compare_pair`, `extract_tissue`, `run_cellchat_from_parts`, `pathway_stats`, `top_diff_pathways`) are similarly module-level because they are CellChat-specific or depend on CellChat-specific data structures.

---

## 8. Punch-List Resolution

### P2-6 — `projectData()` removed

**Status: RESOLVED.**

`projectData(cellchat, PPI.human)` does not appear anywhere in the v2 module's
executable R code. The function was removed in CellChat v2 and the v1 library
erroneously included it in Script 1's inference skeleton (contradicting the same
file's own Critical Constraints table).

The v2 module:
- Does NOT call `projectData()` in any code block
- Includes a comment at the exact location where `projectData()` would have appeared (between `identifyOverExpressedInteractions` and `computeCommunProb`) explicitly documenting that it was removed and why
- Lists `projectData(cellchat, PPI.human)` in the Critical Constraints table under "Do not use" with the explanation: "Function removed in CellChat v2. The v1 library erroneously included this call in its inference skeleton. Any v1 script that included projectData() must have that line removed before running on a v2 install."

---

## 9. Project-Specific Values Staged for Phase 4

### `examples/kidneynew_cellchat.md`

- `LABEL_COL <- "ec_subtype"` (or `"ec_subtype_final"` after post-hoc annotation)
- `ORGANISM <- "human"`
- Cell type vectors: EC-GC, EC-AEA, EC-AVR, EC-DVR, EC-LYM, EC-PTC + non-EC types
- Full `PATHWAY_CATEGORIES` list with 3 categories (Chemokine_Adhesion, Angiocrine, ECM) and their pathway name vectors + chord_paths subsets
- Note on APP exclusion from chord_paths (total_prob ~5x next strongest)
- EC subtype color palette (all 6 subtypes with validated hex values)
- Note that EC-PTC was added post-hoc; the existing CellChat object predates it
- Validated 2026-03-11

### `examples/humanfat_cellchat.md`

- `LABEL_COL <- "mylabel"` (EC subset); `GROUP_COL <- "tissue_type"`
- `ORGANISM <- "human"`
- Source types: AEC, CapEC, CapEC2, VenEC1, VenEC2, VenEC3 (exclude RibHighEC)
- Script 6 comparison pairs: LiposuctionFat vs BreastFat; 4-tissue variant
- EC subtype color palette (all 6 subtypes with validated hex values)
- Tissue comparison color convention (4 tissue types with hex values)
- `CAT_COLORS` for bar plots (7 categories with hex values, validated 2026-04-02)
- Script 3 stacked bubble structure: 3 panels, ECM excluded, separate Mac/ASPC plots
- Script 5 circos: 4 standard plots; two inference objects required (collapsed vs ASPC16)
- Validated 2026-03-30

---

## 10. Self-Check Grep Results

Grep command run against `modules/cellchat.md` with the full Group D term list.
All hits classified:

| Term | Line(s) | Location | Classification |
|---|---|---|---|
| `HumanFat` | 38, 1268, 1273, 1284, 1292, 1308, 1310, 1312 | Line 38: module header validation prose; lines 1268–1312: Phase 4 staging section | Expected — validation documentation (line 38) and Phase 4 staging |
| `#AEC7E8` (matched "AEC") | 396 | `make_circos()` `qual_pal` vector | **False positive** — `#AEC7E8` is a hex color code (ColorBrewer light blue), not the cell type "AEC". Not a leakage issue. |
| `EC-GC`, `EC-AEA`, `EC-AVR`, `EC-DVR`, `EC-LYM`, `EC-PTC` | 1226–1257 | Phase 4 staging section | Expected — Phase 4 staging |
| `mylabel` | 1271 | Phase 4 staging section | Expected — Phase 4 staging |
| `tissue_type` | 1274 | Phase 4 staging section | Expected — Phase 4 staging |
| `AEC`, `CapEC`, `VenEC` | 1276–1290 | Phase 4 staging section | Expected — Phase 4 staging |
| `LiposuctionFat`, `BreastFat` | 1281 | Phase 4 staging section | Expected — Phase 4 staging |

**Zero hits in executable R code** for any prohibited term. The `#AEC7E8` color code false positive is the only non-Phase-4 hit, and it is clearly a hex color, not a cell type reference.

---

## 11. pdf()/dev.off() Exceptions Used

| Script | Usage | CONVENTIONS.md §4 exception |
|---|---|---|
| Script 2 (`netVisual_circle`, `netVisual_chord_gene`) | `pdf(...); ...; dev.off()` | Exception #3 — CellChat circle and chord diagram functions use base R graphics internally |
| Script 5 (`make_circos`) | `cairo_pdf(...); chordDiagram(...); dev.off()` | Exception #3 — circlize/base R graphics; `cairo_pdf()` is the recommended variant for Illustrator-editable text (produces proper text objects rather than paths) |
| Script 6 (comparison heatmap) | `pdf(...); ComplexHeatmap::draw(ht1 + ht2); dev.off()` | Exception #1 — ComplexHeatmap::draw() sends output to the active graphics device; ggsave() cannot intercept it |

Note on `cairo_pdf()` vs `pdf()` for exception #3: `cairo_pdf()` is a pdf-family graphics device that follows the same open/draw/close pattern as `pdf()`. Using `cairo_pdf()` instead of `pdf()` in the circos script produces Illustrator-editable text labels rather than baked-in paths. This is explicitly validated in the v1 source and carries forward as the recommended pattern for circos output.

All other plots in the module use `ggsave(..., device = "pdf", useDingbats = FALSE)` or `ggsave(..., device = cairo_pdf)` as appropriate.

---

## 12. Cross-File Findings

1. **`netVisual_heatmap()` annotation bug is context-dependent.** The bug triggers when filtering to a pathway subset (per-category plots). It does NOT trigger on the merged-object comparison heatmap in Script 6 (no pathway subset applied there). The v2 module documents both cases: `plot_heatmap()` custom ggplot for the per-category case; `netVisual_heatmap()` + `draw()` for the comparison case.

2. **`cairo_pdf()` vs `pdf()` distinction documented.** The v1 used `cairo_pdf()` for circos plots and `pdf()` for other base R outputs. The v2 module preserves this: `cairo_pdf()` in `make_circos()` (for Illustrator editable text), `pdf()` for `netVisual_circle` and `netVisual_chord_gene` in Script 2. This is documented in the exceptions table.

3. **`ggsave(..., device = cairo_pdf)` for patchwork multi-panel outputs.** Script 3 (stacked bubbles) and Script 6 (comparison scatter) use `ggsave(..., device = cairo_pdf)` for patchwork outputs — following the pattern validated in atlas_co_umap.md (PHASE2C_REPORT: "cairo_pdf is passed as a device to ggsave() — this is permitted").

4. **Module explicitly defines `make_bar_plot()`.** The v1 referenced `make_bar_plot()` inside `save_both()` but never defined it explicitly — it was described only in prose. The v2 module provides a complete, parameterized `make_bar_plot()` definition to close this gap.

5. **The 4-tissue same-axis scatter now uses `ggsave()`.** The v1 used `pdf(); print(wrap_plots(...)); dev.off()` for the 4-tissue scatter. The v2 module uses `ggsave(..., device = cairo_pdf)` for patchwork outputs — conformant with CONVENTIONS.md.

6. **CellChat v2 bug list in validated_examples.yaml KidneyNew notes is a partial duplicate.** The v2 module is now the authoritative location for the complete 10-bug list. The KidneyNew notes in validated_examples.yaml list 3 of the 10 bugs and will remain as project-context notes.

---

## 13. Open Items for Phase 3 and Phase 4

### Phase 4 (examples/)

1. **`examples/kidneynew_cellchat.md`**: needs `INPUT_RDS` path, non-EC cell type list with hex colors (stroma/epithelial/MoMF), and confirmation of which script filenames were used in the validated run.

2. **`examples/humanfat_cellchat.md`**: needs `INPUT_RDS` path, full `PATHWAY_CATEGORIES` list extracted from HumanFat project CLAUDE.md, `TARGET_TYPES` vectors (mac_types, aspc_types, aspc16 names), and script filenames from the validated run.

3. **Both examples**: need the complete `CELL_COLORS` named vector for all cell types (not just EC subtypes) to enable full color consistency across all CellChat plots.

### Phase 3 (deployment agent)

4. **Brief schema update**: The current `CONVENTIONS.md` brief schema shows a minimal `cellchat:` block. The v2 module defines a richer schema (`label_col`, `group_col`, `organism`, `n_workers`, `signaling_pathways`, `source_celltypes`, `target_celltypes`, `output_dir`). CONVENTIONS.md should be updated to reflect the richer schema.

5. **SKILL.md routing entry**: SKILL.md currently maps `downstream_analyses.cellchat` to `modules/interactome_cellchat.md`. This should be updated to `modules/cellchat.md`.

6. **Pathway categorization checkpoint pattern**: The module documents a CHECKPOINT between Script 1 and Script 2 where the user reviews pathways before writing plotting scripts. The deployment agent should be instructed to pause at this checkpoint and not auto-generate Scripts 2–6 before the user reviews the pathway list.

---

## 14. Unexpected Findings

1. **The v1 file contradicted itself on `projectData()`.** The v1 Critical Constraints table (line 31) explicitly said to remove `projectData(cellchat, PPI.human)`, but line 71 of the same file included `cellchat <- projectData(cellchat, PPI.human)` in the executable inference pipeline skeleton. The per_file_inventory.md and function_index.md both flagged this as an "INVALID" call. The v2 module removes the call entirely and documents why at the insertion point.

2. **`netVisual_heatmap()` has TWO different behaviors depending on context.** When called with a pathway subset (per-category), it crashes due to an annotation dimension bug → use `plot_heatmap()`. When called on a merged comparison object WITHOUT a pathway subset (Script 6), it works correctly → safe to use with `draw()`. The v1 documented the bug but did not clearly distinguish these two cases. The v2 module documents both explicitly.

3. **`make_bar_plot()` was a phantom function in v1.** The function was called from `save_both()` but never defined — its behavior was described only in prose and in the validated theme block (lines 718–753). This is analogous to the `is_confound()` phantom from de_comprehensive_csv.md. The v2 module provides a complete implementation.

4. **The stacked bubble script had an important `& coord_cartesian(clip="off")` pitfall.** The v1 noted this at the bottom of Script 3: "NOT needed when x-axis is at bottom — it causes the patchwork theme warning." This is subtle and easily missed. The v2 module documents it prominently as a note immediately after the `build_stacked()` function.

5. **Script 5 (`make_circos`) uses `cairo_pdf()` not `pdf()`.** This is the correct form for circos outputs (produces Illustrator-editable text). However, `cairo_pdf()` at the block level is still a pdf()/dev.off()-style pattern that requires documentation under exception #3. The CONVENTIONS.md exception table lists only `pdf()` but the intent clearly extends to `cairo_pdf()` (both are pdf-family graphics devices). The module documents this variant explicitly.
