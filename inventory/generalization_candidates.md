# Generalization Candidates — ~/claude-skills/
*Generated: 2026-05-11 by inventory agent*

For each project-specific file and each section of non-orphan files that contains generalizable patterns, this document describes: what the file is fundamentally doing, how to generalize it, what inputs would need to be parameterized, what stays in an example file, and honest risk assessment.

---

## G1 — `anatomical_de_analysis.md` (IntegratePublicData/methods/)

### What it's fundamentally doing
Running differential expression between N anatomical sites (or any discrete groups), visualizing results as (a) per-direction top-gene dot plots faceted by group and cell type, (b) functional dot plots grouped into biological sections, and (c) GO:BP enrichment. The novel contribution vs `differential_expression.md` is the **N-group comparison framework** with cell-type-resolved visualization that disentangles DE composition effects from expression effects.

### Proposed primitive/module name
`modules/multi_group_de_analysis.md`

### Proposed signature
```
multi_group_de_analysis(
  obj,                  # Seurat object with a `label_col` and a `group_col`
  label_col,            # metadata column holding cell type labels (e.g., "unified_label")
  group_col,            # metadata column holding group variable (e.g., "site", "condition")
  groups,               # character vector of group levels in display order (N ≥ 2)
  comparisons,          # list of list(label, ident1, ident2) — pairs to compare
  functional_gene_sets, # named list of gene vectors for functional sections
  n_each = 12,          # top genes per direction in topgene dotplot
  padj_cut = 0.05,
  lfc_cut  = 0.5,
  output_dir            # output directory path
)
```

### What needs parameterization
| Currently hardcoded | Generalized parameter |
|---|---|
| `c("Vertebrae","Iliac Crest","Femoral Head")` factor levels | `groups` argument |
| `STROMA_ORDER` for cell type ordering | `label_order` argument (optional; default = sorted) |
| `"Femoral Head"` as rightmost panel for label pinning | `pin_to = groups[length(groups)]` (auto-derived) |
| `"site"` as the group column | `group_col` argument |
| `"unified_label"` as the cell type column | `label_col` argument |
| 3 facets hardcoded in `facet_wrap(~ site_group, ncol = 3)` | `ncol = length(groups)` |
| `label_colors` named vector (requires STROMA_ORDER and palette) | `label_colors` argument (optional) |
| `is_ambient()` (simpler version) | Use the full version from `differential_expression.md` |

### What stays in the example file
`examples/bone_marrow_stroma_de.md`:
- The specific `groups = c("Vertebrae","Iliac Crest","Femoral Head")` instantiation
- The `STROMA_ORDER` cell type ordering
- `functional_gene_sets` for bone marrow stromal biology (MSC, osteogenic, adipogenic gene sets)
- The 3-panel output figure description with expected result documentation

### Risk assessment
1. **`run_findmarkers()` dependency is unresolved.** The module calls `run_findmarkers()`, `make_volcano()`, `make_overall_heatmap()`, `make_functional_heatmap()`, `make_pathway_barplot()` which are all defined in `differential_expression.md`. Generalization requires resolving the cross-pipeline injection problem — these functions must be accessible. Either merge the two files into a shared `modules/de_analysis_core.md` or explicitly inject `differential_expression.md` for IntegratePublicData jobs.
2. **N-group assumption in plot layout.** `make_topgene_dotplot()` computes `divider_y` assuming exactly 2 DE directions (up_ident1 vs up_ident2). For 3-way comparisons (A vs B vs C simultaneously), the section divider logic would need extension. The current implementation only handles pairwise comparisons within the N-group object.
3. **`label_colors` coupling.** The axis color feature (`element_text(color = label_colors[ct_present])`) requires a pre-existing named color vector. This is aesthetically valuable but creates a hard dependency on project-specific palettes. Should be optional in the general version.
4. **`AMBIENT_EXPLICIT` absence.** The local `is_ambient()` misses MALAT1 and platelet markers. Any generalized version must use the full dual-filter.

---

## G2 — `co_umap_embedding.md` (IntegratePublicData/methods/)

### What it's fundamentally doing
Producing a co-embedding UMAP for comparing in-house cells with a public atlas, with downsampling to prevent atlas domination of the geometry. Saves coordinates as CSV for fast replotting (separation of computation from visualization). Produces a 4-panel highlight figure.

### Proposed primitive/module name
`modules/atlas_co_umap.md`

### Proposed signature
```
atlas_co_umap(
  so,                    # merged Seurat object (in-house + atlas)
  source_col,            # metadata column identifying dataset source (e.g., "dataset")
  inhouse_label,         # value of source_col for in-house data (e.g., "WCM")
  atlas_label,           # value of source_col for atlas data (e.g., "TS")
  subtype_col,           # cell type / subtype column for in-house cells
  atlas_group_col,       # organ/tissue column for atlas cells
  highlight_inhouse,     # character: which in-house subtype to highlight in Panel D
  highlight_atlas,       # character: which atlas group to highlight in Panel D
  n_per_subtype = 1500,  # cells to sample per in-house subtype before UMAP
  coords_csv,            # path to save/load UMAP coordinate CSV
  output_dir
)
```

### What needs parameterization
| Currently hardcoded | Generalized parameter |
|---|---|
| `"WCM"`, `"TS"` as dataset names | `inhouse_label`, `atlas_label` |
| `"CapEC"` highlighted subtype | `highlight_inhouse` |
| `"Fat EC"` highlighted atlas organ | `highlight_atlas` |
| `CB_CAPEC = "#E69F00"`, `CB_FAT_EC = "#0072B2"` | `color_inhouse`, `color_atlas` (default to Okabe-Ito) |
| `N_PER_SUBTYPE = 1500` | `n_per_subtype` with documented rationale (≥1500 for 80-85% atlas geometry target) |
| Output filenames referencing `ts_harmony` | Derived from `inhouse_label` and `atlas_label` |
| `ec_colors` (HumanFat palette) | `subtype_colors` argument (optional) |
| `output2/` directory | `output_dir` argument |

### What stays in the example file
`examples/tabula_sapiens_co_umap.md`:
- `highlight_inhouse = "CapEC"`, `highlight_atlas = "Fat EC"` from TabulaSapiensComparison
- WCM/TS naming
- The specific `N_PER_SUBTYPE = 1500` validation with rationale (from ratio analysis)
- Script filenames `06f`, `06g` from the TabulaSapiensComparison project

### Risk assessment
1. **`n_per_subtype` may need tuning.** The 1500/subtype recommendation was validated for the WCM/TS size ratio (in-house ~6 subtypes × 38k cells, atlas ~100k cells). For very different size ratios, the optimal number changes. The general version should document this as a tunable parameter with the validation context.
2. **UMAP geometry sensitivity to sampling.** Minor changes to `n_per_subtype` produce visually different UMAPs. The function should warn if the sampling produces <500 cells for any subtype.
3. **Custom fonts.** `co_umap_embedding.md` documents a font error (`"Source Sans 3"`) in its own pitfalls section. The general version must default to system fonts and document custom font as an optional override.
4. **Cairo PDF requirement.** The `device = cairo_pdf` pattern for patchwork is a documented fix. If cairo is not installed, this silently falls back to `pdf()`. Should be noted in the general version.

---

## G3 — `cross_atlas_dotplot.md` (IntegratePublicData/methods/)

### What it's fundamentally doing
Cross-dataset comparison dotplot with depth correction: z-scores expression within each dataset independently before combining, producing a 3-section dotplot (marker genes, TF program, curated biology) that is comparable across datasets despite sequencing depth differences. Also produces a centroid correlation heatmap in Harmony space.

### Proposed primitive/module name
`modules/cross_dataset_dotplot.md`

### Proposed signature
```
cross_dataset_dotplot(
  so,                    # integrated Seurat object
  source_col,            # metadata column identifying dataset source
  subtype_col,           # cell type column
  sources_order,         # character vector: order of dataset sources (for column grouping)
  marker_genes,          # named list of gene sections (Section 1 = known markers, etc.)
  n_variable = 50,       # top variable genes for centroid correlation
  z_cap = 2,             # z-score cap for depth correction
  col_group_col,         # metadata column for column background highlighting
  output_dir
)
```

### What needs parameterization
| Currently hardcoded | Generalized parameter |
|---|---|
| `"WCM Fat CapEC"`, `"WCM Fat VenEC1"` column label convention | `"{source_col} {subtype}" concat pattern` — configurable separator |
| `EC_SUBTYPES` WCM subtype vector | `inhouse_subtypes` argument |
| `"CapEC"` marker gene selection | `marker_genes` argument (named list) |
| `SEC1_FORCE` vector (hardcoded Section 1 TF genes) | `sec1_force` argument (optional override) |
| Background highlight colors (warm yellow WCM, light red TS Fat) | `group_colors` named vector |
| Custom fonts (`"Source Sans 3"`, `"Playfair Display"`) | Default to `""` (system); optional override |
| `allMarkersEC.csv`, `fat_enriched_markers.Rds` input files | Input arguments |

### What stays in the example file
`examples/tabula_sapiens_dotplot.md`:
- WCM/TS column naming convention
- HumanFat-specific marker gene selections
- Specific validated figure dimensions

### Risk assessment
1. **Depth correction assumes datasets are comparable post-z-scoring.** If the in-house and atlas datasets use different gene panels or have very different transcript detection rates, within-dataset z-scoring may not fully correct the batch effect. The general module should document this assumption.
2. **`SEC1_FORCE` pattern is a project hack.** Forcing specific TF genes into Section 1 because they appear in one dataset but not the other is an analysis decision that needs explicit justification in each project. The general version should flag this as project-specific configuration.
3. **Custom fonts will cause errors.** Both fonts (`"Source Sans 3"`, `"Playfair Display"`) are documented as unavailable. The general version must default to no custom font.
4. **Column ordering is fragile.** The `column_order` logic depends on `cor_mat` column names surviving whitespace→dot conversion. Any project with cell type names containing special characters will break the ordering. Must be documented and handled.

---

## G4 — `bulk_concordance.md` (LargeDataset/methods/)

### What it's fundamentally doing
Three-part framework for cross-referencing in-vivo scRNA-seq with a bulk perturbation experiment: (1) annotate existing DE outputs with bulk concordance status, (2) score individual cells for whole-signature resemblance via AddModuleScore, (3) focus on TF program specifically to identify TFs driving concordance. The pattern is applicable to any bulk perturbation × scRNA-seq dataset pair.

### Proposed primitive/module name
`modules/bulk_perturbation_concordance.md`

### What needs parameterization
| Currently hardcoded | Generalized parameter |
|---|---|
| `BULK_CSV` pointing to PPARG OE CSV | `bulk_csv` path |
| `ec_colors`, `mylabel`, `tissue_type`, `tissue_colors` | `subtype_col`, `subtype_colors`, `group_col`, `group_colors` |
| `RibHighEC` exclusion pattern | `exclude_subtypes` vector (optional) |
| `PPARG` as TF of interest in Part 3 | `target_tf` argument |
| `pparg_concordance` as score column name | Derived from `target_tf` + `"_concordance"` |
| `"output3/"` directory | `output_dir` argument |
| Bulk experiment description in comments | `experiment_label` argument for plot titles |

### What stays in the example file
`examples/pparg_concordance.md`:
- PPARG OE experiment context and bulk CSV path
- HumanFat-specific tissue colors and subtype exclusions
- Validated figures from the 2026-03-12 run

### Risk assessment
1. **Part 3 (TF focus) only makes sense for TF overexpression experiments.** For bulk experiments that are not TF OE (e.g., drug treatments, disease vs control), Part 3 has no interpretation. The general module should mark Part 3 as conditional on `target_tf` being specified.
2. **Module score interpretation depends on gene set quality.** The bulk concordance score is only as meaningful as the bulk DE list is specific. The general module should include a validation step (score distribution check, expected cell-type enrichment pattern).
3. **`bulk_lfc_concordance_heatmap.md` overlap.** Both methods address the same scientific question from different angles. The general module should explicitly position itself relative to the parallel-LFC heatmap approach (as the current files already do). In v2, they should be in the same directory with a clear "choose this when" decision tree.

---

## G5 — `differential_expression.md` `functional_gene_sets` (LargeDataset/methods/)

### What it's fundamentally doing
The `functional_gene_sets` block (lines 119–167) is an EC-biology-specific named list of ~80 genes organized into 4 sections. It is injected into every script generated from this file. Currently it is embedded in the file and never configurable at the brief level — any non-EC project gets EC gene sets, which produce meaningless functional heatmaps.

### Proposed primitive
Move to `examples/ec_functional_gene_sets.md` as an example. Add to brief configuration:
```yaml
deg:
  functional_gene_sets: project_specific  # define inline in CLAUDE.md
```
The `differential_expression.md` brief configuration section (line 804) already says `functional_gene_sets: project_specific` — but this placeholder instruction needs a corresponding template in the `primitives/` directory showing how to define a custom set.

### What needs parameterization
The entire `functional_gene_sets` list. In v2, `differential_expression.md` should read:
```r
# functional_gene_sets: defined in project CLAUDE.md context block — see brief config
```
And the brief template should include an example skeleton.

### Risk assessment
Any project that doesn't define `functional_gene_sets` in the CLAUDE.md context block will fail at `make_functional_heatmap()` and `make_functional_dotplot()`. Must provide a clear default or fallback.

---

## G6 — `metabolic_profile.md` AUCell Framework (LargeDataset/methods/)

### What it's fundamentally doing
Using AUCell to score cells for membership in predefined gene programs, then visualizing scores as UMAPs, violin plots, evidence dot plots, and heatmaps. The 5 core metabolic pathways are reasonable defaults for many contexts. The PPARG signature is project-specific.

### Proposed primitive/module name
`primitives/aucell_scoring.md` (for `run_aucell()` + `add_auc_to_seurat()`)
`modules/metabolic_profile.md` (for the full analysis with visualization)

### What needs parameterization
| Currently hardcoded | Generalized parameter |
|---|---|
| PPARG gene set | Remove from shared version; keep in example |
| `tissue_type != "Myelolipoma"` filter | Remove; project-specific exclusion |
| `ec_subtype_colors` | `subtype_colors` argument |
| Script names referencing "HumanFat" | Generic names |
| OXPHOS/FA_Synthesis invalid range syntax (`"NDUFA1"-"NDUFA10"`) | **Must fix**: enumerate actual gene vectors |
| `output/metabolic/` subdirectory | `output_dir` argument |

### What stays in the example file
`examples/humanfat_metabolic_profile.md`:
- PPARG OE signature gene list
- Myelolipoma exclusion rationale
- Validated output structure and figure descriptions

### Risk assessment
1. **Invalid gene set syntax must be fixed before any use.** The `"NDUFA1"-"NDUFA10"` pattern is not valid R and will error. The general version must include explicit gene vectors. This is not a generalization risk — it's a pre-existing bug that affects the current file.
2. **AUCell memory requirements.** The `run_aucell()` implementation builds full cell rankings in memory. For objects >200k cells, this may exceed available memory even with the documented 2500-cell stratified subsampling. The general version should document the memory scaling more explicitly.
3. **Core pathway gene sets are biology-general but require curation.** The 5 metabolic pathways (Glycolysis, Beta_Oxidation, TCA_Cycle, OXPHOS, FA_Synthesis) are widely applicable but the gene membership reflects the file author's curation choices. Another researcher might define OXPHOS differently.

---

## G7 — `celltype_subclustering.md` (LargeDataset/methods/)

### What it's fundamentally doing
A two-phase exploration → endpoint pattern for subclustering: Phase 1 generates unsupervised outputs for label decisions; Phase 2 applies human-provided labels and produces final figures. This framework is entirely general — the phase structure is the contribution.

### What needs parameterization
| Currently hardcoded | Generalized parameter |
|---|---|
| `HARMONY_BY <- "source_file"` | `batch_correction_var` from brief |
| `"tissue_type"` metadata column | `group_col` from brief |
| `"source_file"` metadata column | `sample_col` from brief |
| `adipose_type_colors` | Remove; use `group_colors` from brief |
| Tissue color palette (6 adipose depots) | Remove; replaced by configurable `group_colors` |
| CELLTYPE_LABELS example | Replace with blank template |

### Recommendation
The two-phase framework is a first-class pattern. In v2, `pipelines/LargeDataset/methods/celltype_subclustering.md` can be kept as a module with all HumanFat column names replaced by references to brief parameters. No new file needed — just parameterization.

### Risk assessment
Phase 1 uses `pdf() + dev.off()` (not `ggsave()`). This bypasses the `useDingbats=FALSE` rule. The general version should standardize on `ggsave()` for all outputs. This is a minor inconsistency but worth fixing.

---

## G8 — Harmony Integration Sequence (cross-pipeline)

### What it's fundamentally doing
Both `LargeDataset/pipeline.md` (Stage 3) and `IntegratePublicData/pipeline.md` (Stage 6) describe the identical Seurat+Harmony processing sequence. This is not a method file but a recipe embedded in two pipeline definitions.

### Proposed primitive
`primitives/harmony_integration.md` with a single canonical recipe parameterized on:
- `n_variable_features` (default: 4000)
- `n_pcs` (default: 30)
- `batch_correction_var`
- `assay` (default: "RNA")

Both pipeline.md files would reference this primitive. The specific defaults (n_variable_features: 4000 vs 2850 for subsets) would be documented in each pipeline's stage description.

### Risk assessment
Minimal. The Harmony sequence is stable and well-validated. The only drift risk is if one pipeline updates its defaults without updating the primitive — but this is solved by explicit default documentation in the pipeline stage.

---

## Generalization Priority Matrix

| Candidate | Effort to generalize | Value | Priority |
|---|---|---|---|
| G5 `functional_gene_sets` | Low — parameterize one block | Prevents EC genes in non-EC projects | **High** |
| G8 Harmony sequence | Low — single canonical recipe | Eliminates maintenance duplication | **High** |
| G7 Subclustering framework | Medium — parameterize 4 column names | Applies to every subclustering job | **High** |
| G6 AUCell primitives (`run_aucell`, `add_auc_to_seurat`) | Low — functions are already clean | Reusable for any gene program scoring | **High** |
| G1 `anatomical_de_analysis.md` | High — requires resolving cross-pipeline injection | Only relevant for multi-group site comparisons | **Medium** |
| G4 Bulk concordance | Medium — parameterize project-specific fields | Applicable to any bulk × scRNA-seq project | **Medium** |
| G2 Atlas co-UMAP | Medium — parameterize dataset labels + colors | Applicable to any atlas comparison | **Medium** |
| G3 Cross-atlas dotplot | High — depth correction logic is subtle | Most project-specific of all three atlas methods | **Lower** |
