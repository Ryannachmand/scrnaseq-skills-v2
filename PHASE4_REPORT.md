---
# Phase 4 Report — Examples Directory + Phase 3 Maintenance Items
# Written: 2026-05-19
# Author: Phase 4 migration agent
---

# Phase 4 Report

Covers the 17 examples/ files (15 full builds + 2 skeletons) and the four Phase 3
maintenance items. Includes per-example source documentation, maintenance item diffs,
self-check results, and open items for Phase 5 and Phase 6.

---

## 1. Files Created

### examples/ (17 files)

| File | Lines | Status | Module instantiated |
|---|---|---|---|
| `examples/humanfat_metabolic_profile.md` | 224 | full-build | metabolic_profile.md |
| `examples/humanfat_trajectory.md` | 66 | full-build | trajectory_monocle3.md |
| `examples/humanfat_pyscenic.md` | 59 | full-build | pyscenic_regulons.md |
| `examples/pparg_bulk_concordance.md` | 145 | full-build | bulk_concordance.md (Mode 1) |
| `examples/humanfat_de_comprehensive.md` | 70 | full-build | de_comprehensive_csv.md |
| `examples/humanfat_cellchat.md` | 140 | full-build | cellchat.md |
| `examples/humanfat_mac_subclustering.md` | 87 | full-build | celltype_subclustering.md |
| `examples/kidneynew_cellchat.md` | 123 | full-build | cellchat.md |
| `examples/bone_marrow_stroma_cohort_plots.md` | 108 | full-build | cohort_plots.md |
| `examples/bone_marrow_stroma_label_harmonization.md` | 91 | full-build | label_harmonization.md |
| `examples/load_formats_bonemarrowstroma.md` | 99 | full-build | load_formats.md |
| `examples/bonemarrow_3site_anatomical.md` | 145 | full-build | multi_group_de_analysis.md |
| `examples/tabula_sapiens_co_umap.md` | 116 | full-build | atlas_co_umap.md |
| `examples/tabula_sapiens_dotplot.md` | 118 | full-build | cross_dataset_dotplot.md |
| `examples/ec_functional_gene_sets.md` | 250 | full-build | multi_group_de_analysis.md / de_comprehensive_csv.md |
| `examples/coculture_feature_umap.md` | 69 | skeleton | feature_umap_plot.md |
| `examples/nkxspleen_bulk_lfc.md` | 111 | skeleton | bulk_concordance.md (Mode 2) |
| **Total** | **2,021** | | |

### Modified files

| File | Change type | Lines delta |
|---|---|---|
| `SKILL.md` | Updated library tree + routing table | +106 |
| `modules/feature_umap_plot.md` | P4-M1 brief key path fix | +6 |
| `modules/multi_group_de_analysis.md` | P4-M3 schema documentation fix | +2/-2 |
| `context/known_atlases.md` | P4-M4 validated_date source comment | +5 |

---

## 2. Per-Example Summary

### 2.1 HumanFat examples

#### humanfat_metabolic_profile.md
- **Module:** metabolic_profile.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml HumanFat entry + Phase 2B staging notes
- **Complete:** Five gene set vectors (Glycolysis, Beta_Oxidation, TCA_Cycle, OXPHOS,
  FA_Synthesis) with explicit gene lists replacing v1 invalid range syntax; Myelolipoma
  exclusion documented; ec_subtype_colors from registry; validated 2026-03-05
- **TODOs:** None

#### humanfat_trajectory.md
- **Module:** trajectory_monocle3.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2B staging notes
- **Complete:** LABEL_COL=mylabel, EXCLUDE_PATTERN=^AEC, START_CLUSTER=CapEC2;
  53,487 EC cells / 13 subtypes; validated 2026-03-05
- **TODOs:** None

#### humanfat_pyscenic.md
- **Module:** pyscenic_regulons.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2B staging notes
- **Complete:** 6 EC subtypes, 320 regulons, machine-specific paths documented;
  validated 2026-03-11
- **TODOs:** None

#### pparg_bulk_concordance.md
- **Module:** bulk_concordance.md (Mode 1 — signature_score)
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2B staging notes + Phase 2B U3 note
- **Complete:** Mode 1 context (PPARG1_ROSI vs Control_ROSI HUVECs); adipose tissue palette;
  EC subtype palette; RibHighEC exclusion; EC-biology TF list (~180 TFs, Lambert 2018 basis)
  as explicit R vector; validated 2026-03-12
- **TODOs:** None — TF list is the most complete item in any example file (addresses Phase 2B U3)

#### humanfat_de_comprehensive.md
- **Module:** de_comprehensive_csv.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2B staging notes
- **Complete:** GROUP_COL=tissue_type, 6 tissue groups, LABEL_COL=mylabel, GROUP_ORDER
  for display; ~37,692 EC cells, ~253 significant genes after is_confound();
  validated 2026-04-23
- **TODOs:** None

#### humanfat_cellchat.md
- **Module:** cellchat.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2D staging notes
- **Complete:** LABEL_COL=mylabel, GROUP_COL=tissue_type, full cell_colors vector for
  all 21 cell types; Script 6 LiposuctionFat vs BreastFat comparison; CAT_COLORS for
  bar plots (7 categories); STACKED_CATS ECM exclusion documented; validated 2026-03-30
- **TODOs:** Full PATHWAY_CATEGORIES list not reconstructed — requires extraction from
  HumanFat project CLAUDE.md. Documented as a TODO in the Script 1→2 section.

#### humanfat_mac_subclustering.md
- **Module:** celltype_subclustering.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2B staging notes
- **Complete:** BATCH_CORRECTION_VAR=source_file, GROUP_COL=tissue_type, target_celltypes
  (Macrophage + Monocyte), subset pipeline params; validated 2026-03-24
- **TODOs:** Exact cluster→label mapping from 2026-03-24 run not recovered from library
  records. Documented with TODO placeholder in the example.

---

### 2.2 KidneyNew examples

#### kidneynew_cellchat.md
- **Module:** cellchat.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml KidneyNew entry + Phase 2D staging notes
- **Complete:** LABEL_COL=ec_subtype_final, 6 EC subtypes (EC-GC, EC-AEA, EC-AVR, EC-DVR,
  EC-LYM, EC-PTC), 3 pathway categories (Chemokine_Adhesion, Angiocrine, ECM) with
  APP exclusion note; assay switching pattern (SCT→RNA→CellChat→SCT); plan("sequential")
  critical note; EC-PTC post-hoc annotation context; validated 2026-04-10
- **TODOs:** Non-EC cell type colors (Stroma, Epithelial, MoMF) are placeholder hex values;
  confirm from KidneyNew project CLAUDE.md. Full PATHWAY_CATEGORIES pathway lists are
  approximate; verify from project records.

---

### 2.3 BoneMarrowStroma examples

#### bone_marrow_stroma_cohort_plots.md
- **Module:** cohort_plots.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2A staging notes
- **Complete:** ann_track_cols=c(Dataset, Organ, Condition, Data.Type); validated figure
  dimensions 13.25 × 7 in; ann_colors structure documented with TODO for actual hex values;
  pdf()/dev.off() exceptions documented (CONVENTIONS.md §4 exceptions #2 and #3)
- **TODOs:** Actual color hex values for dataset/organ/condition/data.type palette not
  recovered from library records. Documented as TODO.

#### bone_marrow_stroma_label_harmonization.md
- **Module:** label_harmonization.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2A staging notes
- **Complete:** broad_label="Mesenchymal"; reference_datasets=[Vertebrae, FemoralHead];
  Li transfer result distribution (Adipo-MSC: 1823, Osteo-MSC: 912, Unknown: 75);
  sort-strategy labels (CD14/PE Isolated Stroma); Passaged Stroma; "Cultured Adipo-like MSC"
  derived display label
- **TODOs:** Actual label_col column names per dataset need confirmation from project records

#### load_formats_bonemarrowstroma.md
- **Module:** load_formats.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2A staging notes
- **Complete:** GSM file pattern conventions for all three sources (Leimkuler/Wang/Li);
  FemoralHead UMAP_dim30 non-standard reduction name; Wang nested directory structure;
  Li spliced/unspliced format; RenameFeatures() over rownames() assignment
- **TODOs:** Actual GSM accessions for Li/FemoralHead dataset not recovered from records

#### bonemarrow_3site_anatomical.md
- **Module:** multi_group_de_analysis.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2C staging notes
- **Complete:** GROUPS=c(Vertebrae, Iliac Crest, Femoral Head); 3 pairwise comparisons;
  functional_gene_sets authored for bone marrow stromal biology (4 categories:
  Osteogenic, Adipogenic, Fibroblastic, Hematopoietic_Niche); validated 2026-03-16
- **TODOs:** STROMA_ORDER (label_order) requires human curation — explicitly documented.
  Group direction strip colors are placeholder. These were never defined in v1.

---

### 2.4 TabulaSapiens examples

#### tabula_sapiens_co_umap.md
- **Module:** atlas_co_umap.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2C staging notes
- **Complete:** INHOUSE_LABEL=WCM, ATLAS_LABEL=TS, DOWNSAMPLE_N=1500 with full WCM/TS
  ratio analysis (WCM ~8% / TS ~92% of UMAP); COORDS_CSV path; validated script filenames
  (06f, 06g); validated figure dimensions (21 × 5.5 in); Known-Atlas Convention documented;
  custom font and cairo_pdf issues documented; validated 2026-04-15
- **TODOs:** subtype_col column name ("ec_subtype") needs confirmation

#### tabula_sapiens_dotplot.md
- **Module:** cross_dataset_dotplot.md
- **Status:** full-build
- **Context sources:** validated_examples.yaml + Phase 2C staging notes
- **Complete:** 3-section marker_genes structure; FORCE_GENES_SEC1 (forced TF list with
  rationale); GROUP_COLORS (warm yellow / light red); column label convention;
  sources_order=c(WCM, TS); custom font issues; validated 2026-03-24
- **TODOs:** Specific gene lists in all three sections are placeholder; require project
  records for exact genes used in validated run. FORCE_GENES_SEC1 lists KLF2/ERG as
  examples; full list from project records needed.

---

### 2.5 Canonical reference

#### ec_functional_gene_sets.md
- **Module:** multi_group_de_analysis.md / de_comprehensive_csv.md
- **Status:** full-build
- **Context sources:** Phase 1 F4 note (gene sets removed from differential_expression.md
  primitive) + published EC biology literature
- **Complete:** Four categories (Adhesion_Immune_Trafficking, Signaling_Angiocrine,
  Extracellular_Matrix, Metabolic_Specialized); full gene vectors with ~250 total genes;
  YAML instantiation instructions; usage guidance for both modules
- **TODOs:** Gene lists are authored based on EC biology literature; specific genes
  from validated runs should be cross-checked against actual script outputs when available

---

### 2.6 Skeleton examples

#### coculture_feature_umap.md
- **Module:** feature_umap_plot.md
- **Status:** skeleton
- **Context sources:** validated_examples.yaml CoCulture entry (sparse) + Phase 2A notes
- **Complete:** Module identified, project context documented, TODOs itemized
- **What's marked TODO:** Gene list, dataset path, validation date, output dimensions,
  label column, add to validated_examples.yaml

#### nkxspleen_bulk_lfc.md
- **Module:** bulk_concordance.md (Mode 2 — parallel_lfc)
- **Status:** skeleton
- **Context sources:** Phase 2B staging notes (80,140 EC cells, 10 subtypes, 12 bulk
  samples, 6 Control / 6 NKX2-3 OE; validated 2026-03-13)
- **Complete:** Module identified, Mode 2 context documented, biology_gene_sets structure
  from Phase 2B (6 categories: Chemokines, Adhesion, Angiocrine, ECM, Sinusoidal_Markers,
  TF_Program) with placeholder gene lists, sc_group assignment pattern documented
- **What's marked TODO:** Bulk CSV path, condition_samples/control_samples, subtype_col,
  exact gene lists, NKXSpleen entry in validated_examples.yaml

---

## 3. Maintenance Items Resolved

### P4-M1: feature_umap_plot.md brief key path

**Problem:** Module used `brief$genes` and `brief$output_dir` at the top level. All other
modules use `brief$downstream_analyses.<module_key>.*`.

**Change in `modules/feature_umap_plot.md`:**

Before:
```yaml
brief_keys:
  optional:
    - genes
```
After:
```yaml
brief_keys:
  optional:
    - downstream_analyses.feature_umap_plot.genes
```

Also updated the Configuration Block:
```r
# Before:
GENES <- brief$genes %||% c("GENE1", "GENE2", "GENE3")
OUTPUT_DIR <- file.path(brief$output_dir %||% "output", "feature_umap")

# After:
GENES <- brief$downstream_analyses$feature_umap_plot$genes %||%
         brief$genes %||%    # legacy fallback
         c("GENE1", "GENE2", "GENE3")
OUTPUT_DIR <- file.path(
  brief$downstream_analyses$feature_umap_plot$output_dir %||%
  brief$output_dir %||%    # legacy fallback
  "output",
  "feature_umap"
)
```

Legacy `brief$genes` fallback is preserved for any briefs written before this fix.

**Design note:** The Phase 3 report (P3-1) originally said this fix should be in the
deployment agent's Stage 8 dispatch layer. On reflection, fixing it in the module itself
is cleaner and self-documenting. The deployment agent no longer needs a special case.

---

### P4-M2: SKILL.md routing table

**Changes made:**

1. **Library tree:** Updated to show complete structure — `modules/` (14 files with names),
   `pipelines/` (complete with pipeline.md + brief_template.txt for each), `examples/`
   (all 17 files listed), `context/known_atlases.md` added.

2. **Routing table:** Replaced 8-row "Phase 2 (not yet written)" table with 15-row complete
   table including: all 14 modules, Stage 8 dispatch note (execution order: multi_group_de
   first, bulk_concordance second), atlas-aware modules noted.

3. **Never-auto-inject section:** Updated to document three conditions under which examples/
   files ARE consumed (registry hit, explicit user reference, Phase 6 agent reference).

4. **How to Generate Project Files:** Updated step-by-step instructions to include
   Known-Atlas Convention lookup (step 9), pipeline.md reading (step 7), and correct
   brief_template paths.

5. **Available Pipelines table:** Updated to show brief_template.txt paths and mark Phase 3 complete.

---

### P4-M3: multi_group_de_analysis.md schema fix

**Problem:** `functional_gene_sets` was labeled REQUIRED in both the frontmatter
`brief_keys.required` block and the Brief Schema section comment. But the module code
gates all functional plots conditionally on `!is.null(functional_gene_sets)`.

**Change in `modules/multi_group_de_analysis.md`:**

Frontmatter before:
```yaml
brief_keys:
  required:
    - ...
    - downstream_analyses.multi_group_de.functional_gene_sets
```

Frontmatter after:
```yaml
brief_keys:
  required:
    - ...
    # functional_gene_sets moved to optional
  optional:
    - downstream_analyses.multi_group_de.functional_gene_sets  # null = skip functional plots
    - ...
```

Brief Schema section comment change:
```yaml
# Before:
functional_gene_sets: project_specific  # REQUIRED: named list of gene vectors, biology-specific

# After:
functional_gene_sets: null              # OPTIONAL: named list of gene vectors, biology-specific
                                        # if null: make_functional_dotplot and make_functional_heatmap
                                        # steps are skipped entirely; base DE + dotplot still runs
```

---

### P4-M4: Tabula Sapiens validated_date

**Action taken:** No date change. Added a source comment.

**Evidence found:**
- Phase 2C staging notes: "Validated run: TabulaSapiensComparison, 2026-04-15"
- Phase 3 report (section 4.1): "validated_date: 2026-04-15"
- known_atlases.md: `validated_date: "2026-04-15"` (set in Phase 3)
- validated_examples.yaml: `last_validated: null  # TODO`

**Assessment:** The date is consistent across Phase 2C and Phase 3 sources. The discrepancy
is that `validated_examples.yaml` still shows `null` for TabulaSapiensComparison, while
`known_atlases.md` has 2026-04-15. These could refer to different events (atlas registry
validated_date vs project last_validated). No contradicting date was found.

**Comment added** to known_atlases.md line documenting the source and noting the
validated_examples.yaml discrepancy for Phase 5 to resolve.

---

## 4. Self-Check Results

### 4.1 Cross-project name leakage

Initial grep found three leakage issues, all fixed:

| Issue | File | Violation | Fix |
|---|---|---|---|
| KidneyNew in humanfat_cellchat.md | Line: "KidneyNew CellChat validated separately" | Cross-project reference prose | Removed the sentence |
| KidneyNew in humanfat_pyscenic.md | Line: "pySCENIC also validated separately on KidneyNew" | Cross-project reference prose | Changed to "another project" |
| HumanFat + KidneyNew + BoneMarrowStroma in ec_functional_gene_sets.md | Validation header | Multi-project canonical reference | Changed to "Multiple EC-focused projects" |

After fixes, final leakage check:
- HumanFat: appears only in humanfat_* and pparg_bulk_concordance.md ✓
- KidneyNew: appears only in kidneynew_cellchat.md ✓
- BoneMarrowStroma/BoneMarrow: appears only in bone_marrow_stroma_* and bonemarrow_* and load_formats_bonemarrowstroma.md ✓
- TabulaSapiens: appears only in tabula_sapiens_* ✓
- CoCulture: appears only in coculture_feature_umap.md ✓
- NKXSpleen: appears only in nkxspleen_bulk_lfc.md ✓

**Result: CLEAN.**

### 4.2 @reference resolution

All @modules/ references checked against modules/ directory — all 14 target files exist ✓

@primitives/ references: seurat_v5_rules.md (Rule 1, Rule 5, Rule 6), differential_expression.md —
all present in primitives/ ✓

@examples/pparg_bulk_concordance.md cross-reference in humanfat_metabolic_profile.md — file exists ✓

---

## 5. Cross-File Findings

### 5.1 feature_umap_plot.md P4-M1 design decision

The Phase 3 report (P3-1) proposed fixing this in the deployment agent's Stage 8 dispatch
layer rather than in the module itself. Phase 4 fixed it in the module. Rationale:
- Fixing in the module makes the correct reading path explicit and self-documenting
- The legacy `brief$genes` fallback means existing briefs continue to work
- The deployment agent no longer needs a special case for this one module

This is a design change from what P3-1 proposed. Phase 6 should not implement the
P3-1 dispatch workaround since the module now handles it natively.

### 5.2 STROMA_ORDER remains undefined

STROMA_ORDER (the ordered cell type vector for BoneMarrowStroma x-axis display) was
never defined in any v1 or v2 library file. Phase 4 documents it in
`bonemarrow_3site_anatomical.md` as a `label_order: null` with a prominent TODO comment
and curation instructions. This is the correct approach — it cannot be fabricated
without BoneMarrowStroma project records.

### 5.3 humanfat_mac_subclustering.md cluster→label mapping is incomplete

The exact cluster→label mapping from the 2026-03-24 HumanFat macrophage subclustering
run was not recoverable from library records. The example documents the expected mapping
structure with a TODO placeholder. The validated date and parameters are correct.

### 5.4 kidneynew_cellchat.md non-EC colors are placeholder

The hex colors for Stroma, Epithelial, and MoMF cell types in KidneyNew are `#project_specific`
placeholders. These need to be filled from KidneyNew project CLAUDE.md before running
CellChat scripts. Documented as TODO in the example.

### 5.5 bonemarrow_3site_anatomical.md functional_gene_sets are newly authored

The four functional gene sets (Osteogenic, Adipogenic, Fibroblastic, Hematopoietic_Niche)
are authored based on bone marrow MSC biology literature, not extracted from validated
run outputs. The v1 had EC-biology gene sets; Phase 1 F4 removed them; no bone marrow
stromal gene sets existed anywhere in the library. These are the first stromal gene sets
in the library.

### 5.6 tabula_sapiens_dotplot.md marker_genes are placeholder

The three-section `marker_genes` list in `tabula_sapiens_dotplot.md` shows structure and
a few example genes but the full gene lists from the 2026-03-24 run are not in library
records. Documented as TODO.

---

## 6. Open Items for Phase 5 (Validation) and Phase 6 (Deployment Agent)

### For Phase 5 (validation runs)

1. **Resolve STROMA_ORDER** — Run `sort(unique(so$unified_label))` on the BoneMarrowStroma
   object and curate the label_order vector. Update `bonemarrow_3site_anatomical.md`.

2. **Fill humanfat_mac_subclustering cluster→label mapping** — Retrieve from HumanFat
   2026-03-24 run CLAUDE.md or project notes.

3. **Fill kidneynew_cellchat.md non-EC cell type colors** — Retrieve from KidneyNew project
   CLAUDE.md.

4. **Fill tabula_sapiens_dotplot.md marker_genes** — Retrieve from TabulaSapiensComparison
   project records; should have full gene lists from 06h_ts_dotplot.R or equivalent script.

5. **Fill humanfat_cellchat.md PATHWAY_CATEGORIES** — Retrieve from HumanFat project CLAUDE.md;
   this is the primary outstanding item for the HumanFat CellChat example.

6. **Confirm CoCulture gene list** — The skeleton at `coculture_feature_umap.md` needs the
   gene list from the validated CoCulture run; update to full-build status after.

7. **Validate bonemarrow_3site_anatomical.md functional_gene_sets** — The four stromal gene
   sets are newly authored and not extracted from validated run outputs; cross-check against
   actual BoneMarrowStroma DE outputs when available.

8. **Confirm TabulaSapiensComparison validated_date** — Fill `last_validated:` in
   `validated_examples.yaml` (currently null); determine whether it matches the
   known_atlases.md date of 2026-04-15 or whether it should be different.

### For Phase 6 (deployment agent)

9. **P4-M1 deployment agent note** — Phase 4 fixed the brief key path directly in the module
   (with legacy fallback). Phase 6 does NOT need to implement the P3-1 dispatch workaround.
   The module now reads from `downstream_analyses.feature_umap_plot.genes` natively.

10. **Known-Atlas Convention expansion** — When `brief.atlas: "Tabula Sapiens"` is set,
    inject `source_col`, `atlas_label`, `atlas_group_col`, `recommended_downsample_n` from
    `context/known_atlases.md`. Explicit brief values take precedence. Log which values came
    from registry vs explicit brief.

11. **Registry hit → examples/ injection** — When starting work on a project that matches
    a registry entry in `validated_examples.yaml`, look up whether an examples/ file exists
    for the requested module. If it does, inject it alongside the module for project-specific
    context. Do not auto-inject examples/ for projects that don't match.

12. **NKXSpleen registry entry** — Add NKXSpleen to `validated_examples.yaml` with full
    context_defaults. The skeleton at `nkxspleen_bulk_lfc.md` documents what's needed.

13. **P3-3 (atlas typical_cell_count)** — The `typical_cell_count: 483152` in known_atlases.md
    is an approximation. Confirm against the actual TS download used in TabulaSapiensComparison.

---

## 7. Anything Unexpected

### U1 — PPARG TF list is the most substantive new content in Phase 4

The EC-biology TF list in `pparg_bulk_concordance.md` (~180 genes organized by TF family)
was flagged in Phase 2B (U3) as "the single largest body of biology-specific content in
Group B." Phase 4 authored this as an explicit R vector organized by TF family
(Core EC identity, KLF/SP, FOXO, NFI/NFAT, SOX, GATA, AP-1, bHLH, NR family, YAP/TAZ,
mechanosensing, inflammatory, angiogenic, arterial/venous, metabolic, circadian).
This required composing the list from scratch using published EC TF literature rather than
extracting from a validated run — the v1 list was removed from the module in Phase 1.

### U2 — BoneMarrowStroma functional gene sets are entirely new

No bone marrow stromal gene sets existed anywhere in the v1 or v2 library. The four
categories in `bonemarrow_3site_anatomical.md` (Osteogenic, Adipogenic, Fibroblastic,
Hematopoietic_Niche) were composed from MSC biology literature. These need
validation against actual BoneMarrowStroma DE outputs before being treated as canonical.

### U3 — P4-M1 was resolved more cleanly than P3-1 proposed

P3-1 proposed a deployment agent dispatch workaround for the feature_umap_plot.md
top-level brief key issue. Phase 4's fix (updating the module itself with legacy fallback)
eliminates the need for a deployment agent special case entirely. The fix is cleaner
because: (1) the module is self-documenting about its expected brief key path,
(2) the legacy fallback protects against breaking changes, (3) no deployment agent
logic change is required.

### U4 — ec_functional_gene_sets.md is the largest example file (250 lines)

It was removed from the primitive in Phase 1 (F4) and needed to be re-authored in Phase 4
as a canonical reference. The file is 250 lines because it contains four complete gene
vectors totaling ~250 genes plus usage instructions, attribution, and notes on gene
symbol conventions. This is expected for a canonical reference — it should be dense.

### U5 — tabula_sapiens_dotplot.md marker_genes could not be reconstructed

The three-section marker_genes list from the 2026-03-24 TabulaSapiensComparison run was
not recoverable from library records. The example documents the structure and a few
biologically plausible placeholder genes. Phase 5 must retrieve the actual list from
project records and update the example to full-build status with verified content.
