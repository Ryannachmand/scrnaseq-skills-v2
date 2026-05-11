# Coverage Gaps — ~/claude-skills/
*Generated: 2026-05-11 by inventory agent*

What the library SHOULD have but doesn't: missing files referenced by pipelines, unimplemented functions, uncovered analysis categories, and convention inconsistencies.

---

## 1. References to Methods Files That Don't Exist

### 1.1 No `brief_template.txt` for IntegratePublicData

- **Where referenced:** `SKILL.md` Step 2 instructs: "Read the pipeline's `pipeline.md` and `brief_template.txt`"
- **What exists:** `pipelines/LargeDataset/brief_template.txt` (166 lines) ✅
- **What's missing:** `pipelines/IntegratePublicData/brief_template.txt` — does not exist
- **Impact:** When a user asks to generate project files for an IntegratePublicData job, the agent reads only `pipeline.md` for the brief structure. The inline "Brief Parameters Reference" table in `pipeline.md` lines ~380–416 is a partial substitute but is less explicit and less annotated than the LargeDataset template.
- **v2 action:** Create `IntegratePublicData/brief_template.txt` with REQUIRED/CONFIRM/OPTIONAL annotations.

---

### 1.2 No `methods/aesthetics.md` for IntegratePublicData

- **Where referenced:** `shared/aesthetics.md` line 4 says "For plot-specific implementations, see each pipeline's `methods/aesthetics.md`"
- **What exists:** `pipelines/LargeDataset/methods/aesthetics.md` (207 lines) ✅
- **What's missing:** `pipelines/IntegratePublicData/methods/aesthetics.md` — does not exist
- **Impact:** IntegratePublicData jobs get the shared principles file but no concrete R code for the UMAP, proportion, and cohort plot aesthetics specific to multi-source integration outputs. Plot recipes are scattered inside `cohort_plots.md`, `co_umap_embedding.md`, and `cross_atlas_dotplot.md` without a central coordinator.
- **v2 action:** Create `IntegratePublicData/methods/aesthetics.md` collecting the cohort-specific plot recipes.

---

### 1.3 `anatomical_de_analysis.md` Not in IntegratePublicData Read List

- **Where it's injected:** `SKILL.md` "read ALL files in methods/" causes it to be read for all IntegratePublicData jobs
- **Where it's NOT listed:** `IntegratePublicData/pipeline.md` "Always Read Before Starting" (lines 31–37) does not include it
- **Impact:** The file is read (via SKILL.md blanket instruction) but there is no pipeline stage that tells the agent when and how to apply it. A new job gets `anatomical_de_analysis.md` injected into its context even if there is no anatomical comparison. More critically, the file's phantom function dependencies (`run_findmarkers`, `make_volcano`, etc. from `differential_expression.md`) are never resolved because `differential_expression.md` is not injected.
- **v2 action:** Either (a) add a conditional `methods/anatomical_de_analysis.md` read to the IntegratePublicData pipeline for jobs with multi-site comparisons, or (b) resolve the dependency by making the DE functions available cross-pipeline.

---

### 1.4 `skills_gap_analysis.md` Referenced in Task Description But Doesn't Exist

- **Location searched:** `~/claude-skills/jobs/skills_gap_analysis.md`
- **Status:** File does not exist. The `jobs/` directory itself does not exist.
- **Note:** The inventory task description says to "Read that first to orient." The file was presumably created as a scratch document during an earlier analysis session but was never committed to the library. This gap document cannot reference it.

---

## 2. Unimplemented Functions Referenced in Methods Files

### 2.1 `is_confound()` — Phantom in `de_comprehensive_csv.md`

- **Called at:** `de_comprehensive_csv.md` line 116 in the `sig_genes` filter
- **Defined:** Nowhere in the library
- **What it needs to do:** Based on the "Exclusion filters" section at ~line 215, `is_confound()` should exclude:
  1. Ambient genes (covered by `is_ambient()`)
  2. Sex-linked genes (chromosome X/Y: XIST, UTY, DDX3Y, EIF1AY, KDM5D, USP9Y, etc.)
  3. HLA class II genes (HLA-DRA, HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DPA1, HLA-DPB1, etc.)
  4. Histone genes (HIST1H*, HIST2H*, H2A*, H2B*, H3F3*, etc.)
  5. Unannotated ENSG IDs (genes matching `^ENSG`)
- **v2 action:** Write `is_confound()` as a proper named function that calls `is_ambient()` internally plus adds the four additional exclusion categories. This is a clear primitive candidate.

---

### 2.2 Undefined Color Variables Used Across Multiple Files

The following named color vectors are referenced in methods files without being defined in those files. They are expected to exist in the project CLAUDE.md context block (populated from `lab_context.md`). If the CLAUDE.md was generated without the color blocks, all scripts referencing them will fail with `object 'ec_colors' not found`.

| Variable | Referenced in | Expected source |
|---|---|---|
| `ec_colors` | `bulk_concordance.md`, `cross_atlas_dotplot.md`, `metabolic_profile.md` | `lab_context.md` `ec_subtype_colors_adipose` |
| `tissue_colors` | `differential_expression.md` (heatmap annotation), `bulk_concordance.md` | `lab_context.md` `organ_tissue_colors` |
| `ec_subtype_colors` | `metabolic_profile.md` | `lab_context.md` `ec_subtype_colors_adipose` |
| `adipose_type_colors` | `celltype_subclustering.md` | `lab_context.md` `adipose_type_colors` |
| `label_colors` | `anatomical_de_analysis.md` `make_topgene_dotplot()` | Not defined in `lab_context.md` — would need a BoneMarrowStroma color vector |
| `STROMA_ORDER` | `anatomical_de_analysis.md` | Not defined anywhere in the library |

**`STROMA_ORDER` is a critical gap:** it is referenced in `make_topgene_dotplot()` at line 107 and in the Critical Constraints table but is never defined in any library file. It presumably existed in the BoneMarrowStroma CLAUDE.md but was never brought back into the library.

---

## 3. Brief Parameters Without Corresponding Methods Recipes

### 3.1 `cohort_plots: true` — Gated but Underdocumented

- `IntegratePublicData/pipeline.md` Stage 7 says: "Run cohort_plots.md if `brief.cohort_plots == true`"
- `cohort_plots.md` exists but has no brief configuration block. There's no `brief_template.txt` field for cohort plot options (which of the 3 plot types to run, figure size overrides, etc.)
- **Gap:** A user setting `cohort_plots: true` gets all three plot types unconditionally. No way to request just the proportion heatmap or just the chord diagram.

### 3.2 Trajectory Analysis Configuration

- `LargeDataset/pipeline.md` Stage 8 lists `trajectory` as an optional analysis key
- `trajectory_monocle3.md` has a brief configuration block: `exclude_pattern`, `start_cluster`, `end_cluster`
- **Gap:** `brief_template.txt` does not include the `trajectory:` block under `optional_analyses`. A user who wants trajectory analysis has to find the configuration schema in `trajectory_monocle3.md` directly.
- **v2 action:** Add `trajectory:` block to `brief_template.txt` optional_analyses section.

### 3.3 pySCENIC Analysis Has No Brief Configuration

- `LargeDataset/pipeline.md` Stage 8 lists `pyscenic` as an optional analysis key
- `pyscenic_regulon_analysis.md` contains no brief configuration block
- **Gap:** Database paths, TF list path, Python binary, number of workers, output directory — all are hardcoded in the file. There's no parameterization pattern and no `pyscenic:` key in `brief_template.txt`.
- **v2 action:** Add a `pyscenic:` brief block; parameterize the three machine-specific paths.

---

## 4. Asymmetric Pipeline Coverage

The two pipelines have different capabilities documented. Below is the cross-pipeline capability matrix:

| Capability | LargeDataset | IntegratePublicData | Gap |
|---|---|---|---|
| QC thresholds and filtering | ✅ Stage 2 | ✅ Stage 2 | None |
| Harmony integration | ✅ Stage 3 | ✅ Stage 6 | Duplicated recipe; should share primitive |
| Cell type annotation | ✅ Stages 4–5 | ✅ Stage 3 | Methods differ: marker-based vs label transfer |
| Visualization suite | ✅ methods/aesthetics.md | ❌ No equivalent | **IntegratePublicData gap** |
| Subclustering | ✅ celltype_subclustering.md | ❌ Not documented | **IntegratePublicData gap** — may apply to subset analysis after integration |
| DE between groups | ✅ differential_expression.md | ⚠️ anatomical_de_analysis.md (ORPHAN, phantom deps) | **Critical gap** — functional but broken |
| Multi-group KW statistics | ✅ de_comprehensive_csv.md | ❌ Not documented | **IntegratePublicData gap** |
| Bulk perturbation concordance | ✅ bulk_concordance.md | ❌ Not documented | Arguably LargeDataset-specific; may be intentional |
| Parallel LFC heatmap | ✅ bulk_lfc_concordance_heatmap.md | ❌ Not documented | Same |
| CellChat interactome | ✅ interactome_cellchat.md | ❌ Not documented | Probably LargeDataset-specific |
| Metabolic profile | ✅ metabolic_profile.md | ❌ Not documented | Could apply to integration projects |
| TF regulon (pySCENIC) | ✅ pyscenic_regulon_analysis.md | ❌ Not documented | Could apply to integration projects |
| Trajectory analysis | ✅ trajectory_monocle3.md | ❌ Not documented | Could apply post-integration subset |
| Atlas co-UMAP | ❌ Not documented | ✅ co_umap_embedding.md | LargeDataset datasets rarely compared to atlas |
| Cross-atlas dotplot | ❌ Not documented | ✅ cross_atlas_dotplot.md | Same |
| Cohort plots | ❌ Not documented | ✅ cohort_plots.md | LargeDataset equivalent is proportion plot in aesthetics.md |
| Label harmonization | ❌ Not documented | ✅ label_harmonization.md | LargeDataset uses simpler per-project labeling |

**Most critical asymmetric gap:** DE between groups exists only for LargeDataset in a functional form. The IntegratePublicData version (`anatomical_de_analysis.md`) is an orphan with phantom function dependencies. Any integration project requiring post-integration group comparisons has no working recipe.

---

## 5. Missing Conventions

### 5.1 No Standardized Output Directory Convention

Each methods file specifies output paths differently:
- `differential_expression.md`: `output2/`
- `bulk_concordance.md`: `output3/`
- `metabolic_profile.md`: `output/metabolic/`
- `pyscenic_regulon_analysis.md`: not standardized
- `trajectory_monocle3.md`: `out_dir` variable (no default)

**Missing:** A convention for output directory naming that the brief_template.txt enforces. In v2, consider: `output/{method_name}/` as the default pattern.

---

### 5.2 No Standardized File Naming Convention for Script Files

Script names within methods files are inconsistent:
- Some use numeric prefixes: `06f_ts_umap_presample.R`
- Some use descriptive names: `01_cellchat_inference.R`
- `de_comprehensive_csv.md` does not specify script file names at all

**Missing:** A script naming convention. Suggested: `{nn}_{analysis_name}.R` where `{nn}` is a zero-padded stage number.

---

### 5.3 No Brief Field for Cell Type Column Name

Every methods file has a different assumption about which metadata column holds cell type labels:
- `differential_expression.md`: `ec_sub$mylabel` (hardcoded)
- `celltype_subclustering.md`: configurable via `LABEL_COL` variable
- `anatomical_de_analysis.md`: `unified_label` (IntegratePublicData standard)
- `interactome_cellchat.md`: `ec_subtype` (project-specific)
- `metabolic_profile.md`: not parameterized

**Missing:** A `cell_type_col` brief field that all methods files read from. Currently, the assumption is embedded per-file.

---

### 5.4 No Doublet Detection Recipe

- Neither pipeline includes a doublet detection step
- QC Stage 2 in LargeDataset mentions `percent_mt_max` and feature counts but no doublet filtering
- Common approaches (DoubletFinder, scDblFinder) are not documented
- This is a notable gap for any multi-sample dataset

---

### 5.5 No Multi-Batch Harmony Correction Recipe

- Both pipelines use `RunHarmony(group.by.vars = batch_correction_var)` where `batch_correction_var` is a single column
- No documentation for multi-batch correction (e.g., `group.by.vars = c("dataset", "patient_id")`)
- `NKXSpleen` bulk DESeq2 used `~batch + condition` formula, implying a two-batch structure, but this is not reflected in any scRNA-seq Harmony recipe

---

### 5.6 `conda run` Flag Inconsistency

- User's CLAUDE.md (system rules) mandates: `conda run --no-capture-output -n r-env Rscript`
- `shared/r_environment.md` shows: `conda run -n r-env Rscript` (without `--no-capture-output`)
- Multiple methods files copy the pattern from `r_environment.md`
- **The canonical rule in the library is wrong.** Without `--no-capture-output`, the conda run command hangs after the R script completes.
- **v2 action:** Fix `shared/r_environment.md` to include `--no-capture-output`.

---

### 5.7 No RDS Versioning Convention Enforced

- `shared/r_environment.md` mentions "versioned RDS saves" as a script convention
- Methods files use varying patterns: `saveRDS(obj, "output/obj.Rds")`, `saveRDS(obj, paste0("obj_v", VERSION, ".Rds"))`, or no caching at all
- No library-wide convention for checkpoint RDS file naming
- The two-stage caching architecture in `bulk_lfc_concordance_heatmap.md` is the best-documented example but not extracted as a shared pattern
