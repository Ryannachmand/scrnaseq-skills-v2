# Deletion Candidates — ~/claude-skills/
*Generated: 2026-05-11 by inventory agent*

Content that probably should NOT carry forward into v2. For each: location, why it's a deletion candidate, what (if anything) should replace it.

---

## D1 — `projectData(cellchat, PPI.human)` in Inference Pipeline Skeleton

**Location:** `pipelines/LargeDataset/methods/interactome_cellchat.md` line 71

**Why delete:** The same file's Critical Constraints section explicitly states this function was removed in CellChat v2 and must not be called. Line 71 is inside an "Inference Pipeline" code block that will be copied into generated scripts. Any script generated from this file will include a call to a non-existent function and fail with `Error: could not find function "projectData"`. This is a live bug, not a cosmetic issue.

**What should replace it:** The line should be deleted. The inference pipeline skeleton should be updated to reflect the v2 API without `projectData`. The Critical Constraints section is correct; the code block is wrong.

---

## D2 — `eye-test-01` Stub Entry in `validated_examples.yaml`

**Location:** `lab_context/validated_examples.yaml` lines ~1–40 (the first job entry)

**Why delete:** The entry is entirely unfilled — every value field contains `# TODO` placeholders. It was created as a schema template but never filled in. It provides no reference value (the schema is better documented in the file's own comments). Having an unfilled stub as the first entry misleads any reader or agent that iterates over job records expecting validated data.

**What should replace it:** The entry should be removed. The file's YAML schema (field names + type annotations) can be preserved as a comment block at the top of the file, or promoted to a standalone `examples/validated_job_template.yaml` stub in v2.

---

## D3 — Commented-Out `#- Bone marrow stroma` in `lab_context.md`

**Location:** `lab_context.md` line 22

**Why delete:** The commented-out entry signals that BoneMarrowStroma is no longer an active tissue focus, but it leaves a trace of past work that future agents might misinterpret as a conditional (i.e., "bone marrow stroma is available if uncommented"). It provides no useful information in its current commented state.

**What should replace it:** Remove the comment. If BoneMarrowStroma context is needed for future work, it can be added back when the tissue becomes active again.

---

## D4 — Commented-Out `unified_label_vocabulary` and `known_aliases` in `lab_context.md`

**Location:** `lab_context.md` lines 63–75

**Why delete:** Both sections are entirely commented out with `#` markers and contain no actual content — only placeholder example lines (`#- MSC`, `#"mesenchymal stromal cell": MSC`). These sections are the declared home for the lab's canonical cell type vocabulary, but the vocabulary has never been populated. An empty commented section occupies space and creates a false impression that the vocabulary exists.

**What should replace it:** Two options:
1. Fill these sections with the actual vocabulary (cell type names and aliases actively used in validated runs), then uncomment — this is the right answer but requires human curation.
2. Remove the sections and add a note: "vocabulary defined per-project in label_harmonization.yaml" — honest about the current state.

Recommend option 1 as v2 action: populate from `label_harmonization.md` controlled vocabulary (lines 138–146) and from project run records.

---

## D5 — `🔲 placeholder` Interactome Entry in `LargeDataset/pipeline.md`

**Location:** `pipelines/LargeDataset/pipeline.md` Stage 8 (optional_analyses section), `interactome` entry

**Why delete:** The entry still shows `"🔲 placeholder"` despite `interactome_cellchat.md` being a validated, 1302-line methods file (confirmed validated on KidneyNew run per SKILL.md). The stale placeholder contradicts SKILL.md and creates confusion about whether CellChat is available.

**What should replace it:** Update Stage 8 interactome entry to `"✅ validated KidneyNew + HumanFat runs"` with a reference to `methods/interactome_cellchat.md`. This is a minor edit, not a deletion, but the placeholder marker is the deletion candidate.

---

## D6 — Invalid R Range Syntax in `metabolic_profile.md` Gene Sets

**Location:** `pipelines/LargeDataset/methods/metabolic_profile.md` OXPHOS and FA_Synthesis gene sets

**Why delete:** The gene set definitions use pseudo-notation like `"NDUFA1"-"NDUFA10"` and `"ELOVL1"-"ELOVL7"` to suggest ranges. This is NOT valid R syntax — the `-` operator on character strings does not produce a vector of names. Any script that copies this directly will produce an error or, worse, silently produce `NA`.

**What should replace it:** The gene sets must be replaced with explicit character vectors:
```r
OXPHOS = c("NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10", ...)
```
The actual gene membership should be validated against a reference (e.g., MSigDB HALLMARK_OXIDATIVE_PHOSPHORYLATION or KEGG Oxidative phosphorylation) since the range pseudo-notation does not tell us the exact intended genes.

---

## D7 — `de_comprehensive_csv.md` Reference to Non-Existent `is_confound()`

**Location:** `pipelines/LargeDataset/methods/de_comprehensive_csv.md` line 116 and ~line 215

**Why delete / fix:** The reference to `is_confound()` (defined in `differential_expression.md`) in line 116 is factually wrong on two counts: (1) `differential_expression.md` defines `is_ambient()`, not `is_confound()`; (2) even `is_ambient()` does not cover the additional exclusions (sex-linked, HLA class II, histones, ENSG IDs) that line ~215 says `is_confound()` should cover. The comment at line ~215 is thus a misleading pointer to a function that doesn't exist.

**What should replace it:** Write `is_confound()` as a real function (see coverage_gaps.md §2.1 for the specification) and update line 116 to call it. Until `is_confound()` is written, the file should use `!is_ambient(gene)` as a documented interim fallback with a `# TODO: replace with is_confound()` comment.

---

## D8 — `source_file` as Default `batch_correction_var` in `brief_template.txt`

**Location:** `pipelines/LargeDataset/brief_template.txt` (the `batch_correction_var` default field)

**Why delete:** `source_file` is a HumanFat-specific metadata column name. It is not a general default and it contradicts `lab_context.md`'s `batch_correction_var: sample_id`. Keeping `source_file` as the default means any new project that follows the template without reading the notes will attempt to run Harmony on a non-existent column.

**What should replace it:** Change the default to `sample_id` (matching `lab_context.md`) or remove the default entirely with a comment: "# REQUIRED: set to the metadata column that identifies samples/batches".

---

## D9 — HumanFat-Specific `metadata_display_cols` Default in `brief_template.txt`

**Location:** `pipelines/LargeDataset/brief_template.txt` — `metadata_display_cols` default value and `deg.col: tissue_type`

**Why delete:** These defaults (`tissue_type`, `additional_notes`, `age`, `bmi`, `sex`) are HumanFat metadata column names. A new project without these columns will produce errors when the pipeline tries to visualize metadata that doesn't exist.

**What should replace it:** Clear or generalize the defaults to blank entries with comment annotations like:
```yaml
metadata_display_cols: []  # OPTIONAL: list metadata columns to show on diagnostic UMAPs
```

---

## D10 — Duplicate KidneyNew / fat-scrnaseq-continued Notes in `lab_context.md` and `validated_examples.yaml`

**Location:** 
- `lab_context.md` Per-Project Notes sections for fat-scrnaseq-continued, KidneyNew, EyePublicData
- `lab_context/validated_examples.yaml` corresponding entries

**Why delete (from one location):** Both files contain the same gotchas and run notes for these three projects. The information is maintained in two places with no clear ownership. Any correction to one file must be mirrored to the other manually.

**What should replace it:** In v2, `validated_examples.yaml` should be the single source of per-project run notes. The per-project notes sections in `lab_context.md` should be removed, leaving only the structural defaults (organism, tissues, QC thresholds, color palettes, pipeline conventions) that apply globally.

---

## D11 — `shared/r_environment.md` Stale `conda run` Command

**Location:** `shared/r_environment.md` (the conda run command example)

**Why delete / fix:** Shows `conda run -n r-env Rscript` without `--no-capture-output`. Per the user's CLAUDE.md (system rules), `--no-capture-output` is required because omitting it causes the command to hang after the R script completes. The flag was added after the library file was written. Every methods file that copies this pattern from `r_environment.md` will produce hanging scripts.

**What should replace it:** Fix the canonical command in `r_environment.md` to `conda run --no-capture-output -n r-env Rscript`. All methods files that have hardcoded the old pattern should be updated to match.

---

## Summary Table

| # | Location | Type | Priority | Replacement |
|---|---|---|---|---|
| D1 | interactome_cellchat.md line 71 | Invalid code (live bug) | **CRITICAL** | Delete `projectData()` call from skeleton |
| D6 | metabolic_profile.md gene sets | Invalid code (silent fail) | **CRITICAL** | Replace with explicit gene vectors |
| D7 | de_comprehensive_csv.md | Phantom function reference | **CRITICAL** | Write `is_confound()` or use `is_ambient()` interim |
| D11 | shared/r_environment.md | Stale command (hanging bug) | **HIGH** | Add `--no-capture-output` flag |
| D8 | brief_template.txt batch_correction_var | Wrong default | **HIGH** | Change to `sample_id` or blank |
| D9 | brief_template.txt metadata_display_cols | Wrong defaults | **HIGH** | Clear to generic/blank |
| D2 | validated_examples.yaml eye-test-01 | Empty stub | **Medium** | Remove or fill |
| D5 | LargeDataset/pipeline.md Stage 8 | Stale placeholder | **Medium** | Update to reflect validated interactome |
| D10 | lab_context.md per-project notes | Duplicate content | **Medium** | Remove from lab_context.md; keep in validated_examples.yaml |
| D3 | lab_context.md commented Bone marrow | Stale comment | **Low** | Remove |
| D4 | lab_context.md unified_label_vocabulary | Empty section | **Low** | Fill or remove |
