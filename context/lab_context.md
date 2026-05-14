---
# Lab Context and Conventions — v2
# Migrated from ~/claude-skills/lab_context.md
# Stripped per D3/D4/D10 deletion candidates and Phase 1 scope:
#   - Removed per-project notes (fat-scrnaseq-continued, KidneyNew, EyePublicData)
#     → now live in context/validated_examples.yaml per-project context_defaults
#   - Removed project-specific color palettes (ec_subtype_colors_adipose, ec_subtype_colors_kidney,
#     adipose_type_colors, organ_tissue_colors, celltype_colors) → now in validated_examples.yaml
#   - Removed commented-out bone marrow stroma (D3)
#   - Removed commented-out unified_label_vocabulary and known_aliases sections (D4)
#   - batch_correction_var set to sample_id (lab standard) — conflict with v1 brief_template.txt
#     source_file default was D8 deletion candidate; lab_context.md was always the authority here
# This file describes lab-level conventions only. Per-project overrides go in validated_examples.yaml.
---

# Lab Context and Conventions

Lab-level defaults for all analyses. Read by the deployment agent when generating CLAUDE.md.
Per-project overrides are in @context/validated_examples.yaml — always check that file for
project-specific palettes, column names, and gene sets before generating scripts.

---

## Default Organism

```yaml
default_organism: Homo sapiens
```

---

## Tissues and Systems Studied

Primary tissues and cell types this lab routinely works with.

```yaml
primary_tissues:
  - Adipose (multi-depot: visceral fat, subcutaneous fat, orbital fat, lipoma, breast fat, myelolipoma)
  - Kidney (glomerular + tubular vasculature, whole organ multi-region)
  - Eye (RPE / Choroid / posterior segment — normal, AMD, nAMD, fetal)

cell_compartments_of_interest:
  - Retinal Pigment Epithelium (RPE)
  - Choroidal endothelial cells
  - Choroidal stromal / fibroblasts
  - Pericytes / Mural cells (choroidal)
  - Macrophages / Microglia (choroidal immune)
  - Endothelial cells (EC subtypes)
  - Stromal cells / ASPCs / Fibroblasts
  - Immune cells (Macrophages/MoMF, Monocytes, T cells, B cells, Mast cells, DCs, NK cells)
  - Vascular smooth muscle cells (SMC)
```

---

## Standard QC Thresholds

Leave blank to have the pipeline suggest thresholds from the data distribution.
Override in the analysis brief for specific projects.

```yaml
default_qc_thresholds:
  nFeature_RNA_min:      # data-driven suggestion
  nFeature_RNA_max:      # data-driven suggestion
  nCount_RNA_min:        # data-driven suggestion
  nCount_RNA_max:        # data-driven suggestion
  percent_mt_max:        # data-driven suggestion
```

---

## Label Harmonization

Canonical cell type vocabulary used across projects. Controlled vocabulary
reduces annotation drift between projects. Keys are the canonical label to
use in all downstream analysis; values list known aliases.

```yaml
unified_label_vocabulary:
  # Populate from label_harmonization.md controlled vocabulary as projects are completed.
  # See label_harmonization.md (Phase 2 module) for the full controlled vocabulary rules.
  # TODO: human review needed to populate from validated project runs.
```

---

## Pipeline Conventions

Defaults applied to all jobs unless overridden in the brief.

```yaml
batch_correction_var: sample_id     # lab default — column passed to RunHarmony
                                     # IMPORTANT: this is sample_id at the lab level.
                                     # Individual projects may override this in their
                                     # context_defaults in validated_examples.yaml.
                                     # (HumanFat used source_file — see that project entry)

whole_object_defaults:
  n_variable_features: 4000
  n_pcs: 30
  clustering_resolution: 0.5

subset_defaults:
  n_variable_features: 2850
  n_pcs: 40
  clustering_resolution: 0.39

downsample_ceiling:                  # e.g. 5000 (cells per sample) — blank = no ceiling
downsample_strategy: sample_level    # sample_level or dataset_level
```

---

## Aesthetics and Color Preferences

Canonical palettes are defined in @context/color_palettes.md.
Project-specific palettes (EC subtypes, tissue types, adipose depots) are in
each project's `context_defaults.palettes` block in validated_examples.yaml.

The deployment agent should never hardcode project-specific color vectors in a
shared module or pipeline file. Always reference palettes via context injection.

---

## Notes for the Deployment Agent

1. Always check validated_examples.yaml for the target project before generating CLAUDE.md.
   The per-project `context_defaults` block supersedes lab-level defaults for that project.

2. `batch_correction_var` at the lab level is `sample_id`. If the project brief uses a
   different column name (e.g., `source_file` for HumanFat, `sample` for KidneyNew),
   that override is in validated_examples.yaml `context_defaults.batch_correction_var`.

3. Color palettes for specific tissues/cell types are in validated_examples.yaml under
   each project entry. Do NOT use color palettes from one project for another project.

4. For returning-to-an-old-project scenarios: read the project's validated_examples.yaml entry
   completely before generating any scripts — it contains the key column names, gotchas, and
   analysis decisions from previous runs.
