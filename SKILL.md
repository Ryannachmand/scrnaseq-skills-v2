---
# scRNAseq Skills Library v2 — Master Router
# Authored in v2 Phase 1
---

# scRNAseq Skills Library v2

Entry point for the deployment agent. Read this file first when generating any project CLAUDE.md.

---

## Library Structure

```
~/claude-skills-v2/
├── SKILL.md                      ← you are here
├── CONVENTIONS.md                ← naming rules, reference syntax, brief schema
├── PROJECT_CLAUDE_TEMPLATE.md    ← blank template for new project CLAUDE.md files
├── primitives/                   ← reusable functions and rules (always inject)
│   ├── r_environment.md          ← conda run command, script structure, package list
│   ├── seurat_v5_rules.md        ← permanent Seurat v5 bug fixes
│   ├── file_downloads.md         ← Box.com download URL pattern
│   ├── geo_download.md           ← GEO acquisition priority rules
│   ├── aesthetics.md             ← visual principles (typography, color, layout)
│   ├── visualization.md          ← parameterized R plot recipes
│   ├── differential_expression.md ← is_ambient, is_confound, run_findmarkers, make_* functions
│   ├── harmony_integration.md    ← canonical Seurat+Harmony processing sequence
│   └── aucell_scoring.md         ← run_aucell, add_auc_to_seurat helpers
├── context/                      ← lab-level context (always inject)
│   ├── lab_context.md            ← default organism, tissues, QC defaults, pipeline conventions
│   ├── color_palettes.md         ← canonical color palettes for all contexts
│   └── validated_examples.yaml   ← project registry with per-project overrides
├── modules/                      ← Phase 2 (not yet written)
└── pipelines/                    ← Phase 2 (not yet written)
```

---

## Routing Instructions for the Deployment Agent

### Always inject (every project CLAUDE.md)

1. All files in `primitives/` — these contain universal rules and reusable functions
2. All files in `context/` — these provide lab-level defaults and the project registry

### Conditionally inject (based on brief's `downstream_analyses` keys)

| Brief key | File to inject |
|---|---|
| `downstream_analyses.deg` | (covered by primitives/differential_expression.md) |
| `downstream_analyses.metabolic_profile` | modules/metabolic_profile.md (Phase 2) |
| `downstream_analyses.cellchat` | modules/interactome_cellchat.md (Phase 2) |
| `downstream_analyses.trajectory` | modules/trajectory_monocle3.md (Phase 2) |
| `downstream_analyses.pyscenic` | modules/pyscenic_regulon_analysis.md (Phase 2) |
| `downstream_analyses.bulk_concordance` | modules/bulk_concordance.md (Phase 2) |
| `downstream_analyses.subclustering` | modules/celltype_subclustering.md (Phase 2) |
| `downstream_analyses.cohort_plots` | modules/cohort_plots.md (Phase 2) |

### Never auto-inject

- `examples/` — project-specific instantiation files; only inject when explicitly working on that project
- Modules not referenced in the brief's `downstream_analyses` block

### Per-project overrides

When starting work on a known project, look up its entry in `context/validated_examples.yaml`
and inject the project's `context_defaults` block into the CLAUDE.md. These are the
project-specific palettes, column names, and gene sets the deployment agent should use.

---

## How to Generate Project Files

1. Identify the pipeline the user wants (LargeDataset or IntegratePublicData — ask if unclear)
2. Read this SKILL.md
3. Read CONVENTIONS.md
4. Read all files in primitives/
5. Read all files in context/
6. If the project is in validated_examples.yaml, read its context_defaults
7. Conditionally read modules/ files based on brief's downstream_analyses keys
8. Generate two files in the user's project directory:
   - `CLAUDE.md` — assembled from primitives + context + relevant modules + project overrides
   - `analysis_brief.txt` — from the pipeline's brief_template (Phase 2) or inline schema from CONVENTIONS.md
9. Tell the user which fields in the brief need to be filled

---

## Available Pipelines

| Pipeline | Description | Status |
|---|---|---|
| LargeDataset | Multi-sample scRNAseq from .h5 files: assemble → QC → cluster → annotate → subset → analyze | Active |
| IntegratePublicData | Integrate heterogeneous public + in-house datasets into one harmonized object | Active |

Pipeline definitions live in pipelines/ (Phase 2). Until Phase 2 is complete, reference
the v1 pipeline definitions in ~/claude-skills/ for stage-by-stage instructions.
