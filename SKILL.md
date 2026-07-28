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
│   ├── known_atlases.md          ← known-atlas registry (Tabula Sapiens + placeholders)
│   └── validated_examples.yaml   ← project registry with per-project overrides
├── modules/                      ← v2 modules (14 files; Phase 2 complete)
│   ├── load_formats.md           ← multi-format data loading (h5/mex/rds/geo)
│   ├── cohort_plots.md           ← correlation heatmap + chord + proportion (IntegratePublicData)
│   ├── label_harmonization.md    ← YAML mapping + Seurat label transfer → unified_label
│   ├── feature_umap_plot.md      ← custom FeaturePlot; matches Seurat exactly
│   ├── celltype_subclustering.md ← two-phase subset + re-cluster workflow
│   ├── metabolic_profile.md      ← AUCell scoring of pathway gene sets
│   ├── trajectory_monocle3.md    ← Monocle3 pseudotime trajectory (SeuratWrappers bypass)
│   ├── pyscenic_regulons.md      ← pySCENIC TF regulon analysis (Python + R)
│   ├── bulk_concordance.md       ← bulk OE concordance (Mode 1: score / Mode 2: LFC heatmap)
│   ├── de_comprehensive_csv.md   ← vectorized Kruskal-Wallis comprehensive DE table
│   ├── multi_group_de_analysis.md ← N-group pairwise DE with dotplots + GO functional plots
│   ├── atlas_co_umap.md          ← in-house vs atlas co-embedding UMAP
│   ├── cross_dataset_dotplot.md  ← depth-corrected 3-section cross-dataset dotplot
│   └── cellchat.md               ← CellChat v2 ligand-receptor analysis (6 scripts)
├── pipelines/                    ← pipeline orchestration layer (Phase 3 complete)
│   ├── large_dataset/
│   │   ├── pipeline.md           ← 9-stage LargeDataset orchestration + Stage 8 dispatch
│   │   └── brief_template.txt    ← v2 brief template for LargeDataset projects
│   └── integrate_public_data/
│       ├── pipeline.md           ← 8-stage IntegratePublicData + Known-Atlas Convention
│       └── brief_template.txt    ← v2 brief template for IntegratePublicData projects
└── examples/                     ← project-specific instantiations (public subset)
    ├── cellchat_example_A.md
    ├── cohort_plots_example_A.md
    ├── de_comprehensive_example_A.md
    ├── ec_functional_gene_sets.md
    ├── load_formats_example_A.md
    ├── tabula_sapiens_co_umap.md
    ├── tabula_sapiens_dotplot.md
    └── trajectory_example_A.md
```

---

## Routing Instructions for the Deployment Agent

### Always inject (every project CLAUDE.md)

1. All files in `primitives/` — these contain universal rules and reusable functions
2. All files in `context/` — these provide lab-level defaults and the project registry

### Conditionally inject (based on brief's `downstream_analyses` keys)

Stage 8 dispatch: the deployment agent reads the `downstream_analyses` block of the brief
and injects the corresponding module file. Execution order: `multi_group_de` first,
`bulk_concordance` second (may require DE outputs), all others in parallel.
See pipeline.md Stage 8 dispatch table for the full specification.

| Brief key | Module file | Notes |
|---|---|---|
| `downstream_analyses.deg` | (covered by primitives/differential_expression.md) | Basic DE only |
| `downstream_analyses.multi_group_de` | modules/multi_group_de_analysis.md | N-group pairwise DE; run first |
| `downstream_analyses.bulk_concordance` | modules/bulk_concordance.md | Mode 1 (score) or Mode 2 (LFC); run second |
| `downstream_analyses.de_comprehensive_csv` | modules/de_comprehensive_csv.md | Vectorized KW comprehensive table |
| `downstream_analyses.metabolic_profile` | modules/metabolic_profile.md | AUCell pathway scoring |
| `downstream_analyses.cellchat` | modules/cellchat.md | CellChat v2 ligand-receptor (6 scripts) |
| `downstream_analyses.trajectory` | modules/trajectory_monocle3.md | Monocle3 pseudotime |
| `downstream_analyses.pyscenic` | modules/pyscenic_regulons.md | pySCENIC TF regulon analysis |
| `downstream_analyses.subclustering` | modules/celltype_subclustering.md | Two-phase subset + re-cluster |
| `downstream_analyses.feature_umap_plot` | modules/feature_umap_plot.md | Custom FeaturePlot overlay |
| `downstream_analyses.cohort_plots` | modules/cohort_plots.md | IntegratePublicData only |
| `downstream_analyses.atlas_co_umap` | modules/atlas_co_umap.md | IntegratePublicData + Known-Atlas Convention |
| `downstream_analyses.cross_dataset_dotplot` | modules/cross_dataset_dotplot.md | IntegratePublicData only |
| `label_harmonization.*` | modules/label_harmonization.md | IntegratePublicData Stage 3 |
| `inputs.datasets[*]` | modules/load_formats.md | IntegratePublicData Stage 1 |

### Never auto-inject

- `examples/` — project-specific instantiation files; consumed only when:
  1. The deployment agent's known-project lookup matches (validated_examples.yaml registry hit)
  2. The user explicitly references an example for documentation purposes
  3. As reference material for the Phase 6 deployment agent
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
5. Read all files in context/ (including known_atlases.md)
6. If the project is in validated_examples.yaml, read its context_defaults
7. Read the appropriate pipeline.md (`pipelines/large_dataset/pipeline.md` or
   `pipelines/integrate_public_data/pipeline.md`) for stage-by-stage instructions
8. Conditionally read modules/ files based on brief's downstream_analyses keys
   (see Stage 8 dispatch table in the pipeline.md for execution order)
9. If the brief contains `atlas:`, look up the atlas in context/known_atlases.md
   and apply the Known-Atlas Convention defaults to atlas_co_umap and cross_dataset_dotplot
10. Generate two files in the user's project directory:
    - `CLAUDE.md` — assembled from primitives + context + relevant modules + project overrides
    - `analysis_brief.txt` — from `pipelines/<pipeline>/brief_template.txt`
11. Tell the user which fields in the brief need to be filled

---

## Available Pipelines

| Pipeline | Description | Brief template | Status |
|---|---|---|---|
| LargeDataset | Multi-sample scRNAseq from .h5 files: assemble → QC → cluster → annotate → subset → analyze | pipelines/large_dataset/brief_template.txt | Active |
| IntegratePublicData | Integrate heterogeneous public + in-house datasets into one harmonized object | pipelines/integrate_public_data/brief_template.txt | Active |

Pipeline definitions: `pipelines/large_dataset/pipeline.md` and
`pipelines/integrate_public_data/pipeline.md`. Both are complete (Phase 3).
