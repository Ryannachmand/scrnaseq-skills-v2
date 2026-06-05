---
# Project CLAUDE.md Template — v2
# Place this file in any new analysis project directory as CLAUDE.md.
# Migrated approach from ~/claude-skills/PROJECT_CLAUDE_TEMPLATE.md — updated to reference v2.
---

# Project CLAUDE.md
# This file connects this analysis directory to the shared skills library v2.

@~/claude-skills-v2/SKILL.md

# ── Project-specific overrides ────────────────────────────────────────────────
# The deployment agent populates the block below from validated_examples.yaml
# when working on a known project, or from the analysis_brief.txt for new projects.

# project_name: project_specific
# organism: project_specific
# tissue: project_specific
# output_dir: project_specific
#
# label_col: project_specific     # metadata column holding cell type labels
# group_col: project_specific     # metadata column holding comparison groups
# sample_col: project_specific    # metadata column identifying samples/batches
#
# context_overrides:              # project-specific palettes, gene sets, label order
#   group_colors: project_specific
#   subtype_colors: project_specific
#   functional_gene_sets: project_specific
#   label_order: project_specific

## Universal Rules

At every stage and in every continuation session, re-read all primitives referenced
by the current stage before writing any code or plot. Do not rely on recalled content
from earlier sessions.

1. For any cross-dataset comparison plot (in-house data vs public atlas), use
   @modules/cross_dataset_dotplot.md. Never call DotPlot(group.by="organ") across
   datasets directly. The module applies within-dataset z-scoring (depth correction)
   and enforces SOURCES_ORDER column ordering (in-house left, atlas right). Re-read
   this module before writing any cross-dataset visualization code.

2. Before writing any subclustering script, re-read @modules/celltype_subclustering.md
   in full. The Mode A endpoint section includes required proportion plots
   (CellType_proportion_by_group, CellType_proportion_by_sample) that must be produced
   alongside the labeled UMAP. These are not optional.

## Cross-Dataset Dotplot Context (deployment agent: inject when brief contains cross_dataset_dotplot or IntegratePublicData pipeline)

When generating CLAUDE.md for a project whose brief contains `downstream_analyses.cross_dataset_dotplot`
or whose `pipeline_chain` includes IntegratePublicData, inject the following block:

  ## Cross-Dataset Dotplot Context
  Re-read @modules/cross_dataset_dotplot.md and @examples/tabula_sapiens_dotplot.md
  before writing Step 2 of Stage 8. Key conventions:
  - SOURCES_ORDER: in-house source first (left), atlas source second (right)
  - Depth correction: within-dataset z-scoring per module Part B
  - Marker sections: use @context/lab_context.md EC biological gene categories
  - Column label convention: "{source} {subtype}" for in-house, organ name for atlas

Detect the IntegratePublicData pipeline in pipeline_chain and auto-inject this block.
