---
# Conventions — claude-skills-v2
---

# Conventions

Rules and schemas that govern how the library is structured, how files reference
each other, and how briefs are written. Read this before writing any module or pipeline.

---

## 1. Reference Syntax

Use the `@reference` syntax to tell the deployment agent which library files to read
for a given module or pipeline stage.

```
@primitives/r_environment.md
@primitives/differential_expression.md
@modules/interactome_cellchat.md
@context/color_palettes.md
```

All paths are relative to `~/claude-skills-v2/`. The deployment agent expands
`@reference` paths and injects their content into the generated CLAUDE.md.

---

## 2. Context Dependency Declaration

Every module and pipeline file declares its context dependencies in a YAML frontmatter
`requires_context` block at the top of the file:

```yaml
---
requires_context:
  palettes:
    - group_colors          # required — named vector keyed by group values
    - subtype_colors        # optional — named vector for EC/cell subtype colors
  metadata_columns:
    required:
      - label_col           # cell type label column in Seurat metadata
      - group_col           # comparison group column
    optional:
      - sample_col          # sample/batch column (default: sample_id)
      - batch_col           # batch correction column (default: sample_id)
  brief_keys:
    required:
      - project_name
      - organism
      - output_dir
    optional:
      - downstream_analyses.deg
---
```

The deployment agent checks each `required` field against the brief and the
project's `context_defaults` in `validated_examples.yaml`. If a required field
is missing, the agent must prompt the user before generating scripts.

---

## 3. Function Naming Conventions

| Pattern | Usage |
|---|---|
| `snake_case` | All function and variable names |
| `make_<plot>` | Plot-generating functions that return a ggplot object or ComplexHeatmap |
| `run_<analysis>` | Functions that execute a statistical or computational step |
| `is_<predicate>` | Predicate functions that return a logical vector |
| `get_<accessor>` | Accessor/extraction functions |
| `add_<operation>` | Functions that modify an object by adding content (e.g., metadata) |

---

## 4. Output Directory Conventions

```
output/<stage_name>/     ← for pipeline stages (e.g. output/qc/, output/annotation/)
output/<module_name>/    ← for modules (e.g. output/deg/, output/metabolic/)
```

Every plot function accepts an `output_dir` argument and writes via:
```r
ggsave(file.path(output_dir, filename), device = "pdf", useDingbats = FALSE, units = "in")
```

Use `ggsave()` for all ggplot outputs. The exception list below covers libraries with
no ggsave path. Any other use of `pdf()` / `dev.off()` in a generated script is a violation.

**Documented exceptions (closed list):**

| # | Library / pattern | Reason | How to use |
|---|---|---|---|
| 1 | `ComplexHeatmap::draw()` | `draw()` sends output to the active graphics device; `ggsave()` cannot intercept it | `pdf(path); draw(ht); dev.off()` |
| 2 | `pheatmap` with `grid.text()` title overlay | `grid.text()` draws into the active device mid-render; the overlay cannot be captured by `ggsave()` | `pdf(path); grid.draw(ph$gtable); grid.text(...); dev.off()` |
| 3 | `circlize::chordDiagram()` | circlize uses base R graphics with no ggplot/ggsave path | `pdf(path); chordDiagramFromMatrix(...); dev.off()` |

To add a new exception: update this table in CONVENTIONS.md first, then use the pattern in
the module. Agents must not use `pdf()` / `dev.off()` for any other purpose.

---

## 5. Brief Schema

The analysis brief is a YAML file with the following top-level structure:

```yaml
# ── Identity (REQUIRED) ──────────────────────────────────────────────────────
project_name: "MyProject"        # human-readable name, used in plot titles
organism: "Homo sapiens"         # or "Mus musculus"
tissue: "Kidney"                 # anatomical origin
output_dir: "output"             # root output directory (relative to project dir)

# ── Data Inputs (REQUIRED) ────────────────────────────────────────────────────
inputs:
  type: "h5"                     # h5 | mex | rds | geo
  paths:                         # list of file paths or GEO IDs
    - "data/sample1.h5"
    - "data/sample2.h5"

# ── Pipeline (REQUIRED) ──────────────────────────────────────────────────────
pipeline: LargeDataset           # LargeDataset | IntegratePublicData

# ── Processing Parameters (CONFIRM — review before running) ──────────────────
pipeline_params:
  n_variable_features: 4000      # lab default: 4000 (whole object), 2850 (subset)
  n_pcs: 30                      # lab default: 30 (whole object), 40 (subset)
  clustering_resolution: 0.5     # lab default: 0.5 (whole), 0.39 (subset)
  batch_correction_var: sample_id # REQUIRED: metadata column holding sample/batch IDs

# ── Metadata (CONFIRM) ────────────────────────────────────────────────────────
metadata:
  label_col: project_specific    # metadata column holding cell type labels
  sample_col: sample_id          # metadata column identifying samples
  group_col: project_specific    # metadata column holding comparison groups (for DE)
  display_cols: []               # metadata columns to show on diagnostic UMAPs

# ── QC Thresholds (OPTIONAL — blank = data-driven suggestion) ─────────────────
qc:
  nFeature_RNA_min:
  nFeature_RNA_max:
  percent_mt_max:

# ── Downstream Analyses (OPTIONAL) ────────────────────────────────────────────
downstream_analyses:
  deg:
    enabled: false
    group_var: project_specific          # metadata column holding group labels for DE
    comparisons: project_specific        # list of {label, ident1, ident2, col} dicts
    subsets: project_specific            # list of cell type subsets to run separately
    functional_gene_sets: project_specific  # named list of gene vectors (biology-specific)

  metabolic_profile:
    enabled: false
    gene_sets: project_specific          # named list of gene set vectors (caller-provided)
    label_col: project_specific

  cellchat:
    enabled: false
    cell_types: project_specific         # list of cell types to include

  trajectory:
    enabled: false
    label_col: project_specific
    start_cluster: project_specific
    end_cluster: project_specific
    exclude_pattern: project_specific

  pyscenic:
    enabled: false
    label_col: project_specific
    database_dir: project_specific       # path to pySCENIC database files
    tf_list: project_specific            # path to TF list file

# ── Context Overrides (OPTIONAL) ─────────────────────────────────────────────
context_overrides:
  palettes:
    group_colors: project_specific       # named vector of group → hex color
    subtype_colors: project_specific     # named vector of subtype → hex color
  functional_gene_sets: project_specific # see downstream_analyses.deg.functional_gene_sets
  label_order: project_specific          # ordered character vector of cell type names
  batch_correction_var: project_specific # override lab default (sample_id)
```

---

## 6. The `project_specific` Sentinel

When a primitive or module requires project-specific input that cannot be generalized,
it uses `project_specific` as a sentinel value:

```yaml
label_col: project_specific
functional_gene_sets: project_specific
```

`project_specific` means: **the deployment agent must replace this value** with the
actual project value from either (a) the user's brief, or (b) the project's
`context_defaults` block in `validated_examples.yaml`.

If neither source provides a value, the agent must prompt the user before generating scripts.
Any generated R script that contains `project_specific` as a literal value will fail — this
is intentional and serves as a flag that something was not resolved.

In R code blocks, project-specific inputs appear as comments:
```r
LABEL_COL <- "project_specific"  # REPLACE with actual column name from brief
```

---

## 7. Self-Check Before Committing a File

Before finalizing any module or primitive, verify:
- No project-specific column names hardcoded (every column is a function argument or sentinel)
- No project-specific color vectors hardcoded (plots accept color arguments)
- No project-specific gene sets hardcoded (gene sets are caller-provided)
- Every plot function uses `ggsave(..., device='pdf', useDingbats=FALSE)`
- No `pdf()` / `dev.off()` patterns — unless the use falls under the three documented exceptions in §4
- Function names follow the `make_`, `run_`, `is_`, `get_`, `add_` conventions
