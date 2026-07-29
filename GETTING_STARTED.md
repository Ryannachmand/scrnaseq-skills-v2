# Getting Started

This walkthrough gets you from a fresh clone to a running pipeline against
your own scRNA-seq data. It assumes you are comfortable with scRNA-seq
analysis (Seurat, standard QC/clustering/annotation workflows) but new to
Claude Code and agent-driven pipelines.

Estimated time: 30–60 minutes, most of which is environment setup.

---

## 1. Claude Code, briefly

Claude Code is Anthropic's command-line tool that lets you delegate coding
and analysis tasks to Claude from your terminal. Unlike the Claude web
interface, Claude Code can read your files, write scripts, and run commands
in your working directory — which is what makes it able to execute a full
scRNA-seq pipeline end-to-end.

Install and authenticate per Anthropic's docs:
https://docs.claude.com/claude-code

For the purposes of this library, you'll interact with Claude Code in two
ways:

1. **Interactive sessions** — you `cd` into a project directory, launch
   `claude`, and describe what you want. Claude reads the project's
   `CLAUDE.md` file, which points it at this skills library, and works
   with you interactively.

2. **Non-interactive prompts** — you write a prompt as a plain text file
   and pipe it in with `claude --print`. This is how you kick off long
   pipeline runs without staying attached to a session.

You'll use pattern (1) for setup and pattern (2) for the pipeline run.

---

## 2. Repo setup

### Clone

```bash
git clone https://github.com/Ryannachmand/scrnaseq-skills-v2.git ~/claude-skills-v2
```

Several files in the library reference the path `~/claude-skills-v2/`. Cloning
there is the least-friction option. If you clone elsewhere, update
`PROJECT_CLAUDE_TEMPLATE.md` and any references accordingly.

### R environment (`r-env`)

The library expects a conda environment named `r-env` containing R 4.5+,
Seurat v5, and standard scRNA-seq packages. Create it however you normally
manage R environments; the essentials:

- R 4.5 or later
- Bioconductor 3.20+
- Seurat 5.x
- SeuratObject, SeuratDisk, harmony
- ggplot2, patchwork, dplyr, readr, tibble
- clusterProfiler, msigdbr, fgsea (for pathway analysis modules)
- CellChat, monocle3 (for the corresponding modules — install per those
  packages' upstream instructions, as they have specific dependency chains)

You can verify the environment is set up correctly:

```bash
conda run --no-capture-output -n r-env Rscript -e 'library(Seurat); sessionInfo()'
```

The `--no-capture-output` flag is important — omitting it causes hangs in
some setups. Every R command run by the library uses this pattern.

### Python environment for pySCENIC (`scenicenv`) — optional

Only needed if you plan to use the `pyscenic` module. Create a separate
conda env named `scenicenv` with pySCENIC installed per the aertslab
docs. You'll also need the pySCENIC ranking databases (feather files) and
a TF list, both downloadable from the aertslab pySCENIC resources.

Skip this until you actually need regulon analysis; it's not required for
your first run.

---

## 3. The library's structure, in one paragraph

Every analysis Claude Code runs from this library composes three layers:

- **Primitives** (`primitives/`) — generic reusable functions and
  environment conventions (differential expression, visualization, doublet
  removal, etc.). You almost never edit these.
- **Modules** (`modules/`) — higher-level analyses that call primitives
  (CellChat, subclustering, trajectory, pySCENIC). Each module is a
  self-contained workflow.
- **Pipelines** (`pipelines/`) — end-to-end chains of modules for a
  specific analytical goal. Two are shipped: `large_dataset` (analyze one
  in-house dataset end-to-end) and `integrate_public_data` (integrate
  your data with a public reference atlas).

You interact with the library by filling in a **brief** — a YAML-like
config file that names your inputs, sets pipeline parameters, and
optionally activates downstream analyses. Claude Code reads the brief,
loads the skills it needs, and executes the pipeline.

The `context/` directory holds cross-cutting reference material (color
palettes, known public atlases, lab conventions). The `examples/`
directory contains worked examples showing how each module is invoked in
practice.

---

## 4. Suggested reading order

Before running anything, spend 20 minutes reading these in order:

1. `README.md` (you've probably done this)
2. `SKILL.md` — the entry point Claude reads first. Shows how the library
   is meant to be composed.
3. `CONVENTIONS.md` — naming, paths, config style. Small but important.
4. `pipelines/large_dataset/pipeline.md` — the pipeline you're about to
   run, described stage by stage.
5. One example from `examples/` — pick `cellchat_example_A.md` or
   `de_comprehensive_example_A.md` to see how a module invocation actually
   looks.

Don't try to read every primitive or module upfront. Claude Code loads
them as needed; you can read individual files when you're curious about a
specific analysis.

---

## 5. Your first pipeline run

The goal: run the `large_dataset` pipeline on your own scRNA-seq data
(assumed to be 10x cellranger `filtered_feature_bc_matrix.h5` files).

### Step 5.1 — Set up your project directory

```bash
mkdir -p ~/my_project
cd ~/my_project
mkdir -p data output
```

Copy or symlink your cellranger `.h5` files into `~/my_project/data/`.

### Step 5.2 — Create a project `CLAUDE.md`

Copy the template:

```bash
cp ~/claude-skills-v2/PROJECT_CLAUDE_TEMPLATE.md ~/my_project/CLAUDE.md
```

Open `CLAUDE.md` and confirm the reference at the top points at
`~/claude-skills-v2/SKILL.md` (or the path where you actually cloned the
repo).

### Step 5.3 — Fill in the brief

Copy the template:

```bash
cp ~/claude-skills-v2/pipelines/large_dataset/brief_template.txt ~/my_project/brief.txt
```

Open `~/my_project/brief.txt` and fill in the required fields. The literal
string `project_specific` is a sentinel — any field still containing that
string will cause the pipeline to halt loudly, so you can't accidentally
run with unresolved values.

A minimal filled brief for a first run looks like:

```yaml
project_name: "my_first_run"
organism: "human"
tissue: "peripheral_blood"    # whatever your tissue actually is
output_dir: "output/"

inputs:
  type: "h5"
  paths:
    - "data/sample_01_filtered_feature_bc_matrix.h5"
    - "data/sample_02_filtered_feature_bc_matrix.h5"
    - "data/sample_03_filtered_feature_bc_matrix.h5"
  metadata_csv: null

pipeline: "large_dataset"

pipeline_params:
  n_variable_features: 4000
  n_pcs: 30
  clustering_resolution: 0.8
  batch_correction_var: "sample_id"
  qc_thresholds:
    n_feature_min_percentile: 2
    n_feature_max_percentile: 98
    percent_mt_max: 20

annotation:
  top_n_markers: 25

subset:
  target_celltypes: null
  subset_name: "subset"

downstream_analyses: {}

context_overrides:
  palettes:
    group_colors: null
    subtype_colors: null
  metadata_columns:
    group_col: null
    subtype_col: null
    sample_col: "sample_id"
```

For a first run, leave `downstream_analyses` empty (`{}`) and
`subset.target_celltypes` null. This runs the core pipeline — load, QC,
integration, clustering, annotation — without any downstream modules. You
can re-run with downstream analyses activated once the core pipeline
succeeds.

### Step 5.4 — Kick off the pipeline

From `~/my_project`:

```bash
cd ~/my_project
claude
```

Once the interactive session opens, paste this prompt:

```
I'd like to run the large_dataset pipeline on the samples listed in
brief.txt. Please read the brief, then read
~/claude-skills-v2/pipelines/large_dataset/pipeline.md to understand the
stages. Execute the pipeline stage by stage. Pause at any checkpoint
that requires my input before continuing.
```

Claude Code will read the brief, load the pipeline definition, load the
relevant modules and primitives, and start executing. Expect it to pause
at natural checkpoints (after QC thresholds are proposed, after cell type
annotation is drafted, etc.) for you to confirm or override its choices
before it continues.

### Step 5.5 — What to expect

- The pipeline creates outputs in `~/my_project/output/`, organized by
  stage.
- Log messages appear in the Claude Code session as each stage runs.
- Typical wall-clock time for 3–5 samples of standard 10x data on a
  reasonable workstation: 30 minutes to a few hours, dominated by
  integration and clustering steps.

---

## 6. Where to go next

Once your first run completes:

- **Activate downstream analyses.** Uncomment blocks in `brief.txt`
  (multi_group_de, cellchat, trajectory, etc.) and re-run. Each block's
  fields are documented in the corresponding module file under `modules/`.
- **Try the other pipeline.** `pipelines/integrate_public_data/` covers
  integrating your data with a public reference atlas.
- **Read individual skills.** When you want to understand or customize
  what a module is doing under the hood, open the corresponding file in
  `modules/` and follow its references down into `primitives/`.

---

## 7. When things go wrong

- **"project_specific" appears in an error message.** You left a required
  field unfilled in `brief.txt`. Find it and set a value.
- **R session hangs during a `conda run` call.** You forgot
  `--no-capture-output`. All library-generated R calls include this flag;
  if you're running R commands manually and they hang, add it.
- **Claude Code doesn't seem to find the skills.** Check that your
  project `CLAUDE.md` references the correct path to `SKILL.md` in the
  cloned repo location.
- **A module fails partway through a run.** Check `STATUS.md` for
  known-rough areas. Report reproducible issues on GitHub.

For anything else: read the relevant module's `.md` file — the library is
its own documentation.
