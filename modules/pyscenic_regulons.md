---
requires_context:
  palettes:
    - subtype_colors   # optional — named vector for cell subtype colors in UMAP/violin plots
  metadata_columns:
    required:
      - label_col      # cell subtype label column for AUCell comparison plots
    optional:
      - sample_col     # sample column exported with metadata (for provenance)
  brief_keys:
    required:
      - output_dir
      - downstream_analyses.pyscenic.label_col
      - downstream_analyses.pyscenic.scenic_python_path
      - downstream_analyses.pyscenic.database_dir
      - downstream_analyses.pyscenic.tf_list
    optional:
      - downstream_analyses.pyscenic.n_workers
      - downstream_analyses.pyscenic.min_pct_cells
      - downstream_analyses.pyscenic.min_mean_expr
references:
  - "@primitives/seurat_v5_rules.md"
  - "@primitives/r_environment.md"
---

# Module: TF Regulon Analysis — pySCENIC

pySCENIC infers transcription factor (TF) regulons — TF + co-expressed target genes
supported by motif evidence — and scores each cell for regulon activity (AUC score).

**Three-step pipeline:**
1. **GRN** (arboreto/GRNBoost2): gene regulatory network inference (~7–8 h)
2. **ctx** (cisTarget): motif enrichment to validate TF-target links → regulons (~15–20 min)
3. **AUCell**: score each cell for regulon activity (~5 min)

**Typical use case:** identify subtype-specific TF programs after cell type subset
annotation is finalized.

---

## P2-8 Fix: Machine-Specific Paths

v1 contained three hardcoded absolute paths that break on any other machine.
All three are now `project_specific` sentinel values resolved from the brief.

| v1 hardcoded | v2 brief field | Notes |
|---|---|---|
| `<YOUR_PYSCENIC_DATABASE_DIR>` | `database_dir` | Directory containing the hg38 feather ranking databases and motif annotation table |
| `<YOUR_TF_LIST_PATH>` | `tf_list` | Path to TF list file. The pySCENIC install location is machine-specific. |
| `<YOUR_SCENICENV_PYTHON>` | `scenic_python_path` | Full path to the Python binary in the scenicenv conda environment. Use `conda run -n scenicenv which python` to find it on your machine. |

---

## Critical Bugs & Fixes in pySCENIC 0.12.1

| Bug | Symptom | Fix |
|---|---|---|
| `np.float` removed in NumPy >= 1.20 | `pyscenic aucell` CLI crashes during binarization with `AttributeError: module 'numpy' has no attribute 'float'` | Bypass CLI; use Python API directly via `scenic_01b_aucell.py` (Script 4) |
| `load_motifs()` returns DataFrame, not GeneSignature list | `aucell()` crashes with `ValueError: The truth value of a DataFrame is ambiguous` | Convert with `df2regulons(motifs_df)` before passing to `aucell()` |

**Never use `pyscenic aucell` CLI** on this environment — always use the Python API script.

---

## Resource Requirements

| Step | Tool | Runtime (32k cells, ~9k genes, 16 workers) |
|---|---|---|
| GRN | arboreto/GRNBoost2 | ~7–8 hours |
| ctx | pyscenic ctx | ~15–20 min |
| AUCell | Python API | ~5 min |
| Plots | scenic_02_plots.py | ~2 min |

Runtime scales roughly with n_genes². Filter genes aggressively (see Script 1).

**Always run with `nohup` using the direct Python binary** — `conda run` kills
background processes when the Claude Code session ends.

---

## Brief Schema

```yaml
downstream_analyses:
  pyscenic:
    enabled: true
    label_col: project_specific      # REPLACE: metadata column holding cell subtype labels
    scenic_python_path: project_specific  # REPLACE: full path to scenicenv Python binary
                                           # Find with: conda run -n scenicenv which python
    database_dir: project_specific    # REPLACE: directory with hg38 ranking databases
                                       # Must contain: *.feather files + motifs-*.tbl file
    tf_list: project_specific         # REPLACE: path to allTFs_hg38.txt (or equivalent)
                                       # Typically inside the pySCENIC install directory
    n_workers: 16                     # number of parallel workers for GRN step
    min_pct_cells: 0.05              # gene expressed in >= this fraction of cells (for filtering)
    min_mean_expr: 0.05              # gene mean expression minimum (for filtering)
```

---

## Database Files (hg38 human) — in `database_dir`

| File | Used for |
|---|---|
| `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` | ctx step |
| `hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` | ctx step |
| `motifs-v10-nr.hgnc-m0.00001-o0.0.tbl` | ctx motif annotation |

All gene names must be **HGNC symbols** (not Ensembl IDs). The Seurat RNA assay
already uses HGNC symbols — no conversion needed.

---

## Script Architecture (4 scripts)

### Script 1: `scenic_00_export.R` — R (r-env)

Exports log-normalized RNA counts and cell metadata from Seurat.

```r
library(Seurat)
library(data.table)
set.seed(42)

# ── CONFIG ────────────────────────────────────────────────────────────────────
LABEL_COL         <- "project_specific"   # REPLACE: cell subtype label column
SAMPLE_COL        <- "project_specific"   # REPLACE: sample column to export; or NULL
RDS_IN            <- "project_specific"   # REPLACE: path to Seurat subset object
OUTPUT_SCENIC     <- "output_scenic"
MIN_PCT_CELLS     <- 0.05
MIN_MEAN_EXPR     <- 0.05

dir.create(OUTPUT_SCENIC, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(RDS_IN)

# CRITICAL settings — see @primitives/seurat_v5_rules.md Rule 1
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)

expr_mat <- GetAssayData(obj, layer = "data")

# Gene filter — both thresholds required to keep runtime tractable
min_cells <- ceiling(MIN_PCT_CELLS * ncol(expr_mat))
mean_expr  <- rowMeans(expr_mat)
gene_keep  <- (rowSums(expr_mat > 0) >= min_cells) & (mean_expr >= MIN_MEAN_EXPR)
# Target: ~8,000–10,000 genes. If more, tighten MIN_PCT_CELLS to 0.07–0.10.
message(sprintf("Genes passing filter: %d", sum(gene_keep)))

# Export cells x genes (dense) — use data.table::fwrite for speed
expr_df <- t(as.matrix(expr_mat[gene_keep, ]))
data.table::fwrite(
  as.data.frame(expr_df),
  file      = file.path(OUTPUT_SCENIC, "expr.csv"),
  row.names = TRUE,
  quote     = FALSE
)
# WARNING: Dense CSV for 32k cells x 9k genes is ~1.5 GB.
# Delete it after the loom is created — it is not needed again.

# Metadata: include subtype label column + UMAP coordinates
meta_cols <- c(LABEL_COL)
if (!is.null(SAMPLE_COL) && SAMPLE_COL %in% colnames(obj@meta.data)) {
  meta_cols <- c(meta_cols, SAMPLE_COL)
}
meta <- obj@meta.data[, meta_cols, drop = FALSE]

# Detect UMAP name — see @primitives/seurat_v5_rules.md Rule 5
umap_key <- names(obj@reductions)[grepl("umap", names(obj@reductions), ignore.case = TRUE)][1]
umap_coords <- as.data.frame(Embeddings(obj, umap_key))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
meta <- cbind(meta, umap_coords)

write.csv(meta, file.path(OUTPUT_SCENIC, "metadata.csv"), quote = FALSE)
message("Export complete. Run scenic_00b_prep_loom.py next.")
```

---

### Script 2: `scenic_00b_prep_loom.py` — Python (scenicenv)

Converts cells×genes CSV to genes×cells loom for pySCENIC.

```python
import loompy, numpy as np, pandas as pd
import os

OUTPUT_SCENIC = "output_scenic"

expr_df = pd.read_csv(os.path.join(OUTPUT_SCENIC, "expr.csv"), index_col=0)
# pySCENIC loom: rows=genes, cols=cells
matrix    = expr_df.values.T.astype(np.float32)
row_attrs = {"Gene": np.array(expr_df.columns.tolist())}
col_attrs = {"CellID": np.array(expr_df.index.tolist())}
loompy.create(os.path.join(OUTPUT_SCENIC, "expr.loom"), matrix, row_attrs, col_attrs)
print("Loom created. DELETE expr.csv to free ~1.5 GB:")
print(f"  rm {os.path.join(OUTPUT_SCENIC, 'expr.csv')}")
```

---

### Script 3: `scenic_01_run.sh` — bash (scenicenv)

Run GRN + ctx in sequence. **Do NOT include the aucell step** — it crashes
(see bug table above). AUCell is handled separately in Script 4.

```bash
#!/bin/bash
# ── CONFIG — fill in from brief ──────────────────────────────────────────────
WORKDIR="$(pwd)"               # project_specific: working directory for this run
DBDIR="project_specific"       # REPLACE: database_dir from brief
TF_LIST="project_specific"     # REPLACE: tf_list path from brief
NWORKERS="project_specific"    # REPLACE: n_workers from brief (e.g. 16)

OUTPUT_SCENIC="$WORKDIR/output_scenic"

# Step 1: GRN (arboreto GRNBoost2)
arboreto_with_multiprocessing.py \
    "$OUTPUT_SCENIC/expr.loom" \
    "$TF_LIST" \
    --method grnboost2 \
    --output "$OUTPUT_SCENIC/grn.csv" \
    --num_workers "$NWORKERS" \
    --seed 42

# Step 2: cisTarget (motif enrichment)
pyscenic ctx \
    "$OUTPUT_SCENIC/grn.csv" \
    "$DBDIR/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    "$DBDIR/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    --annotations_fname "$DBDIR/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl" \
    --expression_mtx_fname "$OUTPUT_SCENIC/expr.loom" \
    --output "$OUTPUT_SCENIC/regulons.csv" \
    --mask_dropouts \
    --mode "dask_multiprocessing" \
    --num_workers "$NWORKERS" \
    --min_genes 5
```

**Launch with nohup using the direct Python binary from the brief:**
```bash
# scenic_python_path from brief — fill in from: conda run -n scenicenv which python
SCENIC_PYTHON="project_specific"   # REPLACE: scenic_python_path from brief

nohup "$SCENIC_PYTHON" -c "import subprocess; subprocess.run(['bash', 'scenic_01_run.sh'])" \
  > output_scenic/scenic_run.log 2>&1 &

# Or equivalently — launch the bash script directly with nohup:
nohup bash scenic_01_run.sh > output_scenic/scenic_run.log 2>&1 &
echo "PID: $!"
```

---

### Script 4: `scenic_01b_aucell.py` — Python (scenicenv)

AUCell scoring via Python API, bypassing the CLI binarization bug.

```python
from pyscenic.aucell import aucell
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
import loompy, pandas as pd, numpy as np

OUTPUT_SCENIC = "output_scenic"

# Load expression (cells x genes DataFrame)
with loompy.connect(f"{OUTPUT_SCENIC}/expr.loom", mode="r") as ds:
    expr_df = pd.DataFrame(
        ds[:, :].T,
        index   = ds.ca["CellID"],
        columns = ds.ra["Gene"]
    )

# Load regulons — MUST use df2regulons() to convert DataFrame → GeneSignature list
motifs_df = load_motifs(f"{OUTPUT_SCENIC}/regulons.csv")
regulons  = df2regulons(motifs_df)   # typically ~300 regulons after motif filtering
print(f"Regulons loaded: {len(regulons)}")

# Run AUCell — num_workers from brief
n_workers = 16   # REPLACE: n_workers from brief
auc_mtx = aucell(expr_df, regulons, num_workers=n_workers)
# auc_mtx: cells x regulons, values in [0, 1]

auc_mtx.to_csv(f"{OUTPUT_SCENIC}/auc_scores.csv")
print("AUCell complete. Scores saved to output_scenic/auc_scores.csv")
```

---

## Analysis & Plots (`scenic_02_plots.py`)

Five plots generated from `auc_scores.csv` + `metadata.csv`:

| Output file | Description |
|---|---|
| `01_heatmap_mean_auc.pdf` | Top 30 regulons by variance, z-scored mean AUC, rows = subtypes |
| `02_rss_heatmap.pdf` | Regulon Specificity Score — primary figure for subtype-defining TFs |
| `03_dotplot_rss_regulons.pdf` | Dot size = % cells active (threshold 0.05), color = mean AUC |
| `04_violins_top_regulons.pdf` | AUC distributions for top RSS regulon per subtype |
| `05_umap_tf_activity.pdf` | UMAP panels: subtype colors + one panel per top TF |

```python
import pandas as pd, numpy as np, matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
from scipy.spatial.distance import jensenshannon
import os

OUTPUT_SCENIC = "output_scenic"
LABEL_COL     = "project_specific"   # REPLACE: label column name

auc_df  = pd.read_csv(f"{OUTPUT_SCENIC}/auc_scores.csv", index_col=0)
meta_df = pd.read_csv(f"{OUTPUT_SCENIC}/metadata.csv",   index_col=0)

auc_sub = auc_df.loc[auc_df.index.isin(meta_df.index)].copy()
auc_sub[LABEL_COL] = meta_df.loc[auc_sub.index, LABEL_COL]

regulon_cols = [c for c in auc_sub.columns if c != LABEL_COL]
SUBTYPE_ORDER = sorted(auc_sub[LABEL_COL].unique().tolist())  # REPLACE with desired order

# ── Compute mean AUC and z-score ──────────────────────────────────────────────
mean_auc = auc_sub.groupby(LABEL_COL)[regulon_cols].mean()
z_data   = mean_auc.apply(zscore, axis=0)

# ── Plot 1: Top 30 by variance, grouped by peak subtype ──────────────────────
variance    = mean_auc.var(axis=0).sort_values(ascending=False)
top30       = variance.head(30).index.tolist()
peak_subtype = z_data[top30].idxmax(axis=0)

ordered_cols = []
for st in SUBTYPE_ORDER:
    group = peak_subtype[peak_subtype == st].index.tolist()
    group_sorted = sorted(group, key=lambda r: z_data.loc[st, r], reverse=True)
    ordered_cols.extend(group_sorted)

# Column labels always display with (+ ) suffix
col_labels = [f"{r}(+)" if "(+)" not in r else r for r in ordered_cols]

fig, ax = plt.subplots(figsize=(max(10, len(ordered_cols) * 0.35 + 2),
                                 max(4, len(SUBTYPE_ORDER) * 0.6 + 2)))
sns.heatmap(z_data[ordered_cols].loc[SUBTYPE_ORDER],
            xticklabels=col_labels, yticklabels=SUBTYPE_ORDER,
            cmap="RdBu_r", center=0, ax=ax,
            linewidths=0.3, linecolor="white")
ax.tick_params(axis="x", labelsize=12, rotation=45)
ax.tick_params(axis="y", labelsize=13)

# White dividers between groups
counts   = [sum(peak_subtype[ordered_cols] == st) for st in SUBTYPE_ORDER]
boundary = 0
for c in counts[:-1]:
    boundary += c
    ax.axvline(boundary - 0.5, color="white", linewidth=1.2)

plt.tight_layout()
plt.savefig(f"{OUTPUT_SCENIC}/plots/01_heatmap_mean_auc.pdf", dpi=200)
plt.close()

# ── RSS computation ───────────────────────────────────────────────────────────
def compute_rss(auc_sub, subtypes, regulon_cols, label_col):
    rss = {}
    for st in subtypes:
        mask      = (auc_sub[label_col] == st).astype(float).values
        mask     /= mask.sum()
        scores    = {}
        for reg in regulon_cols:
            auc_vals = auc_sub[reg].values
            auc_prob = auc_vals / (auc_vals.sum() + 1e-10)
            scores[reg] = 1 - jensenshannon(auc_prob, mask)
        rss[st] = scores
    return pd.DataFrame(rss)

rss_df = compute_rss(auc_sub, SUBTYPE_ORDER, regulon_cols, LABEL_COL)
rss_df.to_csv(f"{OUTPUT_SCENIC}/rss_scores.csv")
# Plots 02–05 follow same architecture — see pattern above; omitted for brevity
```

### Font sizes (validated settings)
```python
# Heatmap x-axis (regulon names): fontsize=12
# Heatmap y-axis (subtype names):  fontsize=13
# Violin x-axis (subtype names):   fontsize=11, rotation=45, ha="right"
# Violin panel titles:              fontsize=13
# UMAP panel titles:                fontsize=13
# Use 30 regulons in heatmap (not 50 — 50 is too crowded)
```

---

## TF Candidate Selection CSVs

Always generate two candidate CSVs after AUCell so the user can curate which
regulons appear in each plot.

```python
# Variance-based candidates — for heatmap selection
variance  = mean_auc.var(axis=0).sort_values(ascending=False)
top200    = variance.head(200).index.tolist()
peak_st   = z_data[top200].idxmax(axis=0)
peak_z    = z_data[top200].max(axis=0)

rows = []
for rank, reg in enumerate(top200, 1):
    row = {"rank": rank, "regulon": reg, "variance": variance[reg],
           "peak_subtype": peak_st[reg], "peak_zscore": peak_z[reg]}
    for st in SUBTYPE_ORDER:
        row[f"mean_auc_{st}"] = mean_auc.loc[st, reg] if st in mean_auc.index else None
        row[f"zscore_{st}"]   = z_data.loc[st, reg]   if st in z_data.index   else None
    rows.append(row)
pd.DataFrame(rows).to_csv(f"{OUTPUT_SCENIC}/regulon_candidates.csv", index=False)

# RSS-based candidates — for RSS heatmap selection
rows = []
for st in SUBTYPE_ORDER:
    for rank, (reg, rss_val) in enumerate(rss_df[st].nlargest(20).items(), 1):
        rows.append({"subtype": st, "rank": rank, "regulon": reg,
                     "rss": rss_val,
                     "zscore": float(z_data.loc[st, reg]) if st in z_data.index else None,
                     "mean_auc": mean_auc.loc[st, reg] if st in mean_auc.index else None})
pd.DataFrame(rows).to_csv(f"{OUTPUT_SCENIC}/rss_candidates.csv", index=False)
```

**Applying manual exclusions/inclusions:**
```python
EXCLUDE = {"GENE1", "GENE2"}   # regulons to remove
INCLUDE = ["GENE3", "GENE4"]   # regulons to force-add

regs = [r for r in base_list if r not in EXCLUDE]
for r in INCLUDE:
    if r in mean_auc.columns and r not in regs:
        regs.append(r)
```

---

## TF-Specific Network Plot (`scenic_03_network.py`)

Radial hub-and-spoke network for a single TF's regulon. Uses `regulons.csv` (motif-filtered),
NOT `grn.csv` (correlation-only). See GRN vs Regulons note below.

```python
import matplotlib.pyplot as plt, numpy as np

# ── CONFIG ────────────────────────────────────────────────────────────────────
TARGET_TF = "project_specific"   # REPLACE: TF of interest (e.g. the top RSS regulon)

# Validated figure design parameters
fig, ax = plt.subplots(figsize=(10, 10))
r_scaled = 0.55 + 0.1 * abs(correlations[gene]["r"]) / max_r  # node radius
width    = 1.2 + 3.5 * abs(r_val)                              # edge width
size     = 40  + 120 * abs(r_val)                              # target gene node size
# TF hub node: s=700

# Edge coloring: RdBu_r, normalized to [-0.3, 0.3]
# Red = positively correlated targets; Blue = negatively correlated
# Node coloring: YlOrRd by mean expression in peak subtype
```

**GRN vs Regulons:**

| File | Method | Typical targets | Evidence |
|---|---|---|---|
| `grn.csv` | GRNBoost2 co-expression only | ~3,000+ | Correlation — no causality |
| `regulons.csv` | GRN + cisTarget motif filtering | ~30–50 | Co-expression AND binding motif |

Always use `regulons.csv` for the network plot.

---

## Output Directory Structure

```
output_scenic/
├── expr.loom                # pySCENIC input (keep — needed for AUCell)
├── metadata.csv             # cell metadata with subtype labels + UMAP coords
├── grn.csv                  # GRN output (~700k TF-target pairs)
├── regulons.csv             # cisTarget output (~3500 rows before df2regulons)
├── auc_scores.csv           # AUCell output: cells x regulons [0, 1]
├── mean_auc_per_subtype.csv # summary table
├── rss_scores.csv           # RSS scores per regulon per subtype
├── regulon_candidates.csv   # top 200 by variance — for user curation
├── rss_candidates.csv       # top 20/subtype by RSS — for user curation
└── plots/
    ├── 01_heatmap_mean_auc.pdf
    ├── 02_rss_heatmap.pdf
    ├── 03_dotplot_rss_regulons.pdf
    ├── 04_violins_top_regulons.pdf
    ├── 05_umap_tf_activity.pdf     # rasterized=False (vector)
    └── 06_{TF}_network.pdf         # one file per TF of interest
```

**Do not keep `expr.csv`** — it is a 1.5 GB intermediate that is fully
redundant once the loom file exists.

---

## UMAP Plot Notes

- Always use `rasterized=False` in scatter calls for publication-quality vector PDFs
- Point size `s=0.8` works well for ~30k cells
- UMAPs must use coordinates exported from Seurat (in `metadata.csv`) —
  do not re-compute UMAP in Python

---

## Project-Specific Values (Stage for Phase 4 examples/)

- `examples/nkxspleen_pyscenic.md`:
  - NKXSpleen project paths and validated configuration
  - SUBTYPE_ORDER for spleen EC subtypes
  - Validated 317 regulons from KidneyNew run context
- `examples/humanfat_pyscenic.md`:
  - HumanFat EC subtype configuration, 6 subtypes, 320 regulons
  - ec_subtype column mapping
