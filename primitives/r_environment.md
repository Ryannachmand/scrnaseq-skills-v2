---
# R Environment Rules — v2
# Migrated from ~/claude-skills/shared/r_environment.md with F1 bug fix applied.
# F1 FIX: Added --no-capture-output flag. Without this flag, conda run hangs
#         after the R script completes. This was missing in v1.
---

# R Environment Rules

Applies to all pipelines on this machine.

---

## Always Use r-env

Run all R scripts with:
```bash
conda run --no-capture-output -n r-env Rscript <script.R> > <logfile> 2>&1
```

**CRITICAL:** The `--no-capture-output` flag is required. Without it, the conda run
command hangs indefinitely after the R script completes. Never omit this flag.

Never use the system `Rscript`.

**Why:** System R (4.3.1) has ABI incompatibilities with the `sp` package,
causing Seurat and Harmony to fail to load. The `r-env` conda environment
has R 4.5.1 with all required packages installed.

---

## Packages Available in r-env

- Seurat (v5), harmony, ggplot2, dplyr, tidyr, patchwork, scales
- monocle3, SingleCellExperiment
- cowplot, ggrepel
- AUCell (for gene set scoring)
- ComplexHeatmap, circlize
- CellChat (v2)
- clusterProfiler, org.Hs.eg.db
- matrixStats, data.table, readxl
- pheatmap
- Note: SeuratWrappers is NOT installed — see modules/trajectory_monocle3.md for workaround
- Note: tibble is available but avoid `tibble::column_to_rownames()` — see seurat_v5_rules.md Rule 6

A separate `scenicenv` Python environment handles pySCENIC (not R-based).

---

## Script Structure Convention

Every generated R script must be:

- **Self-contained:** load all libraries at the top; define all paths as variables at the top
- **Verbose:** use `cat()` and `message()` liberally for progress visibility
- **Memory-safe:** call `gc()` after memory-intensive steps (ScaleData, RunPCA, RunHarmony)
- **Checkpointed:** save intermediate RDS at every major stage — never overwrite; use versioned
  suffixes (`_raw`, `_filtered`, `_harmony`, `_annotated`)
- **PDF-safe:** always save plots with `ggsave(..., device = "pdf", useDingbats = FALSE)`;
  never use `pdf()` / `dev.off()` patterns
