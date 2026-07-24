---
# Example: HumanFat pySCENIC TF Regulon Analysis
# Status: full-build
# Validated: 2026-03-11
---

# Example: HumanFat pySCENIC TF Regulon Analysis

## Instantiates
- @modules/pyscenic_regulons.md

## Project context
- Project: HumanFat
- Validated: 2026-03-11 (HumanFat run; pySCENIC also validated separately on another project)
- Status: full-build

## Brief block

```yaml
downstream_analyses:
  pyscenic:
    enabled: true
    label_col: mylabel              # EC subset label column
    scenic_python_path: "<YOUR_SCENICENV_PYTHON>"  # find with: conda run -n scenicenv which python
                                    # machine-specific path — confirm before running
    database_dir: "<YOUR_PYSCENIC_DATABASE_DIR>"  # feather databases downloadable from aertslab GitHub
                                    # path to pySCENIC database files (rankings + motifs)
                                    # machine-specific — confirm before running
    tf_list: "<YOUR_TF_LIST_PATH>"  # allTFs_hg38.txt available from aertslab pySCENIC resources
                                    # human TF list for grn inference
                                    # machine-specific — confirm before running
    n_workers: 16                   # parallelization workers for pyscenic grn step
    min_pct_cells: 0.05             # minimum fraction of cells expressing a TF
    min_mean_expr: 0.1              # minimum mean expression cutoff
```

## Validation notes

- Dataset: HumanFat EC subset (6 EC subtypes: AEC, CapEC, CapEC2, VenEC1, VenEC2, VenEC3)
- RibHighEC excluded before running pySCENIC
- 320 regulons identified in the validated HumanFat run
- Label column: `mylabel` (the EC subset label column)
- Outputs written to `output/pyscenic/`
- Three-step CLI pipeline:
  1. `scenic_01a_grn.py` — GRN inference (most compute-intensive; use nohup)
  2. `scenic_01b_aucell.py` — AUCell regulon scoring (Python API; CLI is broken due to np.float removal)
  3. `scenic_01c_regulons.R` — AUCell → Seurat metadata, UMAP overlay, heatmap
- `pyscenic aucell` CLI is broken (numpy `np.float` removal); always use Python API via `scenic_01b_aucell.py`

## Known issues / quirks

- All three paths (scenic_python_path, database_dir, tf_list) are machine-specific; confirm
  existence before running — the module will fail with a clear error if any path is missing
- n_workers = 16 is validated on the lab machine; adjust for available CPU cores
- HumanFat used RNA assay (not SCT) for SCENIC input — JoinLayers must be applied before
  `GetAssayData()` extraction
- pySCENIC outputs use the subtype ordering from SUBTYPE_ORDER in the Python scripts;
  for HumanFat, ordered as: AEC, CapEC, CapEC2, VenEC1, VenEC2, VenEC3 (excluding RibHighEC)
- For the comparable run on a kidney dataset, see the sibling pySCENIC example
