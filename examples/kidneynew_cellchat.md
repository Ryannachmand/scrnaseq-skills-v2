---
# Example: KidneyNew CellChat Ligand-Receptor Analysis
# Status: full-build
# Validated: 2026-04-10
---

# Example: KidneyNew CellChat Ligand-Receptor Analysis

## Instantiates
- @modules/cellchat.md

## Project context
- Project: KidneyNew
- Validated: 2026-04-10 (TFEC Expression Atlas run)
- Status: full-build

## Brief block

```yaml
downstream_analyses:
  cellchat:
    label_col: ec_subtype_final     # EC subset label column (post-hoc annotation)
                                    # use ec_subtype_final after EC-PTC was added;
                                    # early-phase scripts may use ec_subtype (pre-hoc)
    group_col: null                 # no tissue comparison in KidneyNew; Script 6 skipped
    organism: "human"
    n_workers: 4                    # plan("sequential") is MANDATORY — see critical note
    source_celltypes:               # kidney EC subtypes as signal senders
      - "EC-GC"                     # glomerular capillary
      - "EC-AEA"                    # ascending thin limb
      - "EC-AVR"                    # ascending vasa recta
      - "EC-DVR"                    # descending vasa recta
      - "EC-LYM"                    # lymphatic
      - "EC-PTC"                    # peritubular capillary (post-hoc annotation)
    target_celltypes:               # EC subtypes + stromal/epithelial/MoMF as receivers
      - "EC-GC"
      - "EC-AEA"
      - "EC-AVR"
      - "EC-DVR"
      - "EC-LYM"
      - "EC-PTC"
      - "Stroma"                    # stromal cell types (confirm exact label from object)
      - "Epithelial"                # epithelial cell types (confirm exact label from object)
      - "MoMF"                      # monocyte-derived macrophage/macrophage-like cells

context_overrides:
  palettes:
    cell_colors:
      # Kidney EC subtypes
      "EC-GC":   "#10B8F5"
      "EC-AEA":  "#ECA7A0"
      "EC-AVR":  "#A9B100"
      "EC-DVR":  "#36E2A6"
      "EC-LYM":  "#EFADFF"
      "EC-PTC":  "#4B6584"
      # Non-EC cell types — confirm hex values from project records
      "Stroma":     "#8FBC8F"       # TODO: confirm from KidneyNew project CLAUDE.md
      "Epithelial": "#DDA0DD"       # TODO: confirm from KidneyNew project CLAUDE.md
      "MoMF":       "#F0E68C"       # TODO: confirm from KidneyNew project CLAUDE.md
```

## Assay switching (critical)

KidneyNew uses SCT-normalized objects. CellChat requires RNA assay:

```r
# Before Script 1 — switch to RNA for CellChat
DefaultAssay(so) <- "RNA"
so <- JoinLayers(so)  # MANDATORY before GetAssayData — see @primitives/seurat_v5_rules.md Rule 1
# ... run CellChat ...
# After Script 1 — restore SCT for any visualization
DefaultAssay(so) <- "SCT"
```

## Pathway categorization (Script 1 → Script 2)

KidneyNew CellChat identified three primary pathway categories:

```r
PATHWAY_CATEGORIES <- list(
  Chemokine_Adhesion = list(
    pathways   = c("CXCL", "CCL", "ICAM", "PECAM1", "VCAM"),
    chord_paths = c("CXCL", "CCL", "ICAM"),  # APP excluded from chord (see note)
    grad_high  = "#7B3FA0"
  ),
  Angiocrine = list(
    pathways   = c("VEGF", "ANGPT", "NOTCH", "WNT", "BMP"),
    chord_paths = c("VEGF", "ANGPT", "NOTCH"),
    grad_high  = "#E53935"
  ),
  ECM = list(
    pathways   = c("COLLAGEN", "FN1", "LAMININ", "MK", "HSPG"),
    chord_paths = c("COLLAGEN", "FN1", "LAMININ"),
    grad_high  = "#3D70BE"
  )
)
# APP pathway note: APP (amyloid precursor protein) typically dominates total interaction
# probability (~5x next strongest pathway). Excluded from chord_paths to prevent visual
# compression of biologically informative pathways. Still included in pathway list for
# statistical summaries.
```

## Validation notes

- Validated on KidneyNew EC subset (7 subtypes including post-hoc EC-PTC)
- EC-PTC was added post-hoc via subclustering (seurat_cluster 3, res=0.2, subcluster 3_1, 991 cells)
- The existing CellChat object (pre-EC-PTC annotation) does NOT include EC-PTC as a cell type;
  if re-running CellChat after EC-PTC annotation, use `ec_subtype_final` (includes EC-PTC)
- KidneyNew NA handling: NA values in mylabel stored as string "NA", not R native NA;
  check both `is.na(x)` and `x == "NA"` when filtering
- Outputs written to `output/cellchat/`

## Known issues / quirks

- CRITICAL: `plan("sequential")` at the very top of every CellChat script, before library loads.
  `plan("multisession")` crashes chord diagram rendering. This is a CellChat v2 known bug.
- KidneyNew uses SCT for visualization and RNA for CellChat — ALWAYS switch assay before
  running any CellChat script; ALWAYS restore SCT after
- EC-Prolif (proliferating EC subtype) shares soft purple (`#EFADFF`) with EC-LYM in some plots;
  distinguish by label if both are present in the same visualization
- Script 6 is skipped (`group_col: null`) because KidneyNew has no tissue comparison groups
- The non-EC cell types (Stroma, Epithelial, MoMF) hex colors above are placeholders;
  confirm exact values from KidneyNew project CLAUDE.md before running Scripts 2-5
