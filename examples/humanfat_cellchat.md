---
# Example: HumanFat CellChat Ligand-Receptor Analysis
# Status: full-build
# Validated: 2026-03-30 (full analysis); 2026-04-02 (bar plots CAT_COLORS)
---

# Example: HumanFat CellChat Ligand-Receptor Analysis

## Instantiates
- @modules/cellchat.md

## Project context
- Project: HumanFat
- Validated: 2026-03-30
- Status: full-build

## Brief block

```yaml
downstream_analyses:
  cellchat:
    label_col: mylabel              # EC subset label column
    group_col: tissue_type          # tissue comparison column (enables Script 6)
                                    # Script 6 compares LiposuctionFat vs BreastFat
    organism: "human"
    n_workers: 4                    # plan("sequential") is MANDATORY — see critical note
    source_celltypes:               # EC subtypes as signal senders
      - "AEC"
      - "CapEC"
      - "CapEC2"
      - "VenEC1"
      - "VenEC2"
      - "VenEC3"
    target_celltypes:               # all cell types in the object as potential receivers
      - "AEC"
      - "CapEC"
      - "CapEC2"
      - "VenEC1"
      - "VenEC2"
      - "VenEC3"
      - "Stromal Cell / ASPC"
      - "Macrophage"
      - "Naive Monocyte"
      - "Polarized Monocyte"
      - "Prolif. Monocyte"
      - "Dendritic Cell"
      - "pDC"
      - "Mast Cell"
      - "Neutrophil"
      - "B cell"
      - "Plasma Cell"
      - "CD3+ T cell"
      - "CD4+ T cell"
      - "CD8+ T cell"
      - "Epithelial"

context_overrides:
  palettes:
    cell_colors:
      # EC subtypes
      "AEC":                        "#F4433C"
      "CapEC":                      "#FF9800"
      "CapEC2":                     "#FFEB3B"
      "VenEC1":                     "#4CAF50"
      "VenEC2":                     "#2196F3"
      "VenEC3":                     "#9C27B0"
      # Non-EC cell types
      "Stromal Cell / ASPC":        "#2D6A4F"
      "Macrophage":                 "#3D405B"
      "Naive Monocyte":             "#4A7C59"
      "Polarized Monocyte":         "#6B8F71"
      "Prolif. Monocyte":           "#7B6888"
      "Dendritic Cell":             "#B8860B"
      "pDC":                        "#DAA520"
      "Mast Cell":                  "#81B29A"
      "Neutrophil":                 "#9E9E9E"
      "B cell":                     "#B0C4DE"
      "Plasma Cell":                "#D4A5A5"
      "CD3+ T cell":                "#4ECDC4"
      "CD4+ T cell":                "#45B7D1"
      "CD8+ T cell":                "#2196A8"
      "Epithelial":                 "#A0785A"
```

## Script 6 tissue comparison configuration

Script 6 runs only when `group_col` is set. For HumanFat, the primary comparison is
LiposuctionFat vs BreastFat. The 4-tissue variant includes Subcutaneous Fat and Visceral Fat.

```r
# Primary Script 6 comparison
cc_lipo <- extract_tissue(so, group_col = "tissue_type", group_val = "Liposuction Fat")
cc_breast <- extract_tissue(so, group_col = "tissue_type", group_val = "Breast Fat")
compare_pair(cc_lipo, cc_breast, label1 = "LiposuctionFat", label2 = "BreastFat")
```

## Pathway categorization (Script 1 → Script 2)

After Script 1 runs inference, review the pathway list before authoring PATHWAY_CATEGORIES.
The HumanFat run identified categories in this structure (validated 2026-03-30):

```r
PATHWAY_CATEGORIES <- list(
  "CAT_COLORS assignment validated 2026-04-02" = list(
    # Category names and colors validated in the bar plot run
  )
)
# CAT_COLORS:
CAT_COLORS <- c(
  "Angiocrine"     = "#FF9800",
  "Metabolic"      = "#4CAF50",
  "Immune_Traffic" = "#E53935",
  "ECM"            = "#3D70BE",
  "Neuronal"       = "#9C27B0",
  "Adhesion"       = "#00BCD4",
  "Other"          = "#9E9E9E"
)
# Note: ECM category excluded from Script 3 stacked bubble (STACKED_CATS omits ECM)
# Rationale: ECM pathways have very high chord probability; their visual weight in
# the stacked bubble overwhelms the more biologically interesting angiocrine signals
```

## Validation notes

- Validated on HumanFat EC subset (6 EC subtypes + all non-EC cell types)
- Script 6 primary comparison: LiposuctionFat vs BreastFat (adipose tissue depots)
- Script 5 circos: 4 standard plots; two inference objects required (full vs collapsed ASPC16 sub-population)
- Outputs written to `output/cellchat/`

## Known issues / quirks

- CRITICAL: `plan("sequential")` at the very top of every CellChat script, before library loads.
  `plan("multisession")` crashes chord diagram rendering. This is a CellChat v2 known bug.
- RibHighEC: excluded from SOURCE_TYPES (not biologically interpretable as sender);
  may be present in target types via the whole-object label column — verify before Script 1
- Label column for CellChat is `mylabel` (EC subset) — this is different from `celltype`
  (whole-object label). Make sure the Seurat object loaded for CellChat has `mylabel` populated
- HumanFat CellChat object was validated on the EC subset object, NOT the whole-object RDS
- Script 3 stacked bubbles: separate ASPC16 sub-population plot required; `mac_types`,
  `aspc_types`, `aspc16` are subsets of target_celltypes defined per-run
