---
# CellChat Ligand-Receptor Analysis — worked example (generic)
#
# This example shows how to invoke the cellchat module on a typical
# dataset. Placeholders in angle brackets and values marked with
# "# REPLACE:" must be customized for your project.
#
# For the module reference itself, see: modules/cellchat.md
# Status: full-build
---

# CellChat Ligand-Receptor Analysis — worked example (generic)

## Instantiates
- @modules/cellchat.md

## Brief block

```yaml
downstream_analyses:
  cellchat:
    label_col: mylabel              # REPLACE: label column for your cell type annotations
    group_col: tissue_type          # tissue comparison column (enables Script 6)
                                    # Script 6 compares two groups within this column
    organism: "human"
    n_workers: 4                    # plan("sequential") is MANDATORY — see critical note
    source_celltypes:
      # REPLACE: cell type names from your dataset; the values below are illustrative examples
      - "AEC"
      - "CapEC"
      - "CapEC2"
      - "VenEC1"
      - "VenEC2"
      - "VenEC3"
    target_celltypes:
      # REPLACE: cell type names from your dataset; the values below are illustrative examples
      - "AEC"
      - "CapEC"
      - "CapEC2"
      - "VenEC1"
      - "VenEC2"
      - "VenEC3"
      - "Stromal Cell"
      - "Macrophage"
      - "Monocyte"
      - "Dendritic Cell"
      - "T cell"
      - "B cell"
      - "Epithelial"

context_overrides:
  palettes:
    cell_colors:
      # REPLACE: assign colors for your cell types; see context/color_palettes.md
      "AEC":        "#REPLACE_HEX"
      "CapEC":      "#REPLACE_HEX"
      "CapEC2":     "#REPLACE_HEX"
      "VenEC1":     "#REPLACE_HEX"
      "VenEC2":     "#REPLACE_HEX"
      "VenEC3":     "#REPLACE_HEX"
      "Stromal Cell": "#REPLACE_HEX"
      "Macrophage": "#REPLACE_HEX"
      "T cell":     "#REPLACE_HEX"
      "B cell":     "#REPLACE_HEX"
```

## Script 6 tissue comparison configuration

Script 6 runs only when `group_col` is set. Set the `group_val` arguments to the two
tissue or condition groups you want to compare (REPLACE with your actual group labels).

```r
# REPLACE: substitute your group labels from the group_col column
cc_groupA <- extract_tissue(so, group_col = "tissue_type", group_val = "TissueTypeA")
cc_groupB <- extract_tissue(so, group_col = "tissue_type", group_val = "TissueTypeB")
compare_pair(cc_groupA, cc_groupB, label1 = "TissueTypeA", label2 = "TissueTypeB")
```

## Pathway categorization (Script 1 → Script 2)

After Script 1 runs inference, review the pathway list before authoring PATHWAY_CATEGORIES.
REPLACE the category names and pathway assignments with those appropriate for your biology.

```r
# REPLACE: define categories relevant to your biological question
# Example category structure for vascular/stromal biology:
PATHWAY_CATEGORIES <- list(
  "Angiocrine"     = c("VEGF", "ANGPT", "NOTCH"),    # example pathways — REPLACE
  "Metabolic"      = c("ADIPONECTIN"),               # example pathways — REPLACE
  "Immune_Traffic" = c("CXCL", "CCL"),               # example pathways — REPLACE
  "ECM"            = c("COLLAGEN", "FN1"),            # example pathways — REPLACE
  "Other"          = c()                              # catch-all for unlisted pathways
)
# REPLACE: assign colors from context/color_palettes.md for each category
CAT_COLORS <- c(
  "Angiocrine"     = "#REPLACE_HEX",
  "Metabolic"      = "#REPLACE_HEX",
  "Immune_Traffic" = "#REPLACE_HEX",
  "ECM"            = "#REPLACE_HEX",
  "Other"          = "#REPLACE_HEX"
)
# Note: consider whether any high-probability categories should be excluded from
# Script 3 stacked bubble to avoid overwhelming lower-signal categories visually
```

## Validation notes

- Example validated on a single-cell-type subset object (one cell type as senders, all cell types as receivers)
- Script 6: two-group comparison using the tissue_type column (TissueTypeA vs TissueTypeB)
- Script 5 circos: 4 standard plots; two CellChat inference objects may be required if a
  sub-population needs separate treatment
- Outputs written to `<YOUR_OUTPUT_DIR>/cellchat/`

## Known issues / quirks

- CRITICAL: `plan("sequential")` at the very top of every CellChat script, before library loads.
  `plan("multisession")` crashes chord diagram rendering. This is a CellChat v2 known bug.
- Confirm which label column to use before Script 1: the subset label column may differ from
  the whole-object annotation column — verify the correct column is populated in your object
- Build the CellChat object from the cell-type-subset object, not the whole-object RDS,
  unless all cell types are included in both source and target lists
- Script 3 stacked bubbles: sub-population plots may require separate CellChat objects;
  define sub-population membership vectors from target_celltypes per-run
