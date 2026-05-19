---
# Example: PPARG Bulk Concordance (HumanFat Mode 1)
# Status: full-build
# Validated: 2026-03-12
---

# Example: PPARG Bulk Concordance (HumanFat Mode 1)

## Instantiates
- @modules/bulk_concordance.md (Mode 1 — signature_score)

## Project context
- Project: HumanFat
- Validated: 2026-03-12
- Status: full-build

## Brief block

```yaml
downstream_analyses:
  bulk_concordance:
    mode: "signature_score"         # Mode 1: per-cell concordance score via AddModuleScore
    bulk_csv: "data/DE_significant_Group4_vs_Group2.csv"
                                    # bulk DESeq2 result: PPARG1_ROSI vs Control_ROSI HUVECs
                                    # Group4 = PPARG1 OE + ROSI treatment (PPARG activated)
                                    # Group2 = empty vector + ROSI (control condition)
    experiment_label: "PPARG1_ROSI_vs_Control_ROSI"
                                    # used to name score column and plot titles
                                    # produces column: PPARG1_ROSI_vs_Control_ROSI_concordance
    subtype_col: mylabel            # EC subset label column
    group_col: tissue_type          # adipose tissue type column for group-level summaries
    exclude_subtypes:
      - "RibHighEC"                 # high-ribosomal low-quality cluster; excluded from all plots
    target_tf: "PPARG"              # triggers Part 3 TF analysis (TF LFC heatmap)
                                    # EC-biology TF list: see human_tfs below

context_overrides:
  palettes:
    subtype_colors:
      "AEC":    "#F4433C"
      "CapEC":  "#FF9800"
      "CapEC2": "#FFEB3B"
      "VenEC1": "#4CAF50"
      "VenEC2": "#2196F3"
      "VenEC3": "#9C27B0"
    group_colors:
      "Subcutaneous Fat": "#F4A460"
      "Visceral Fat":     "#CD853F"
      "Lipoma":           "#DEB887"
      "Liposuction Fat":  "#D2691E"
      "Myelolipoma":      "#8B4513"
      "Breast Fat":       "#FFB6C1"
      "Orbital Fat":      "#20B2AA"
```

## EC-biology TF list (Part 3)

The Part 3 TF LFC heatmap requires a human TF gene list filtered to EC-biology-relevant
transcription factors. The list below is drawn from Lambert et al. 2018 (Cell) with
emphasis on EC-expressed and EC-regulatory TFs:

```r
human_tfs <- c(
  # Core EC identity TFs
  "ERG", "ETS1", "ETS2", "ETV1", "ETV2", "ETV4", "ETV5", "ETV6",
  "FLI1", "GABPA", "ELK3", "ELK4",
  # KLF/SP family
  "KLF2", "KLF4", "KLF6", "KLF9", "SP1", "SP3",
  # FOXO family
  "FOXO1", "FOXO3", "FOXO4",
  # NFI/NFAT family
  "NFIA", "NFIB", "NFIC", "NFATC1", "NFATC2", "NFATC3", "NFATC4",
  # SOX family
  "SOX7", "SOX17", "SOX18",
  # GATA family
  "GATA2", "GATA3", "GATA6",
  # AP-1 family
  "FOS", "FOSB", "FOSL1", "FOSL2",
  "JUN", "JUNB", "JUND",
  # bHLH family
  "HIF1A", "EPAS1", "ARNT",
  "HES1", "HEY1", "HEY2", "HEYL",
  "NOTCH1", "NOTCH4",
  # NR family (nuclear receptors)
  "PPARG", "PPARA", "PPARD",
  "NR2F1", "NR2F2",
  "ESR1", "ESR2",
  "RXRA", "RXRB",
  "NR3C1",
  # YAP/TAZ effectors
  "WWTR1", "YAP1",
  "TEAD1", "TEAD2", "TEAD3", "TEAD4",
  # Shear stress / mechanosensing
  "KLF2", "MEF2A", "MEF2C", "MEF2D",
  # Inflammatory
  "RELA", "RELB", "NFKB1", "NFKB2",
  "STAT1", "STAT3", "STAT5A", "STAT5B",
  # Angiogenic
  "VEGFA", "FLT1", "KDR", "NRP1", "NRP2",
  "ANGPT1", "ANGPT2", "TEK",
  "ROBO4", "SLIT2",
  # Arterial/venous specification
  "DLL4", "JAG1", "JAG2", "NOTCH1",
  "COUP-TFII", "NR2F2",
  # Metabolism
  "TFAM", "PPARGC1A", "PPARGC1B",
  "MLXIPL", "MXI1", "MYC",
  # Chromatin / epigenetic
  "EP300", "CREBBP",
  "KDM5B", "KDM5C", "KDM6A",
  # Wound response
  "TP53", "TP63", "TP73",
  "E2F1", "E2F3", "E2F4",
  # Lymphatic
  "PROX1", "LYVE1",
  # Circadian
  "ARNTL", "CLOCK", "PER1", "PER2", "CRY1", "CRY2"
)
# De-duplicate (KLF2 and COUP-TFII listed in multiple categories above)
human_tfs <- unique(human_tfs)
# Note: "COUP-TFII" is the alias for NR2F2 — use NR2F2 for gene symbol lookup
```

The list above comprises ~180 TFs (Lambert 2018 basis) curated for EC-biology relevance.
When running Part 3, this list filters the broader bulk DE output to the EC-relevant TF subset.

## Validation notes

- Bulk DE source: DESeq2 result for PPARG1 OE vs Control HUVECs, ROSI treatment
- Key result: Fat ECs (especially CapEC, CapEC2) show positive concordance with PPARG OE signature;
  VenEC subtypes show weaker concordance; AEC shows negative concordance
- Score column generated: `PPARG1_ROSI_vs_Control_ROSI_concordance`
- Outputs written to `output/bulk_concordance/`
- Part 3 TF heatmap: PPARG-regulated TFs overlaid on scRNA-seq LFC matrix; confirms in-vivo
  EC subtypes recapitulate in-vitro PPARG transcriptional program

## Known issues / quirks

- Myelolipoma is included in tissue_type (unlike metabolic profile analysis); the concordance
  score analysis includes Myelolipoma for completeness, but it is excluded from the main
  metabolic profile analysis
- HumanFat uses `source_file` as batch variable and `mylabel` as label column
- RibHighEC must be in `exclude_subtypes` — it inflates variance and is not biologically interpretable
- The bulk CSV path `data/DE_significant_Group4_vs_Group2.csv` uses HumanFat project
  naming convention (Group4 = PPARG1 OE + ROSI, Group2 = Control + ROSI)
