---
# Example: BoneMarrowStroma Label Harmonization
# Status: full-build
# Validated: BoneMarrowStroma run
---

# Example: BoneMarrowStroma Label Harmonization

## Instantiates
- @modules/label_harmonization.md

## Project context
- Project: BoneMarrowStroma
- Validated: BoneMarrowStroma run
- Status: full-build

## Brief block

```yaml
label_harmonization:
  label_col:                        # REQUIRED: per-dataset source label columns
    Vertebrae:  "cell_type"         # label column in Leimkuler (Vertebrae) dataset
    IliacCrest: "cell_type"         # label column in Wang (IliacCrest) dataset
    FemoralHead: "cell_type"        # label column in Li (FemoralHead) dataset
                                    # TODO: confirm actual column names from each dataset

  broad_label: "Mesenchymal"        # REQUIRED: fallback label for low-confidence cells
                                    # Mesenchymal is the appropriate broad category for
                                    # stromal cells that do not pass score_high threshold
                                    # Rationale: low-confidence stromal cells are still
                                    # mesenchymal lineage; "Unknown" would lose biological info

  label_transfer_reference:         # datasets used as reference in Seurat label transfer
    - "Vertebrae"
    - "FemoralHead"
    # IliacCrest (Wang) used as query; Vertebrae + FemoralHead as reference

  label_transfer_score_high: 0.75   # lab default; cells above this score get reference label
  label_transfer_score_low: 0.40    # cells below this score get broad_label ("Mesenchymal")

  known_label_mappings:             # manual YAML mappings for well-characterized labels
    # Sort-strategy labels (cells sorted by surface markers before sequencing)
    "CD14 Isolated Stroma":  "CD14 Isolated Stroma"   # preserve as-is; sort strategy
    "PE Isolated Stroma":    "PE Isolated Stroma"      # preserve as-is; sort strategy
    # Passaged/cultured cells
    "Passaged Stroma":       "Passaged Stroma"         # passaged stromal cells (in vitro)
    # Derived display labels
    "Adipo-MSC":             "Cultured Adipo-like MSC" # Li dataset "Adipo-MSC" → display as
                                                        # "Cultured Adipo-like MSC" in plots
                                                        # Rationale: these are cultured cells
                                                        # from the FemoralHead dataset; the
                                                        # "Cultured" prefix clarifies origin
```

## Label transfer result distribution (validated)

From the validated BoneMarrowStroma Li-dataset transfer:

```r
# Li dataset label transfer result distribution (Adipo-MSC example):
# Adipo-MSC: 1823 cells transferred with score >= 0.75
# Osteo-MSC: 912  cells transferred with score >= 0.75
# Unknown:   75   cells with score < 0.40 → assigned broad_label "Mesenchymal"
# (cells between 0.40 and 0.75 assigned broad_label for safety)

# The derived display label mapping:
# "Adipo-MSC" (reference label) → "Cultured Adipo-like MSC" (display label in plots)
# This transformation is applied in the known_label_mappings block above
```

## Validation notes

- Reference datasets: Vertebrae (Leimkuler) + FemoralHead (Li) used as reference
- Query dataset: IliacCrest (Wang) transferred labels against the reference
- Three source datasets produce a unified unified_label column
- Sort-strategy labels (CD14 Isolated Stroma, PE Isolated Stroma) are kept as-is because
  they represent a cell isolation protocol that is biologically meaningful
- Passaged Stroma cells are explicitly labeled as cultured/passaged to distinguish from
  fresh isolates
- Broad label "Mesenchymal" applied to cells below score_low threshold (< 0.40)

## Known issues / quirks

- The `stop()` in label_harmonization.md fires if `broad_label` is not set; always provide
  `broad_label: "Mesenchymal"` in the brief for BoneMarrowStroma
- FemoralHead UMAP: uses non-standard UMAP reduction name "UMAP_dim30" — see
  @primitives/seurat_v5_rules.md Rule 5 for dynamic UMAP detection; label transfer
  itself does not use the UMAP reduction, so this only affects visualization steps
- JoinLayers must be called on the reference merged object before TransferData (Rule 1)
- The `%||%` operator used in the module is from rlang (loaded via Seurat); confirm
  Seurat is loaded before running the module script
