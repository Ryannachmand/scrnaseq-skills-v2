---
# Example: HumanFat Metabolic Profile
# Status: full-build
# Validated: 2026-03-05
---

# Example: HumanFat Metabolic Profile

## Instantiates
- @modules/metabolic_profile.md

## Project context
- Project: HumanFat
- Validated: 2026-03-05
- Status: full-build

## Brief block

```yaml
downstream_analyses:
  metabolic_profile:
    enabled: true
    label_col: mylabel            # EC subset label column (AEC, CapEC, CapEC2, VenEC1, VenEC2, VenEC3, RibHighEC)
    group_col: tissue_type        # tissue comparison column
    gene_sets:
      Glycolysis:
        - HK1
        - HK2
        - GPI
        - PFKL
        - PFKM
        - PFKP
        - ALDOA
        - ALDOB
        - ALDOC
        - TPI1
        - GAPDH
        - PGK1
        - PGAM1
        - ENO1
        - ENO2
        - PKM
        - PKLR
        - LDHA
        - LDHB
      Beta_Oxidation:
        - CPT1A
        - CPT1B
        - CPT2
        - ACADVL
        - ACADL
        - ACADM
        - ACADS
        - ACADSB
        - HADHA
        - HADHB
        - ECHS1
        - HADH
        - ACAA2
        - ACAA1
        - ACSL1
        - ACSL4
        - ACSL5
        - ACSL6
        - SLC25A20
      TCA_Cycle:
        - CS
        - ACO1
        - ACO2
        - IDH1
        - IDH2
        - IDH3A
        - IDH3B
        - IDH3G
        - OGDH
        - OGDHL
        - DLST
        - DLD
        - SUCLA2
        - SUCLG1
        - SUCLG2
        - SDHA
        - SDHB
        - SDHC
        - SDHD
        - FH
        - MDH1
        - MDH2
        - PC
      OXPHOS:
        - NDUFA1
        - NDUFA2
        - NDUFA3
        - NDUFA4
        - NDUFA5
        - NDUFA6
        - NDUFA7
        - NDUFA8
        - NDUFA9
        - NDUFA10
        - NDUFA11
        - NDUFA12
        - NDUFA13
        - NDUFB1
        - NDUFB2
        - NDUFB3
        - NDUFB4
        - NDUFB5
        - NDUFB6
        - NDUFB7
        - NDUFB8
        - NDUFB9
        - NDUFB10
        - NDUFB11
        - NDUFC1
        - NDUFC2
        - NDUFS1
        - NDUFS2
        - NDUFS3
        - NDUFS4
        - NDUFS5
        - NDUFS6
        - NDUFS7
        - NDUFS8
        - NDUFV1
        - NDUFV2
        - NDUFV3
        - UQCRC1
        - UQCRC2
        - UQCRFS1
        - UQCRH
        - UQCRQ
        - COX4I1
        - COX5A
        - COX5B
        - COX6A1
        - COX6B1
        - COX7A1
        - COX7A2
        - COX7B
        - COX8A
        - ATP5A1
        - ATP5B
        - ATP5C1
        - ATP5D
        - ATP5E
        - ATP5F1
        - ATP5G1
        - ATP5G2
        - ATP5G3
        - ATP5H
        - ATP5I
        - ATP5J
        - ATP5J2
        - ATP5L
        - ATP5O
      FA_Synthesis:
        - ACLY
        - ACACA
        - ACACB
        - FASN
        - ELOVL1
        - ELOVL2
        - ELOVL3
        - ELOVL4
        - ELOVL5
        - ELOVL6
        - ELOVL7
        - SCD
        - SCD5
        - FADS1
        - FADS2
        - GPAM
        - GPAT3
        - GPAT4
        - AGPAT1
        - AGPAT2
        - AGPAT3

context_overrides:
  palettes:
    subtype_colors:
      "AEC":    "#F4433C"
      "CapEC":  "#FF9800"
      "CapEC2": "#FFEB3B"
      "VenEC1": "#4CAF50"
      "VenEC2": "#2196F3"
      "VenEC3": "#9C27B0"
```

## Pre-processing note

Before running the metabolic profile module, apply the Myelolipoma exclusion:

```r
# Exclude Myelolipoma from metabolic analysis
# Rationale: Myelolipoma is a rare benign tumor (hematopoietic + adipose tissue);
# its metabolism is not representative of adipose EC biology and inflates variance.
ec <- ec[, ec$tissue_type != "Myelolipoma"]
```

This filter is applied before passing the object to `run_aucell()` / `metabolic_profile.md`.
The module itself does not filter — this is caller responsibility.

## Validation notes

- Validated on the HumanFat EC subset from the HumanFat_Yang run
- Batch correction variable: source_file (HumanFat convention — not sample_id)
- Active assay at scoring time: RNA (JoinLayers applied)
- Gene sets authored as explicit vectors to replace v1 invalid range syntax:
  - v1 used `"NDUFA1"-"NDUFA10"` (invalid R) — replaced with explicit OXPHOS gene list above
  - v1 used `"ELOVL1"-"ELOVL7"` — replaced with explicit FA_Synthesis list above
- AUCell outputs: per-cell AUC scores added as metadata columns; summary heatmap and UMAP per pathway
- Outputs written to `output/metabolic_profile/`

## Known issues / quirks

- RibHighEC should be excluded from any analysis plot but is present in the input object;
  the module's visualization code will show it unless explicitly filtered via SUBSET_VAL
- HumanFat uses `source_file` (not `sample_id`) as the batch variable; pass as `BATCH_CORRECTION_VAR`
- The PPARG_Targets gene set (transcriptional targets of PPARG OE) is a separate analysis
  and does not belong in the five metabolic pathways above; see @examples/pparg_bulk_concordance.md
- tissue_type column contains tissue categories; additional_notes contains higher-resolution
  organ location; remove empty strings from additional_notes before grouping
