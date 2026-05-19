---
# Example: NKXSpleen Bulk LFC Concordance (Skeleton)
# Status: skeleton
# Validated: TODO
---

# Example: NKXSpleen Bulk LFC Concordance

## Instantiates
- @modules/bulk_concordance.md (Mode 2 — parallel_lfc)

## Project context
- Project: NKXSpleen
- Validated: 2026-03-13 (per Phase 2B staging notes; formal registry entry not yet created)
- Status: skeleton

## Brief block

```yaml
downstream_analyses:
  bulk_concordance:
    mode: "parallel_lfc"           # Mode 2: gene-level LFC concordance heatmap
    bulk_csv: "project_specific"   # TODO: path to bulk DESeq2 result CSV
                                    # NKX2-3 OE vs Control HUVECs (12 bulk samples total)
    experiment_label: "NKX2-3_OE_vs_Control"
                                    # TODO: confirm exact label used in validated run
    subtype_col: "project_specific"  # TODO: confirm EC subtype column name for NKXSpleen
    group_col: "project_specific"    # TODO: confirm group column
    condition_samples:               # bulk samples in condition (NKX2-3 OE arm)
      # TODO: fill in from DESeq2 sample metadata
      # 6 NKX2-3 OE samples
    control_samples:                 # bulk samples in control arm
      # TODO: fill in from DESeq2 sample metadata
      # 6 Control samples
    biology_gene_sets:               # gene categories for Mode 2 heatmap row annotation
      Chemokines:
        # TODO: fill in from NKXSpleen project records (v1 had NKX2-3-specific list)
        - "CXCL9"
        - "CXCL10"
        - "CXCL16"
        - "CX3CL1"
      Adhesion:
        - "ICAM1"
        - "VCAM1"
        - "SELE"
        - "SELP"
      Angiocrine:
        - "ANGPT1"
        - "ANGPT2"
        - "VEGFA"
        - "KITLG"
      ECM:
        - "COL4A1"
        - "COL4A2"
        - "FN1"
        - "LAMA4"
      Sinusoidal_Markers:
        - "CLEC4G"       # liver sinusoidal EC marker; present in spleen sinusoids
        - "LYVE1"
        - "STAB1"
        - "STAB2"
        - "FCGR2B"
      TF_Program:
        # TFs regulated by NKX2-3 OE — TODO: fill in from project records
        - "NKX2-3"
        - "project_specific"
```

## Dataset context (from Phase 2B staging notes)

- 80,140 EC cells in the NKXSpleen scRNA-seq dataset
- 10 EC subtypes (spleen sinusoidal EC subtypes — names TODO)
- 12 bulk RNA-seq samples: 6 Control + 6 NKX2-3 OE (HUVECs or spleen ECs — confirm)
- 2-batch DESeq2 design (batches affect sample selection — see sc_group assignment pattern)
- Validated 2026-03-13 (per Phase 2B bulk_concordance staging notes)

## sc_group assignment pattern (Mode 2 specific)

Mode 2 requires defining a `sc_group` column that maps sc cells to "condition" vs "control":

```r
# sc_group assignment for NKXSpleen
# TODO: fill in from project records — which sc cells correspond to each bulk condition?
# Pattern from bulk_concordance.md Mode 2:
so$sc_group <- dplyr::case_when(
  so[[SUBTYPE_COL]] %in% c("project_specific") ~ "condition",   # target subtypes
  so[[SUBTYPE_COL]] %in% c("project_specific") ~ "control",     # reference subtypes
  TRUE ~ NA_character_
)
```

## Validation notes

- Phase 2B confirmed Mode 2 was validated on NKXSpleen (2026-03-13)
- NKXSpleen is NOT yet in validated_examples.yaml; adding the registry entry is a Phase 5/6 task
- The biology_gene_sets above are partially reconstructed from Phase 2B staging notes
  (chemokines, adhesion, angiocrine, ECM, sinusoidal markers, TF program);
  specific gene lists are placeholders — TODO to fill from project records

## TODOs for future wrap-up agent

1. Retrieve exact bulk CSV path and confirm experiment_label from NKXSpleen run records
2. Retrieve condition_samples / control_samples from DESeq2 sample metadata
3. Confirm `subtype_col` and `group_col` column names for NKXSpleen scRNA-seq object
4. Fill in biology_gene_sets with NKX2-3-specific gene lists from project records
5. Fill in sc_group assignment (which sc subtypes map to condition vs control)
6. Add NKXSpleen to validated_examples.yaml with full context_defaults
7. Retrieve EC subtype names (10 subtypes)
8. Confirm validated figure dimensions from 2026-03-13 run
9. Mark this example as full-build once the above are resolved
10. Consider adding nkxspleen_pyscenic.md skeleton if SCENIC was also run on this project
