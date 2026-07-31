---
# Example: Cross-Dataset Dotplot Example (in-house vs. Tabula Sapiens)
# Status: full-build
# Example structure for an in-house vs. public atlas dotplot comparison
---

# Example: Cross-Dataset Dotplot Example (in-house vs. Tabula Sapiens)

## Instantiates
- @modules/cross_dataset_dotplot.md

## Project context
- Validated: 2026-03-24
- Status: full-build

## Brief block

```yaml
atlas: "Tabula Sapiens"             # Known-Atlas Convention; auto-populates source_col and
                                    # atlas_group_col from context/known_atlases.md

downstream_analyses:
  cross_dataset_dotplot:
    enabled: true
    source_col: "dataset"           # metadata column identifying WCM vs TS
    subtype_col: "ec_subtype"       # EC subtype column
                                    # TODO: confirm exact column name
    sources_order:
      - "WCM"                       # in-house dataset (first column group)
      - "TS"                        # Tabula Sapiens (second column group)
    col_group_col: "dataset"        # column used for column grouping (same as source_col here)
    marker_genes:
      # Three-section structure:
      # Section 1: shared EC markers (present in both WCM and TS)
      # Section 2: WCM-specific markers (in-house EC biology not represented in TS)
      # Section 3: TS-specific markers (atlas markers below in-house detection threshold)
      # REPLACE: fill in gene lists from your own project data
      Section_1_Shared:
        - KDR
        - FLT1
        - PECAM1
        - CDH5
        - VWF
        - CLDN5
        - ESAM
        - ERG
        - ETS1
        - PROX1             # lymphatic marker (shared)
      Section_2_ProjectSpecific:
        - "<GENE_1>"  # REPLACE: markers identified in your in-house data
        - "<GENE_2>"
        - "<GENE_3>"
      Section_3_TS_Specific:
        - ACE2               # kidney EC marker in TS
        # TODO: add TS-specific markers from project records
    sec1_force:
      # TFs present in in-house data but below atlas detection threshold due to sequencing depth.
      # These are forced into Section 1 to show in-house expression even when they
      # fall below the automatic variable gene selection threshold in the atlas.
      # FORCE_GENES_SEC1 is a project-specific analytical decision that must be justified per-gene.
      - "<GENE_A>"  # REPLACE: gene confirmed in in-house data but below atlas detection threshold
      - "<GENE_B>"
    group_colors:
      "WCM": "#FFF8E1"          # warm yellow background for in-house WCM column group
      "TS":  "#FFEBEE"          # light red background for primary TS atlas column group
```

## Column label convention

The v1 used the following column naming pattern for the dotplot x-axis:

```r
# In-house columns: "{inhouse_label} {subtype}"
# e.g., "WCM CapEC", "WCM VenEC1", "WCM AEC"
# 
# Atlas columns: "{organ_tissue} EC" (toTitleCase)
# e.g., "Fat EC", "Kidney EC", "Liver EC"
# 
# This is constructed via: paste0(SOURCE_COL_value, " ", SUBTYPE_COL_value)
# for in-house, and the atlas cell_ontology_class values for TS columns
```

## Figure dimensions and output directory

```r
# Validated figure dimensions from 2026-03-24 run:
# TODO: confirm exact width/height from project records
# Typical range: width = 18-22 in, height = 10-14 in (varies with gene count × subtype count)
# Output dir: output2/ (this project used output2/ — set OUTPUT_DIR <- "output2/cross_dataset_dotplot")
```

## Validation notes

- Illustrative example of the cross-dataset dotplot pattern
- Three-section dotplot: shared markers / WCM-specific / TS-specific
- TF diamond variant was also run (Part C); uses `shape = 23` patch on the dotplot
- Depth correction: within-dataset z-scoring applied separately to WCM and TS
  (each dataset z-scored independently per gene before combining)
- Output dir: `output2/` (not `output/`)

## Known issues / quirks

- CRITICAL: Custom fonts "Source Sans 3" (x-axis labels) and "Playfair Display" (y-axis
  gene labels) are NOT available in r-env. The v1 used both; they silently fall back to
  system fonts. The module defaults to system fonts — do not pass `family` arguments.
- The `sec1_force` list is a project-specific analytical decision: genes are forced into
  Section 1 because their in-house expression is biologically confirmed but falls below the
  automatic variable gene threshold in the deeper atlas dataset. Each gene requires justification.
- Column whitespace→dot conversion: if loading `cor_mat` from a CSV round-trip, apply
  `gsub("\\.", " ", colnames(cor_mat))` to restore spaces in column names
- `cairo_pdf` required for patchwork outputs (same issue as co_umap — `useDingbats = FALSE` omitted)
- `SOURCES_ORDER = c("WCM", "TS")` — WCM first ensures in-house columns appear on the left;
  this is the conventional display order for "our data vs atlas" comparisons
