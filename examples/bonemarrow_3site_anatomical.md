---
# Example: BoneMarrowStroma 3-Site Anatomical DE Analysis
# Status: full-build
# Validated: BoneMarrowStroma run (2026-03-16 per Phase 2C staging notes)
---

# Example: BoneMarrowStroma 3-Site Anatomical DE Analysis

## Instantiates
- @modules/multi_group_de_analysis.md

## Project context
- Project: BoneMarrowStroma
- Validated: 2026-03-16 (BoneMarrowStroma run per Phase 2C staging notes)
- Status: full-build

## Brief block

```yaml
downstream_analyses:
  multi_group_de:
    enabled: true
    group_col: site                 # anatomical site column
    label_col: unified_label        # harmonized cell type label column
    groups:                         # 3 anatomical sites in display order
      - "Vertebrae"
      - "Iliac Crest"
      - "Femoral Head"
    comparisons:
      - label: "Vertebrae_vs_IliacCrest"
        ident1: "Vertebrae"         # positive LFC (expressed more in Vertebrae)
        ident2: "Iliac Crest"
      - label: "Vertebrae_vs_FemoralHead"
        ident1: "Vertebrae"
        ident2: "Femoral Head"
      - label: "IliacCrest_vs_FemoralHead"
        ident1: "Iliac Crest"
        ident2: "Femoral Head"
    functional_gene_sets:
      Osteogenic:
        - RUNX2
        - SP7
        - COL1A1
        - COL1A2
        - ALPL
        - BGLAP
        - SPP1
        - IBSP
        - SPARC
        - FN1
        - MMP13
        - CTSK
      Adipogenic:
        - PPARG
        - CEBPA
        - CEBPB
        - ADIPOQ
        - FABP4
        - LPL
        - FASN
        - SCD
        - PLIN1
        - RETN
        - LEP
        - ADIPSIN
      Fibroblastic:
        - COL3A1
        - COL6A1
        - COL6A2
        - COL6A3
        - VIM
        - S100A4
        - PDGFRA
        - PDGFRB
        - FAP
        - THY1
        - ACTA2
      Hematopoietic_Niche:
        - CXCL12
        - KITLG
        - ANGPT1
        - LEPR
        - IGF1
        - VCAM1
        - IL7
        - FLT3L
        - THPO
        - EPO
    label_order: null               # TODO: STROMA_ORDER — requires human curation from
                                    # BoneMarrowStroma project records. STROMA_ORDER is the
                                    # ordered vector of stromal cell type names for x-axis
                                    # display in make_anatomical_dotplot. It was referenced
                                    # in v1 but never formally defined. Must be established
                                    # before running this module.
    label_colors: null              # TODO: per-cell-type axis label colors; keyed by
                                    # STROMA_ORDER values; null = all grey30

context_overrides:
  palettes:
    group_colors:
      "Vertebrae":    "#project_specific"   # TODO: confirm direction strip colors
      "Iliac Crest":  "#project_specific"   # TODO: confirm direction strip colors
      "Femoral Head": "#project_specific"   # TODO: confirm direction strip colors
```

## STROMA_ORDER curation note

`STROMA_ORDER` was referenced in the v1 `anatomical_de_analysis.md` but was never
formally defined in any library file. Before running this module, the user must
curate the ordered vector of stromal cell type names from BoneMarrowStroma project
records. Suggested curation process:

1. Run `sort(unique(so$unified_label))` to list all cell types
2. Order them biologically (e.g., progenitor → adipogenic → osteogenic → fibroblastic → other)
3. Set `label_order` in the brief to this ordered vector
4. Optionally set `label_colors` as a named vector of hex colors keyed by the ordered names

## Functional gene sets note

The `functional_gene_sets` above are authored for bone marrow stromal biology. They represent
the four major functional programs of bone marrow stromal cells (MSCs):
- Osteogenic: bone formation / mineralization
- Adipogenic: fat cell differentiation
- Fibroblastic: connective tissue / matrix production
- Hematopoietic_Niche: support of hematopoietic progenitor niches

## Validation notes

- Validated BoneMarrowStroma run (2026-03-16 per Phase 2C staging notes)
- 3 anatomical sites: Vertebrae, Iliac Crest, Femoral Head
- 3 pairwise comparisons (all combinations)
- `functional_gene_sets` above are new — v1 had EC-biology gene sets that were removed
  in Phase 1; the bone marrow stromal gene sets above are Phase 4 additions
- Outputs written to `output/multi_group_de/`

## Known issues / quirks

- `STROMA_ORDER` must be curated before running — use `label_order: null` for alphabetical
  ordering as a first pass, then refine after reviewing the output
- `label_order` is the v2 equivalent of the v1 `STROMA_ORDER` constant
- The comparisons use ASCII dash not em-dash in labels (mbcsToSbcs rendering error with em-dash)
- The `is_ambient()` function must be loaded from @primitives/differential_expression.md,
  not redefined locally — the module's YAML frontmatter ensures this via `references:` block
- `make_anatomical_dotplot()` requires ALL groups (not just the comparison pair) to be present
  in the Seurat object; do NOT subset to ident1/ident2 before calling this function
