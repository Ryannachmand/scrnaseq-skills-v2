---
# Example: Functional Gene Sets (Reference)
# Status: full-build
---

# Example: Functional Gene Sets (Reference)

This is an example `functional_gene_sets` reference for use in
`@modules/multi_group_de_analysis.md` and `@modules/de_comprehensive_csv.md`.

The four named categories below were hardcoded in the v1
`differential_expression.md` and removed in Phase 1 (F4 fix) to decouple
biology-specific gene sets from the primitive. They are documented here as the
authoritative source for functional annotations.

**Attribution:**
Gene list organization reflects biological pathway consensus from published literature
(Aird 2007, Goveia 2019, Potente & Carmeliet 2017).

---

## How to instantiate in a brief

```yaml
downstream_analyses:
  multi_group_de:
    functional_gene_sets:
      Adhesion_Immune_Trafficking:
        - ICAM1
        - ICAM2
        - VCAM1
        - SELE
        - SELP
        - PECAM1
        - ESAM
        - CDH5
        - CLDN5
        - OCLN
        - TJP1
        - JAM2
        - JAM3
        - MADCAM1
        - CX3CL1
        - CXCL12
        - CCL2
        - CCL5
        - CCL14
        - CXCL1
        - CXCL5
        - CXCL8
        - IL6
        - IL8
        - ACKR1
        - DARC
        - MCAM
        - CD34
        - CD44
        - CD47
        - CD99
        - CD93
      Signaling_Receptors_Ligands:
        - VEGFA
        - VEGFB
        - VEGFC
        - VEGFD
        - FLT1
        - KDR
        - NRP1
        - NRP2
        - ANGPT1
        - ANGPT2
        - ANGPT4
        - TEK
        - TIE1
        - NOTCH1
        - NOTCH4
        - DLL4
        - JAG1
        - JAG2
        - HES1
        - HEY1
        - HEY2
        - HEYL
        - EFNB2
        - EPHB4
        - BMP4
        - BMP9
        - BMPR2
        - ROBO4
        - SLIT2
        - WNT2
        - WNT5A
        - FGF1
        - FGF2
        - FGF8
        - FGFR1
        - PDGFB
        - PDGFC
        - IGFBP3
        - IGF1
        - KIT
        - KITLG
        - THPO
        - EPO
        - EPOR
        - IL7
        - FLT3L
      Extracellular_Matrix:
        - COL1A1
        - COL1A2
        - COL3A1
        - COL4A1
        - COL4A2
        - COL15A1
        - COL18A1
        - COL8A1
        - COL8A2
        - FN1
        - LAMA4
        - LAMA5
        - LAMB1
        - LAMB2
        - LAMC1
        - HSPG2
        - SDC1
        - SDC2
        - SDC4
        - GPC1
        - GPC4
        - ELN
        - FBLN1
        - FBLN2
        - FBLN5
        - LTBP1
        - LTBP2
        - LTBP3
        - LTBP4
        - MMP1
        - MMP2
        - MMP9
        - MMP14
        - TIMP1
        - TIMP2
        - TIMP3
        - ADAM10
        - ADAM17
        - ADAMTS1
        - ADAMTS5
        - SPARC
        - SPARCL1
        - OGN
        - MFAP2
        - MFAP4
        - MFAP5
        - EMILIN1
        - EMILIN2
        - EMCN
      Metabolic_Specialized:
        - GLUT1
        - SLC2A1
        - SLC2A3
        - HK1
        - HK2
        - LDHA
        - LDHB
        - PDHA1
        - PDHB
        - CS
        - IDH1
        - IDH2
        - SDHA
        - SDHB
        - COX4I1
        - ATP5A1
        - PPARG
        - PPARA
        - PPARGC1A
        - FABP4
        - FABP5
        - LPL
        - APOE
        - ABCA1
        - NOS3
        - NOS1
        - NOS2
        - PTGIS
        - TBXA2R
        - HMOX1
        - HMOX2
        - SOD1
        - SOD2
        - CAT
        - GPX1
        - GPX4
        - TXNIP
        - BNIP3L
        - BNIP3
        - LC3B
        - BECN1
        - ATG5
        - ATG7
        - TFAM
        - MFN1
        - MFN2
        - OPA1
        - DRP1
```

---

## Category descriptions

| Category | Biological function | Typical use |
|---|---|---|
| Adhesion_Immune_Trafficking | Leukocyte adhesion, junctional proteins, chemokine expression | Inflammation response; cell activation state |
| Signaling_Receptors_Ligands | Paracrine niche signals; growth factors, Notch/Wnt/BMP/VEGF | Tissue support function; identity in specialized niches |
| Extracellular_Matrix | Basement membrane, fibrillar collagens, matrix remodeling | Structural identity; tissue-specific ECM niche |
| Metabolic_Specialized | Glycolysis, OXPHOS, lipid metabolism, redox balance, autophagy | Metabolic state; relevant in vascularized tissues |

---

## Usage in multi_group_de_analysis brief

Reference the categories above directly in the `functional_gene_sets` block.
The module runs `make_functional_dotplot()` and `make_functional_heatmap()` from
`@primitives/differential_expression.md` conditioned on these gene sets being non-null.

```yaml
# Minimal example — use all four categories:
downstream_analyses:
  multi_group_de:
    functional_gene_sets: !include functional_gene_sets_example    # conceptual shorthand
    # In practice: copy the four categories from this file into the brief directly
```

## Usage in de_comprehensive_csv brief

`de_comprehensive_csv.md` does not use `functional_gene_sets` directly — that module
produces a ranked CSV output. But the functional categories above can be used post-hoc
to annotate the CSV output (e.g., labeling genes in specific functional categories).

## Notes on gene list curation

- All genes verified as human HGNC symbols
- GLUT1 = SLC2A1 (same gene; use SLC2A1 for Seurat row matching since symbols are canonical)
- LC3B = MAP1LC3B (same gene; use MAP1LC3B if the Seurat object uses HGNC canonical symbols)
- KDR = VEGFR2 (same gene; KDR is the HGNC symbol)
- For mouse: convert to mouse gene symbols (typically lowercase with initial cap:
  Icam1, Vcam1, etc.) if `organism: "Mus musculus"` is set in the brief
