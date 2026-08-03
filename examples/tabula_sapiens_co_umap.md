---
# Example: Cross-Dataset UMAP Example (in-house vs. Tabula Sapiens)
# Status: full-build
# Example structure for cross-dataset UMAP
---

# Example: Cross-Dataset UMAP Example (in-house vs. Tabula Sapiens)

## Instantiates
- @modules/atlas_co_umap.md

## Project context
- Validated: 2026-04-15
- Status: full-build

## Brief block

```yaml
atlas: "Tabula Sapiens"             # Known-Atlas Convention — auto-populates several params
                                    # from context/known_atlases.md registry:
                                    #   source_col: "dataset"
                                    #   atlas_label: "TS"
                                    #   atlas_group_col: "cell_ontology_class"
                                    #   recommended_downsample_n: 1500

downstream_analyses:
  atlas_co_umap:
    enabled: true
    source_col: "dataset"           # metadata column identifying dataset source (WCM vs TS)
    inhouse_label: "WCM"            # in-house data label in source_col
    atlas_label: "TS"               # Tabula Sapiens label in source_col
    subtype_col: "cell_subtype"     # cell subtype column (in-house data)
                                    # TODO: confirm exact column name in this project's object
    atlas_group_col: "cell_ontology_class"   # TS cell type column
    highlight_inhouse: "CapEC"      # in-house subtype to highlight in comparison panels
    highlight_atlas: "<ATLAS_CELL_TYPE>"  # atlas cell group to highlight in comparison panels
    n_per_subtype: 1500             # validated downsampling rate
                                    # rationale: see WCM/TS ratio analysis below
    coords_csv: "output2/ts_harmony_umap_presample_coords.csv"
                                    # path to save/load UMAP coordinates CSV

context_overrides:
  palettes:
    subtype_colors:
      "AEC":    "#REPLACE_HEX"
      "CapEC":  "#REPLACE_HEX"
      "CapEC2": "#REPLACE_HEX"
      "VenEC1": "#REPLACE_HEX"
      "VenEC2": "#REPLACE_HEX"
      "VenEC3": "#REPLACE_HEX"
```

## In-house/atlas ratio analysis (downsampling rationale)

```
Dataset composition before downsampling:
  In-house: N in-house cells, K subtypes
  Atlas: ~N cells across multiple tissues

After downsampling at n_per_subtype = 1500:
  in-house contribution: n_per_subtype × K subtypes = ~M cells
  Atlas contribution: ~100,000 cells (no downsampling applied to atlas)

UMAP input composition:
  In-house: ~M / (M + 100,000) = minority of UMAP input
  Atlas:    ~100,000 / (M + 100,000) = majority of UMAP input
  → Atlas dominates UMAP geometry

Interpretation: at 1500/subtype, in-house cells embed into atlas topology rather than
forcing atlas cells into in-house topology. This is the desired behavior for atlas
comparison — it reveals where in-house subtypes fall within the atlas landscape.
Increasing DOWNSAMPLE_N above ~2500 risks in-house cells distorting atlas topology.
```

## Validated script filenames and figure dimensions

```r
# Example script filenames for a cross-dataset UMAP run:
# 06f_ts_umap_presample.R   — pre-downsample embedding + coordinate saving
# 06g_ts_umap_replot.R      — replot from saved coordinates with annotation panels

# Validated figure dimensions:
fig_w <- 21    # inches (4-panel layout: WCM only, TS only, co-UMAP, highlight panel)
fig_h <- 5.5   # inches
# Output: ggsave(..., device = cairo_pdf, width = 21, height = 5.5, units = "in")
# Note: cairo_pdf required for patchwork multi-panel PDFs (standard pdf() only renders
# the first panel in patchwork; useDingbats = FALSE must be OMITTED with cairo_pdf)
```

## Known-Atlas Convention

When `brief.atlas: "Tabula Sapiens"` is set, the deployment agent auto-populates:
- `source_col: "dataset"` (IntegratePublicData standard column)
- `atlas_label: "TS"`
- `atlas_group_col: "cell_ontology_class"` (confirmed TS cell type column)
- `recommended_downsample_n: 1500`

Explicit values in the `atlas_co_umap` config block take precedence over registry defaults.

## Validation notes

- Illustrative example of the cross-dataset UMAP pattern
- Output dir: `output2/` (this project used output2/ — set OUTPUT_DIR accordingly)
- 4-panel co-UMAP output: WCM cells only, TS cells only, co-embedding, highlight comparison

## Known issues / quirks

- CRITICAL: Custom font "Source Sans 3" is NOT available in r-env. The v1 used
  `base_family = "Source Sans 3"` which silently failed and produced unstyled text.
  The module defaults to system fonts — do not pass a `base_family` argument.
- `cairo_pdf` required for patchwork multi-panel output; standard `pdf()` only renders
  the first panel in a patchwork layout
- `useDingbats = FALSE` must be OMITTED when using `cairo_pdf` (it is not a valid argument)
- The UMAP coordinates CSV (`coords_csv`) is saved/loaded to avoid re-running the embedding;
  if the Seurat object or downsampling changes, delete the CSV to force re-embedding
