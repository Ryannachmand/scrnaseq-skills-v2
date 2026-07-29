---
# Known Atlases Registry
# Part of the Known-Atlas Convention — see CONVENTIONS.md
---

# Known Atlases Registry

This file is the data layer for the Known-Atlas Convention. It stores validated
parameters for recognized public atlases so that project briefs can refer to an
atlas by name rather than specifying every parameter explicitly.

**What this file is:**
A curated registry of atlas names, their metadata column conventions, and validated
parameters for the `atlas_co_umap` and `cross_dataset_dotplot` modules.

**What this file is NOT:**
Logic for expanding atlas names into brief parameters — that belongs to the deployment
agent (Phase 6). This file is pure data.

---

## Known-Atlas Convention Summary

When a project brief sets `atlas: "Tabula Sapiens"` (or any recognized alias), the
deployment agent looks up the atlas here and injects defaults into the
`downstream_analyses.atlas_co_umap` and `downstream_analyses.cross_dataset_dotplot`
config blocks before generating scripts. The following fields are auto-populated:

- `source_col` — the metadata column identifying dataset source
- `atlas_label` — short atlas identifier for plot titles
- `atlas_group_col` — the column in atlas metadata holding cell type labels
- `n_per_subtype` / `recommended_downsample_n` — validated downsampling rate

Explicitly set values in the module config block take precedence over registry defaults.
The following parameters are always project-specific and must be set explicitly in the
brief even when `brief.atlas` is recognized:

- `atlas_co_umap`: `inhouse_label`, `subtype_col`, `highlight_inhouse`,
  `highlight_atlas`, `coords_csv`
- `cross_dataset_dotplot`: `sources_order`, `marker_genes`, `col_group_col`,
  and the in-house source_col value (= `inhouse_label` from atlas_co_umap)

Cross-references:
- @modules/atlas_co_umap.md — uses `source_col`, `atlas_label`, `atlas_group_col`,
  `n_per_subtype` from this registry
- @modules/cross_dataset_dotplot.md — uses `source_col`, `atlas_group_col` from
  this registry

---

## How to Add a New Atlas

1. Choose a canonical `registry_key` (lowercase, underscores; this is the Python/YAML dict key)
2. Add at least two recognized `display_names` (the full name users would type in `brief.atlas`)
3. Fill all required fields: `citation`, `source_col`, `atlas_label`,
   `atlas_group_col`, `recommended_downsample_n`
4. Run a real project comparing in-house data against the atlas before setting
   `recommended_downsample_n` — the default of 1500 is a starting point, not a universal constant
5. Set `validated_date` after the first successful atlas_co_umap run
6. Add a note to the `notes` block documenting any quirks (normalization status,
   unusual column names, download format, font restrictions)

---

## Registry

```yaml
known_atlases:

  tabula_sapiens:
    display_names:
      - "Tabula Sapiens"
      - "tabula sapiens"
      - "TS"
      - "Tabula_Sapiens"
      - "TabulaSapiens"
    citation: "Tabula Sapiens Consortium, Science 2022 (doi:10.1126/science.abl4896)"
    source_col: "dataset"              # the standard IntegratePublicData metadata column
    atlas_label: "TS"                  # short name used in plot titles and source_col values
    atlas_group_col: "cell_ontology_class"   # cell type column in TS metadata
    typical_cell_count: 483152         # approximate; varies by tissue subset downloaded
    recommended_downsample_n: 1500     # in-house cells per subtype before co-UMAP
                                       # validated: TabulaSapiensComparison project
                                       # produces ≥75% atlas cells in UMAP (80-85% observed)
                                       # with ~38k in-house cells across 6 subtypes
    validated_date: "2026-04-15"
    notes: |
      Download format: h5ad (anndata); convert to Seurat with SeuratDisk or zellkonverter.
      Cell type column: cell_ontology_class (Cell Ontology terms, not free text).
      Normalization status: raw counts in 'counts' layer; pre-normalized in 'X'.
      Always use raw counts layer for integration with in-house data.
      Font restriction: custom fonts ("Source Sans 3", "Playfair Display") are not
      available in r-env — generated scripts must use system fonts. Use cairo_pdf for
      patchwork-based multi-panel PDFs.
      DOWNSAMPLE_N formula: target_n = floor(0.25 * n_atlas_cells / n_subtypes)
      provides a starting estimate; 1500 was validated for the WCM-EC use case.
      Adjust upward for atlases with fewer cells or in-house datasets with more subtypes.
      The atlas_group_col "cell_ontology_class" is coarse-grained — it maps to broad
      cell types (e.g., "endothelial cell of venule"), not subtypes. For cross_dataset_dotplot,
      col_group_col should be set to a project-specific grouping that matches the
      biological comparison of interest.
```

---

## Placeholder Entries (not yet validated)

The following atlases appear in the lab's project history but have not yet been
formally registered with validated parameters. Add them after completing an
atlas_co_umap or cross_dataset_dotplot run with real data.

```yaml
# human_cell_atlas:
#   display_names: ["Human Cell Atlas", "HCA"]
#   citation: "TODO"
#   source_col: "dataset"
#   atlas_label: "HCA"
#   atlas_group_col: null             # TODO: confirm column name from HCA download
#   recommended_downsample_n: 1500    # starting point; validate with real data
#   validated_date: null
#   notes: |
#     TODO: Multiple HCA portal releases with different column naming conventions.
#     Confirm download format and cell type column before registering.

# mouse_cell_atlas:
#   display_names: ["Mouse Cell Atlas", "MCA"]
#   citation: "TODO"
#   source_col: "dataset"
#   atlas_label: "MCA"
#   atlas_group_col: null             # TODO: confirm column name from MCA download
#   recommended_downsample_n: 1500    # starting point; validate with real data
#   validated_date: null
#   notes: |
#     TODO: For mouse projects. Validate organism field matches before integration.
```

---

## Registry Schema Reference

Each entry has the following structure:

```yaml
<registry_key>:
  display_names: [list of recognized strings for brief.atlas matching]
  citation: "Author, Journal Year (doi:...)"
  source_col: "metadata column in the merged object identifying dataset source"
  atlas_label: "short name; value in source_col for atlas cells"
  atlas_group_col: "column in atlas metadata holding cell type labels"
  typical_cell_count: integer  # approximate; for reference
  recommended_downsample_n: integer  # validated in-house cells per subtype for co-UMAP
  validated_date: "YYYY-MM-DD"  # date of first successful atlas_co_umap validation
  notes: |
    Free text: download format, normalization status, known quirks, validation context.
```

Required fields: `display_names`, `citation`, `source_col`, `atlas_label`,
`atlas_group_col`, `recommended_downsample_n`.

Optional fields: `typical_cell_count`, `validated_date`, `notes`.

---

## Deployment Agent Integration

The deployment agent consumes this registry by:
1. Parsing `brief.atlas` string against all `display_names` lists (case-insensitive)
2. On match: reading the registry entry for the matched atlas
3. Injecting the registry values into `downstream_analyses.atlas_co_umap` and
   `downstream_analyses.cross_dataset_dotplot` config blocks as defaults
4. Logging which values came from the registry vs. the brief for traceability
5. On no match: raising a warning and requiring explicit parameters in the module blocks

The registry should be kept as static YAML data. Logic for expansion belongs in the
deployment agent code, not in this file.
