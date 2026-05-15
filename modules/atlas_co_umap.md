---
requires_context:
  palettes:
    - subtype_colors  # optional: named vector for in-house cell subtype labels (Panel B)
  metadata_columns:
    required:
      - source_col        # metadata column identifying dataset source (in-house vs atlas)
      - subtype_col       # cell type/subtype column for in-house cells
      - atlas_group_col   # organ/tissue column for atlas cells
    optional: []
  brief_keys:
    required:
      - output_dir
      - downstream_analyses.atlas_co_umap.source_col
      - downstream_analyses.atlas_co_umap.inhouse_label
      - downstream_analyses.atlas_co_umap.atlas_label
      - downstream_analyses.atlas_co_umap.subtype_col
      - downstream_analyses.atlas_co_umap.atlas_group_col
      - downstream_analyses.atlas_co_umap.highlight_inhouse
      - downstream_analyses.atlas_co_umap.highlight_atlas
      - downstream_analyses.atlas_co_umap.coords_csv
    optional:
      - downstream_analyses.atlas_co_umap.n_per_subtype
      - downstream_analyses.atlas_co_umap.color_inhouse
      - downstream_analyses.atlas_co_umap.color_atlas
      - downstream_analyses.atlas_co_umap.subtype_colors
references:
  - "@primitives/harmony_integration.md"
  - "@primitives/visualization.md"
  - "@primitives/aesthetics.md"
  - "@primitives/seurat_v5_rules.md"
---

# Module: Atlas Co-UMAP Embedding

Produces a co-embedding UMAP comparing in-house cells with a public atlas after Harmony
batch correction. Downsamples in-house cells before UMAP computation so the atlas drives
the 2D geometry rather than being distorted by in-house cell count. Saves coordinates as
CSV for fast replotting. Produces a 4-panel highlight figure.

**Known-Atlas Convention (Phase 3 concern):** When a brief names a recognized atlas (e.g.
"Tabula Sapiens"), the deployment agent can auto-populate `source_col`, `atlas_label`,
`atlas_group_col`, and `n_per_subtype` from the context registry
(`context/known_atlases` — proposed, not yet implemented). See PHASE2C_REPORT.md for the
design. This module is atlas-agnostic; all atlas identity flows through parameters.

---

## Critical Design Decision: Downsample IN-HOUSE Cells Before UMAP

**Wrong approach:** compute UMAP on all cells, then drop in-house cells at plot time.
- If in-house cell count is comparable to or exceeds atlas cell count, the in-house data
  dominates UMAP geometry, pulling structure toward in-house clusters.

**Correct approach:** downsample in-house cells first, then compute UMAP on the smaller
merged set.
- Atlas cells dominate at 75–85% of the input, so atlas topology drives the 2D layout
  and in-house cells fall into it rather than distorting it.
- Target ratio: atlas should comprise ≥75% of UMAP input cells.

**Tuning the downsampling rate:**
The default `n_per_subtype = 1500` was validated for a specific in-house/atlas size ratio
(see Phase 4 examples/ for the specific ratio analysis). For different ratios:
- If (n_subtypes × n_per_subtype) / n_atlas_cells > 0.25, reduce n_per_subtype
- If any subtype has < 500 cells after sampling, the module warns — that subtype may be
  underrepresented in the UMAP

---

## Risks and Known Issues

**Risk 1 — downsampling rate tuning:** The default 1500 cells/subtype was validated for
a particular in-house/atlas size ratio. For different ratios, recalculate:
`target_n = floor(0.25 * n_atlas_cells / n_subtypes)`. Document this in the examples/ file.

**Risk 2 — UMAP geometry sensitivity:** Minor changes to the downsampling rate produce
visually different UMAPs. A warning fires if any subtype samples < 500 cells.
This is expected for rare subtypes and should be noted in figure legends.

**Risk 3 — Custom fonts:** The v1 source referenced `"Source Sans 3"` which is NOT
available in r-env and causes a font error. This module defaults to system fonts
(`base_family = ""`). Custom font override is possible if the font is confirmed installed.

**Risk 4 — cairo_pdf requirement:** `ggsave(..., device = cairo_pdf)` is required for
multi-panel patchwork PDFs. The standard `pdf()` device only renders the first panel in
some environments. If cairo is unavailable, ggsave falls back silently — check the output.
Note: `useDingbats = FALSE` is NOT compatible with `cairo_pdf` — omit it when using cairo.

---

## Brief Schema

```yaml
downstream_analyses:
  atlas_co_umap:
    enabled: true
    source_col: project_specific        # REQUIRED: metadata column identifying dataset source
    inhouse_label: project_specific     # REQUIRED: source_col value for in-house data
    atlas_label: project_specific       # REQUIRED: source_col value for atlas data
    subtype_col: project_specific       # REQUIRED: cell type/subtype column for in-house cells
    atlas_group_col: project_specific   # REQUIRED: organ/tissue column for atlas cells
    highlight_inhouse: project_specific # REQUIRED: in-house subtype to highlight in Panel D
    highlight_atlas: project_specific   # REQUIRED: atlas group to highlight in Panel D
    n_per_subtype: 1500                 # cells sampled per in-house subtype before UMAP
                                        # validate against in-house/atlas size ratio (see above)
    color_inhouse: "#E69F00"            # Okabe-Ito orange (colorblind-safe default)
    color_atlas: "#0072B2"              # Okabe-Ito blue (colorblind-safe default)
    subtype_colors: null                # optional named color vector for Panel B in-house subtypes
    coords_csv: project_specific        # REQUIRED: path to save/load UMAP coordinate CSV
```

---

## R Packages Required

```r
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
```

---

## Configuration Block

```r
# ── CONFIG ────────────────────────────────────────────────────────────────────
SOURCE_COL       <- "project_specific"  # REPLACE: metadata column identifying dataset source
INHOUSE_LABEL    <- "project_specific"  # REPLACE: source_col value for in-house data
ATLAS_LABEL      <- "project_specific"  # REPLACE: source_col value for atlas data
SUBTYPE_COL      <- "project_specific"  # REPLACE: cell type column for in-house cells
ATLAS_GROUP_COL  <- "project_specific"  # REPLACE: organ/tissue column for atlas cells

HIGHLIGHT_INHOUSE <- "project_specific" # REPLACE: in-house subtype to highlight in Panel D
HIGHLIGHT_ATLAS   <- "project_specific" # REPLACE: atlas group to highlight in Panel D

DOWNSAMPLE_N <- 1500   # cells per in-house subtype before UMAP
                        # tune based on in-house/atlas size ratio — see Risk 1 above

# Colorblind-safe highlight colors (Okabe-Ito palette defaults)
COLOR_INHOUSE <- "#E69F00"   # orange for in-house highlight
COLOR_ATLAS   <- "#0072B2"   # blue for atlas highlight

# Optional: named color vector for Panel B in-house subtype labels
# subtype_colors <- c("project_specific" = "#project_specific")  # REPLACE if needed

PT <- 0.4   # uniform point size for all panels

RDS_IN     <- "project_specific"  # REPLACE: path to merged Harmony-corrected Seurat object
                                   # Harmony reduction must already be present in the object
COORDS_CSV <- "project_specific"  # REPLACE: path to save/load UMAP coordinate CSV
OUTPUT_DIR <- file.path("output", "atlas_co_umap")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
```

---

## Step 1: Pre-UMAP Downsampling

Downsample in-house cells to `DOWNSAMPLE_N` per subtype. Atlas cells are kept in full.
This ensures atlas topology dominates the UMAP geometry.

```r
merged <- readRDS(RDS_IN)
merged <- JoinLayers(merged)   # @primitives/seurat_v5_rules.md Rule 1

meta <- merged@meta.data

# Identify atlas vs in-house cells
atlas_cells  <- rownames(meta)[meta[[SOURCE_COL]] == ATLAS_LABEL]
inhouse_meta <- meta[meta[[SOURCE_COL]] == INHOUSE_LABEL, ]

# Downsample in-house cells: DOWNSAMPLE_N cells per subtype
inhouse_sampled <- inhouse_meta %>%
  tibble::rownames_to_column("cell") %>%
  group_by(!!rlang::sym(SUBTYPE_COL)) %>%
  slice_sample(n = DOWNSAMPLE_N, replace = FALSE) %>%
  ungroup() %>%
  pull(cell)

# Warn if any subtype has < 500 cells after sampling
subtype_counts <- table(inhouse_meta[inhouse_sampled, SUBTYPE_COL])
small_subtypes <- names(subtype_counts)[subtype_counts < 500]
if (length(small_subtypes) > 0) {
  warning(sprintf(
    "DOWNSAMPLE_N = %d produces < 500 cells for subtype(s): %s. These may be underrepresented in UMAP.",
    DOWNSAMPLE_N, paste(small_subtypes, collapse = ", ")
  ))
}

keep_cells <- c(atlas_cells, inhouse_sampled)
message(sprintf("UMAP input: %d total (%d atlas = %.0f%%, %d in-house)",
  length(keep_cells), length(atlas_cells),
  100 * length(atlas_cells) / length(keep_cells), length(inhouse_sampled)))

sub <- subset(merged, cells = keep_cells)
```

---

## Step 2: UMAP Computation and Coordinate CSV

```r
# Detect Harmony reduction name — @primitives/seurat_v5_rules.md Rule 5
harm_key <- names(sub@reductions)[grepl("harmony", names(sub@reductions), ignore.case = TRUE)][1]
if (is.na(harm_key)) stop("No Harmony reduction found. Run @primitives/harmony_integration.md first.")

sub <- RunUMAP(sub, reduction = harm_key, dims = 1:30,
               umap.method = "uwot", metric = "cosine", verbose = FALSE)

# Extract and save coordinates — all downstream plotting loads this CSV, no Rds reload needed
umap_key <- names(sub@reductions)[grepl("umap", names(sub@reductions), ignore.case = TRUE)][1]
coords <- as.data.frame(Embeddings(sub, umap_key))
colnames(coords) <- c("UMAP1", "UMAP2")
coords[[SOURCE_COL]]     <- sub@meta.data[[SOURCE_COL]]
coords[[SUBTYPE_COL]]    <- ifelse(sub@meta.data[[SOURCE_COL]] == INHOUSE_LABEL,
                                    sub@meta.data[[SUBTYPE_COL]], NA_character_)
coords[[ATLAS_GROUP_COL]] <- ifelse(sub@meta.data[[SOURCE_COL]] == ATLAS_LABEL,
                                     sub@meta.data[[ATLAS_GROUP_COL]], NA_character_)
write.csv(coords, COORDS_CSV, row.names = TRUE, quote = FALSE)
message("Coordinates saved to: ", COORDS_CSV)
```

---

## Step 3: Load Coordinates and Build 4-Panel Figure

Load from CSV for fast replotting — no Seurat object needed after this point.

```r
coords <- read.csv(COORDS_CSV, row.names = 1, stringsAsFactors = FALSE)

inhouse_df <- coords[!is.na(coords[[SUBTYPE_COL]]),  ]
atlas_df   <- coords[!is.na(coords[[ATLAS_GROUP_COL]]), ]

# ── Panel A: Dataset source ───────────────────────────────────────────────────
# Colors derived from parameters — no hardcoded dataset labels
source_colors <- c(
  setNames("#C62828", INHOUSE_LABEL),
  setNames("#1565C0", ATLAS_LABEL)
)
p_source <- ggplot(coords, aes(UMAP1, UMAP2, color = !!rlang::sym(SOURCE_COL))) +
  geom_point(size = PT, alpha = 0.5, stroke = 0) +
  scale_color_manual(values = source_colors, name = "Dataset") +
  theme_classic(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 10), legend.position = "right") +
  labs(title = "Dataset source")

# ── Panel B: In-house subtypes ────────────────────────────────────────────────
# Atlas cells grey; in-house cells colored by subtype
sub_cols <- if (exists("subtype_colors")) subtype_colors else NULL

p_inhouse <- ggplot() +
  geom_point(data = atlas_df,   aes(UMAP1, UMAP2),
             color = "grey80", size = PT, alpha = 0.3, stroke = 0) +
  geom_point(data = inhouse_df, aes(UMAP1, UMAP2,
             color = !!rlang::sym(SUBTYPE_COL)),
             size = PT, alpha = 0.9, stroke = 0) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 10), legend.position = "right") +
  labs(title = "In-house subtypes")

if (!is.null(sub_cols)) {
  p_inhouse <- p_inhouse + scale_color_manual(values = sub_cols, name = "Subtype")
}

# ── Panel C: Atlas groups ─────────────────────────────────────────────────────
# In-house cells grey; atlas cells colored by group
p_atlas <- ggplot() +
  geom_point(data = inhouse_df, aes(UMAP1, UMAP2),
             color = "grey80", size = PT, alpha = 0.3, stroke = 0) +
  geom_point(data = atlas_df,   aes(UMAP1, UMAP2,
             color = !!rlang::sym(ATLAS_GROUP_COL)),
             size = PT, alpha = 0.5, stroke = 0) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 10), legend.position = "right") +
  labs(title = "Atlas groups")

# ── Panel D: Colorblind-safe highlight ────────────────────────────────────────
# Isolates primary comparison of interest. Everything else grey.
# COLOR_INHOUSE and COLOR_ATLAS default to Okabe-Ito orange + blue (colorblind-safe).
highlight_label_inhouse <- HIGHLIGHT_INHOUSE
highlight_label_atlas   <- HIGHLIGHT_ATLAS

p_highlight <- ggplot() +
  geom_point(data = atlas_df[atlas_df[[ATLAS_GROUP_COL]] != highlight_label_atlas, ],
             aes(UMAP1, UMAP2), color = "grey85", size = PT, alpha = 0.25, stroke = 0) +
  geom_point(data = inhouse_df[inhouse_df[[SUBTYPE_COL]] != highlight_label_inhouse, ],
             aes(UMAP1, UMAP2), color = "grey60", size = PT, alpha = 0.30, stroke = 0) +
  geom_point(data = atlas_df[atlas_df[[ATLAS_GROUP_COL]] == highlight_label_atlas, ],
             aes(UMAP1, UMAP2), color = COLOR_ATLAS, size = PT, alpha = 0.85, stroke = 0) +
  geom_point(data = inhouse_df[inhouse_df[[SUBTYPE_COL]] == highlight_label_inhouse, ],
             aes(UMAP1, UMAP2), color = COLOR_INHOUSE, size = PT, alpha = 0.95, stroke = 0) +
  geom_point(data = data.frame(UMAP1 = NA_real_, UMAP2 = NA_real_,
                                grp = c(highlight_label_inhouse, highlight_label_atlas)),
             aes(UMAP1, UMAP2, color = grp), na.rm = TRUE, size = 2) +
  scale_color_manual(
    values = setNames(c(COLOR_INHOUSE, COLOR_ATLAS),
                      c(highlight_label_inhouse, highlight_label_atlas)),
    name = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 10), legend.position = "right") +
  labs(title = sprintf("%s vs %s", highlight_label_inhouse, highlight_label_atlas))
```

---

## Step 4: Combine and Save

Use `cairo_pdf` for patchwork multi-panel output. Do NOT use `useDingbats = FALSE`
with `cairo_pdf` — that argument is only valid for the base `pdf()` device.

```r
combined <- (p_source | p_inhouse | p_atlas | p_highlight) +
  plot_annotation(title = sprintf("Co-embedding UMAP: %s + %s", INHOUSE_LABEL, ATLAS_LABEL))

out_base <- file.path(OUTPUT_DIR, "co_umap_4panel")

ggsave(paste0(out_base, ".pdf"),
       combined, width = 21, height = 5.5, device = cairo_pdf)
ggsave(paste0(out_base, ".png"),
       combined, width = 21, height = 5.5, dpi = 200, bg = "white")

message("4-panel figure saved to: ", out_base)
```

---

## Common Pitfalls

| Problem | Cause | Fix |
|---------|-------|-----|
| Only first panel renders in PDF | `pdf()` + `print(combined)` used instead of ggsave | Use `ggsave(..., device = cairo_pdf)` |
| Error: `useDingbats` not recognized | `useDingbats = FALSE` passed with `cairo_pdf` | Remove `useDingbats` when `device = cairo_pdf` |
| Font error in PDF | Custom `base_family` font not installed in r-env | Remove custom `family` argument; use default system font |
| In-house cells dominate UMAP shape | UMAP computed before downsampling | Always downsample in-house cells before `RunUMAP` |
| Rare subtype invisible in UMAP | `n_per_subtype` exceeds actual subtype cell count | Watch for `< 500 cells` warning; note in figure legend |
| atlas_df missing rows | `ATLAS_GROUP_COL` has NAs in in-house rows | The `!is.na(coords[[ATLAS_GROUP_COL]])` filter handles this correctly |

---

## Project-Specific Values (Stage for Phase 4 examples/)

`examples/tabula_sapiens_co_umap.md` must define:

- `INHOUSE_LABEL = "WCM"`, `ATLAS_LABEL = "TS"`, `SOURCE_COL = "dataset"`
- `SUBTYPE_COL` = EC subtype column name (confirm from TabulaSapiensComparison project)
- `ATLAS_GROUP_COL` = atlas organ column name
- `HIGHLIGHT_INHOUSE = "CapEC"`, `HIGHLIGHT_ATLAS = "Fat EC"`
- `N_PER_SUBTYPE = 1500` with ratio analysis context:
  - In-house: ~38k cells across 6 subtypes; atlas: ~100k cells
  - 1500/subtype × 6 subtypes = ~9k in-house vs ~100k atlas → atlas = 91% of UMAP input
  - Result: atlas topology dominates while all 6 in-house subtypes are visible
- `COORDS_CSV = "output2/ts_harmony_umap_presample_coords.csv"`
- Validated script filenames: `06f_ts_umap_presample.R` (embed), `06g_ts_umap_replot.R` (replot)
- Validated figure dimensions: `width = 21, height = 5.5` for 4-panel
- Validated run: TabulaSapiensComparison, 2026-04-15
- **Known font issue:** v1 used `base_family = "Source Sans 3"` which is not available in
  r-env and causes a font error (documented in v1 pitfalls section). System fonts used instead.
