---
# Shared Aesthetics Principles — v2
# Migrated from ~/claude-skills/shared/aesthetics.md
# These are visual rules and philosophy — not plot-specific code.
# Plot-specific R code lives in primitives/visualization.md.
# Last validated: HumanFat_Yang run, 2026-03-05
---

# Shared Aesthetics Principles

Generalizable visual standards that apply across ALL pipelines.
These are rules and philosophy — not plot-specific code.
For R code implementations of these principles, see @primitives/visualization.md.

---

## Typography

### Text Sizes
| Element | Size | Notes |
|---|---|---|
| Plot title | 16–18pt | 16pt standard; 18pt for standalone figure titles |
| Axis title | 15pt | Both x and y |
| Axis tick labels | 13–14pt | 14pt preferred; 13pt when space is tight |
| Legend title | 12.5–14pt | Match axis tick labels |
| Legend text | 12.5pt | |
| Annotation text (inside plot) | 4.5–7pt | Section labels: 4.5pt; significance stars: 7pt |
| Small tile labels (metadata tracks) | 1.8pt | Values printed inside small tiles |
| Gene names (y-axis on violin/dot plots) | 13pt | Always italic |
| Cell type labels (x-axis) | 13.5–14.5pt | 14.5pt on dot plots; 13.5pt on violins |
| Subtitle / caption | 12.5pt | |

### Font Styles
- Gene names: always **italic** regardless of plot type
- Plot titles: bold for differential abundance and standalone figures; plain elsewhere
- Cell type x-axis labels on dot plots: **bold**, colored by cell type
- Section labels on grouped plots: **bold.italic**
- Annotation text inside plots: plain or bold.italic depending on context

---

## Color Philosophy

### General Rules
- Never rely on ggplot2 default color palettes — always specify explicitly
- Colors assigned to cell types must be **consistent across all plots in a project**
  Define them once in the project context and reference everywhere
- Use **colorblind-friendly** palettes where possible for categorical data
- Continuous data: use sequential palettes with meaningful directionality
  (e.g., warm = high expression, cool = low)
- Diverging data (e.g., log2 odds ratio, fold change): always anchor midpoint at 0
  using `midpoint = 0` or `rescale(c(...), to = c(0,1))`

### Canonical Palette References
All canonical color palette definitions live in @context/color_palettes.md.
Reference that file rather than hardcoding color values inline.

### Continuous Scales (quick reference)
| Data type | Palette | Notes |
|---|---|---|
| Gene expression (dot/violin) | `c("#F5F5F5","#FFF9C4","#FFB300","#E53935")` | white→yellow→orange→red |
| TF diverging | `c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#D6604D","#B2182B")` | 6-stop blue→red |
| Differential abundance | `gradient2(low="steelblue", mid="white", high="firebrick")` | anchored at 0 |
| Pseudotime | `scale_color_viridis_c(option="magma")` | never inferno |
| Age metadata track | `c("#FFF5EB","#FD8D3C","#7F2704")` | warm orange |
| BMI metadata track | `c("#EFF3FF","#6BAED6","#084594")` | blue |

### Categorical / Annotation Colors
- Sex: `F = "#E07A9F"`, `M = "#45B7D1"`
- Background section bands on grouped plots: alternate `#F8F8F8` / `#FFFFFF` / `#EEF4FF`
- Separator lines between groups: always white (`color = "white"`)
- Grid lines when needed: `color = "grey92", linewidth = 0.3`
- Panel borders when needed: `color = "grey70", fill = NA, linewidth = 0.5`

---

## Layout and Spacing

### Themes
- Default base theme: `theme_classic()` — use for most plots
- `theme_classic(base_size = 14)` for stacked violins
- `theme_classic(base_size = 16)` for trajectory plots
- `theme_minimal()` only when grid lines are genuinely informative (rarely)
- Never use default `theme_gray()`

### Axes
- Remove x-axis line and ticks on stacked/grouped plots where x is categorical
  (`axis.line.x = element_blank(), axis.ticks.x = element_blank()`)
- Y-axis lines and ticks: `color = "grey60", linewidth = 0.3` (not full black)
- Axis text that needs rotation: always `angle = 45, hjust = 1`
  Exception: dot plot top axis uses `angle = 45, hjust = 0, vjust = 0`

### Legends
- Default position: `"right"`
- For dot plots and TF plots: **always save legend separately** using `cowplot::get_legend()`
  The main plot has `legend.position = "none"`
- Legend key sizes: `unit(0.35–0.5, "cm")` depending on plot
- Multi-panel figures: use `plot_layout(guides = "collect") & theme(legend.position = "right")`

### Margins
- Standard plot margin: `margin(t=10, b=10, l=10, r=10)`
- Wide right margin when section labels extend outside panel: `margin(t=10, r=90, b=10, l=10)`
  Always pair with `coord_cartesian(clip = "off")`
- Metadata track plots: `margin(t=1, b=1, l=5, r=5)`

### Multi-panel Assembly
- Use `patchwork` for all multi-panel figures
- Use `wrap_plots(list, ncol=1)` for stacking variable-length plot lists
  **Never use the `/` operator for more than ~8 panels** — it fails with long lists
- Panel heights for proportion plot + metadata tracks: `c(10, 0.6, 0.6, 0.6, 0.6, 0.6)`

---

## Figure Sizing

### Principles
- Always use **dynamic sizing formulas** based on data dimensions — never hardcode
- Build formulas around `n_genes`, `n_celltypes`, `n_samples` as appropriate
- Target readable label sizes: roughly 0.28–0.51 inches per gene row, 0.35–0.75 inches per cell type column

### Standard Formulas
| Plot type | Height | Width |
|---|---|---|
| Stacked violin | `n_genes * 0.51 + 1.5` | `n_ctypes * 0.35 + 1.5` |
| Dot plot | `n_genes * 0.28 + 2` | `n_ctypes * 0.75 + 3.0` |
| TF diamond plot | `n_genes * 0.29 + 2` | `n_ctypes * 0.75 + 3.0` |
| Proportion plot | `5` (fixed) | `max(10/1.65, n_samples * 0.38/1.65 + 3.5)` |

---

## File Output

- Always use `ggsave(..., device = "pdf", useDingbats = FALSE)` — never `pdf()` / `dev.off()`
- Units: always `units = "in"`
- One file per plot type — never bundle unrelated plots into one PDF
- Legend files: saved separately as `{plotname}_legend.pdf`
- All outputs go to `output/<module_name>/` under the project directory

---

## UMAP Feature Plots (gene expression)

### Color scale
- Always use `c("lightgrey", "blue")` — matches Seurat FeaturePlot default
- **Always anchor limits explicitly:** `limits = c(0, max(obj$gene))` using the full object max
  so zero always maps to lightgrey regardless of which cells are plotted

### Cell plotting
- **Plot ALL cells**, sorted ascending by expression (`arrange(gene)`) so high-expressing cells land on top
- **Do NOT filter to expressing cells only** (e.g. `gene > 0`) — this shifts the color minimum
  off zero and inflates apparent expression in low-expressing groups

### Point sizes
- UMAP point size: `pt.size = 0.3` (cell type / categorical)
- For expression plots on large datasets: use `Seurat:::AutoPointSize()` pattern
  (see feature_umap_plot.md in v1 for the validated implementation)
- **Always `raster = FALSE`** — rasterized points at small sizes render nearly invisible

---

## Points and Scatter

- Dot plot dot size range: `c(0.3, 6)` for circles; `c(0.8, 5.5)` for diamonds
- Dot plot stroke: `0.32` (circles) / `0.40` (diamonds), `color = "grey30"` / `"grey25"`
- Significance indicator dots (e.g., pseudotime median): `size = 2.5, color = "#E53935"`
