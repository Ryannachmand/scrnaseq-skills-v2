---
# Canonical Color Palettes — v2
---

# Canonical Color Palettes

Single authoritative source for all non-project-specific color palettes.
Reference this file rather than hardcoding color values in any primitive or module.

**Project-specific palettes** (EC subtypes, adipose depot types, organ tissue types)
are in `context/validated_examples.yaml` under each project's `context_defaults.palettes` block.
Do NOT define project-specific palettes here.

---

## 1. Diverging Palette (Correlation, log2FC, Fold Change)

6-stop blue→white→red. Use for: correlation heatmaps, z-score heatmaps, log2FC values,
differential abundance log2OR. Always anchor midpoint at 0.

```r
PALETTE_DIVERGING_6 <- c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B")
# Usage: scale_fill_gradientn(colors = PALETTE_DIVERGING_6, values = scales::rescale(c(-3,-1,0,1,2,3)))

# For 3-stop (simpler): steelblue → white → firebrick
PALETTE_DIVERGING_3 <- c("steelblue", "white", "firebrick")
# Usage: scale_fill_gradient2(low="steelblue", mid="white", high="firebrick", midpoint=0)

# Volcano direction colors (endpoint only)
COLOR_UP_IDENT1 <- "#B2182B"    # red — higher in ident1 (positive log2FC)
COLOR_UP_IDENT2 <- "#2166AC"    # blue — higher in ident2 (negative log2FC)
COLOR_NS        <- "grey78"     # non-significant
```

---

## 2. Gene Expression Continuous Scale

For dot plots and violin plots showing scaled average expression.

```r
PALETTE_EXPRESSION <- c("#F5F5F5", "#FFF9C4", "#FFB300", "#E53935")
# white → pale yellow → amber → red
# Usage: scale_fill_gradientn(colors = PALETTE_EXPRESSION)
# Or:   scale_color_gradientn(colors = PALETTE_EXPRESSION)
```

For UMAP feature plots (matches Seurat FeaturePlot default):
```r
PALETTE_FEATURE_UMAP <- c("lightgrey", "blue")
# Usage: scale_color_gradientn(colors = PALETTE_FEATURE_UMAP, limits = c(0, max_val))
# CRITICAL: always anchor limits at c(0, max_val) so zero maps to lightgrey
```

---

## 3. Differential Abundance Heatmap

```r
PALETTE_DIFF_ABUNDANCE <- c("steelblue", "white", "firebrick")
# Usage: scale_fill_gradient2(low="steelblue", mid="white", high="firebrick", midpoint=0)
# Significance threshold: FDR < 0.01 & |log2OR| > 0.58 (OR > 1.5)
# Star label: ifelse(p_adj < 0.01 & abs(log2OR) > 0.58, "*", "")
```

---

## 4. Pseudotime Scale

```r
PALETTE_PSEUDOTIME <- "magma"
# Usage: scale_color_viridis_c(option = "magma")
# NEVER use "inferno" — magma only. This is validated convention.
```

---

## 5. Categorical Default — Okabe-Ito Colorblind-Safe

8-color palette for categorical data where project-specific colors are not defined.

```r
PALETTE_OKABE_ITO <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermilion
  "#CC79A7",  # purple-pink
  "#000000"   # black
)
# Usage: scale_fill_manual(values = PALETTE_OKABE_ITO)
# Or:   scale_color_manual(values = PALETTE_OKABE_ITO)
```

---

## 6. Metadata Track Palettes

For proportion plot metadata annotation strips.

```r
# Age gradient (warm orange)
PALETTE_AGE <- c("#FFF5EB", "#FD8D3C", "#7F2704")

# BMI gradient (blue)
PALETTE_BMI <- c("#EFF3FF", "#6BAED6", "#084594")

# Sex (categorical)
COLOR_SEX_F <- "#E07A9F"
COLOR_SEX_M <- "#45B7D1"
```

---

## 7. Project-Specific Palettes

Project-specific palettes are NOT defined here. They are in `validated_examples.yaml`
under each project's `context_defaults.palettes` block.

Current projects with defined palettes in validated_examples.yaml:
- `HumanFat` — celltype_colors, adipose_type_colors, ec_subtype_colors, organ_tissue_colors
- `KidneyNew` — ec_subtype_colors_kidney
- `BoneMarrowStroma` — TODO: label_colors for STROMA_ORDER (not yet in library)

To add a new project palette: add an entry to `context_defaults.palettes` in
`validated_examples.yaml`. Never add project palettes to this file.

---

## Usage Notes

1. Import these into generated scripts as R variables at the top of the script.
2. When a module or primitive accepts a `group_colors` or `subtype_colors` argument,
   the caller provides values from this palette file OR from the project's context_defaults.
3. The `is_confound()` / `is_ambient()` filter functions in differential_expression.md
   use `COLOR_UP_IDENT1` / `COLOR_UP_IDENT2` for significance coloring — these endpoint
   colors derive from PALETTE_DIVERGING_6 and should not be changed independently.
