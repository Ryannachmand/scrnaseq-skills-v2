---
# Lab Context and Conventions — v2
---

# Lab Context and Conventions

Lab-level defaults for all analyses. Read by the deployment agent when generating CLAUDE.md.
Per-project overrides are in @context/validated_examples.yaml — always check that file for
project-specific palettes, column names, and gene sets before generating scripts.

---

## Default Organism

```yaml
default_organism: Homo sapiens
```

---

## Tissues and Systems Studied

Primary tissues and cell types this lab routinely works with.

```yaml
primary_tissues:
  # Replace with tissues relevant to your lab's work.
  # These drive default marker selection and palette expectations.
  - "<TISSUE_1>"    # e.g., a primary organ, tissue type, or disease context
  - "<TISSUE_2>"

cell_compartments_of_interest:
  # Replace with cell types routinely studied in your lab.
  # These inform module defaults for subclustering, DE, and downstream analyses.
  - "<CELL_TYPE_1>"
  - "<CELL_TYPE_2>"
  # Add as many entries as needed.
```

---

## Standard QC Thresholds

Leave blank to have the pipeline suggest thresholds from the data distribution.
Override in the analysis brief for specific projects.

```yaml
default_qc_thresholds:
  nFeature_RNA_min:      # data-driven suggestion
  nFeature_RNA_max:      # data-driven suggestion
  nCount_RNA_min:        # data-driven suggestion
  nCount_RNA_max:        # data-driven suggestion
  percent_mt_max:        # data-driven suggestion
```

---

## Label Harmonization

Canonical cell type vocabulary used across projects. Controlled vocabulary
reduces annotation drift between projects. Keys are the canonical label to
use in all downstream analysis; values list known aliases.

```yaml
unified_label_vocabulary:
  # Populate from label_harmonization.md controlled vocabulary as projects are completed.
  # See label_harmonization.md (Phase 2 module) for the full controlled vocabulary rules.
  # TODO: human review needed to populate from validated project runs.
```

---

## Pipeline Conventions

Defaults applied to all jobs unless overridden in the brief.

```yaml
batch_correction_var: sample_id     # lab default — column passed to RunHarmony
                                     # IMPORTANT: this is sample_id at the lab level.
                                     # Individual projects may override this in their
                                     # context_defaults in validated_examples.yaml.
                                     # (some projects use source_file — see per-project registry)

whole_object_defaults:
  n_variable_features: 4000
  n_pcs: 30
  clustering_resolution: 0.5

subset_defaults:
  n_variable_features: 2850
  n_pcs: 40
  clustering_resolution: 0.39

downsample_ceiling:                  # e.g. 5000 (cells per sample) — blank = no ceiling
downsample_strategy: sample_level    # sample_level or dataset_level
```

---

## Subclustering Conventions

### Minimum cluster size

Before running Mode A autonomous annotation, drop clusters with fewer than 50 cells:
  cells_keep <- names(obj$seurat_clusters)[
    obj$seurat_clusters %in% names(table(obj$seurat_clusters))[table(obj$seurat_clusters) >= 50]
  ]
  obj <- obj[, cells_keep]

Clusters under 50 cells are likely noise or doublet aggregates. Annotating them adds
spurious labels and inflates the subtype count. Log dropped clusters to decision_log.txt.

Expected subtype count varies by cell type; typical range is 4-12 for a multi-sample dataset. If FindClusters
returns > 12 clusters at the requested resolution, evaluate whether resolution should
be lowered before proceeding with Mode A.

### Subset embedding recompute

Subclustering always recomputes variable features, PCA, Harmony, and UMAP on the
cell-type subset. After `whole[, subset_cells]`, the resulting object retains
parent-object pca/harmony/umap embeddings for those cells; the Step 1 re-cluster
block in @modules/celltype_subclustering.md overwrites all of them.

Doublet removal via run_doublet_removal() (@primitives/doublet_removal.md) runs
on the subset BEFORE normalization, catching doublets that were not resolved in
whole-object space.

Plotting or caching a subsetted object before RunUMAP has completed on the subset
will silently display parent-object geometry.

### Sample-cluster proportion barplot -- always required

The labeled UMAP and the sample-cluster proportion barplot are ALWAYS produced as a pair.
The proportion barplot (CellType_proportion_by_sample.pdf) is not optional -- it is the
standard diagnostic that reveals batch effects, sample imbalance, and over-clustering.

Enforcement: After generating the labeled UMAP endpoint file, always call
make_proportion_plot() for both:
  - CellType_proportion_by_group.pdf (x = GROUP_COL, fill = subtype)
  - CellType_proportion_by_sample.pdf (x = SAMPLE_COL, grouped by GROUP_COL, fill = subtype)

These are already specified in @modules/celltype_subclustering.md Mode A endpoint section
(lines 427-459). This rule exists to prevent them from being silently skipped.

---

## Aesthetics and Color Preferences

Canonical palettes are defined in @context/color_palettes.md.
Project-specific palettes (cell subtypes, tissue types) are in
each project's `context_defaults.palettes` block in validated_examples.yaml.

The deployment agent should never hardcode project-specific color vectors in a
shared module or pipeline file. Always reference palettes via context injection.

---

## Cross-Dataset Comparison Conventions

When comparing in-house data against a public atlas (e.g. Tabula Sapiens):

### Column ordering in cross-dataset dotplots

**Rule:** In-house data ALWAYS appears on the LEFT, atlas data on the RIGHT.

Implement by setting SOURCES_ORDER in @modules/cross_dataset_dotplot.md:
  SOURCES_ORDER <- c("inhouse_label", "atlas_label")

The in-house label (e.g. "inhouse") must be the first element.
This rule applies regardless of cell count -- do not let the atlas dominate column
order simply because it has more cells.

### Depth correction requirement

Cross-dataset dotplots MUST use within-dataset z-scoring (depth correction) per
@modules/cross_dataset_dotplot.md Part B. Never combine raw avg_exp across datasets
without z-scoring first -- atlas cells typically have higher absolute expression
values due to sequencing depth differences.

### Required module

All cross-dataset dotplots use @modules/cross_dataset_dotplot.md.
Never call DotPlot(group.by="organ") directly across datasets -- this skips depth
correction and destroys column ordering.

---

## Functional Gene Categories

These are example functional gene sets that can be referenced across the library.
Used for: cross-context dotplots, subcluster dotplots split by Condition or Group,
module scoring.

These are REFERENCE UNIVERSES for filtering and scoring, NOT plot-ready panels.
For plotting, apply the "Dotplot gene selection rule" below to subset to genes
showing statistical change relevant to the plot's stratification axis.

Genes may appear in multiple categories (e.g. THBS1-4 in Angiocrine and ECM;
SPARC/SPARCL1 in both; BMP2/4/6/7 in both). This is biologically correct --
these molecules are pleiotropic. Do not deduplicate across categories.
If module scoring multiple categories together is desired, note that genes
in multiple categories will contribute to each module's score.

### Vascular Signaling Factors

```r
angiocrine_factors <- c(
  "VEGFA", "VEGFB", "VEGFC", "VEGFD", "PGF",
  "ANGPT1", "ANGPT2", "ANGPT4", "ANGPTL1", "ANGPTL2", "ANGPTL3", "ANGPTL4",
  "FGF1", "FGF2", "FGF7", "FGF9", "FGF10", "FGF18", "FGF21", "FGF23",
  "DLL1", "DLL3", "DLL4", "JAG1", "JAG2",
  "WNT1", "WNT2", "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A", "WNT5B",
  "WNT7A", "WNT7B", "WNT9A", "WNT9B", "WNT10A", "WNT10B", "WNT11", "WNT16",
  "BMP2", "BMP4", "BMP6", "BMP7", "GDF11", "GDF15",
  "PDGFA", "PDGFB", "PDGFC", "PDGFD",
  "HGF", "IGF1", "IGF2",
  "CTGF", "CYR61", "NOV",
  "SEMA3A", "SEMA3C", "SEMA3E", "SEMA3F", "SEMA4A", "SEMA4D", "SEMA4F",
  "EFNA1", "EFNA2", "EFNA3", "EFNA4", "EFNA5",
  "EFNB1", "EFNB2", "EFNB3",
  "NTN1", "NTN4", "SLIT2", "SLIT3",
  "APLN", "APELA",
  "EDN1", "EDN2", "EDN3",
  "ADM", "ADM2",
  "NRG1", "NRG2", "NRG3", "NRG4",
  "THBS1", "THBS2", "THBS3", "THBS4",
  "SPARC", "SPARCL1",
  "RSPO1", "RSPO2", "RSPO3", "RSPO4",
  "KITLG", "FLT3LG"
)
```

### Cytokines & Immune Modulatory

Union of surface receptors and cytokines/chemokines. Single category per lab
convention. Captures cell-surface signaling receptors AND secreted ligands
in one biology-oriented set.

```r
cytokines_immune_modulatory <- c(
  # --- Surface receptors ---
  "IL1R1", "IL1R2", "IL1RL1", "IL1RL2", "IL2RA", "IL2RB", "IL2RG",
  "IL3RA", "IL4R", "IL5RA", "IL6R", "IL6ST", "IL7R", "IL9R",
  "IL10RA", "IL10RB", "IL11RA", "IL12RB1", "IL12RB2", "IL13RA1", "IL13RA2",
  "IL15RA", "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE",
  "IL18R1", "IL18RAP", "IL20RA", "IL20RB", "IL21R", "IL22RA1", "IL22RA2",
  "IL23R", "IL27RA", "IL31RA",
  "ICAM1", "ICAM2", "ICAM3", "ICAM4", "ICAM5",
  "VCAM1", "PECAM1", "MCAM", "ALCAM", "NCAM1", "NCAM2",
  "CADM1", "CADM2", "CADM3", "CADM4",
  "ESAM", "JAM2", "JAM3", "JAML",
  "SELE", "SELP", "SELL",
  "ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "ITGA7", "ITGA8", "ITGA9",
  "ITGAV", "ITGAE", "ITGAL", "ITGAM", "ITGAX",
  "ITGB1", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8",
  "TGFBR1", "TGFBR2", "TGFBR3", "ACVR1", "ACVR1B", "ACVR1C", "ACVR2A", "ACVR2B",
  "BMPR1A", "BMPR1B", "BMPR2",
  "EDNRA", "EDNRB", "EGFR", "FGFR1", "FGFR2", "FGFR3", "FGFR4",
  "KDR", "FLT1", "FLT4", "TEK", "TIE1",
  "PDGFRA", "PDGFRB", "MET", "AXL", "MERTK", "TYRO3",
  "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
  "EPHA1", "EPHA2", "EPHA3", "EPHA4", "EPHB1", "EPHB2", "EPHB3", "EPHB4",
  "CXCR1", "CXCR2", "CXCR3", "CXCR4", "CXCR5", "CXCR6",
  "CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCR8", "CCR9", "CCR10",
  "ACKR1", "ACKR2", "ACKR3", "ACKR4",
  "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2",
  "TNFRSF1A", "TNFRSF1B", "NGFR", "LIFR", "OSMR", "CNTFR",
  "CSF1R", "CSF2RA", "CSF2RB", "CSF3R",
  # --- Cytokines, chemokines, growth factors ---
  "IL1A", "IL1B", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL8", "IL9",
  "IL10", "IL11", "IL12A", "IL12B", "IL13", "IL14", "IL15", "IL16", "IL17A", "IL17B",
  "IL17C", "IL17D", "IL17F", "IL18", "IL19", "IL20", "IL21", "IL22", "IL23A", "IL24",
  "IL25", "IL26", "IL27", "IL31", "IL32", "IL33", "IL34", "IL36A", "IL36B", "IL36G",
  "IL37", "EBI3",
  "IFNA1", "IFNA2", "IFNB1", "IFNG", "IFNL1", "IFNL2", "IFNL3",
  "TNF", "LTA", "LTB", "TNFSF4", "TNFSF8", "TNFSF9", "TNFSF10", "TNFSF11",
  "TNFSF12", "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", "TNFSF18",
  "TGFB1", "TGFB2", "TGFB3", "BMP2", "BMP4", "BMP6", "BMP7",
  "GDF15", "MSTN", "INHBA", "INHBB", "NODAL",
  "CSF1", "CSF2", "CSF3",
  "LIF", "OSM", "CTF1", "CNTF", "KITLG", "FLT3LG", "THPO", "EPO",
  "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXCL9", "CXCL10",
  "CXCL11", "CXCL12", "CXCL13", "CXCL14", "CXCL16", "CXCL17",
  "CCL1", "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8", "CCL11", "CCL13",
  "CCL14", "CCL15", "CCL16", "CCL17", "CCL18", "CCL19", "CCL20", "CCL21",
  "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28",
  "CX3CL1", "XCL1", "XCL2"
)
```

### Extracellular Matrix (ECM)

```r
extracellular_matrix <- c(
  # --- Collagens ---
  "COL1A1", "COL1A2", "COL2A1", "COL3A1",
  "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6",
  "COL5A1", "COL5A2", "COL5A3",
  "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6",
  "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3",
  "COL10A1", "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1",
  "COL16A1", "COL17A1", "COL18A1", "COL19A1", "COL20A1", "COL21A1", "COL22A1",
  "COL23A1", "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1",
  # --- Laminins ---
  "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5",
  "LAMB1", "LAMB2", "LAMB3", "LAMB4",
  "LAMC1", "LAMC2", "LAMC3",
  # --- Fibronectin and elastin ---
  "FN1", "FNDC1", "FNDC3A", "FNDC3B", "FNDC4", "FNDC5",
  "ELN", "EMILIN1", "EMILIN2", "EMILIN3",
  # --- Proteoglycans ---
  "DCN", "BGN", "LUM", "FMOD", "PRELP", "OMD", "OGN", "ASPN", "ECM2", "KERA",
  "VCAN", "ACAN", "NCAN", "BCAN",
  "HSPG2",
  "SDC1", "SDC2", "SDC3", "SDC4",
  "GPC1", "GPC2", "GPC3", "GPC4", "GPC5", "GPC6",
  # --- Matricellular ---
  "THBS1", "THBS2", "THBS3", "THBS4", "THBS5",
  "SPARC", "SPARCL1", "SPOCK1", "SPOCK2", "SPOCK3",
  "TNC", "TNN", "TNR", "TNXB",
  "POSTN", "TGFBI",
  "FBLN1", "FBLN2", "FBLN5", "FBLN7",
  "MFAP1", "MFAP2", "MFAP3", "MFAP4", "MFAP5",
  "LTBP1", "LTBP2", "LTBP3", "LTBP4",
  "NID1", "NID2",
  # --- MMPs / TIMPs / ADAMTSs ---
  "MMP1", "MMP2", "MMP3", "MMP7", "MMP8", "MMP9", "MMP10", "MMP11", "MMP12", "MMP13",
  "MMP14", "MMP15", "MMP16", "MMP17", "MMP19", "MMP20", "MMP21", "MMP23B", "MMP24",
  "MMP25", "MMP26", "MMP27", "MMP28",
  "TIMP1", "TIMP2", "TIMP3", "TIMP4",
  "ADAMTS1", "ADAMTS2", "ADAMTS3", "ADAMTS4", "ADAMTS5", "ADAMTS6", "ADAMTS7",
  "ADAMTS8", "ADAMTS9", "ADAMTS10", "ADAMTS12", "ADAMTS13", "ADAMTS14",
  "ADAMTS15", "ADAMTS16", "ADAMTS17", "ADAMTS18", "ADAMTS19", "ADAMTS20",
  # --- Crosslinking and other ---
  "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4",
  "VTN", "PCOLCE", "PCOLCE2",
  "SERPINE1", "SERPINE2", "SERPINH1",
  "CTSK", "CTSB", "CTSL", "CTSS"
)
```

---

## Dotplot Gene Selection Rule

When generating a dotplot from any functional gene category (or any other large
reference gene list), DO NOT plot the full reference list. Subset to genes
showing statistical change along the plot's stratification axis.

### Filter

For each gene in the reference list, include in the plot if:

```
p_val_adj < 0.1  AND  abs(avg_log2FC) > 0.25
```

These thresholds correspond to "significant or near-significant" by lab
convention. The logFC floor excludes biologically trivial changes that
happen to be statistically significant in well-powered comparisons.

### DE source -- matches the plot's stratification axis

The DE results used for filtering must match the biological question the
plot is asking:

- Plot stratified by Condition (Tumor vs Normal) -> filter by Condition DE
- Plot stratified by Group (Normal, ConditionA, ConditionB, ConditionC) -> filter by
  per-Group DE vs Normal (union: include a gene if it passes the filter
  for ANY non-Normal group vs Normal)
- Plot stratified by cell type or subcluster -> filter by per-cell-type DE
  for the relevant contrast
- Cross-organ plot -> filter by organ-vs-rest DE for at least one organ

If multiple plots are produced from the same gene category (e.g. one by
Condition and one by Group), each uses its own gene selection -- they will
typically share some genes and differ on others. This is correct.

### Floor

If fewer than 4 genes pass the filter, show the top 10 genes from the
reference list by abs(avg_log2FC), regardless of significance. Annotate
the plot title or caption noting that the filter floor was triggered (e.g.
"showing top 10 by |logFC|; <4 genes passed p_adj<0.1, |logFC|>0.25").

### Ceiling

If more than 30 genes pass the filter, show the top 30 by abs(avg_log2FC).
Annotate the plot title or caption noting truncation (e.g. "top 30 of N
genes passing filter").

### Single-direction signature lists

If the input is already a pre-filtered signature gene list (e.g. signature
derivation output, custom curated list), the filter does not apply --
plot all genes in the signature list directly, up to the 30-gene ceiling.

### What to do with the full reference list

The full reference lists ARE used for:
- AddModuleScore (use the full list -- breadth is the goal)
- Pathway enrichment universes
- DE result filtering (when reporting which genes in a category changed)
- Heatmaps with > 30 genes that include only filtered-significant genes
  but where individual gene labels are not the focus

The filter rule applies primarily to DOTPLOTS where each gene is
individually labeled and the goal is interpretability of specific genes.

---

## Notes for the Deployment Agent

1. Always check validated_examples.yaml for the target project before generating CLAUDE.md.
   The per-project `context_defaults` block supersedes lab-level defaults for that project.

2. `batch_correction_var` at the lab level is `sample_id`. If the project brief uses a
   different column name (e.g., `source_file` or `sample` depending on the project),
   that override is in validated_examples.yaml `context_defaults.batch_correction_var`.

3. Color palettes for specific tissues/cell types are in validated_examples.yaml under
   each project entry. Do NOT use color palettes from one project for another project.

4. For returning-to-an-old-project scenarios: read the project's validated_examples.yaml entry
   completely before generating any scripts — it contains the key column names, gotchas, and
   analysis decisions from previous runs.
