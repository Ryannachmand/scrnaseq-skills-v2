---
# Module: Loading Heterogeneous Data Formats
# Pipeline: IntegratePublicData
# Migrated from ~/claude-skills/pipelines/IntegratePublicData/methods/load_formats.md
# Phase 2 Group A — parameterized; project-specific values removed
requires_context:
  palettes: []
  metadata_columns:
    required: []
    optional:
      - sample_col      # column to assign sample identity (default: sample_id)
  brief_keys:
    required:
      - inputs.paths    # list of file paths or GEO accessions
      - inputs.type     # h5 | mex | rds | geo
    optional:
      - project_name    # used in CreateSeuratObject project= argument
references:
  - "@primitives/seurat_v5_rules.md"
  - "@primitives/r_environment.md"
---

# Module: Loading Heterogeneous Data Formats

Recipes for every supported input format in the IntegratePublicData pipeline.
Run this stage before any QC or merge operations.

See @primitives/seurat_v5_rules.md for JoinLayers (Rule 1) and dynamic UMAP
reduction name detection (Rule 5) which both appear in recipes below.

---

## Format Decision Tree

Before writing any loading code, identify the format:

```
Is it an RDS file?
  └─ Yes → Load with readRDS(); go to "RDS Objects" section
  └─ No → Is it an HDF5 (.h5) file?
      └─ Yes → Does filename contain "molecule_info"?
          └─ Yes  → SKIP — this is Cell Ranger molecule info, not a count matrix
          └─ No   → Use Read10X_h5()
      └─ No → Is it a directory with barcodes/features/matrix files?
          └─ Yes → Use Read10X(dir)
          └─ No → Are the barcodes/features/matrix files loose (not in a subdir)?
              └─ Yes → Read manually (see "Loose MEX files" section)
              └─ No → Is it a .mtx.gz with separate barcodes/features files?
                  └─ Yes, with spliced/unspliced → Use spliced only (see "Spliced/unspliced")
```

---

## 10x HDF5 Files

```r
# Standard filtered count matrix
counts <- Read10X_h5(filepath)

# Multi-modal h5 (returns a list)
if (is.list(counts)) counts <- counts[["Gene Expression"]]

obj <- CreateSeuratObject(counts, project = dataset_name, min.cells = 3, min.features = 200)
```

**Which h5 files to use from GEO:**
- `*_filtered_feature_bc_matrix.h5` → preferred (Cell Ranger filtered)
- `*_raw_gene_bc_matrices_h5.h5` → use when filtered not available (apply QC after)
- `*_molecule_info.h5` → **SKIP** — not a count matrix
- `*.html.gz` → **SKIP** — QC report

**When a directory mixes file types:** filter by pattern to select the correct h5 files.
Inspect the download directory first; the pattern is dataset-specific. Example:
```r
h5_files <- list.files(dir, pattern = "_filtered_feature_bc_matrix\\.h5$", full.names = TRUE)
# Adjust pattern to match your dataset's actual filenames
```

---

## 10x MEX Directory (Standard)

```r
# Directory contains: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
counts <- Read10X(data.dir = mex_directory)
obj <- CreateSeuratObject(counts, project = dataset_name, ...)
```

**Nested subdirectory structure:** Some GEO deposits nest the MEX files inside a
sample-specific subdirectory. Inspect the download to find the correct path:
```r
# Example layout: {sample_dir}/{version_dir}/barcodes.tsv.gz
mex_dir <- file.path(base_dir, sample_id, "filtered_feature_bc_matrix")
counts   <- Read10X(data.dir = mex_dir)
```

---

## Loose MEX Files (Not in Subdirectory)

When barcodes/features/matrix files share a directory with other files,
`Read10X()` cannot be used directly. Read files manually:

```r
# Inspect downloaded files to get exact filenames before filling these paths
barcodes_path <- file.path(dir, paste0(accession, "_barcodes.tsv.gz"))  # REPLACE with actual filename
features_path <- file.path(dir, paste0(accession, "_features.tsv.gz"))  # REPLACE with actual filename
matrix_path   <- file.path(dir, paste0(accession, "_matrix.mtx.gz"))    # REPLACE with actual filename

barcodes <- read.table(gzcon(file(barcodes_path, "rb")), header = FALSE)$V1
features <- read.table(gzcon(file(features_path, "rb")), header = FALSE)$V2  # col 2 = gene symbols
mat      <- readMM(gzcon(file(matrix_path, "rb")))
rownames(mat) <- features
colnames(mat) <- barcodes

obj <- CreateSeuratObject(mat, project = dataset_name, ...)
```

Alternative: copy the three files into a temporary subdirectory with standard names,
then use `Read10X()`:
```r
tmp <- tempdir()
file.copy(c(barcodes_path, features_path, matrix_path), tmp)
# rename to standard names if needed: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
counts <- Read10X(tmp)
```

---

## Spliced/Unspliced Matrix (scVelo format)

Some GEO deposits provide RNA velocity matrices instead of standard count matrices.
**Always use the spliced matrix only.** The unspliced matrix is not appropriate for
standard Seurat analysis.

```r
library(Matrix)

# Inspect downloaded files to get exact filenames
barcodes_gz <- file.path(dir, paste0(accession, "_barcodes.tsv.gz"))  # REPLACE
features_gz <- file.path(dir, paste0(accession, "_features.tsv.gz"))  # REPLACE
spliced_gz  <- file.path(dir, paste0(accession, "_spliced.mtx.gz"))   # REPLACE

barcodes <- read.table(gzcon(file(barcodes_gz, "rb")), header = FALSE)$V1
features <- read.table(gzcon(file(features_gz, "rb")), header = FALSE)$V2

# CRITICAL: readMM returns cells×genes for some deposits — must check and transpose
mat <- readMM(gzcon(file(spliced_gz, "rb")))
cat("Matrix dimensions before transpose check:", dim(mat), "\n")

# If nrow << ncol, matrix is cells×genes → transpose to genes×cells
if (nrow(mat) < ncol(mat)) {
  mat <- t(mat)
  cat("Transposed to", dim(mat), "\n")
}

rownames(mat) <- features
colnames(mat) <- barcodes

obj <- CreateSeuratObject(mat, project = dataset_name, min.cells = 3, min.features = 200)
```

**Note on reading binary files:** `.mtx.gz` requires `gzcon(file(..., "rb"))`.
Text files (barcodes, features) can use either `gzcon(file(..., "rb"))` or
`read.table(..., header=FALSE)` with R's built-in gz decompression.

---

## RDS Objects — In-House

```r
so <- readRDS(filepath)
# Immediately inspect:
cat("Dims:", dim(so), "\n")
cat("Assays:", Assays(so), "\n")
cat("Reductions:", names(so@reductions), "\n")
cat("Layers:", Layers(so), "\n")
print(head(so@meta.data))
```

### Multi-layer Seurat v5 objects (JoinLayers)

If object has >3 layers (e.g., `counts.1, counts.2, data.1, data.2, ...`),
it was created by merging in Seurat v5 without joining. Fix immediately:
```r
so <- JoinLayers(so, assay = "RNA")
```
See @primitives/seurat_v5_rules.md Rule 1 for full explanation.

### Ensembl ID rownames

Some objects have Ensembl IDs (ENSG...) as gene names instead of symbols.
Detect and convert before merge:

```r
# Detect
if (any(grepl("^ENSG", rownames(so)))) {
  message("Ensembl IDs detected — converting to gene symbols via biomaRt")
  library(biomaRt)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters    = "ensembl_gene_id",
                   values     = rownames(so),
                   mart       = mart)
  # Remove empty symbols; keep first symbol per ENSG
  mapping <- mapping[mapping$hgnc_symbol != "", ]
  mapping <- mapping[!duplicated(mapping$ensembl_gene_id), ]
  # Subset to mappable genes
  keep <- rownames(so) %in% mapping$ensembl_gene_id
  so   <- so[keep, ]
  sym  <- mapping$hgnc_symbol[match(rownames(so), mapping$ensembl_gene_id)]

  # Seurat v5: rename features using RenameFeatures() across all layers
  # TODO: verify RenameFeatures() signature in your installed Seurat v5 version.
  # Preferred (Seurat >= 5.0): RenameFeatures(so, new.names = sym, assay = "RNA")
  so <- RenameFeatures(so, new.names = sym, assay = "RNA")

  # Fallback if RenameFeatures() is unavailable:
  # counts_mat <- GetAssayData(so, assay="RNA", layer="counts"); rownames(counts_mat) <- sym
  # so <- SetAssayData(so, layer="counts", new.data=counts_mat)
  # data_mat <- GetAssayData(so, assay="RNA", layer="data"); rownames(data_mat) <- sym
  # so <- SetAssayData(so, layer="data", new.data=data_mat)

  # VERIFY:
  cat("First 5 rownames after rename:", rownames(so)[1:5], "\n")
}
```
Note: biomaRt requires internet access. Expect ~5 min for large gene sets.
Unmapped ENSG IDs are dropped at the gene intersection step in Stage 5.

### Non-standard reduction names

Objects from external labs or public datasets may store UMAP under a non-standard
reduction name. Never hardcode `reduction = "umap"`.
See @primitives/seurat_v5_rules.md Rule 5 for the dynamic detection pattern.

```r
red_name <- grep("umap", names(so@reductions), ignore.case = TRUE, value = TRUE)[1]
if (is.na(red_name)) stop("No UMAP found. Available: ", paste(names(so@reductions), collapse = ", "))
cat("Using reduction:", red_name, "\n")
# For pre-integration individual plots only; after re-integration, standard name applies
DimPlot(so, reduction = red_name)
```

### CDS (Monocle3) objects

If `readRDS()` returns a `cell_data_set` object, convert to Seurat:
```r
library(SingleCellExperiment)
sce       <- as(cds, "SingleCellExperiment")
counts_mat <- counts(sce)
so        <- CreateSeuratObject(counts_mat, meta.data = as.data.frame(colData(sce)))
```

---

## Sample Demultiplexing from Barcodes

When a single matrix file contains multiple samples (common in GEO deposits):
```r
# Check for sample-encoding prefixes
prefixes <- unique(gsub("_[ACGT]{16}.*", "", colnames(obj)))
cat("Unique prefixes:", length(prefixes), "\n")
print(head(prefixes, 20))

# If meaningful prefixes found: split
if (length(prefixes) > 1 && !all(prefixes == "")) {
  obj$sample_id <- gsub("_[ACGT]{16}.*", "", colnames(obj))
  sample_list   <- SplitObject(obj, split.by = "sample_id")
} else {
  # Treat as single batch
  obj$sample_id <- paste0(dataset_name, "_pooled")
}
```

---

## GEO Bulk Supplementary Files

### GSM vs GSE Accessions

- **GSM** = single sample. Supplementary file is often a Cufflinks tmap or FPKM file, **not** a count matrix.
- **GSE** = series. May be an aggregated expression matrix (xlsx, csv) with samples as columns.
- Always inspect the downloaded filename and column names before assuming format.
  Use `getGEOSuppFiles(accession, baseDir = cache_dir, makeDirectory = TRUE)` — it caches locally,
  so subsequent runs use the cached file without re-downloading.

### Cufflinks tmap Files (`*cuffcmp.transcripts.gtf.tmap.txt.gz`)

These appear when a GEO depositor ran Cufflinks/Cuffcompare. They are **transcript-level**, not gene-level.

**Column names:**
`ref_gene_id, ref_id, class_code, cuff_gene_id, cuff_id, FMI, FPKM, FPKM_conf_lo, FPKM_conf_hi, cov, len, major_iso_id, ref_match_len`

**class_code meanings (critical for interpretation):**
- `=` — assembled transcript exactly matches reference → full-length, highest confidence FPKM
- `c` — assembled transcript is *contained within* reference → partial coverage, typically underestimated FPKM
  A `c` code on a lowly expressed gene does **not** imply novel isoforms — it reflects insufficient
  read depth to assemble the full transcript.

**Gene extraction rule:**
```r
# Prefer the "=" transcript if available; fall back to the row with highest FPKM
get_gene_fpkm <- function(tbl, gene) {
  rows  <- tbl[tbl$ref_gene_id == gene, ]
  if (nrow(rows) == 0) return(NA_real_)
  exact <- rows[rows$class_code == "=", ]
  if (nrow(exact) > 0) return(exact$FPKM[1])   # use full-length transcript
  rows$FPKM[which.max(rows$FPKM)]               # fallback: highest isoform
}
```

**FPKM validation — always check these before trusting the values:**
```r
# Expected ranges in a typical poly-A bulk RNA-seq sample
# GAPDH:  500–2000 FPKM
# ACTB:   1000–5000 FPKM
# B2M:    500–3000 FPKM
# If these are off by orders of magnitude, the FPKM column is not standard.
check_genes <- c("GAPDH", "ACTB", "B2M")
tbl[tbl$ref_gene_id %in% check_genes, c("ref_gene_id", "class_code", "FPKM")]
```
Also check tissue-specific marker genes for your sample type to confirm sample identity.
Example: for liver sinusoidal endothelial cells, markers like CLEC4G and LYVE1 should
be detectably expressed (>100 FPKM). Substitute appropriate markers for your tissue.

**Note on human reference gene names in mouse-aligned samples:**
Cufflinks tmap files always use the reference genome's gene names, regardless of the sample species.
A mouse sample aligned to the human genome will report human gene symbols (not mouse lowercase).
Always verify sample species from the GEO accession description, not from gene name case in the file.

---

### Excel Expression Matrices (`.xlsx`)

When a GSE supplementary file is `.xlsx`, `read.delim()` will fail silently with null-byte errors.

```r
# readxl is required — verify it is installed in r-env
tbl <- readxl::read_xlsx(fpath)
tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
```

Column names frequently encode biology rather than sample IDs — inspect before assuming:
```r
cat("Columns:", paste(colnames(tbl), collapse = ", "), "\n")
# e.g. "SEC.one", "SEC.three", "high.four", "low.one" → cell type groups, not just replicates
```
Group columns by prefix (grep) and summarize within each group rather than averaging all numeric columns.

---

## RDS Objects — HTO-Demultiplexed (No Cell Type Labels)

When a Seurat RDS object has **HTO barcodes as Idents** (e.g., `B0301-ACCCACC...`, `B0302-GGTCGAG...`)
and no cell type annotation column in `meta.data`, the object was hashtag-demultiplexed and
never annotated beyond sample identity.

**Detection:**
```r
id_tab <- head(sort(table(Idents(so)), decreasing = TRUE))
# If Ident names look like "B0301-ACCCACCAGTAAGAC" → HTO barcodes, not cell types
```

**Handling:**
- If the object is described in the brief as a pre-filtered, cell-type-specific object
  (e.g., "BM ECs"), treat all cells as that cell type — assign `so$cell_type <- "project_specific"`.
- Do NOT attempt to re-cluster and re-annotate unless the brief explicitly requires it.
- Check `HTO_classification.global` column: keep only `"Singlet"` cells if cell counts allow.

---

## Post-loading Checklist

After loading any object, verify:
- [ ] `dim(so)` looks reasonable (not 0 cells, not millions of genes)
- [ ] `rownames(so)[1:5]` are gene symbols (not ENSG IDs, not numeric)
- [ ] `so$sample_id` is set and unique values look correct
- [ ] `so$condition` is set correctly from brief
- [ ] `so$data_type` is set correctly from brief
- [ ] No duplicate cell barcodes: `any(duplicated(colnames(so)))`
