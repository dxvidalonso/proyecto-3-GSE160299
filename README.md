# GSE160299 — DESeq2 Differential Expression Analysis

Reproducible RNA-seq analysis pipeline for the GEO dataset [GSE160299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160299), comparing **Parkinson's Disease (PD)** vs **healthy control (NC)** samples using DESeq2.

---

## Repository structure

```
GSE160299_scripts/
├── 00_configuration.R               # Global parameters and output directories
├── 01_packages.R                    # Package installation and loading
├── 02_plot_helpers.R                # ggplot2 theme and helper functions
├── 03_geo_download.R                # Download GEO metadata and supplementary files
├── 04_metadata_parsing.R            # Parse sample condition and detect batch
├── 05_count_matrix_parsing.R        # Read, clean, and align the count matrix
├── 06_metadata_library_qc.R         # Library size and metadata QC plots
├── 07_deseq2_model.R                # Fit DESeq2 model and compute VST
├── 08_geo2r_qc_plots.R              # QC plots: boxplot, PCA, UMAP, heatmap
├── 09_differential_expression.R     # DE results, volcano, MA, p-value diagnostics
├── 10_heatmaps_gene_profiles.R      # Top-gene heatmap and per-sample profiles
├── 11_run_manifest_session_info.R   # Output manifest, run summary, sessionInfo
└── run_all.R                        # ▶ Master script — runs all steps in order
```

---

## Quickstart

```r
# 1. Clone the repository
# git clone https://github.com/<your-user>/GSE160299-DESeq2.git
# cd GSE160299-DESeq2

# 2. Open R (≥ 4.2) and set the working directory to the scripts folder
setwd("GSE160299_scripts")

# 3. Run the full pipeline
source("run_all.R")
```

All packages are installed automatically on first run if `INSTALL_MISSING <- TRUE` in `00_configuration.R`.

---

## Key parameters (`00_configuration.R`)

| Parameter | Default | Description |
|---|---|---|
| `GSE_ID` | `"GSE160299"` | GEO accession to download |
| `TEST_LEVEL` | `"PD"` | Test group label |
| `REF_LEVEL` | `"NC"` | Reference group label |
| `ALPHA` | `0.05` | FDR significance threshold |
| `LFC_CUTOFF` | `1.0` | Minimum absolute log2 fold change |
| `MIN_COUNT` | `10` | Minimum count for the expression filter |
| `MIN_SAMPLES` | `2` | Minimum samples passing `MIN_COUNT` |
| `ADJUST_FOR_BATCH_IF_POSSIBLE` | `TRUE` | Include batch in design if non-confounded |
| `TOP_HEATMAP_GENES` | `50` | Genes shown in the top-DE heatmap |
| `TOP_PROFILE_GENES` | `12` | Genes shown in profile plots |

---

## Outputs

All outputs are written to `GSE160299_DESeq2_analysis/` (created automatically):

```
GSE160299_DESeq2_analysis/
├── data_raw/          # Downloaded GEO supplementary files
├── results/           # CSV tables (DE results, metadata, run summary)
└── figures/           # PDF and PNG figures
```

### Main result tables

| File | Description |
|---|---|
| `*_all_genes.csv` | Full DESeq2 results for every tested gene |
| `*_significant_*.csv` | Significant genes (FDR < α and \|log2FC\| ≥ cutoff) |
| `*_DEG_summary.csv` | Count of up / down / not-significant genes |
| `*_normalized_counts.csv` | DESeq2 size-factor normalized counts |
| `*_sample_metadata_final.csv` | Sample metadata used in the model |
| `run_summary.csv` | Single-row run parameters and result counts |
| `output_manifest.csv` | Description of every figure produced |
| `sessionInfo.txt` | R and package versions for reproducibility |

### Figures (PDF + PNG)

| # | Figure | Description |
|---|---|---|
| 01 | `metadata_group_counts` | Sample count per group |
| 02 | `metadata_overview_tiles` | Metadata field overview |
| 03 | `qc_library_sizes` | Raw library sizes |
| 04 | `qc_detected_genes` | Detected genes per sample |
| 05 | `geo2r_style_boxplot` | Expression distribution boxplot |
| 06 | `geo2r_style_density` | Expression density curves |
| 07 | `pca_by_condition` | PCA colored by condition |
| 08–09 | `pca_by_batch / pca_after_removal` | PCA batch diagnostics (if batch detected) |
| 10 | `geo2r_style_umap` | UMAP sample relationship |
| 11 | `sample_distance_heatmap` | Sample-to-sample distances |
| 12 | `deseq2_dispersion_estimates` | DESeq2 dispersion fit |
| 13 | `mean_variance_trend` | Mean-variance trend |
| 14 | `deg_counts_summary` | Up / down / NS gene counts |
| 15 | `volcano` | Volcano plot |
| 16 | `MD_MA_plot` | MA / mean-difference plot |
| 17 | `raw_pvalue_histogram` | Raw p-value histogram |
| 18 | `adjusted_pvalue_histogram` | Adjusted p-value histogram |
| 19 | `qq_plot_pvalues` | Q-Q plot of p-values |
| 20 | `top_de_genes_heatmap` | Top DE gene heatmap (z-scored VST) |
| 21 | `gene_profile_top_genes` | Per-sample expression profiles |

---

## Dependencies

| Package | Source | Purpose |
|---|---|---|
| `DESeq2` | Bioconductor | Differential expression |
| `limma` | Bioconductor | Batch effect removal (visualization only) |
| `GEOquery` | Bioconductor | GEO data download |
| `org.Hs.eg.db` | Bioconductor | Gene symbol annotation |
| `AnnotationDbi` | Bioconductor | Annotation database interface |
| `apeglm` | Bioconductor | LFC shrinkage (optional but recommended) |
| `tidyverse` | CRAN | Data wrangling and plotting |
| `ggrepel` | CRAN | Non-overlapping plot labels |
| `pheatmap` | CRAN | Heatmaps |
| `uwot` | CRAN | UMAP embedding |
| `matrixStats` | CRAN | Fast row/column statistics |
| `scales` | CRAN | Axis label formatting |
| `janitor` | CRAN | Column name cleaning |
| `ggvenn` | CRAN | Venn diagrams (optional) |

Tested with **R ≥ 4.2** and **Bioconductor ≥ 3.16**.

---

## Citation

If you use this pipeline, please cite:

- **DESeq2:** Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. https://doi.org/10.1186/s13059-014-0550-8
- **GEO dataset:** [GSE160299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160299)

---

## License

MIT
