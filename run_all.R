#!/usr/bin/env Rscript
# =============================================================================
# GSE160299 DESeq2 Workflow — MASTER SCRIPT
# Runs all steps in order by sourcing each individual script.
#
# Steps:
#   00_configuration.R          — Global parameters and output directories
#   01_packages.R               — Package installation and loading
#   02_plot_helpers.R           — Plot theme and helper functions
#   03_geo_download.R           — Download GEO metadata and supplementary files
#   04_metadata_parsing.R       — Parse and infer condition/batch from metadata
#   05_count_matrix_parsing.R   — Read, clean, and align the count matrix
#   06_metadata_library_qc.R    — Library size and metadata QC plots
#   07_deseq2_model.R           — Fit DESeq2 model and compute VST
#   08_geo2r_qc_plots.R         — GEO2R-style QC plots (boxplot, PCA, UMAP, etc.)
#   09_differential_expression.R — DE results, volcano, MA, p-value diagnostics
#   10_heatmaps_gene_profiles.R — Heatmap of top genes and per-sample profiles
#   11_run_manifest_session_info.R — Output manifest, run summary, sessionInfo
# =============================================================================

# Ensure the working directory is the folder containing these scripts.
# If running from a different location, uncomment and adjust:
# setwd("/path/to/GSE160299_scripts")

cat("========================================================\n")
cat(" GSE160299 DESeq2 Workflow — MASTER SCRIPT\n")
cat(" Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

steps <- c(
  "00_configuration.R",
  "01_packages.R",
  "02_plot_helpers.R",
  "03_geo_download.R",
  "04_metadata_parsing.R",
  "05_count_matrix_parsing.R",
  "06_metadata_library_qc.R",
  "07_deseq2_model.R",
  "08_geo2r_qc_plots.R",
  "09_differential_expression.R",
  "10_heatmaps_gene_profiles.R",
  "11_run_manifest_session_info.R"
)

for (step in steps) {
  cat(sprintf("\n--- Running: %s ---\n", step))
  source(step)
}

cat("\n========================================================\n")
cat(" All steps completed successfully.\n")
cat(" Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n")
