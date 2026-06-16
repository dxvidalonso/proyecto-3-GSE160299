#!/usr/bin/env Rscript
# =============================================================================
# GSE160299 DESeq2 Workflow — Step 7: DESeq2 model
# =============================================================================

source("00_configuration.R")
source("01_packages.R")
source("02_plot_helpers.R")
source("03_geo_download.R")
source("04_metadata_parsing.R")
source("05_count_matrix_parsing.R")
source("06_metadata_library_qc.R")

message_step("Constructing DESeq2 object...")

dds0 <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ condition
)

keep_genes <- rowSums(DESeq2::counts(dds0) >= MIN_COUNT) >= MIN_SAMPLES
dds <- dds0[keep_genes, ]
message_step("Genes retained after expression filter: ", nrow(dds))

# Use batch in the design only if it exists and the model matrix is full-rank.
has_real_batch <- nlevels(colData(dds)$batch) > 1 && !all(colData(dds)$batch == "not_available")
design_formula <- ~ condition
batch_adjusted <- FALSE

if (ADJUST_FOR_BATCH_IF_POSSIBLE && has_real_batch) {
  mm_try <- try(model.matrix(~ batch + condition, data = as.data.frame(colData(dds))), silent = TRUE)
  if (!inherits(mm_try, "try-error") && qr(mm_try)$rank == ncol(mm_try)) {
    design_formula <- ~ batch + condition
    batch_adjusted <- TRUE
  } else {
    warning("A batch-like column was detected but is confounded or rank-deficient. Using ~ condition only.")
  }
}

design(dds) <- design_formula
message_step("DESeq2 design formula: ", paste(deparse(design_formula), collapse = ""))

dds <- DESeq2::DESeq(dds)
vsd <- DESeq2::vst(dds, blind = FALSE)

norm_counts <- DESeq2::counts(dds, normalized = TRUE)
readr::write_csv(as.data.frame(norm_counts) %>% tibble::rownames_to_column("gene_id"),
                 file.path(RESDIR, paste0(GSE_ID, "_DESeq2_normalized_counts.csv")))

message("[07_deseq2_model] DESeq2 model fitted and normalized counts saved.")
