#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(copykat)
})

# -----------------------------
# Command-line options
# -----------------------------
option_list <- list(
  make_option("--input_h5ad",type = "character",help = "Path to input h5ad object "),
  make_option("--prefix",type = "character",help = "Prefix for all output files"),
  make_option("--genome",type = "character",default = "hg20",help = "Reference genome: hg20 (hg38) or mm10"),
  make_option("--n_cores",type = "integer",help = "Number of CPU cores for CopyKAT")
)

opts <- parse_args(OptionParser(option_list = option_list))

# Check inputs
stopifnot(file.exists(opts\$input_h5ad))
stopifnot(!is.null(opts\$prefix))

# -----------------------------
# Load h5ad object
# -----------------------------
adata <- read_h5ad(opts\$input_h5ad)
seu <- adata\$as_Seurat()


# -----------------------------
# Extract raw counts
# -----------------------------
exp.rawdata <- as.matrix(GetAssayData(seu, layer="counts"))

# -----------------------------
# Run CopyKAT
# -----------------------------
copykat_results <- copykat(
  rawmat = exp.rawdata,
  id.type = "symbol",
  sam.name = opts\$prefix,  # Use prefix for all outputs
  genome = opts\$genome,
  n.cores = opts\$n_cores
)

# -----------------------------
# Save CopyKAT results 
# -----------------------------
# if needed uncomment
#copykat_rds_file <- paste0(opts\$prefix, "_copykat_results.rds")
#saveRDS(copykat_results, file = copykat_rds_file)

# -----------------------------
# Add predictions back to Seurat object
# -----------------------------
pred_df <- data.frame(
  #cell_id = copykat_results\$prediction\$cell.names,
  copykat_pred = copykat_results\$prediction\$copykat.pred
)
rownames(pred_df) <- copykat_results\$prediction\$cell.names

seurat_obj <- AddMetaData(seu, metadata = pred_df)


# -----------------------------
# Save updated Seurat object
# -----------------------------
#seurat_out_file <- paste0(opts\$prefix, "_add_copykat_to_seurat.rds")
#saveRDS(seurat_obj, file = seurat_out_file)

adata_processed <- as_AnnData(seurat_obj)
adata_file_name <- paste0(opts\$prefix, ".h5ad")
write_h5ad(adata_processed, adata_file_name)

# Write version information
writeLines(
    c(
        '"${task.process}":',
        paste('    r:', paste(R.Version()\$major, R.Version()\$minor, sep = ".")),
        paste('    copykat:', as.character(packageVersion('copykat'))),
        paste('    optparse:', as.character(packageVersion('infercnv'))),
        paste('    Seurat:', as.character(packageVersion('Seurat')))
    ),
    'versions.yml'
)

