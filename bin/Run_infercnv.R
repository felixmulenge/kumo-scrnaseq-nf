#!/usr/bin/env Rscript


# This script runs inferCNV on a Seurat object and export result files:
# - inferCNV heatmap png
# - inferCNV RDS results file
# - wide metadata table with CNVs

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(Matrix)
  library(infercnv)
  library(anndataR)
  library(Seurat)
  library(fs)
})

# set option per inferCNV:
# """Please use "options(scipen = 100)" before running infercnv if you are using
# the analysis_mode="subclusters" option or you may encounter an error while the
# hclust is being generated."""
# analysis_mode="subclusters" is the default, which we do use
options(scipen = 100)

# Command-line arguments
option_list <- list(
  make_option(opt_str = "--input_h5ad",type = "character",default = "",help = "Path to the h5ad file to run inferCNV on"),
  #make_option(opt_str = "--annotation",type = "character",default = "",help = "Path to annotation file"),
  make_option(c("--ref_groups"), type="character", default="", help="Comma-separated reference group names"),
  make_option(opt_str = "--gene_order_file",type = "character",default = "",help = "Path to gene order file (tab-delimited, no header)"),
  make_option(opt_str = "--temp_dir",type = "character",default = "infercnv_tmp",help = "Temporary directory to save inferCNV output"),
  make_option(opt_str = c("-t", "--threads"),type = "integer",help = "Number of threads to use"),
  make_option(opt_str = "--random_seed",type = "integer",default = 2025,help = "Random seed for reproducibility"),
  make_option(opt_str = "--prefix",type = "character",default = NULL,help = "Prefix for all output files")
)


# Parse options
opts <- parse_args(OptionParser(option_list = option_list))

# Check required files
stopifnot(
  "input_h5ad does not exist" = file.exists(opts\$input_h5ad),
  "annotation file does not exist" = file.exists(opts\$annotation),
  "gene_order_file does not exist" = file.exists(opts\$gene_order_file),
  "prefix was not provided" = !is.null(opts\$prefix)
)

# Set seed
set.seed(opts\$random_seed)

# Define output file paths using prefix
opts\$output_rds        <- paste0(opts\$prefix, "_run_final_infercnv_obj.rds")
opts\$output_table      <- paste0(opts\$prefix, "_map_metadata_from_infercnv.tsv")
opts\$output_heatmap    <- paste0(opts\$prefix, "_infercnv.png")
opts\$output_seurat_rds <- paste0(opts\$prefix, "_add_infercnv_to_seurat.rds")


# Read the AnnData objects and convert directly to Seurat
adata <- read_h5ad(opts\$input_h5ad)
seu <- adata\$as_Seurat()


# ---- Extract raw counts matrix (genes x cells) ----
raw_counts_matrix <- GetAssayData(seu, slot = "counts")     #seurat v4


# define relevant infercnv output files for later use/checks
# infercnv will automatically create these files at these hardcoded paths
fs::dir_create(opts\$temp_dir) # ensure output directory exists, to be safe
scratch_infercnv_rds <- file.path(opts\$temp_dir, "run_infercnv_obj")
scratch_metadata_file <- file.path(opts\$temp_dir, "map_metadata_from_infercnv.txt")
scratch_png_file <- file.path(opts\$temp_dir, "infercnv.png")

# create annotations_df table which requires:
# - rownames as cell barcodes
# - a single column with annotation labels; we use "predicted_labels" from celltypist
# see  df_celltypist = predictions_adata.obs.loc[adata.obs.index, ["predicted_labels", "conf_score"] for context

annotations_df <- data.frame(
  annotation = seu@meta.data["predicted_labels"],
  row.names = rownames(seu@meta.data)
)

# ---- Parse reference groups ----
ref_groups <- character(0)
if (nchar(opt\$ref_groups) > 0) {
  ref_groups <- strsplit(opt\$ref_groups, ",", fixed=TRUE)[[1]]
  ref_groups <- trimws(ref_groups)
}
message("Reference groups: ", ifelse(length(ref_groups) == 0, " ", paste(ref_groups, collapse=", ")))     # if ref_group_name not specified append blank!


# Run inferCNV
infercnv_result <- tryCatch({
    # create the inferCNV object and pipe into run
    # so errors from either step are caught
  infercnv::CreateInfercnvObject(
    raw_counts_matrix = raw_counts_matrix,
    annotations_file = opts\$annotation,
    delim = "\t",
    gene_order_file = opts\$gene_order_file,
    ref_group_name = ref_groups
  ) |>
    infercnv::run(
      cutoff = 0.1,       # 10x Genomics
      out_dir = opts\$temp_dir,    # save all intermediate files here
      denoise = TRUE,
      HMM = TRUE,
      HMM_type = "i6",
      save_rds = FALSE,     # don't save intermediate RDS files
      num_threads = opts\$threads
    )
}, error = function(e) {
  message("inferCNV failed; creating empty result files")
  file.create(opts\$output_rds, opts\$output_table, opts\$output_heatmap, opts\$output_seurat_rds)
  NULL
})


# If inferCNV succeeded
if (!is.null(infercnv_result)) {
  stopifnot("inferCNV did not write expected output" = file.exists(scratch_infercnv_rds))
  
  # Add inferCNV results to Seurat object
  seurat_obj <- infercnv::add_to_seurat(
    seurat_obj = NULL,
    infercnv_output_path = opts\$temp_dir
  )
  
  # Move results to final output files
  fs::file_move(scratch_infercnv_rds, opts\$output_rds)
  fs::file_move(scratch_metadata_file, opts\$output_table)
  fs::file_move(scratch_png_file, opts\$output_heatmap)
  
  # Save Seurat object
  saveRDS(seurat_obj, opts\$output_seurat_rds)
  message("Saved Seurat object with inferCNV results to: ", opts\$output_seurat_rds)
}

# Confirm all files exist
stopifnot(
  "inferCNV results file not created" = file.exists(opts\$output_rds),
  "inferCNV metadata table file not created" = file.exists(opts\$output_table),
  "PNG file not created" = file.exists(opts\$output_heatmap),
  "Seurat object file not created" = file.exists(opts\$output_seurat_rds)
)


# Write version information
writeLines(
    c(
        '"${task.process}":',
        paste('    r:', paste(R.Version()\$major, R.Version()\$minor, sep = ".")),
        paste('    anndataR:', as.character(packageVersion('anndataR'))),
        paste('    infercnv:', as.character(packageVersion('infercnv'))),
        paste('    Seurat:', as.character(packageVersion('Seurat')))
    ),
    'versions.yml'
)
