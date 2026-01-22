#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(infercnv)
  library(anndataR)
  library(Seurat)
})

option_list <- list(
  make_option(c("--h5ad"), type="character", help="Input .h5ad file"),
  make_option(c("--annotations"), type="character", help="Cell annotations file (2-column TSV: cell_id<TAB>group)"),
  make_option(c("--gene_order"), type="character", help="Gene order file (gene<TAB>chr<TAB>start<TAB>end)"),
  make_option(c("--ref_groups"), type="character", default="", help="Comma-separated reference group names"),
  make_option(c("--outdir"), type="character", default="infercnv_out", help="Output directory"),

  # infercnv::run parameters
  make_option(c("--cutoff"), type="double", default=0.10, help="infercnv cutoff"),
  make_option(c("--cluster_by_groups"), type="logical", default=TRUE, help="cluster_by_groups"),
  make_option(c("--denoise"), type="logical", default=TRUE, help="denoise"),
  make_option(c("--HMM"), type="logical", default=TRUE, help="HMM (must be TRUE if extracting features)"),
  make_option(c("--analysis_mode"), type="character", default="subclusters", help="analysis_mode (e.g. 'subclusters')"),
  make_option(c("--num_threads"), type="integer", default=4, help="Number of threads"),
  make_option(c("--tumor_subcluster_partition_method"), type="character", default="leiden", help="tumor_subcluster_partition_method"),

  # Feature extraction (inferCNV wiki: Extracting features)
  make_option(c("--extract_features"), type="logical", default=TRUE,
              help="Extract per-cell CNV feature vectors after inferCNV run (requires HMM outputs)."),
  make_option(c("--feature_top_n"), type="integer", default=10,
              help="top_n passed to infercnv::add_to_seurat()"),
  make_option(c("--write_feature_table"), type="logical", default=TRUE,
              help="Write extracted feature metadata table to outdir as map_metadata_from_infercnv.txt")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$h5ad) || is.null(opt$annotations) || is.null(opt$gene_order)) {
  stop("Missing required arguments: --h5ad, --annotations, --gene_order")
}

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)


# Read the AnnData objects and convert directly to Seurat
adata <- read_h5ad("${h5ad}")
seu <- adata\$as_Seurat()


# ---- Extract raw counts matrix (genes x cells) ----
raw_counts_matrix <- GetAssayData(seu, slot = "counts")     #seurat v4


#assay_names <- SummarizedExperiment::assayNames(sce)
#message("Available assays: ", paste(assay_names, collapse=", "))

#get_counts <- function(sce) {
#  an <- SummarizedExperiment::assayNames(sce)
#  if ("counts" %in% an) return(SummarizedExperiment::assay(sce, "counts"))
#  if ("raw" %in% an) return(SummarizedExperiment::assay(sce, "raw"))
#  if ("raw_counts" %in% an) return(SummarizedExperiment::assay(sce, "raw_counts"))
#  if ("X" %in% an) {
#    warning("Falling back to assay 'X' (may be normalized/log-transformed; not ideal for inferCNV).")
#    return(SummarizedExperiment::assay(sce, "X"))
#  }
#  stop("No suitable assay found for raw counts. Expected one of: counts/raw/raw_counts (or at least X).")
#}

#raw_counts_matrix <- get_counts(sce)

# Ensure sparse dgCMatrix if possible (helps memory)
if (!inherits(raw_counts_matrix, "dgCMatrix")) {
  raw_counts_matrix <- as(raw_counts_matrix, "dgCMatrix")
}

# Ensure dimnames exist
if (is.null(rownames(raw_counts_matrix)) || is.null(colnames(raw_counts_matrix))) {
  stop("Counts matrix must have rownames (genes) and colnames (cells).")
}

#message("Counts matrix dims: genes=", nrow(raw_counts_matrix), " cells=", ncol(raw_counts_matrix))

# ---- Parse reference groups ----
ref_groups <- character(0)
if (nchar(opt$ref_groups) > 0) {
  ref_groups <- strsplit(opt$ref_groups, ",", fixed=TRUE)[[1]]
  ref_groups <- trimws(ref_groups)
}
message("Reference groups: ", ifelse(length(ref_groups) == 0, "(none specified)", paste(ref_groups, collapse=", ")))

# ---- Create infercnv object ----
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = raw_counts_matrix,
  annotations_file  = opt$annotations,
  gene_order_file   = opt$gene_order,
  ref_group_names   = ref_groups
)

# ---- Run infercnv ----
if (isTRUE(opt$extract_features) && !isTRUE(opt$HMM)) {
  stop("--extract_features=TRUE requires --HMM=TRUE so that HMM outputs exist for feature extraction.")
}

message("Running inferCNV -> outdir: ", opt$outdir)
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = opt$cutoff,
  out_dir = opt$outdir,
  cluster_by_groups = opt$cluster_by_groups,
  denoise = opt$denoise,
  HMM = opt$HMM,
  analysis_mode = opt$analysis_mode,
  num_threads = opt$num_threads,
  tumor_subcluster_partition_method = opt$tumor_subcluster_partition_method
)

# ---- Extract features (inferCNV wiki: Extracting features) ----
# ALWAYS:
#  - instantiate Seurat object
#  - call infercnv::add_to_seurat()
#  - write Seurat RDS to fixed filename in current working directory
if (isTRUE(opt$extract_features)) {
  message("Extracting inferCNV features from: ", opt$outdir)


  # Add inferCNV-derived metadata/features to the seurat object
  seurat_obj <- infercnv::add_to_seurat(
    infercnv_output_path = opt$outdir,
    seurat_obj = seu,
    top_n = opt$feature_top_n
  )

  # Write feature table to outdir (wiki naming)
  if (isTRUE(opt$write_feature_table)) {
    feature_tsv <- file.path(opt$outdir, "map_metadata_from_infercnv.txt")
    utils::write.table(
      seurat_obj@meta.data,
      file = feature_tsv,
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
    message("Wrote feature table: ", feature_tsv)
  }

  # ALWAYS write Seurat RDS to a fixed output name (no longer optional)
  #fixed_rds <- "seurat_with_infercnv_features.rds"
  #saveRDS(seurat_obj, fixed_rds)
  #message("Wrote Seurat object with inferCNV features: ", file.path(getwd(), fixed_rds))
  write_h5ad(seurat_obj, "${prefix}.h5ad")
}

# Save the output


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





#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(infercnv)
  library(Seurat)
})

# -----------------------------
# Nextflow-injected variables
# -----------------------------
h5ad_in            <- "!{h5ad}"
annotations_file   <- "!{annotations_tsv}"
gene_order_file    <- "!{gene_order_tsv}"

ref_groups_str     <- "!{ref_groups}"
cutoff             <- as.numeric("!{cutoff}")

cluster_by_groups  <- tolower("!{cluster_by_groups}") %in% c("true","t","1","yes","y")
denoise            <- tolower("!{denoise}") %in% c("true","t","1","yes","y")
analysis_mode      <- "!{analysis_mode}"

feature_top_n      <- as.integer("!{feature_top_n}")
write_feature_tbl  <- tolower("!{write_feature_table}") %in% c("true","t","1","yes","y")

prefix             <- "!{prefix}"

# -----------------------------
# Output conventions (SOUPX-style prefix)
# -----------------------------
outdir         <- paste0(prefix, ".infercnv_out")
seurat_rds_out <- paste0(prefix, ".seurat_with_infercnv_features.rds")
versions_yml   <- "versions.yml"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading h5ad: ", h5ad_in)
sce <- zellkonverter::readH5AD(h5ad_in)

# ---- Extract raw counts matrix (genes x cells) ----
get_counts <- function(sce) {
  an <- SummarizedExperiment::assayNames(sce)
  if ("counts" %in% an) return(SummarizedExperiment::assay(sce, "counts"))
  if ("raw" %in% an) return(SummarizedExperiment::assay(sce, "raw"))
  if ("raw_counts" %in% an) return(SummarizedExperiment::assay(sce, "raw_counts"))
  if ("X" %in% an) {
    warning("Falling back to assay 'X' (may be normalized/log-transformed; not ideal for inferCNV).")
    return(SummarizedExperiment::assay(sce, "X"))
  }
  stop("No suitable assay found for raw counts. Expected: counts/raw/raw_counts (or X as last resort).")
}

raw_counts_matrix <- get_counts(sce)

if (!inherits(raw_counts_matrix, "dgCMatrix")) {
  raw_counts_matrix <- as(raw_counts_matrix, "dgCMatrix")
}
if (is.null(rownames(raw_counts_matrix)) || is.null(colnames(raw_counts_matrix))) {
  stop("Counts matrix must have rownames (genes) and colnames (cells).")
}

message("Counts matrix dims: genes=", nrow(raw_counts_matrix), " cells=", ncol(raw_counts_matrix))

# ---- Parse reference groups ----
ref_groups <- character(0)
if (nchar(ref_groups_str) > 0) {
  ref_groups <- trimws(strsplit(ref_groups_str, ",", fixed = TRUE)[[1]])
}
message("Reference groups: ", ifelse(length(ref_groups) == 0, "(none)", paste(ref_groups, collapse = ", ")))

# ---- Create inferCNV object ----
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = raw_counts_matrix,
  annotations_file  = annotations_file,
  gene_order_file   = gene_order_file,
  ref_group_names   = ref_groups
)

# ---- Run inferCNV (force HMM because features are required) ----
num_threads <- as.integer(Sys.getenv("NXF_TASK_CPUS", "4"))

message("Running inferCNV -> outdir: ", outdir)
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = cutoff,
  out_dir = outdir,
  cluster_by_groups = cluster_by_groups,
  denoise = denoise,
  HMM = TRUE,
  analysis_mode = analysis_mode,
  num_threads = num_threads
)

# ---- Extract features (always) via add_to_seurat ----
message("Extracting features (top_n=", feature_top_n, ")")

seurat_obj <- Seurat::CreateSeuratObject(counts = raw_counts_matrix)

seurat_obj <- infercnv::add_to_seurat(
  infercnv_output_path = outdir,
  seurat_obj = seurat_obj,
  top_n = feature_top_n
)

if (isTRUE(write_feature_tbl)) {
  feature_tsv <- file.path(outdir, "map_metadata_from_infercnv.txt")
  utils::write.table(
    seurat_obj@meta.data,
    file = feature_tsv,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  message("Wrote feature table: ", feature_tsv)
}

saveRDS(seurat_obj, seurat_rds_out)
message("Wrote Seurat RDS: ", seurat_rds_out)

# ---- versions.yml (simple, no extra deps) ----
# Keep YAML minimal to avoid needing the 'yaml' package in the container.
lines <- c(
  "INFERCNV_FROM_H5AD:",
  paste0("  r-base: \"", R.version.string, "\""),
  paste0("  infercnv: \"", as.character(utils::packageVersion("infercnv")), "\""),
  paste0("  seurat: \"", as.character(utils::packageVersion("Seurat")), "\""),
  paste0("  zellkonverter: \"", as.character(utils::packageVersion("zellkonverter")), "\""),
  paste0("  SingleCellExperiment: \"", as.character(utils::packageVersion("SingleCellExperiment")), "\"")
)
writeLines(lines, con = versions_yml)

message("Done.")

