library(optparse)
library(Seurat)
library(tidyverse)

sessionInfo()
set.seed(42)


# Parse command line arguments --------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='Path to sample Seurat object (.rds)')
parser <- add_option(parser, c('-e', '--exclude'), default=NA, type='character',
                     help='Clusters to exclude in .yml file')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for Seurat object with module scores')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

PROJECT_DIR <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium'
ANALYSIS_DIR <- file.path(PROJECT_DIR, 'A07_seurat_cluster_QC')
INPUT_DIR <- file.path(PROJECT_DIR, 'A06_seurat_clustering_per_patient')

# sample_name <- 'A1_24320.rds'
# opt$sample <- file.path(INPUT_DIR, sample_name)
# opt$exclude <- file.path(ANALYSIS_DIR, 'exclude.yml')
# opt$output_dir <- file.path(file.path(ANALYSIS_DIR, 'temp'))
# 
# opt$sample
# 
# sample_name <- 'D1_24319.rds'

# Load data ---------------------------------------------------------------


st <- readRDS(opt$sample)
fn <- tools::file_path_sans_ext(basename(opt$sample)) 
prefix <- file.path(opt$output_dir, fn)

# Cropping spatial feature plot distorts image
# https://github.com/satijalab/seurat/issues/5141

img_coords = st@images$slice1@coordinates

get_aspect_ratio <- function(x){
  (max(x$imagerow) - min(x$imagerow)) / (max(x$imagecol) - min(x$imagecol))
}

aspect_ratio <- get_aspect_ratio(img_coords)

# Load clusters to exlude
Sys.setenv(R_CONFIG_ACTIVE = 'samples')
exclude <- config::get(file = opt$exclude, use_parent = FALSE)

exclude_cluster <- exclude[[fn]]

# is.null(exclude_cluster)
# is.null(exclude$D1_24319)
# 
# is.null(exclude[['D1_24319']])

# Filter out clusters -----------------------------------------------------

if (!is.null(exclude_cluster)) {
  print(paste0('Removing cluster ', exclude_cluster, ' from ', fn))
  st <- subset(st, subset = seurat_clusters != exclude_cluster)
}

# Spatial plots -----------------------------------------------------------

jpeg(paste(prefix, 'clusters_spatial.jpeg', sep='_'), width=5, height=5, units='in', res=300)
SpatialDimPlot(st, label = TRUE, label.size = 3) + theme(aspect.ratio = aspect_ratio)
dev.off()

# Export seurat object ----------------------------------------------------

saveRDS(st, file = paste(prefix, '.rds', sep=''))
print(paste('Exported:', fn, sep=' '))
