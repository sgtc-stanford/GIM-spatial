library(optparse)
library(Seurat)
library(spacexr)
library(tidyverse)

sessionInfo()

set.seed(42)

# Parse command line arguments  -------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='RCTD object (.rds)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/E01_seurat_to_RCTD/B1_24321_RCTD.rds'
# opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/E02_RCTD_deconvolution/temp'

# Load data ---------------------------------------------------------------

myRCTD <- readRDS(opt$sample)

fn <- tools::file_path_sans_ext(basename(opt$sample)) 
prefix <- paste(opt$output_dir, fn, sep='/')

# Run RCTD ----------------------------------------------------------------

# Assign 1-2 cell types per spot
myRCTD_db <- run.RCTD(myRCTD, doublet_mode = 'doublet')

# Assign up to 4 (default) cell types per spot
myRCTD_multi <-  run.RCTD(myRCTD, doublet_mode = 'multi')

# Assign any no. of cell types per spot
myRCTD_full <- run.RCTD(myRCTD, doublet_mode = 'full')

# Export outputs ----------------------------------------------------------

saveRDS(myRCTD_db, paste(prefix, 'doublet.rds', sep='_'))
saveRDS(myRCTD_full, paste(prefix, 'full.rds', sep='_'))
saveRDS(myRCTD_multi, paste(prefix, 'multi.rds', sep='_'))
print(paste('Exported:', fn, sep=' '))
