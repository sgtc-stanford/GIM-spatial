library(optparse)
library(Seurat)
library(tidyverse)

sessionInfo()


# Parse command line arguments  -------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='Path to sample Seurat object (.rds)')
parser <- add_option(parser, c('-b', '--barcodes'), default=NA, type='character',
                     help='Path to .csv file matching barcodes to sample')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for Seurat objects with replicate metadata')

opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A02_seurat_demux/seurat_objs/B1_24321_24322.rds'
# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A02_seurat_demux/seurat_objs_demuxed/B1_24322.rds'
# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A02_seurat_demux/seurat_objs_demuxed/A1_24320.rds'
# opt$barcodes <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A03_seurat_label_replicates/B1_replicate_barcodes.csv'
# opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A03_seurat_label_replicates/temp'


# Load data ---------------------------------------------------------------

st <- readRDS(opt$sample)

if (is.na(opt$barcodes)) {
  replicates <- tibble(
    Barcode = rownames(st@meta.data),
    replicate = factor(1)
    )
} else {
  replicates <- read_csv(opt$barcodes, col_types = cols(
    replicate = col_factor()))
}

fn <- tools::file_path_sans_ext(basename(opt$sample)) 
prefix <- paste(opt$output_dir, fn, sep='/')

# Merge replicate values into metadata ------------------------------------

meta.data_tbl <- as_tibble(st@meta.data, rownames='Barcode') %>%
  left_join(replicates, by='Barcode')

meta.data_df <- as.data.frame(meta.data_tbl %>% select(-(Barcode))) 
rownames(meta.data_df) <- meta.data_tbl %>% pull(Barcode)

st@meta.data <- meta.data_df

# Print levels to check for NA's
print(levels(st$replicate))


# Plot replicates ---------------------------------------------------------

# Cropping spatial feature plot distorts image
# https://github.com/satijalab/seurat/issues/5141

img_coords = st@images$slice1@coordinates
aspect_ratio <- (max(img_coords$imagerow) - min(img_coords$imagerow)) / (max(img_coords$imagecol) - min(img_coords$imagecol))

SpatialDimPlot(st, group.by = 'replicate') + theme(legend.position = "right", aspect.ratio = aspect_ratio)
ggsave(paste(prefix, 'replicates.jpeg', sep='_'), dpi = 300)


# Export seurat object ----------------------------------------------------

saveRDS(st, file = paste(prefix, '.rds', sep=''))
print(paste('Exported:', fn, sep=' '))
