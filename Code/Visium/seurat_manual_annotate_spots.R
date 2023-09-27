library(optparse)
library(Seurat)
library(tidyverse)

sessionInfo()
set.seed(42)

# Parse command line arguments  -------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='Path to sample Seurat object (.rds)')
parser <- add_option(parser, c('-b', '--barcodes'), default=NA, type='character',
                     help='Path to .csv file matching barcodes to sample')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for Seurat object with annotated spots')

opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

print(opt$sample)
# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A07_cluster_annotation/path_annotation/A1_24320.rds'
# opt$barcodes <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A08_manual_spot_annotation/A1_24320_Base_Pit_Metaplasia_barcodes.csv'
# opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A08_manual_spot_annotation/temp'

# Load data ---------------------------------------------------------------

st <- readRDS(opt$sample)
fn <- tools::file_path_sans_ext(basename(opt$sample)) 
prefix <- paste(opt$output_dir, fn, sep='/')

ann <-  read_csv(opt$barcodes) %>%
  rename(Annotation = 2) %>%
  mutate(Annotation = as.factor(Annotation))

# Cropping spatial feature plot distorts image
# https://github.com/satijalab/seurat/issues/5141

img_coords = st@images$slice1@coordinates
aspect_ratio <- (max(img_coords$imagerow) - min(img_coords$imagerow)) / (max(img_coords$imagecol) - min(img_coords$imagecol))


# Annotate spots ----------------------------------------------------------

# Subset only annotated spots
cell_barcodes <- ann %>% pull(Barcode)
st <- subset(st, cells = cell_barcodes)

meta.data_tbl <- as_tibble(st@meta.data, rownames='Barcode') %>%
  left_join(ann, by='Barcode')

meta.data_df <- as.data.frame(meta.data_tbl %>% select(-(Barcode))) 
rownames(meta.data_df) <- meta.data_tbl %>% pull(Barcode)

st@meta.data <- meta.data_df

Idents(st) <- 'Annotation'

# Print levels to check for NA's
print(levels(st$Annotation))
# Idents(st, cells = cell_barcodes) <- ann %>% pull(Annotation)


# Plot --------------------------------------------------------------------

DimPlot(st, reduction = "umap", label = F)
ggsave(paste(prefix, 'annotations_UMAP.jpeg', sep='_'), dpi = 300)

SpatialDimPlot(st, label = FALSE, repel = TRUE, pt.size.factor = 3) + theme(aspect.ratio = aspect_ratio) +
  labs(fill='Region')
ggsave(paste(prefix, 'annotations_spatial_cropped.jpeg', sep='_'), dpi = 300)

SpatialDimPlot(st, label = FALSE, repel = TRUE, pt.size.factor = 1, crop=FALSE) + 
  labs(fill='Region')
ggsave(paste(prefix, 'annotations_spatial_full_size.jpeg', sep='_'), dpi = 300)

# Export Seurat object ----------------------------------------------------

saveRDS(st, file = paste(prefix, '.rds', sep=''))
print(paste('Exported:', fn, sep=' '))
