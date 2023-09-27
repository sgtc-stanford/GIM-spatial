library(optparse)
library(spacexr)
library(Seurat)
library(tidyverse)

sessionInfo()

set.seed(42)

# Parse command line arguments  -------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='Seurat object (.rds)')
parser <- add_option(parser, c('-r', '--reference'), default=NA, type='character',
                     help='Path to scRNAseq reference (.rds)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A04_seurat_QC/B1_24322.rds'
# opt$reference <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/D01_scRNA_reference_QC/aggr.final_celltype.rds'
# opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A08_seurat_to_RCTD/temp'

# Load data ---------------------------------------------------------------

# ST
st <- readRDS(opt$sample)

fn <- tools::file_path_sans_ext(basename(opt$sample)) 
prefix <- paste(opt$output_dir, fn, sep='/')

# Cropping spatial feature plot distorts image
# https://github.com/satijalab/seurat/issues/5141

img_coords = st@images$slice1@coordinates
aspect_ratio <- (max(img_coords$imagerow) - min(img_coords$imagerow)) / (max(img_coords$imagecol) - min(img_coords$imagecol))

# scRNA-seq reference
ref <- readRDS(opt$reference)


# Prepare reference -------------------------------------------------------

DefaultAssay(ref) <- 'SCT'

counts_ref <- as.matrix(ref@assays$RNA@counts)
cell_types_ref <- as.factor(ref$final_celltype)
names(cell_types_ref) <- colnames(ref)
reference <- spacexr::Reference(counts = counts_ref, cell_types = cell_types_ref)

# Prepare ST --------------------------------------------------------------

# Load and transform coordinates
coords <-st@images$slice1@coordinates[c(5,4)]
colnames(coords) <- c('x', 'y')
# Mirror/flip along x-axis
coords$y <- coords$y * -1
coords <- as.data.frame(coords)

ggplot(coords, aes(x=x, y=y)) +
  geom_point() +
  coord_fixed()
ggsave(paste(prefix, 'spot_scatterplot.jpeg', sep='_'), dpi = 300)

# Get counts 
counts <- as.matrix(st@assays$Spatial@counts)

# Create RCTD spatial object
puck <- SpatialRNA(coords = coords, counts = counts)
barcodes <- colnames(puck@counts)

plot_puck_continuous(puck, barcodes, puck@nUMI, 
                     ylimit = c(0, round(quantile(puck@nUMI, 0.9))),
                     title = 'plot of nUMI', 
                     size = 2)
ggsave(paste(prefix, 'spatial_nUMI.jpeg', sep='_'), dpi = 300)


# Create RCTD object ------------------------------------------------------

myRCTD <- create.RCTD(puck, reference, max_cores = 10)

# Export RCTD object ------------------------------------------------------

saveRDS(myRCTD, file = paste(prefix, 'RCTD.rds', sep='_'))
print(paste('Exported:', fn, sep=' '))

