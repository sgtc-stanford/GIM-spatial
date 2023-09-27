library('optparse')
library(Seurat)

# Parse command line arguments

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='Path to sample Seurat object (.rds)')
parser <- add_option(parser, c('-b', '--barcodes'), default=NA, type='character',
                     help='Path to .csv file matching barcodes to sample')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for demultiplexed Seurat objects')

opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# Load data
#opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A02_seurat_demux/seurat_objs/B1_24321_24322.rds'
#opt$barcodes <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A02_seurat_demux/B1_dissection_barcodes.csv'
#opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A02_seurat_demux/temp'

st <- readRDS(opt$sample)
barcodes <- read.csv(opt$barcodes)

fn <- tools::file_path_sans_ext(basename(opt$sample)) 
fn

st_array <- strsplit(fn, '_')[[1]][1]

# Demultiplex samples
st[['dissection']] <- barcodes$dissection

plot_name <- paste(fn, 'demux_spatial_plot.jpeg', sep='_')
jpeg(paste(opt$output_dir, plot_name, sep='/'))
SpatialDimPlot(st, group.by = 'dissection')
dev.off()

# Save each dissection/sample individually

for (dissection_value in unique(barcodes$dissection)) {
  st_subset <- subset(st, subset = dissection == dissection_value )
  subset_name <- paste(st_array, dissection_value, sep='_')
  
  # Plot individual dissections/samples
#  img_coords = st_subset@images$slice1@coordinates
#  aspect_ratio <- (max(img_coords$imagerow) - min(img_coords$imagerow)) / (max(img_coords$imagecol) - min(img_coords$imagecol))
  
#  jpeg(paste(opt$output_dir, '/', subset_name, '.jpeg', sep=''))
#  SpatialDimPlot(st_subset, alpha=0.3) + theme(aspect.ratio = aspect_ratio)
#  dev.off()
  saveRDS(st_subset, file = paste(opt$output_dir, '/', subset_name, '.rds', sep=''))
  print(paste('Exported:', subset_name, sep=' '))
}
