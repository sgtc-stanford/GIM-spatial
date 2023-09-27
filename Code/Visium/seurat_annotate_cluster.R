library(dplyr)
library(ggplot2)
library('optparse')
library(patchwork)
library(Seurat)

#################################
# Parse command line arguments #
#################################

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='Path to sample Seurat object (.rds)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for annotated Seurat objects')
parser <- add_option(parser, c('-a', '--annotation'), default=NA, type='character',
                     help='Cluster annotation .yml file')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A03_seurat_clustering/A1_24320.rds'
# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A06_seurat_clustering_per_patient/D1_24319.rds'
# opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A07_cluster_annotation/temp'
# opt$annotation <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A07_cluster_annotation/path_ann_for_infercnv.yml'

#############
# Load data #
#############

st <- readRDS(opt$sample)

# Extract metadata from file name
fn <- tools::file_path_sans_ext(basename(opt$sample)) 
st_array <- strsplit(fn, '_')[[1]][1]

# Cropping spatial feature plot distorts image
# https://github.com/satijalab/seurat/issues/5141

img_coords = st@images$slice1@coordinates
aspect_ratio <- (max(img_coords$imagerow) - min(img_coords$imagerow)) / (max(img_coords$imagecol) - min(img_coords$imagecol))

prefix <- paste(opt$output_dir, fn, sep='/')

# Select annotation to use # 
Sys.setenv(R_CONFIG_ACTIVE = fn)
ann <- config::get(file = opt$annotation, use_parent = FALSE)

#####################
# Annotate clusters #
#####################

st <- RenameIdents(st, ann)

# Subset sample to exclude certain clusters
if ('Exclude' %in% ann) {
  st <- subset(st, idents = 'Exclude', invert = TRUE)
}

# Plot UMAP
jpeg(paste(prefix, 'cell_types_UMAP.jpeg', sep='_'), width=5, height = 5, units='in', res=300)
DimPlot(st, reduction = "umap", label = TRUE)
dev.off()

# Plot on tissue
jpeg(paste(prefix, 'cell_types_spatial.jpeg', sep='_'), width=5, height=5, units='in', res=300)
SpatialDimPlot(st, label = TRUE, label.size = 3) + theme(aspect.ratio = aspect_ratio)
dev.off()

###################
# Cluster markers #
###################

st.markers <- FindAllMarkers(st, min.pct = 0.1, logfc.threshold = 0.25)
st.markers %>% 
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)

# Plot markers for each cluster
st.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

jpeg(paste(prefix, 'cell_type_top_10_markers_heatmap.jpeg', sep='_'), width=10, height=10, units='in', res=300)
DoHeatmap(st, features = top10$gene) + NoLegend()
dev.off()

#########################
# Export seurat object #
#########################

saveRDS(st, file = paste(prefix, '.rds', sep=''))
print(paste('Exported:', fn, sep=' '))
