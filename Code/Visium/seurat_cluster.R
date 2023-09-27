library(dplyr)
library(ggplot2)
library('optparse')
library(patchwork)
library(Seurat)

set.seed(42)

# Parse command line arguments --------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='Path to sample Seurat object (.rds)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for clustered Seurat objects')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A02_seurat_demux/seurat_objs_demuxed/B1_24321.rds'
#. opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A03_seurat_clustering/temp'

#############
# Load data #
#############

st <- readRDS(opt$sample)

fn <- tools::file_path_sans_ext(basename(opt$sample)) 
st_array <- strsplit(fn, '_')[[1]][1]

# Cropping spatial feature plot distorts image
# https://github.com/satijalab/seurat/issues/5141

img_coords = st@images$slice1@coordinates
aspect_ratio <- (max(img_coords$imagerow) - min(img_coords$imagerow)) / (max(img_coords$imagecol) - min(img_coords$imagecol))

prefix <- paste(opt$output_dir, fn, sep='/')

######
# QC #
######

# Counts per spot pre-filtering
jpeg(paste(prefix, 'counts_pre-filter.jpeg', sep='_'), width=10, height=5, units='in', res = 300)
plot1 <- VlnPlot(st, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(st, features = "nCount_Spatial") + theme(legend.position = "right", aspect.ratio = aspect_ratio)
wrap_plots(plot1, plot2)
dev.off()

# Features per spot pre-filtering
jpeg(paste(prefix, 'features_pre-filter.jpeg', sep='_'), width=10, height=5, units='in', res = 300)
plot1 <- VlnPlot(st, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(st, features = "nFeature_Spatial") + theme(legend.position = "right", aspect.ratio = aspect_ratio)
wrap_plots(plot1, plot2)
dev.off()

### Filter spots by counts ###
spot_num_pre_filter <- dim(st)[2]

st = st[, st$nFeature_Spatial > 500]

print(c("No. of spots pre-filter: ", spot_num_pre_filter, "No. of spots post-filter: ", dim(st)[2]))

### Filter genes by spots ###
gene_num_pre_filter <- dim(st)[1]

min_spots = 3
counts_binary = Matrix::rowSums(st@assays$Spatial@counts > 0)
selected_features <- rownames(st)[counts_binary > min_spots]
st <- subset(st, features = selected_features)

print(c("No. of genes pre-filter:", gene_num_pre_filter, "No. of genes post-filter:", dim(st)[1]))

# Counts per spot post-filtering
jpeg(paste(prefix, 'counts_post-filter.jpeg', sep='_'), width=10, height=5, units='in', res = 300)
plot1 <- VlnPlot(st, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(st, features = "nCount_Spatial") + theme(legend.position = "right", aspect.ratio = aspect_ratio)
wrap_plots(plot1, plot2)
dev.off()

# Features per spot post-filtering
jpeg(paste(prefix, 'features_post-filter.jpeg', sep='_'), width=10, height=5, units='in', res = 300)
plot1 <- VlnPlot(st, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(st, features = "nFeature_Spatial") + theme(legend.position = "right", aspect.ratio = aspect_ratio)
wrap_plots(plot1, plot2)
dev.off()

# Top expressed genes
counts = st@assays$Spatial@counts
counts@x = counts@x/rep.int(colSums(counts), diff(counts@p))
most_expressed <- order(Matrix::rowSums(counts), decreasing = T)[20:1]

jpeg(paste(prefix, 'top_20_genes.jpeg', sep='_'),  width=8, height=8, units='in', res = 300)
boxplot(as.matrix(t(as.matrix(counts[most_expressed, ]))), cex = 0.1, las = 1, xlab = "% total count per spot",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

dev.off()

# Top variable genes
st_var_features <- FindVariableFeatures(st, selection.method = "vst", nfeatures = 2000)
top15 <- head(VariableFeatures(st_var_features), 15)

jpeg(paste(prefix, 'top_15_var_genes.jpg', sep='_'))
plot1 <- VariableFeaturePlot(st_var_features)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
plot2
dev.off()

# Normalise gene expression
st <- SCTransform(st, assay = "Spatial", verbose = FALSE)

# Top variable genes
st_var_features <- FindVariableFeatures(st, selection.method = "vst", nfeatures = 2000)
top15 <- head(VariableFeatures(st_var_features), 15)

jpeg(paste(prefix, 'top_15_var_genes.jpeg', sep='_'), width=10, height = 5, units='in', res=300)
plot1 <- VariableFeaturePlot(st_var_features)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
plot2
dev.off()

############################
# Dimentionality Reduction #
############################

st <- RunPCA(st, assay = "SCT", verbose = FALSE)

jpeg(paste(prefix, 'PCA_Elbow_plot.jpeg', sep='_'), width=10, height = 5, units='in', res=300)
ElbowPlot(st, ndims=50)
dev.off()

##############
# Clustering #
##############

pca_dims = 20
st <- FindNeighbors(st, reduction = "pca", dims = 1:pca_dims)
st <- FindClusters(st, verbose = FALSE, resolution = 1.4)
st <- RunUMAP(st, reduction = "pca", dims = 1:pca_dims)

# Plot UMAP
jpeg(paste(prefix, 'clusters_UMAP.jpeg', sep='_'), width=5, height = 5, units='in', res=300)
DimPlot(st, reduction = "umap", label = TRUE)
dev.off()

# Plot clusters on tissue
jpeg(paste(prefix, 'clusters_spatial.jpeg', sep='_'), width=5, height=5, units='in', res=300)
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

jpeg(paste(prefix, 'clusters_top_10_markers_heatmap.jpeg', sep='_'), width=10, height=10, units='in', res=300)
DoHeatmap(st, features = top10$gene) + NoLegend()
dev.off()

#########################
# Export seurat object #
#########################

saveRDS(st, file = paste(prefix, '.rds', sep=''))
print(paste('Exported:', fn, sep=' '))
