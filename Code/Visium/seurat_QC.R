library(dplyr)
library(ggplot2)
# library(ggpubr)
library('optparse')
library(patchwork)
library(tidyverse)
library(Seurat)
library(viridis)

sessionInfo()

#################################
# Parse command line arguments #
#################################

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='Path to sample Seurat object (.rds)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for clustered Seurat objects')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

metrics_path <- file.path(opt$output_dir, 'QC_metrics.csv')

# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A02_seurat_demux/seurat_objs_demuxed/A1_24320.rds'
# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A02_seurat_demux/seurat_objs_demuxed/B1_24321.rds'
# opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A03_seurat_QC/temp'

#############
# Load data #
#############

st <- readRDS(opt$sample)

fn <- tools::file_path_sans_ext(basename(opt$sample)) 
st_array <- strsplit(fn, '_')[[1]][1]

if (file.exists(metrics_path)){
  metrics <- read_csv(metrics_path) %>% add_row(fn = fn)
} else {
  metrics <- tibble(fn = fn)
}

# Cropping spatial feature plot distorts image
# https://github.com/satijalab/seurat/issues/5141

img_coords = st@images$slice1@coordinates
aspect_ratio <- (max(img_coords$imagerow) - min(img_coords$imagerow)) / (max(img_coords$imagecol) - min(img_coords$imagecol))

prefix <- paste(opt$output_dir, fn, sep='/')

Idents(st) <- fn

################
# Filter Spots #
################

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
metrics[metrics$fn == fn, 'spots_pre_filter'] <- dim(st)[2]

st = st[, st$nFeature_Spatial > 500]

metrics[metrics$fn == fn, 'spots_post_filter'] <- dim(st)[2]

################
# Filter Genes #
################

metrics[metrics$fn == fn, 'genes_pre_filter'] <- dim(st)[1]

min_spots = 3
counts_binary = Matrix::rowSums(st@assays$Spatial@counts > 0)
selected_features <- rownames(st)[counts_binary > min_spots]
st <- subset(st, features = selected_features)

metrics[metrics$fn == fn, 'genes_post_filter'] <- dim(st)[1]

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

####################################
# Top expressed and variable genes #
####################################

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

jpeg(paste(prefix, 'top_15_var_genes.jpeg', sep='_'))
plot1 <- VariableFeaturePlot(st_var_features)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
plot2
dev.off()

#############################
# Normalise gene expression #
#############################

st <- SCTransform(st, assay = "Spatial", verbose = FALSE)

######################
# Cell Cycle scoring #
######################

# Assume data has already been normalised (i.e. SCTransform)
# https://github.com/satijalab/seurat/issues/1679

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

st <- CellCycleScoring(st, 
                       s.features = s.genes,
                       g2m.features = g2m.genes)

SpatialFeaturePlot(st, features = 'S.Score') + 
  scale_fill_viridis() +
  theme(legend.position = "right", aspect.ratio = aspect_ratio)
ggsave(paste(prefix, 'cell_cycle_S.jpeg', sep='_'), dpi = 300)

SpatialFeaturePlot(st, features = 'G2M.Score') +
  scale_fill_viridis() +
  theme(legend.position = "right", aspect.ratio = aspect_ratio)
ggsave(paste(prefix, 'cell_cycle_G2M.jpeg', sep='_'), dpi = 300)

SpatialDimPlot(st, group.by = 'Phase') + theme(legend.position = "right", aspect.ratio = aspect_ratio)
ggsave(paste(prefix, 'cell_cycle_phase.jpeg', sep='_'), dpi = 300)


#####################################
# Export Seurat object  and metrics #
#####################################

saveRDS(st, file = paste(prefix, '.rds', sep=''))
write_csv(metrics, metrics_path)
print(paste('Exported:', fn, sep=' '))
