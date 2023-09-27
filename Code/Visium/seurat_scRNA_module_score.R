library(ggpubr)
library(optparse)
library(Seurat)
library(tidyverse)
library(viridis)

sessionInfo()

set.seed(42)

# Parse command line arguments  -------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='Path to sample Seurat object (.rds)')
parser <- add_option(parser, c('-m', '--markers'), default=NA, type='character',
                     help='Path to .csv file containing markers for modules')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for Seurat objects')

opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# opt$sample <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A04_seurat_QC/B1_24322.rds'
# opt$markers <- '/mnt/ix1/Projects/M039_170911_Gastric_scRNA/A00_interpatient_merge/P5846_5866_5931_6207_6342_6649_6592_6709_N_T_PBMC/A00_seurat_regressnUMI_GEMbatch/subset_epithelial/A00_seurat_20PC_res0.8/filter_subset_markers_class.csv'
# opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A05_seurat_scRNA_module_score/temp'

# Load data ---------------------------------------------------------------

st <- readRDS(opt$sample)

markers <- read_csv(opt$markers, col_types = cols(
  cluster = col_factor(NULL),
  gene = col_factor(NULL)
))

fn <- tools::file_path_sans_ext(basename(opt$sample)) 
prefix <- paste(opt$output_dir, fn, sep='/')

# Cropping spatial feature plot distorts image
# https://github.com/satijalab/seurat/issues/5141

img_coords = st@images$slice1@coordinates
aspect_ratio <- (max(img_coords$imagerow) - min(img_coords$imagerow)) / (max(img_coords$imagecol) - min(img_coords$imagecol))


# Add module scores -------------------------------------------------------

normal_epi <- markers %>% filter(cluster=="normal epithelium", p_val_adj <= 0.05, avg_logFC >= 0.25)
gc_type1 <- markers %>% filter(cluster=="tumor epithelium", p_val_adj <= 0.05, avg_logFC >= 0.25)
gc_type2 <- markers %>% filter(cluster=="tumor like epithelium", p_val_adj <= 0.05, avg_logFC >= 0.25)

st <- AddModuleScore(st, features = list(normal_epi$gene), assay = "Spatial", name = "scRNA_normal")
st <- AddModuleScore(st, features = list(gc_type1$gene), assay = "Spatial", name = "scRNA_GC_type1")
st <- AddModuleScore(st, features = list(gc_type2$gene), assay = "Spatial", name = "scRNA_GC_type2")

st@meta.data <- st@meta.data %>% rename(
  scRNA_normal = scRNA_normal1,
  scRNA_GC_type1 = scRNA_GC_type11,
  scRNA_GC_type2 = scRNA_GC_type21)

# Plot module score -------------------------------------------------------

SpatialFeaturePlot(st, features = 'scRNA_normal') + 
  scale_fill_viridis() +
  theme(legend.position = "right", aspect.ratio = aspect_ratio)
ggsave(paste(prefix, 'scRNA_normal_score.jpeg', sep='_'), dpi = 300)

SpatialFeaturePlot(st, features = 'scRNA_GC_type1') +
  scale_fill_viridis() +
  theme(legend.position = "right", aspect.ratio = aspect_ratio)
ggsave(paste(prefix, 'scRNA_GC_type1_score.jpeg', sep='_'), dpi = 300)

SpatialFeaturePlot(st, features = 'scRNA_GC_type2') +
  scale_fill_viridis() +
  theme(legend.position = "right", aspect.ratio = aspect_ratio)
ggsave(paste(prefix, 'scRNA_GC_type2_score.jpeg', sep='_'), dpi = 300)

VlnPlot(st, features = c('scRNA_normal', 'scRNA_GC_type1', 'scRNA_GC_type2'),
        same.y.lims = T, log = T)
ggsave(paste(prefix, 'scRNA_scores_violin_plot.jpeg', sep='_'), dpi = 300)

# Export seurat object ----------------------------------------------------

saveRDS(st, file = paste(prefix, '.rds', sep=''))
print(paste('Exported:', fn, sep=' '))

