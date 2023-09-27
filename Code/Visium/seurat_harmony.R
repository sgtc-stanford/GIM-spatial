library(optparse)
library(Seurat)
library(tidyverse)
library(harmony)

sessionInfo()

set.seed(42)

# Parse command line arguments  -------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample_dir'), default=NA, type='character', 
                     help='Directory containing Seurat objects (.rds)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

opt$sample_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A04_seurat_QC'
opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A05_seurat_harmony/testing'

# Load Data ---------------------------------------------------------------

metadata <- tibble(
  files = list.files(path=opt$sample_dir, pattern = '*.rds'),
  sample = tools::file_path_sans_ext(files),
)

# Load each Seurat object
seu_objs <- c()
for (file in metadata$files){
  st <- file.path(opt$sample_dir, file) %>%
    readRDS()
  st$orig.ident <- filter(metadata, files == file)$sample
  seu_objs <- append(seu_objs, st)
}

### Merge seurat objects ### 

st_all <- merge(x = seu_objs[[1]], y = seu_objs[-1], add.cell.ids = metadata$sample )

prefix <- paste(opt$output_dir, 'ABCD_24318-24322', sep='/')

# Preprocess --------------------------------------------------------------

# Unclear if ScTransform should be run before Harmony
# https://github.com/immunogenomics/harmony/issues/41

# Select top features across all patients
st_all_features <- SelectIntegrationFeatures(object.list = seu_objs, nfeatures = 2000)
VariableFeatures(st_all) <- st_all_features

st_all <- RunPCA(st_all, assay = "SCT", verbose = FALSE)

ElbowPlot(st_all, n=30)
ggsave(paste(prefix, 'pre_harmony_PCA_elbow.jpeg', sep='_'), dpi = 300)

DimPlot(st_all, dims=c(1,2), reduction = "pca", group.by = 'orig.ident')
ggsave(paste(prefix, 'pre_harmony_PCA_patients.jpeg', sep='_'), dpi = 300)


pca_dims = 20
st_all <- FindNeighbors(st_all, reduction = "pca", dims = 1:pca_dims)
st_all <- RunUMAP(st_all, reduction = "pca", dims = 1:pca_dims)

DimPlot(st_all, reduction = "umap", label = FALSE, group.by = 'orig.ident')
ggsave(paste(prefix, 'pre_harmony_UMAP_patients.jpeg', sep='_'), dpi = 300)

# Harmony -----------------------------------------------------------------

st_all <- RunHarmony(st_all, assay.use = 'SCT', reduction = 'pca', dims.use = 1:20, group.by.vars = 'orig.ident')

# Vis batch correction ----------------------------------------------------

DimPlot(st_all, dims=c(1,2), reduction = "harmony", group.by = 'orig.ident')
ggsave(paste(prefix, 'post_harmony_PCA_patients.jpeg', sep='_'), dpi = 300)

st_all <- FindNeighbors(st_all, reduction = "harmony", dims = 1:pca_dims)
st_all <- RunUMAP(st_all, reduction = "harmony", dims = 1:pca_dims)

DimPlot(st_all, reduction = "umap", label = FALSE, group.by = 'orig.ident')
ggsave(paste(prefix, 'post_harmony_UMAP_patients.jpeg', sep='_'), dpi = 300)


# Export Seurat object ----------------------------------------------------

saveRDS(st, file = paste(prefix, '.rds', sep=''))
print(paste('Exported:', fn, sep=' '))
