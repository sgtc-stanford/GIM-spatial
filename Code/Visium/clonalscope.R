library(Clonalscope)
library(gplots)
library(Seurat)
library(tidyverse)

set.seed(42)
sessionInfo()

# Parse command line arguments  -------------------------------------------

PROJECT_DIR <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/G01_clonalscope'
OUTPUT_DIR <- file.path(PROJECT_DIR, 'temp')

GITHUB_DIR <- file.path(PROJECT_DIR, 'Clonalscope')

CELLRANGER_DIR <- '/mnt/ix1/Seq_Runs/20220405_NV1_1188'
# FEATURES_DIR <- file.path(CELLRANGER_DIR, '24320_intramucosal_ca_visium/outs/filtered_feature_bc_matrix')
FEATURES_DIR <- file.path(CELLRANGER_DIR, '24319_GIM_visium/outs/filtered_feature_bc_matrix')

SEU_DIR <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A07_cluster_annotation/path_ann_for_infercnv'
# SEU_PATH <- file.path(SEU_DIR, 'A1_24320.rds')
# SEU_PATH <- file.path(SEU_DIR, 'B1_24321.rds')
SEU_PATH <- file.path(SEU_DIR, 'D1_24319.rds')

# Select whether to use counts, barcodes and features from Cell ranger or Seurat
USE_SEU <- TRUE

# Load data ---------------------------------------------------------------

# Size of each chromosome (hg19 and GRCh38 are provided.)
size=read.table(file.path(GITHUB_DIR, "data-raw/sizes.cellranger-GRCh38-1.0.0.txt"), 
                stringsAsFactors = F)

# List of cyclegenes retrieved from the "CopyKAT"package (https://github.com/navinlabcode/copykat)
cyclegenes=readRDS(file.path(GITHUB_DIR, "data-raw/cyclegenes.rds"))

# bed file indicating gene positions (hg19 and GRCh38 are provided.)
bed=read.table(file.path(GITHUB_DIR, "data-raw/hg38.genes.bed"), 
               sep='\t', header = T)

# Chromosome arm position files
chrarm=read.table(file.path(GITHUB_DIR, "data-raw/cytoarm_table_hg38.txt"), stringsAsFactors = F, sep='\t', header=T)
chrarm=chrarm[order(as.numeric(chrarm[,1]),as.numeric(chrarm[,3])),]
bin_bed=chrarm[,-2]

# Load spatial data -------------------------------------------------------
st <- readRDS(SEU_PATH)

# Cropping spatial feature plot distorts image
# https://github.com/satijalab/seurat/issues/5141

img_coords = st@images$slice1@coordinates
aspect_ratio <- (max(img_coords$imagerow) - min(img_coords$imagerow)) / (max(img_coords$imagecol) - min(img_coords$imagecol))

# Load: 
# (1) A (sparse) matrix with each row being a feature and each column being a cell barcode. (Output from cellranger)
# (2) A matrix/data.frame with barcodes for each cell in the first column.
# (3) A matrix/data.frame with the columns- 1st: gendID; 2nd: gene name. Each row is a feature whose order should correspond to the rows of "mtx"

features <- read.table(file.path(FEATURES_DIR, 'features.tsv.gz'), stringsAsFactors = F, sep = '\t', header = F)

if (USE_SEU == TRUE){
  print('Using seurat counts')
  mtx <- st@assays$Spatial@counts
  barcodes <- as.data.frame(rownames(st@meta.data))
  features <- features[match(rownames(mtx), features[,2]),]
}else{
  print('Using spaceranger counts')
  mtx <- readMM(file.path(FEATURES_DIR, 'matrix.mtx.gz'))
  barcodes <- read.table(file.path(FEATURES_DIR, 'barcodes.tsv.gz'), stringsAsFactors = F, sep = '\t', header = F)
}

# Cell type annotations
# celltype0=readRDS(file.path(GITHUB_DIR, "data-raw/BC_ductal2/ST/celltype_all.rds"))



# Classify as tumor or normal
# Assign cells to be used as baseline for coverage normalization as "normal"
colnames(celltype0)
cell_types <- as_tibble(Idents(st), rownames = 'barcode') %>%
  rename(celltype = value) %>%
  mutate(celltype = str_replace(celltype, 'Pit', 'normal'))

celltype = str_replace(celltype, 'Base', 'normal')
cell_types <- as_tibble(st$replicate, rownames = 'barcode') %>%
  rename(celltype = value)

cell_types$celltype <- recode_factor(cell_types$celltype, `3` = 'normal')

cell_types

SpatialDimPlot(st, group.by = 'replicate', label = F, label.size = 3, stroke = 0, pt.size.factor = 1) + theme(legend.position = "right", aspect.ratio = aspect_ratio)
# celltype = str_replace(celltype, 'Cancer1', 'tumor'),
# celltype = str_replace(celltype, 'Cancer3', 'tumor'),
# celltype = str_replace(celltype, 'Cancer4', 'tumor')
  # mutate(celltype = replace(celltype, celltype != 'normal', 'tumor'))
cell_types
cell_types_df <- as.data.frame(cell_types)

# Generate segmentation table for each chromosome arm ---------------------

seg_table_filtered=data.frame("chr"=bin_bed[,1], 'start'=as.numeric(bin_bed[,2]),
                              'end'=as.numeric(bin_bed[,3]), 'states'=1, 'length'=as.numeric(bin_bed[,3])-as.numeric(bin_bed[,2]),
                              'mean'=0, 'var'=0, 'Var1'=1:nrow(bin_bed),'Freq'=50000,
                              'chrr'=paste0(bin_bed[,1],":", bin_bed[,2]), stringsAsFactors = F)


# Filter input features ---------------------------------------------------

Input_filtered=FilterFeatures(mtx=mtx, barcodes=barcodes, features=features, cyclegenes=cyclegenes)

# Remove raw inputs
# rm(mtx); rm(barcodes); rm(features)

Input_filtered


# Subclone detection ------------------------------------------------------

Cov_obj=RunCovCluster(mtx=Input_filtered$mtx, barcodes=Input_filtered$barcodes, 
                      features=Input_filtered$features, bed=bed, 
                      celltype0=cell_types_df, var_pt=0.99, var_pt_ctrl=0.99, include='all',
                      alpha_source='all', ctrl_region=NULL, 
                      seg_table_filtered=seg_table_filtered, size=size,
                      dir_path=OUTPUT_DIR, breaks=50, prep_mode = 'intersect', seed=200) 
# save the object
# saveRDS(Cov_obj,paste0(dir_path,"/Cov_obj.rds"))


# Viz ---------------------------------------------------------------------

clustering= Cov_obj$result_final$clustering
clustering2= Cov_obj$result_final$clustering2
result=Cov_obj$result_final$result
Zest=result$Zest
table(result$Zest)

# result=AssignCluster(clustering2, mincell = 0)
# Zest=result$Zest
# table(result$Zest)

PlotClusters(df = clustering$data, celltype = cell_types_df, Assign_obj =result, mode = "genome",  fontsize = 7, lab_mode='annot')

emb=PlotUMAP(df = clustering$data, celltype = cell_types_df, Assign_obj =result, mode = "Zest")


# Spatial plotting --------------------------------------------------------

cols=col2hex(colors()[c(1,609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])

seu2 <- st

clusters=Zest[match(names(seu2$orig.ident), names(Zest))]; clusters[is.na(clusters)]=0
Idents(seu2)=factor(clusters, levels=sort(as.numeric(unique(clusters))))
SpatialDimPlot(seu2, label = F, label.size = 3, stroke = 0, pt.size.factor = 1) + theme(legend.position = "right", aspect.ratio = aspect_ratio)
