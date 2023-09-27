# File: seurat_umap.R
# Author: Sue Grimes
# Desc: Script reads saved Seurat R object (with PCA calculated) and runs UMAP
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parse input parameters                                                      #                                                 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
suppressPackageStartupMessages(require(optparse))
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf

option_list = list(
  make_option(c("-r", "--seurat"), action="store", default=NA, type='character',
              help="Seurat/PCA R object"),
  make_option(c("-c", "--config_ini"), action="store", default=NA, type='character',
              help="Configuration parameters file, default=NA: pipeline defaults will be used")			  
)
opt = parse_args(OptionParser(option_list=option_list))

cat("User input >>","\n")
cat("seurat:", opt$seurat,"\n")
cat("config.ini:", opt$config_ini,"\n")

fn_in  = opt$seurat

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load packages, and custom helper functions, parse config.ini if given       #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)
library(pracma)
library(ini)
library(stringr)

source('/mnt/ix1/Resources/LabSoftware/software/scRNA-v3/pipeline/_scRNA_methods.R')

# Set standard default parameters (overwritten by config.ini if given)
params = list(dims_reduce=20, cluster_res=0.8, do_reorder='Y', algorithm=1, vst_features=2000)

if ( !is.null(opt$config_ini) ) {
  config_ini = read.ini(opt$config_ini)

  if ( 'umap' %in% names(config_ini) ) {
    config_params = type.convert(config_ini$umap, as.is=TRUE)
    params = modifyList(params, config_params)
	}
  }
  
cat("Seurat parameters: >>","\n")
print(params)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load Seurat object, run UMAP                                                #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
if ( str_sub(fn_in, start=-4) == 'Robj' ) {
  seurat_pca <- loadRObj(fn_in)
} else {
  seurat_pca <- readRDS(fn_in)
}  

#UMAP dimension reduction and associated plots
seurat_umap <- FindNeighbors(seurat_pca, dims=1:params$dims_reduce)
seurat_umap <- FindClusters(seurat_umap, resolution=params$cluster_res, algorithm=params$algorithm)
if ( params$do_reorder == 'Y' ) {
  seurat_umap <- BuildClusterTree(seurat_umap, reorder=T, reorder.numeric=T, dims=1:params$dims_reduce)
}	
seurat_umap <- RunUMAP(seurat_umap, dims=1:params$dims_reduce, do.label=T)

pdf(paste0(seurat_umap@project.name, '_umap_clusters.pdf'))
  DimPlot(seurat_umap, reduction='umap', label=T)
dev.off()

cellnumbers <- table(Idents(seurat_umap), seurat_umap@meta.data$orig.ident)
write.csv(cellnumbers, file=paste0(seurat_umap@project.name,'.cellnumbers.csv'))

# Normalize, scale RNA slot in preparation for find markers
seurat_umap <- NormalizeData(seurat_umap, assay="RNA")
seurat_umap <- FindVariableFeatures(seurat_umap, assay="RNA", selection.method="vst", nfeatures=params$vst_features)
seurat_umap <- ScaleData(seurat_umap, assay="RNA", features=rownames(seurat_umap), vars.to.regress="nCount_RNA")

saveRDS(seurat_umap, file=paste0(seurat_umap@project.name,'.seurat_umap.rds'))

