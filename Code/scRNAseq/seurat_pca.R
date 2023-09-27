# File: seurat_pca.R
# Author: Sue Grimes
# Desc: Script reads saved Seurat R object and performs PCA dimension reduction #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parse input parameters                                                      #                                                 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
suppressPackageStartupMessages(require(optparse))
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf

option_list = list(
  make_option(c("-r", "--seurat"), action="store", default=NA, type='character',
              help="Seurat R object"),
  make_option(c("-c", "--config_ini"), action="store", type='character',
              help="Configuration parameters file, default=NA: pipeline defaults will be used")
)
opt = parse_args(OptionParser(option_list=option_list))

cat("User input >>","\n")
cat("seurat:", opt$seurat,"\n")
cat("config_ini:", opt$config_ini,"\n")

if ( is.null(opt$seurat) ) {
  print(paste0("Required parameter: seurat is missing"))
  quit(status=10)
} else {
  fn_in   = opt$seurat
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load packages and custom helper functions; parse config.ini parameters if provided      #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)
library(pracma)
library(stringr)
library(ini) 

source('/mnt/ix1/Resources/LabSoftware/software/scRNA-v3/pipeline/_scRNA_methods.R')

# Set standard default parameters (overwritten by config.ini if given)
params = list(nr_pcs=50, do_approx=TRUE)

if ( !is.null(opt$config_ini) ) {
  config_ini = read.ini(opt$config_ini)

  if ( 'pca' %in% names(config_ini) ) {
    config_params = type.convert(config_ini$pca, as.is=TRUE)
    params = modifyList(params, config_params)
  }
}

cat("Seurat PCA parameters: >>","\n")
print(params)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load Seurat object, perform PCA                                               #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
if ( str_sub(fn_in, start=-4) == 'Robj' ) {
  seurat_obj <- loadRObj(fn_in)
} else {
  seurat_obj <- readRDS(fn_in)
}  

# Run principal component analysis using variable genes
# defaults: npcs: 50; do_approx=T
#seurat_pca <- RunPCA(object=seurat_obj, npcs=params$nr_pcs, approx=params$do_approx, ndims.print=1:5, nfeatures.print=min(params$nr_pcs,30))
seurat_pca <- RunPCA(seurat_obj, npcs=params$nr_pcs, approx=params$do_approx)
saveRDS(seurat_pca, file=paste0(seurat_pca@project.name, '.seurat_pca.rds'))


