# File: seurat_metadata.R
# Author: Sue Grimes
# Desc: Script metadata information from supplied config file and saves to Seurat object
#
library(Seurat)
library(dplyr)
library(ini)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parse command line parameters                                               #                                                 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
suppressPackageStartupMessages(require(optparse))
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf

option_list = list(
  make_option(c("-r", "--seurat"), action="store", default=NA, type='character',
              help="Seurat R object"),
  make_option(c("-c", "--config_ini"), action="store", default=NA, type='character',
              help="Configuration parameters: metadata fields/values")	   
  )
  
required_args=c('seurat', 'config_ini')  
opt = parse_args(OptionParser(option_list=option_list))

cat("User input >>","\n")
cat("seurat:", opt$seurat,"\n")
cat("config_ini:", opt$config_ini,"\n")

for (arg in required_args) {
  if ( is.null(opt[[arg]]) ) {
    print(paste0("Required parameter: ", arg, " is missing"))
	quit(status=10)
  }
}  
  
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Read metadata config file, and save metadata fields/values to Seurat object             #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
seurat_obj <- readRDS(opt$seurat)
params <- read.ini(opt$config_ini)

for (fld in names(params$metadata)) { 
  seurat_obj@meta.data[[fld]] = params$metadata[[fld]] 
  }

saveRDS(seurat_obj, file=opt$seurat)

