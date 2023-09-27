# File: remove_doublets.R
# Author: Sue Grimes
# Desc: Remove marked doublets from Seurat object
#
library(Seurat)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parse command line parameters                                               #                                                 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
suppressPackageStartupMessages(require(optparse))
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf

option_list = list(
  make_option(c("-r", "--seurat"), action="store", default=NA, type='character',
              help="Seurat R object")	   
  )
  
required_args=c('seurat')  
opt = parse_args(OptionParser(option_list=option_list))

cat("User input >>","\n")
cat("seurat:", opt$seurat,"\n")

for (arg in required_args) {
  if ( is.null(opt[[arg]]) ) {
    print(paste0("Required parameter: ", arg, " is missing"))
	quit(status=10)
  }
}  
  
seurat_obj <- readRDS(opt$seurat)
seurat_obj <- subset(seurat_obj, subset= DF!='Doublet')
saveRDS(seurat_obj, file=paste0(seurat_obj@project.name,'.seurat_dfx.rds'))
