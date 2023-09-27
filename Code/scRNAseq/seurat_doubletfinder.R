#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parse input parameters                                                      #                                                 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
suppressPackageStartupMessages(require(optparse))
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf

option_list = list(
  make_option(c("-r", "--seurat"), action="store", default=NA, type='character',
              help="Seurat R umap object"),
  make_option(c("-c", "--config_ini"), action="store", default=NA, type='character',
              help="Configuration parameters file, default=NA: pipeline defaults will be used")
)

opt = parse_args(OptionParser(option_list=option_list))

cat("User input >>","\n")
cat("seurat:", opt$seurat,"\n")
cat("config.ini:", opt$config_ini, "\n")

if ( is.na(opt$seurat) ) {
  print(paste0("Required parameter: seurat is missing"))
  quit(status=10)
}  
fn_in = opt$seurat

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load packages, and custom helper functions                                  #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
suppressPackageStartupMessages(require(ggplot2))
# manual: https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf
suppressPackageStartupMessages(require(Matrix))
# manual: https://cran.r-project.org/web/packages/Matrix/Matrix.pdf
library(Seurat)
library(dplyr)
library(DoubletFinder)
library(ini)
library(stringr)

source('/mnt/ix1/Resources/LabSoftware/software/scRNA-v3/pipeline/_scRNA_methods.R')

# Set standard default parameters (overwritten by config.ini if given)
params = list(remove_doublets = 'Y')

if ( !is.na(opt$config_ini) ) {
  config_ini = read.ini(opt$config_ini)
  if ( 'doublets' %in% names(config_ini) ) {
    config_params = type.convert(config_ini$doublets, as.is=TRUE)
    params = modifyList(params, config_params)
  }
}

cat("DoubletFinder parameters: >>","\n")
print(params)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Read Seurat object, annotate doublets, save for further analysis            #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
if ( str_sub(fn_in, start=-4) == 'Robj' ) {
  seurat_obj <- loadRObj(fn_in)
} else {
  seurat_obj <- readRDS(fn_in)
}  
#get pK value
sweep.res.list <- paramSweep_v3(seurat_obj, PCs = 1:20, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT=F)
bcmvn <- find.pK(sweep.stats)
tmp <- bcmvn[with(bcmvn,order(-BCmetric)),]
tmp$pK <- as.character(tmp$pK)
pk_val= as.numeric(tmp[1,2])

# pN (proportion of doublets vs merged real-artificial data) defaults to 0.25; results are largely pN invariant
# pK (PC neighborhood size as a proportion of merged real-artificial data); from find.pK above
#two thresholds for doublets based on cell loading and homotypic doublets. First with cell loading
#only cell loading used
cell_loading <- read.csv('/mnt/ix1/Resources/LabSoftware/software/scRNA-v3/10x_multiplet_rate.csv')
seurat_nr_cells = length(Idents(seurat_obj))
doublet_rate_for_seurat_cells <- cell_loading[which(abs(cell_loading$cells_recovered-seurat_nr_cells)==min(abs(cell_loading$cells_recovered-seurat_nr_cells))),]
droplet_doublet_rate <- doublet_rate_for_seurat_cells[1,1]
nExp_poi <- round(droplet_doublet_rate*seurat_nr_cells)
seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN=0.25, pK=pk_val, nExp=nExp_poi, reuse.pANN=F, sct = TRUE)
names(seurat_obj@meta.data) <- gsub("^.*DF.*$","DF", names(seurat_obj@meta.data))
names(seurat_obj@meta.data) <- gsub("^.*pANN.*$","pANN", names(seurat_obj@meta.data))

## With homotypic adjustment
#homotypic.prop <- modelHomotypic(seurat_obj@meta.data$cell_label)
#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN=0.25, pK=pk_val, nExp=nExp_poi.adj, reuse.pANN="pANN")
#names(seurat_obj@meta.data) <- gsub("^.*DF.classifications.*$","DF_homotypic", names(seurat_obj@meta.data))

#seurat_obj@meta.data[,"DF_hi.lo"] <- seurat_obj@meta.data$DF
#seurat_obj@meta.data$DF_hi.lo[which(seurat_obj@meta.data$DF_hi.lo == "Doublet" & seurat_obj@meta.data$DF_homotypic == "Singlet")] <- "Doublet_lo"
#seurat_obj@meta.data$DF_hi.lo[which(seurat_obj@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"

pdf(paste0(seurat_obj@project.name, '_df.pdf'))
DimPlot(seurat_obj, group.by="DF")
dev.off()
 
doublets <- table(Idents(seurat_obj), seurat_obj@meta.data$DF)
write.csv(doublets, file=paste0(seurat_obj@project.name,'.doublets.csv'))

if ( params$remove_doublets == 'Y' ) {
  seurat_obj <- subset(seurat_obj, subset= DF!='Doublet')
  length(Idents(seurat_obj))
  saveRDS(seurat_obj, file=paste0(seurat_obj@project.name,'.seurat_dfx.rds'))
} else {
  saveRDS(seurat_obj, file=paste0(seurat_obj@project.name,'.seurat_df.rds'))
}
