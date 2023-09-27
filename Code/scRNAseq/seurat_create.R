# File: seurat_create.R
# Author: Sue Grimes
# Desc: Script reads cellranger sparse matrix format files, and creates Seurat object
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parse command line parameters                                               #                                                 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
suppressPackageStartupMessages(require(optparse))
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf

option_list = list(
  make_option(c("-m", "--mprefix"), action="store", default=NA, type='character',
              help="Prefix for cellranger matrix file (<mprefix>.matrix.mtx)"),
  make_option(c("-p", "--project"), action="store", default=NA, type='character',
              help="Project name, default=mprefix value"),
  make_option(c("-d", "--dir_name"), action="store", default='.', type='character',
              help="Directory for input files, default='.'"),
  make_option(c("-c", "--config_ini"), action="store", default=NA, type='character',
              help="Configuration parameters file, default=NA: pipeline defaults will be used"),
  make_option(c("-b", "--hq_barcodes"), action="store", default='None', type='character',
              help="High quality barcodes file, default='None'")		  
  
  )
opt = parse_args(OptionParser(option_list=option_list))

cat("User input >>","\n")
cat("dir_name:", opt$dir_name,"\n")
cat("mprefix:", opt$mprefix,"\n")
cat("project:", opt$project,"\n")
cat("config_ini:", opt$config_ini,"\n")
cat("hq barcodes:", opt$hq_barcodes,"\n")

dir_name  = opt$dir_name
if (is.na(opt$project)) { project = opt$mprefix } else { project = opt$project } 

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load packages and custom helper functions; parse config.ini parameters if provided      #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(ini)

source('/mnt/ix1/Resources/LabSoftware/software/scRNA-v3/pipeline/_scRNA_methods.R')

# Set standard default parameters (overwritten by config.ini if given)
params = list(min_cells=3, min_genes=200, max_genes=5000, max_mito=30, cell_cycle='Y', 
              sct_features=2000, sct_regress=NULL, only_var_genes=FALSE)
metadata = FALSE

if ( !is.null(opt$config_ini) ) {
  config_ini = read.ini(opt$config_ini)

  if ( 'create' %in% names(config_ini) ) {
    config_params = type.convert(config_ini$create, as.is=TRUE)
    params = modifyList(params, config_params)
	params$sct_regress <- null_or_vector(params$sct_regress)
	
    if ( 'metadata' %in% names(config_ini) ) { metadata = TRUE }
  }
}

cat("Seurat parameters: >>","\n")
print(params)

# Add metadata from sample_config.ini in current directory, if file exists
if ( file.exists('sample_config.ini') ) {
  sparams = read.ini('sample_config.ini')
  if ( 'metadata' %in% names(sparams) ) {
    if ( metadata == TRUE ) {
      params$metadata <- c(config_ini$metadata, sparams$metadata)
    } else {
      params$metadata <- sparams$metadata
    }
    metadata = TRUE
  }
}  

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# General variables, methods                                                  #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
filesep = .Platform$file.sep
geneID  = 'hugo'
cr_genome = Sys.getenv('GENOME')
organism = ifelse(substr(cr_genome,1,2)=='mm', 'Mus musculus', 'Homo sapiens')

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Initialize Seurat object, and write as R object for subsequent analysis     #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
expr_mtx <- read10Xtrio(dirname=dir_name, prefix10x=opt$mprefix)

seurat_obj <- CreateSeuratObject(counts=expr_mtx, min.cells=params$min_cells, min.features=params$min_genes, project=project)
length(Idents(seurat_obj))
seurat_obj <- annotateMito(seurat_obj, geneID=geneID)

#Initial QC plots
pdf(paste0(seurat_obj@project.name, '_qc_plots.pdf'))
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
  plot1 <- FeatureScatter(seurat_obj, feature1="nCount_RNA", feature2="percent.mito")
  plot2 <- FeatureScatter(seurat_obj, feature1="nCount_RNA", feature2="nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
dev.off()

# Filter out low quality cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > params$min_genes & nFeature_RNA < params$max_genes & percent.mito < params$max_mito)
length(Idents(seurat_obj))

# Log normalize, find variable features, scale data
#seurat_obj <- NormalizeData(seurat_obj, normalization.method="LogNormalize", scale.factor=10000)
#seurat_obj <- FindVariableFeatures(seurat_obj, selection.method="vst", nfeatures=vst_features)
#seurat_obj <- ScaleData(seurat_obj, features = all.genes)

# Alternatively, use SCTransform which replaces: NormalizeData, FindVariableFeatures, ScaleData
# OLD: seurat_obj <- SCTransform(seurat_obj, vars.to.regress=c('nCount_RNA', 'percent.mito'), variable.features.n=nr_features)

# if (params$max_mito > 50) { params$sct_regress <- union(params$sct_regress, c('percent.mito')) } 
seurat_obj <- SCTransform(seurat_obj, vars.to.regress=params$sct_regress, variable.features.n=params$sct_features, 
                          return.only.var.genes=params$only_var_genes)

# Annotate cell cycle
if ( params$cell_cycle == 'Y' ) {
  if (organism == 'Mus musculus') {
    cc.genes <- readRDS('/mnt/ix1/Resources/scRNA_Ref/CellCycle/mouse_cell_cycle_genes.rds')
  }	
  seurat_obj <- CellCycleScoring(seurat_obj, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=F)
  seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score
  cell_cycle_cts <- as.data.frame(table(seurat_obj@meta.data$Phase))
  colnames(cell_cycle_cts) <- c('Phase', 'Cells')
  write.csv(cell_cycle_cts, file=paste0(seurat_obj@project.name,'.cell_cycle.csv'), row.names=F, quote=F)
}

#add metadata if provided in config file (eg patient, condition) 
if ( metadata ) {
  for (fld in names(params$metadata)) { 
    seurat_obj@meta.data[[fld]] = params$metadata[[fld]] 
  }
} 
saveRDS(seurat_obj, file=paste0(seurat_obj@project.name, '.seurat_obj.rds'))

# Variable features plot
pdf(paste0(seurat_obj@project.name, '_vfeatures_plot.pdf'))
  topN <- head(VariableFeatures(seurat_obj), 12)
  plot1 <- VariableFeaturePlot(seurat_obj)
  LabelPoints(plot=plot1, points=topN, repel=TRUE)
dev.off()

