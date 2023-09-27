# File: seurat_markers.R
# Author: Anuja Sathe
# Desc: Script reads saved Seurat R object (with PCA calculated) and runs UMAP
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parse input parameters                                                      #                                                 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
suppressPackageStartupMessages(require(optparse))
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf

option_list = list(
  make_option(c("-r", "--seurat"), action="store", default=NA, type='character',
              help='Seurat/UMAP R object')
  )
opt = parse_args(OptionParser(option_list=option_list))

cat("User input >>","\n")
cat("seurat:", opt$seurat,"\n")
fn_in   = opt$seurat

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load packages, and custom helper functions                                  #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)
library(pracma)  
library(viridis)
library(stringr)

source('/mnt/ix1/Resources/LabSoftware/software/scRNA-v3/pipeline/_scRNA_methods.R')

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load Seurat object and find differentially expressed genes between clusters # 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
if ( str_sub(fn_in, start=-4) == 'Robj' ) {
  seurat_umap <- loadRObj(fn_in)
} else {
  seurat_umap <- readRDS(fn_in)
}  

#print("Variable Genes found are below:")
#print(seurat_umap@var.genes)

#Find cluster markers from RNA slot
DefaultAssay(seurat_umap) <- "RNA"
all.markers <- FindAllMarkers(seurat_umap, assay="RNA", min.pct=0.25, logfc.threshold=0.25)
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)

write.csv(all.markers, file=paste0(seurat_umap@project.name, '.all_markers.csv'))
#write.csv(top20.markers, file=paste0(seurat_umap@project.name, '.top20_markers.csv'))

top5.markers <- top20.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
pdf(paste0(seurat_umap@project.name, '.top5_markers.pdf'))
  DoHeatmap(seurat_umap, features=as.vector(top5.markers$gene)) + theme(axis.text.y = element_text(size=5))
dev.off()

cell_markers<- c("EPCAM", "TFF3", "DCN", "ACTA2", "VWF", "CD14", "CD68", "CD3D", "CD8A","NKG7","CD79A","CD19")
immune_markers <- c("CD3D", "CD4","IL7R", "FOXP3", "IL2RA","CD8A", "CD8B","GNLY", "NKG7","TRDC", "SELL", "CCR7", "CCR6", "KLRB1", "CD19", "MS4A1", "SDC1")

pdf("cell_markers.pdf")
DoHeatmap(seurat_umap, features=cell_markers) + scale_fill_viridis(option="D")
dev.off()

pdf("immune_markers.pdf")
DoHeatmap(seurat_umap, features=immune_markers) + scale_fill_viridis(option="D")
dev.off()
