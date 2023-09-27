#running in /home/asathe/shared/seu_4/ virtual environment that has Harmony, seurat 4
#only batch i.e. patient ID used to integrate in harmony. condition NOT included as a variable 
library(harmony)
library(Seurat)
library(viridis)
library(ggplot2)
library(dplyr)
library(future)

dims.reduce = 20
cluster.res = 0.6
ram_gb = 48
options(future.globals.maxSize = ram_gb * 1024^3)
plan('multiprocess', workers=12)

#loads merged object and continues downstream Seurat on it
load('merged.rds')
aggr <- seu_merged
aggr <- SCTransform(aggr) %>% RunPCA()
aggr <- RunHarmony(aggr, group.by.vars = "patientID", assay.use = "SCT")
aggr <- RunUMAP(aggr, reduction = "harmony", dims = 1:dims.reduce)
aggr <- FindNeighbors(aggr, reduction = "harmony", dims = 1:dims.reduce) %>% FindClusters(resolution=cluster.res)
saveRDS(aggr, file = "aggr.rds", version = 2)

plan('sequential')
png("umap1.png")
DimPlot(aggr, label =TRUE)
dev.off()

png("umap2.png")
DimPlot(aggr, group.by = "condition",label =TRUE)
dev.off()

png("umap3.png")
DimPlot(aggr, group.by = "patientID",label =TRUE)
dev.off()

cellnumbers <- table(Idents(aggr), aggr@meta.data$orig.ident)
write.csv(cellnumbers, file="cellnumbers.csv")

DefaultAssay(aggr) <- "RNA"
plan('multiprocess', workers=12)
# Normalize, scale RNA slot in preparation for find markers
aggr <- NormalizeData(aggr, assay="RNA")
aggr <- FindVariableFeatures(aggr, assay="RNA")
aggr <- ScaleData(aggr, assay="RNA", features=rownames(aggr), vars.to.regress="nCount_RNA")
saveRDS(aggr, file = "aggr.rds", version = 2)

markers_RNA <- FindAllMarkers(aggr, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_RNA, file = "markers_RNA.csv")

cell_markers<- c("EPCAM", "TFF3", "DCN", "ACTA2", "VWF", "CD14", "CD68", "CD3D", "CD8A","NKG7","CD79A","CD19")
immune_markers <- c("CD3D", "CD4","IL7R", "FOXP3", "IL2RA","CD8A", "CD8B","GNLY", "NKG7","TRDC", "SELL", "CCR7", "CCR6", "KLRB1", "CD19", "MS4A1", "SDC1")

pdf("cell_markers.pdf")
DoHeatmap(aggr, features=cell_markers) + scale_fill_viridis(option="D")
dev.off()

pdf("immune_markers.pdf")
DoHeatmap(aggr, features=immune_markers) + scale_fill_viridis(option="D")
dev.off()

top5 <- markers_RNA %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

pdf("top5_cluster.pdf")
DoHeatmap(aggr, features = top5$gene) + scale_fill_viridis(option = "D")
dev.off()