library(ggplot2)
library(ggpubr)
library(optparse)
library(Seurat)
library(tibble)
library(tidyverse)

sessionInfo()

# Parse command line arguments --------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='Path to sample Seurat object (.rds)')
parser <- add_option(parser, c('-g', '--gene_pos'), default=NA, type='character',
                     help='Gene position file (.bed)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for InferCNV inputs')

opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))


# opt$sample <- list.files(path = '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A07_cluster_annotation/path_ann_for_infercnv',
#                          full.names = TRUE,
#                          pattern = '*.rds')[[1]]
# opt$gene_pos <- '/mnt/ix1/Projects/M070_200622_GI_multiomics/scRNA/00_Code/xiangqi/gene_locs.sorted.bed'
# opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/F01_seurat_to_InferCNV/temp'

# Load data ---------------------------------------------------------------

st <- readRDS(opt$sample)

fn <- tools::file_path_sans_ext(basename(opt$sample)) 
prefix <- paste(opt$output_dir, fn, sep='/')

# Cropping spatial feature plot distorts image
# https://github.com/satijalab/seurat/issues/5141

img_coords = st@images$slice1@coordinates
aspect_ratio <- (max(img_coords$imagerow) - min(img_coords$imagerow)) / (max(img_coords$imagecol) - min(img_coords$imagecol))

gene_pos <- read_tsv(opt$gene_pos, col_names = c('Gene', 'chrom', 'chromStart', 'chromEnd')) 

# Vis clusters ------------------------------------------------------------

SpatialDimPlot(st, label = TRUE, label.size = 2) + theme(aspect.ratio = aspect_ratio)
ggsave(paste(prefix, 'spatial_plot.jpeg', sep='_'), dpi = 300)


# Export counts -----------------------------------------------------------
# InferCNV File format: https://github.com/broadinstitute/inferCNV/wiki/File-Definitions

counts_raw <- st@assays$Spatial@counts
saveRDS(counts_raw, paste(prefix, 'counts_raw.rds', sep='_'))

counts_sct <- st@assays$SCT@counts
saveRDS(counts_sct, paste(prefix, 'counts_SCT.rds', sep='_'))


# Export annotations ------------------------------------------------------
# Apparently InferCNV doesn't like numerical cluster numbers: https://github.com/broadinstitute/infercnv/issues/223

annotations <- Idents(st)
colnames(annotations) <- NULL
write.table(annotations, paste(prefix, 'annotations.tsv', sep='_') , sep = '\t',
            quote = FALSE, row.names = TRUE, col.names = FALSE)

# Export gene ordering file -----------------------------------------------

gene_pos_raw <- gene_pos %>% 
  filter(Gene %in% rownames(counts_raw))
write_tsv(gene_pos_raw, file = paste(prefix, 'gene_pos_raw.tsv', sep = '_'), col_names = FALSE)

gene_pos_sct <- gene_pos %>%
  filter(Gene %in% rownames(counts_sct))
write_tsv(gene_pos_sct, file = paste(prefix, 'gene_pos_SCT.tsv', sep = '_'), col_names = FALSE)

print(paste('Exported:', fn, sep=' '))
