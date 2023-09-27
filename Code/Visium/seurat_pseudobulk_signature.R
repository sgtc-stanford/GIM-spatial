
# Setup environment -------------------------------------------------------

library(ComplexHeatmap)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(optparse)
library(RColorBrewer)
library(Seurat)
library(tidyverse)

set.seed(42)
sessionInfo()


# Parse command line arguments --------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-c', '--counts'), default=NA, type='character', 
                     help='Path to counts matrix (genes x samples)')
parser <- add_option(parser, c('-s', '--samples'), default=NA, type='character',
                     help='Path to sample metadata table(.csv)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

opt$counts <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/C01_pseudobulk/pseudobulk_A08_manual_spot_annotation_raw.txt'
opt$samples <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/samples.csv'
opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/C04_pseudobulk_signature'

bulk_signature_path <- '/mnt/ix1/Projects/M075_201130_GIM_P01/bulkRNA/C01_20230608_Intersection_discovery_validation/temp/validated_genes_c5_for_spatial_mapping.csv'
metaplasia_overexpressed_path <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/C03_pseudobulk_DE/metaplasia_overexpressed_genes.csv'

# Load data ---------------------------------------------------------------

counts <- read_tsv(opt$counts)
samples <- read_csv(opt$samples)

bulk_signature_genes <- read.csv(bulk_signature_path) %>% pull(x)
bulk_signature_genes

metaplasia_overexpressed_genes <- read_csv(metaplasia_overexpressed_path) %>% pull(x)
metaplasia_overexpressed_genes

fn <- tools::file_path_sans_ext(basename(opt$counts)) 
prefix <- paste(opt$output_dir, fn, sep='/')


# Normalise gene expression -----------------------------------------------

counts_df <- counts %>% 
  select(-(Gene)) %>% 
  as.data.frame()
rownames(counts_df) <- counts %>% pull(Gene)

pseudobulk_metadata <- tibble(
  pseudobulk = colnames(counts_df)
)

pseudobulk_metadata <- pseudobulk_metadata %>% 
  separate(pseudobulk, c('array', 'LIMS_dissection', 'group'), sep='_', remove=FALSE) %>%
  mutate(LIMS_dissection = as.double(LIMS_dissection)) %>%
  left_join(samples, by='LIMS_dissection')

group <- factor(pseudobulk_metadata$group,
                levels=c('Base', 'Pit', 'Metaplasia'),
                ordered=TRUE)

y <- DGEList(counts=counts_df, group = group, samples = pseudobulk_metadata %>% 
               mutate(dissection = factor(LIMS_dissection)) %>% as.data.frame())
y$samples$array <- factor(pseudobulk_metadata$array.x)
y$samples$dissection <- factor(pseudobulk_metadata$LIMS_dissection)
y$samples$region <- factor(pseudobulk_metadata$region)
y$samples$metaplasia <- factor(pseudobulk_metadata$metaplasia)
y$samples$Stage <- factor(pseudobulk_metadata$OLGIM, levels = unique(pseudobulk_metadata$OLGIM))
y$samples

counts_lcpm <- cpm(y, log=TRUE)

counts_lcpm

# Z-score along genes
counts_lcpm_z <- t(scale(t(counts_lcpm)))

# Pseudobulk signature ----------------------------------------------------

signature <- intersect(bulk_signature_genes, metaplasia_overexpressed_genes)

length(signature)


# Plot pseudobulk signature heatmap ---------------------------------------


ha <- y$samples %>%
  mutate(Region = recode_factor(group, 
                                'Base' = 'Normal Base',
                                'Pit' = 'Normal Pit'))
# ha$Stage <- addNA(ha$Stage)

ha$Stage <- recode_factor(ha$Stage, `2`='OLGIM II', `3`='OLGIM III')
levels(ha$Stage) <- c(levels(ha$Stage), 'Gastric Cancer')
ha$Stage <- replace_na(ha$Stage, 'Gastric Cancer')
ha$Stage
ha

column_ha <- HeatmapAnnotation(df = as.data.frame(select(ha, Stage, Region)),
                               col = list(Stage = c('OLGIM II' = brewer.pal(11, "RdYlBu")[9],
                                                            'OLGIM III' = brewer.pal(11, "RdYlBu")[10],
                                                            'Gastric Cancer' = brewer.pal(9, 'Reds')[6]),
                                          Region = c('Normal Base' = "#4DAF4A",
                                                     'Normal Pit' = "#377EB8",
                                                     'Metaplasia' = "#FF7F00")
                               ),
                               simple_anno_size = unit(3, 'mm'),
                               annotation_name_gp = gpar(fontsize = 10)
)

# labels <- RNAscope_genes
# labels_indices <- which(rownames(counts_lcpm_z[common.degs,]) %in% labels)
# 
# ha_labels = rowAnnotation(foo = anno_mark(at=labels_indices, labels = labels))

sample_names <- tibble(old_names = colnames(counts_lcpm_z)) %>%
  separate(old_names, c('array', 'LIMS_dissection', 'group'), sep='_', remove=FALSE) %>%
  mutate(LIMS_dissection = as.double(LIMS_dissection)) %>%
  left_join(samples, by='LIMS_dissection') %>%
  rowwise() %>%
  mutate(LIMS_patient_long = paste0('P0', LIMS_patient)) %>%
  unite('new_names', c(LIMS_patient_long,group), sep = ' - ')

signature_mat <- counts_lcpm_z[signature,]
colnames(signature_mat) <- sample_names$new_names

col <- as.vector(recode_factor(y$samples$group, 
                               Base = "#4DAF4A", 
                               Pit = "#377EB8", 
                               Metaplasia = "#FF7F00"))

hmap <- Heatmap(signature_mat, name='z-score',
        row_names_gp = gpar(fontsize=9),
        column_title = 'Pseudobulk Samples', column_title_side = 'top',
        column_names_gp = gpar(fontsize = 11, col=col),
        clustering_method_rows = 'ward.D2', clustering_distance_rows = 'manhattan',
        clustering_method_columns = 'ward.D2', clustering_distance_columns = 'manhattan',
        column_split = 2,
        top_annotation = column_ha)

pdf(file = file.path(opt$output_dir, 'Heatmap_36_genes_7x8.pdf'), height=7, width=8)
draw(hmap)
dev.off()


# # Find and export signature genes significantly overexpressed in metaplasia ---
# length(signature_genes)
# genes_35 <- intersect(signature_genes, common.degs)
# genes_35_lst <- list(metaplasia.vs.base = data.frame(symbol = genes_35,  metaplasia.vs.base[genes_35,]),
#                      metaplasia.vs.pit = data.frame(symbol = genes_35, metaplasia.vs.pit[genes_35,]),
#                      base.vs.pit = data.frame(symbol = genes_35, base.vs.pit[genes_35,]))
# openxlsx::write.xlsx(genes_35_lst, file.path(opt$output_dir, '35_genes.xlsx'))
# 
# 
# # Heatmap of 35 genes -----------------------------------------------------
# 
# hmap_35 <- Heatmap(counts_lcpm_z[genes_35,], name='z-score', 
#                    row_title = 'Genes', row_names_gp = gpar(fontsize=8), show_row_names = T,
#                    column_title = 'Pseudobulk Samples', column_title_side = 'top', column_names_gp = gpar(fontsize=11, col= col),
#                    clustering_method_rows = 'ward.D2', clustering_distance_rows = 'manhattan', row_dend_width = unit(1.5, 'cm'),
#                    clustering_method_columns = 'ward.D2', clustering_distance_columns = 'manhattan',
#                    cluster_columns= T,
#                    column_split = 2,
#                    top_annotation = column_ha)
# 
# draw(hmap_35)

# pdf(file = file.path(opt$output_dir, 'Heatmap_35_overexpressed_signature_genes.pdf'), height=6.25, width=10.14)
# draw(hmap_35)
# dev.off()
