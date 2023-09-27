library(optparse)
library(RColorBrewer)
library(Seurat)
library(tidyverse)

sessionInfo()

set.seed(42)

# Parse command line arguments  -------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-r', '--reference'), default=NA, type='character', 
                     help='Directory containing scRNAseq reference dataset as Seurat object')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# Load data ---------------------------------------------------------------

opt$reference <- '/mnt/ix1/Projects/M075_201130_GIM_P01/scRNA/C01_harmony/aggr.final_celltype.rds'
opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/D01_scRNA_reference_QC'

ref <- readRDS(opt$reference)

fn <- tools::file_path_sans_ext(basename(opt$reference)) 
prefix <- paste(opt$output_dir, fn, sep='/')

# Vis reference -----------------------------------------------------------

DefaultAssay(ref) <- 'SCT'

DimPlot(ref, group.by = 'final_celltype')
ggsave(paste(prefix, 'pre_QC_umap.jpeg', sep='_'), dpi = 300)

cell_type_counts <- table(ref$final_celltype)
cell_type_counts_tbl <- tibble(cell_type = names(cell_type_counts),
                               count = cell_type_counts) %>% 
  arrange(desc(count)) %>%
  mutate(cell_type = factor(cell_type, levels = unique(cell_type)))

ggplot(data = cell_type_counts_tbl, mapping = aes(x= cell_type, y = count)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_text(aes(label=count), vjust=-0.3, size=3.5) +
  NoLegend()
ggsave(paste(prefix, 'pre_QC_cell_type_counts.jpeg', sep='_'), dpi = 300)

# Filter scRNA-seq reference ----------------------------------------------

# Remove cell type with less than 25 counts (RCTD default min cells)
cell_type_subset <- names(cell_type_counts[cell_type_counts > 25])

ref_filtered <- subset(x = ref, subset = final_celltype %in% cell_type_subset)


# Vis post QC -------------------------------------------------------------

DimPlot(ref_filtered, group.by = 'final_celltype')
ggsave(paste(prefix, 'post_QC_umap.jpeg', sep='_'), dpi = 300)

cell_type_counts_filtered <- table(ref_filtered$final_celltype)
cell_type_counts_filtered_tbl <- tibble(cell_type = names(cell_type_counts_filtered),
                               count = cell_type_counts_filtered) %>% 
  arrange(desc(count)) %>%
  mutate(cell_type = factor(cell_type, levels = unique(cell_type)))

ggplot(data = cell_type_counts_filtered_tbl, mapping = aes(x= cell_type, y = count)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_text(aes(label=count), vjust=-0.3, size=3.5) +
  NoLegend()
ggsave(paste(prefix, 'post_QC_cell_type_counts.jpeg', sep='_'), dpi = 300)


# Export reference --------------------------------------------------------

saveRDS(ref_filtered, file = paste(prefix, '.rds', sep=''))
print(paste('Exported:', fn, sep=' '))
