library(edgeR)
library('optparse')
library(tidyverse)
library(Seurat)

#################################
# Parse command line arguments #
#################################

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample_dir'), default=NA, type='character', 
                     help='Directory containing Seurat objects (.rds)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for pseudobulk matrix')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# opt$sample_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A04_cluster_annotation/path_annotation'
# opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A05_pseudobulk/temp'

#############
# Load data #
#############

metadata <- tibble(
  files = list.files(path=opt$sample_dir, pattern = '*.rds'),
  sample = tools::file_path_sans_ext(files),
)

# Load each Seurat object
seu_objs <- c()
for (file in metadata$files){
  st <- file.path(opt$sample_dir, file) %>%
    readRDS()
  st$orig.ident <- filter(metadata, files == file)$sample
  seu_objs <- append(seu_objs, st)
}

### Merge seurat objects ### 

st_all <- merge(x = seu_objs[[1]], y = seu_objs[-1], add.cell.ids = metadata$sample )

##############
# Pseudobulk #
##############

counts_list <- list()

for (i in 1:nrow(metadata)) {
  sample_name <- pull(metadata, sample)[i]
  st <- subset(st_all, subset = orig.ident == sample_name)
  for (cell_type in levels(st)) {
    st_subset <- subset(st, idents = cell_type)
    pseudobulk <- apply(st_subset@assays$Spatial@counts, 1, sum)
    experiment <- paste(sample_name, cell_type, sep='_')
    #     print(experiment)
    counts_list[[experiment]] <- pseudobulk
  }
}

counts_raw <- as.data.frame(counts_list) %>% as_tibble(rownames = 'Gene')

counts_cpm <- counts_raw %>% select(-(Gene)) %>% as.matrix()
rownames(counts_cpm) <- counts_raw %>% pull(Gene)
counts_cpm <- counts_cpm %>% cpm() %>% as_tibble(rownames = 'Gene')

####################################
# Export Pseudobulk count matrices #
####################################

fn_prefix = paste('pseudobulk', basename(opt$sample_dir), sep='_')

fn_raw = paste(fn_prefix, 'raw.txt', sep='_')
fn_cpm = paste(fn_prefix, 'cpm.txt', sep='_')

write.table(counts_raw, file = file.path(opt$output_dir, fn_raw), sep='\t', row.names = FALSE)
write.table(counts_cpm, file = file.path(opt$output_dir, fn_cpm), sep='\t', row.names = FALSE)

print(paste('Exported:', fn_raw, sep=' '))
print(paste('Exported:', fn_cpm, sep=' '))
