library(Seurat)
library(infercnv)
library(readxl)
library(optparse)
library(future)


# Configure environment ---------------------------------------------------

options(future.globals.maxSize = 10000 * 1024^2)
plan('multiprocess', workers=24)

sessionInfo()

# Parse command line arguments --------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-c', '--counts'), default=NA, type='character', 
                     help='Path to counts matrix (.rds)')
parser <- add_option(parser, c('-a', '--annotations'), default=NA, type='character',
                     help='Path to annotations file (.tsv')
parser <- add_option(parser, c('-g', '--gene_pos'), default=NA, type='character',
                     help='Path to Gene position file (.tsv)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for InferCNV')

opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# Load data ---------------------------------------------------------------

fn <- tools::file_path_sans_ext(basename(opt$counts)) 
counts <- readRDS(opt$counts)
print(paste('Loaded:', fn, sep=' '))

# Create InferCNV object -------------------------------------------------------

infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix= counts, 
                                            gene_order_file= opt$gene_pos,
                                            annotations_file= opt$annotations,
                                            delim="\t",
                                            ref_group_names=c('Normal Reference')
                                            ) 

# Run InferCNV ------------------------------------------------------------

infercnv_obj <- infercnv::run(infercnv_obj, 
                              cutoff=0.1, 
                              out_dir= opt$output_dir, 
                              cluster_by_groups= TRUE, 
                              scale_data = TRUE,
                              num_threads=40, 
                              denoise=T,
                              noise_filter = 0.2,
                              HMM=T)

print(paste('Completed InferCNV on:', fn, sep=' '))



