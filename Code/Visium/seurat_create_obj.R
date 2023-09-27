library('optparse')
library(Seurat)

# Parse command line arguments

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample_dir'), default=NA, type='character', 
                     help='cell ranger output directory')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for Seurat objects')
parser <- add_option(parser, c('-n', '--output_name'), default=NA, type='character',
                     help='Name of output .rds file containing seurat object')

opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# Create Seurat Object

st <- Load10X_Spatial(data.dir = opt$sample_dir)

# Save Seurat Object
fn = paste(opt$output_name, '.rds', sep='')
saveRDS(st, file=paste(opt$output_dir, fn, sep='/'))

print(paste('Exported:', opt$output_name))