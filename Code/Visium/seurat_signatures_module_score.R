library(ggbeeswarm)
library(ggpubr)
library(optparse)
library(rstatix)
library(Seurat)
library(tidyverse)
library(viridis)

sessionInfo()
set.seed(42)

# Parse command line arguments  -------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-s', '--sample'), default=NA, type='character', 
                     help='Path to sample Seurat object (.rds)')
parser <- add_option(parser, c('-a', '--annotation'), default=NA, type='character',
                     help='Path to sample Seurat object (.rds) with annotated spots')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for Seurat object with module scores')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

# PROJECT_DIR <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium'
# ANALYSIS_DIR <- file.path(PROJECT_DIR, 'A10_seurat_signatures_module_scores')
# ANN_DIR <- file.path(PROJECT_DIR, 'A09_manual_spot_annotation')
# INPUT_DIR <- file.path(PROJECT_DIR, 'A07_seurat_cluster_QC')
# 
# sample_name <- 'B1_24321.rds'
# # sample_name <- 'A1_24320.rds'
# opt$sample <- file.path(INPUT_DIR, sample_name)
# opt$annotation <- file.path(ANN_DIR, sample_name)
# opt$output_dir <- file.path(file.path(ANALYSIS_DIR, 'temp'))

# Load Visium  ---------------------------------------------------------------

st <- readRDS(opt$sample)
fn <- tools::file_path_sans_ext(basename(opt$sample)) 
prefix <- file.path(opt$output_dir, fn)

# Cropping spatial feature plot distorts image
# https://github.com/satijalab/seurat/issues/5141

img_coords = st@images$slice1@coordinates

get_aspect_ratio <- function(x){
  (max(x$imagerow) - min(x$imagerow)) / (max(x$imagecol) - min(x$imagecol))
}

aspect_ratio <- get_aspect_ratio(img_coords)

st_ann <- readRDS(opt$annotation)

# Load gene signatures ----------------------------------------------------

RNAscope_genes <- c('TFF3', 'ANXA13', 'HKDC1', 'DMBT1', 'OLFM4', 'CPS1', 'ANPEP', 'SLC39A5',
                    'ONECUT2', 'CLDN3', 'CDH17', 'CDX1')

# genes_25 <- read_csv('/mnt/ix1/Projects/M075_201130_GIM_P01/bulkRNA/A04_221222_DEA_WGCNA_validation/spatial_signature_25_genes.csv') %>% pull(...1)
genes_26 <- read_csv('/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/A10_seurat_signatures_module_scores/inputs/26_genes.csv') %>% pull(Symbol)

genes_26

genes_stem <- c('CPS1', 'DMBT1','ONECUT2', 'OLFM4', 'CDX1')
genes_differentiated <- setdiff(genes_26, genes_stem)

filter_genes <- function(x) {
  # Find genes are common with the Visium genes
  intersect(rownames(st@assays$Spatial@counts), x)
}

RNAscope_genes_filtered <- filter_genes(RNAscope_genes)
# genes_25_filtered <- filter_genes(genes_25)
genes_26_filtered <- filter_genes(genes_26)

length(RNAscope_genes_filtered)
length(genes_26_filtered)

# Add module scores -------------------------------------------------------

rm_last_char <- function(x){
  # Remove last character from a string
  str_sub(x, end=-2)
}
DefaultAssay(st) <- 'SCT'

# score_name <- 'RNAscope_12_genes'
# score_name <- 'genes_25'
score_name <- 'genes_26'

gene_set <- genes_26_filtered

st <- AddModuleScore(st, features = list(gene_set), name = score_name)

st@meta.data <- st@meta.data %>%
  rename_with(.fn = rm_last_char, .cols=any_of(paste0(score_name, '1')))

SpatialFeaturePlot(st, features = score_name) + theme(legend.position = "right", aspect.ratio = aspect_ratio)
ggsave(paste(prefix, score_name, 'score_spatial.jpeg', sep='_'), dpi = 300)

SpatialFeaturePlot(st, features = score_name, pt.size.factor = 1, crop = FALSE) + 
  theme(legend.position = "right")
ggsave(paste(prefix, score_name, 'score_spatial_full.jpeg', sep='_'), dpi = 300)

st_subset <- subset(st, cells = colnames(st_ann))
st_subset_img_coords <- st_subset@images$slice1@coordinates
st_subset_aspect_ratio <- get_aspect_ratio(st_subset_img_coords)

SpatialFeaturePlot(st_subset, features = score_name, pt.size.factor = 3, crop= TRUE) + 
  theme(legend.position = "right", aspect.ratio = st_subset_aspect_ratio)
ggsave(paste(prefix, score_name, 'score_subset_spatial_cropped.jpeg', sep='_'), dpi = 300)

SpatialFeaturePlot(st_subset, features = score_name, pt.size.factor = 1, crop= FALSE) + 
  theme(legend.position = "right")
ggsave(paste(prefix, score_name, 'score_subset_spatial_full.jpeg', sep='_'), dpi = 300)

# Boxplot module score ----------------------------------------------------

st_subset <- AddMetaData(st_subset, data.frame(Idents(st_ann)) %>%
              rename('Region'='Idents.st_ann.'))

Idents(st_subset) <- st_subset$Region
st_subset$Region <- factor(st_subset$Region, levels=c('Base', 'Pit', 'Metaplasia'))

SpatialDimPlot(st_subset, label = FALSE, repel = TRUE, pt.size.factor = 1, crop=FALSE)

y_axis_txt = '26-gene signature module score (per spot)'

ModScore <- function(module_score, seu_) {
  cat("Running unpaired statistical test.\n")
  tmp <- FetchData(seu_, vars = module_score)
  tmp <- tibble::rownames_to_column(tmp, var="cell_barcode")
  tmp$ident <- plyr::mapvalues(tmp$cell_barcode, from=colnames(seu_), to=as.character(Idents(seu_)))
  tmp <- setNames(tmp, c("cell_barcode", "module_score", "ident"))
  x <- length(levels(Idents(seu_)))
  if (x == 1) {
    cat("Inadequate groups for statistical testing\n")
  } else if (x == 2) {
    cat("Running Bartlett test for homogeneity of variance\n")
    # subtitle <- paste0('Bartlett test, p = ', p_format(bartlett.test(module_score ~ ident, data = tmp)$p.value, 3) )
    if (bartlett.test(module_score ~ ident, data = tmp)$p.value <= 0.05) {
      cat("P-value ≤0.05. Significant difference in variances. Running unpaired Welch Two Sample t-test\n")
      subtitle <- "Welch's t-test"
      tmp_compare <- tmp  %>%
        t_test(module_score ~ ident, detailed = T) %>%
        add_significance()
      return(list(tmp_compare, subtitle))
    } else {
      cat("P-value > 0.05. Equal variance assumed. Running unpaired Student's Two Sample t-test\n")
      subtitle <- "Student's t-test"
      tmp_compare <- tmp  %>%
        t_test(module_score ~ ident, var.equal = T, detailed = T) %>%
        add_significance()
      return(list(tmp_compare, subtitle))
    }
  } else {
    cat("Running Bartlett test for homogeneity of variance\n")
    if (bartlett.test(module_score ~ ident, data = tmp)$p.value <= 0.05) {
      cat("P-value ≤0.05. Variance is not equal between the groups. Running Kruskal-Wallis test.\n")
      subtitle <- get_test_label(kruskal_test(module_score ~ ident, data = tmp), type = 'text')
      if(kruskal_test(module_score ~ ident, data = tmp)$p<0.05) {
        cat("Significant differences between the groups. Running Dunn's test.\n")
        tmp_compare <- tmp  %>%
          dunn_test(module_score ~ ident, p.adjust.method = "bonferroni", detailed = T) %>%
          add_significance()
        return(list(tmp_compare, subtitle))
      } else {
        cat("Kruskal-Wallis P-val > 0.05. No significant differences between the groups.\n")
      }
    } else {
      cat("P-value > 0.05. Equal variances assumed. Running ANOVA test.\n")
      if (summary(aov(module_score ~ ident, data = tmp))[[1]][[5]][1]<=0.05) {
        cat("Significant differences between the groups. Running Tukey HSD\n")
        tmp_compare <- aov(module_score ~ ident, data = tmp) %>% tukey_hsd()
        return(list(tmp_compare , subtitle))
      } else {
        cat("Kruskal-Wallis P-val > 0.05. No significant differences between the groups.\n")
      }
    }
  }
}

# tmp <- FetchData(st_subset, vars = score_name)
# tmp <- tibble::rownames_to_column(tmp, var="cell_barcode")
# tmp$ident <- plyr::mapvalues(tmp$cell_barcode, from=colnames(st_subset), to=as.character(Idents(st_subset)))
# tmp <- setNames(tmp, c("cell_barcode", "module_score", "ident"))
# 
# p_format(bartlett.test(module_score ~ ident, data = tmp)$p.value, 3)

modscore <- ModScore('genes_26', st_subset)[[1]]
subtitle <- ModScore('genes_26', st_subset)[[2]]
subtitle

ModScore('genes_26', st_subset)
# Convert score_name string to symbol
# https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot
score_name_sym <- sym(score_name)
st_subset@meta.data %>%
  ggplot(aes(x=`Region`, y=!!score_name_sym)) +
    geom_boxplot(aes(col=`Region`), width=0.75, alpha=0.6, fill='white') +
    geom_quasirandom(aes(fill=`Region`), dodge.width = 0.75, pch=21, col='black') +
    geom_hline(yintercept = 0, linetype='dashed') +
    scale_color_manual(breaks = c('Base', 'Pit', 'Metaplasia'),
                       values = c('#4DAF4A', '#377EB8', '#FF7F00')) +
    scale_fill_manual(breaks = c('Base', 'Pit', 'Metaplasia'),
                      values = c('#4DAF4A', '#377EB8', '#FF7F00')) +
  labs(y=y_axis_txt,
       subtitle = subtitle) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14)) +
  stat_pvalue_manual(data = modscore %>% add_y_position())
  # stat_pvalue_manual(data = modscore, y.position = c(1.6, 1.4, 1.5))
ggsave(paste(prefix, score_name, 'score_boxplot.pdf', sep='_'), dpi = 300, width=5, height=5)

# Export spatial plot of each gene in module ------------------------------

pdf(file = paste(prefix, paste0(score_name, '.pdf'), sep='_'))
SpatialFeaturePlot(st, features = score_name) + theme(legend.position = "right", aspect.ratio = aspect_ratio)
for (gene in gene_set){
  p <- SpatialFeaturePlot(st, features = gene) + theme(legend.position = "right", aspect.ratio = aspect_ratio)
  print(p)
}
dev.off()

# Export seurat object ----------------------------------------------------

saveRDS(st_subset, file = paste(prefix, '.rds', sep=''))
print(paste('Exported:', fn, sep=' '))


