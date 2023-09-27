library(car)
library(dunn.test)
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
parser <- add_option(parser, c('-s', '--sample_dir'), default=NA, type='character', 
                     help='Directory containing Seurat objects (.rds)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory for plots')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

PROJECT_DIR <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium'

opt$sample_dir <- file.path(PROJECT_DIR, 'A10_seurat_signatures_module_scores')
opt$output_dir <- file.path(PROJECT_DIR, 'A10_seurat_signatures_module_scores')


# Load data ---------------------------------------------------------------

metadata <- tibble(
  files = list.files(path=opt$sample_dir, pattern = '*.rds'),
  sample = tools::file_path_sans_ext(files),
)

prefix <- file.path(opt$output_dir, 'All_patients')

# Load each Seurat object
seu_objs <- c()
for (file in metadata$files){
  st <- file.path(opt$sample_dir, file) %>%
    readRDS()
  st$orig.ident <- filter(metadata, files == file)$sample
  seu_objs <- append(seu_objs, st)
}

# Merge seurat objects

st_all <- merge(x = seu_objs[[1]], y = seu_objs[-1], add.cell.ids = metadata$sample )


# Plot module score across patients ---------------------------------------

# score_name = 'RNAscope_12_genes'
score_name = 'genes_26'
y_axis_txt = '26-gene signature module score (per spot)'

st_all$Region <- factor(st_all$Region, levels = c('Base', 'Pit', 'Metaplasia'))

# Convert score_name string to symbol
# https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot
score_name_sym <- sym(score_name)
st_all@meta.data %>%
  ggplot(aes(x=`Region`, y=!!score_name_sym, col=orig.ident, fill=orig.ident)) +
  geom_boxplot(width=0.75, alpha=0.6, fill='white') +
  geom_quasirandom(dodge.width = 0.75, pch=21, col='black') +
  geom_hline(yintercept = 0, linetype='dashed') +
  labs(col='Patient Sample', fill='Patient Sample', y=y_axis_txt) +
  theme_classic()
# ggsave(paste(prefix, score_name, 'score_boxplot_across_patients.pdf', sep='_'), dpi = 300, width=7, height=7)

# Plot module score across regions ----------------------------------------
st_all@meta.data %>%
  ggplot(aes(x=`Region`, y=!!score_name_sym, col=`Region`, fill=`Region`)) +
  geom_boxplot(width=0.75, alpha=0.6, fill='white') +
  geom_quasirandom(dodge.width = 0.75, pch=21, col='black') +
  geom_hline(yintercept = 0, linetype='dashed') +
  scale_color_manual(breaks = c('Base', 'Pit', 'Metaplasia'),
                     values = c('#4DAF4A', '#377EB8', '#FF7F00')) +
  scale_fill_manual(breaks = c('Base', 'Pit', 'Metaplasia'),
                     values = c('#4DAF4A', '#377EB8', '#FF7F00')) +
  labs(y=y_axis_txt) +
  theme_classic()
  # stat_compare_means(comparisons = list(c('Pit', 'Metaplasia'), c('Base', 'Metaplasia')),
  #                    method = 'wilcox.test')
# ggsave(paste(prefix, score_name, 'score_boxplot_across_regions.pdf', sep='_'), dpi = 300, width=7, height=7)


# Compare module score between regions ------------------------------------

tmp <- FetchData(st_all, vars = score_name)
tmp <- tibble::rownames_to_column(tmp, var="cell_barcode")
tmp$ident <- plyr::mapvalues(tmp$cell_barcode, from=colnames(st_all), to=as.character(Idents(st_all)))
tmp <- setNames(tmp, c("cell_barcode", "module_score", "ident"))

## Test for homogeneity of variance ----------------------------------------
# Test if variance is equal across all samples

bartlett.test(module_score ~ ident, data = tmp)
leveneTest(module_score ~ ident, data = tmp)


## Compare groups ----------------------------------------------------------

kruskal.test(module_score ~ ident, data=tmp)
kruskal_test(module_score ~ ident, data=tmp)
dunn.test(tmp$module_score, tmp$ident, 
          method = 'bonferroni', 
          list = TRUE)

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
    if (bartlett.test(module_score ~ ident, data = tmp)$p.value <= 0.05) {
      cat("P-value ≤0.05. Significant difference in variances. Running unpaired Welch Two Sample t-test\n")
      tmp_compare <- tmp  %>%
        t_test(module_score ~ ident, detailed = T) %>%
        add_significance()
      return(tmp_compare)
    } else {
      cat("P-value > 0.05. Equal variance assumed. Running unpaired Student's Two Sample t-test\n")
      tmp_compare <- tmp  %>%
        t_test(module_score ~ ident, var.equal = T, detailed = T) %>%
        add_significance()
      return(tmp_compare)
    }
  } else {
    cat("Running Bartlett test for homogeneity of variance\n")
    if (bartlett.test(module_score ~ ident, data = tmp)$p.value <= 0.05) {
      cat("P-value ≤0.05. Variance is not equal between the groups. Running Kruskal-Wallis test.\n")
      if(kruskal_test(module_score ~ ident, data = tmp)$p.value<0.05) {
        cat("Significant differences between the groups. Running Dunn's test.\n")
        tmp_compare <- tmp  %>%
          dunn_test(module_score ~ ident, p.adjust.method = "bonferroni", detailed = T) %>%
          add_significance()
        return(tmp_compare)
      } else {
        cat("Kruskal-Wallis P-val > 0.05. No significant differences between the groups.\n")
      }
    } else {
      cat("P-value > 0.05. Equal variances assumed. Running ANOVA test.\n")
      if (summary(aov(module_score ~ ident, data = tmp))[[1]][[5]][1]<=0.05) {
        cat("Significant differences between the groups. Running Tukey HSD\n")
        tmp_compare <- aov(module_score ~ ident, data = tmp) %>% tukey_hsd()
        return(tmp_compare)
      } else {
        cat("Kruskal-Wallis P-val > 0.05. No significant differences between the groups.\n")
      }
    }
  }
}

modscore_stats <- ModScore(score_name, st_all)

# modscore_stats <- modscore_stats %>%
#   mutate(p = p.adj,
#          p.signif = p.adj.signif,
#          p.format = p.adj)
# 
# modscore_stats$.y. <- score_name
# 
# modscore_stats
# Plot module scores across regions with stats ----------------------------

st_all@meta.data %>%
  ggplot(aes(x=`Region`, y=!!score_name_sym)) +
  geom_boxplot(aes(col=`Region`), width=0.75, alpha=0.6, fill='white') +
  geom_quasirandom(aes(fill=`Region`), dodge.width = 0.75, pch=21, col='black') +
  geom_hline(yintercept = 0, linetype='dashed') +
  scale_color_manual(breaks = c('Base', 'Pit', 'Metaplasia'),
                     values = c('#4DAF4A', '#377EB8', '#FF7F00')) +
  scale_fill_manual(breaks = c('Base', 'Pit', 'Metaplasia'),
                    values = c('#4DAF4A', '#377EB8', '#FF7F00')) +
  labs(y=y_axis_txt, 
       subtitle = get_test_label(kruskal_test(module_score ~ ident, data = tmp))) +
  stat_pvalue_manual(data = modscore_stats, y.position = c(1.6, 1.4, 1.5)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14))
# stat_compare_means(comparisons = list(c('Pit', 'Metaplasia'), c('Base', 'Metaplasia')),
#                    method = 'wilcox.test')
ggsave(paste(prefix, score_name, 'score_boxplot_across_regions_stats.pdf', sep='_'), dpi = 300, width=5, height=5)

