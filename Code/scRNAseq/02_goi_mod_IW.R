# 00 - Load required R packages -------------------------------------------

library(optparse)
library(Seurat)
library(dplyr)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)


# 01 - Set up working directory -------------------------------------------

setwd("./F02_IW_signature_GC")
if(dir.exists("./temp")) {
  cat("Output files will be saved in existing temp folder.\n")
} else {
  cat("Output files will be saved in new temp folder.\n")
  dir.create("./temp")
}
if(dir.exists("./plots")) {
  cat("Plots will be saved in existing plots folder.\n")
} else {
  cat("Plots will be saved in new plots folder.\n")
  dir.create("./plots")
}


# 02 - Parse command line arguments ---------------------------------------

parser <- OptionParser()
parser <- add_option(parser,
                     c("-i", "--input_scRNAseq"),
                     default=NA,
                     type="character", 
                     help="Path to scRNAseq files")
parser <- add_option(parser,
                     c("-s", "--gene_signature"),
                     default=NA,
                     type="character", 
                     help="Path to scRNAseq files")
parser <- add_option(parser,
                     c("-o", "--output_dir"),
                     default=NA,
                     type="character",
                     help="Output directory for processed count files")
parser <- add_option(parser,
                     c("-p", "--plot_dir"),
                     default=NA,
                     type="character",
                     help="Output directory for plots")
opt <- parse_args(parser,
                  args = commandArgs(trailingOnly = TRUE))

opt$input_scRNAseq <- "../../../M039_170911_Gastric_scRNA/A00_interpatient_merge/P5846_5866_5931_6207_6342_6649_6592_6709_N_T_PBMC/A00_seurat_regressnUMI_GEMbatch/subset_epithelial/A00_seurat_20PC_res0.8"
opt$gene_signature <- "../../bulkRNA/D01_20230720_TCGA_validation/temp"
opt$output_dir <- "./temp"
opt$plot_dir <- "./plots"



# 03 - Load scRNAseq object and gene signature ----------------------------

load(file.path(opt$input_scRNAseq, "subset_filter.Robj"))


## 3.1 - Inspect the data -------------------------------------------------

class(subset_filter)
seu <- UpdateSeuratObject(subset_filter)
head(seu@meta.data)
table(seu@meta.data$orig.ident, seu@meta.data$condition)
#sanity check
DimPlot(seu, group.by = "condition")
table(seu@meta.data$class, seu@meta.data$condition)


## 3.2 - Remove metaplasia P6649 ------------------------------------------

length(Idents(seu))
seu <- subset(seu, subset = patientID == "6649", invert=TRUE)
seu <- SCTransform(seu) %>% RunPCA()
length(Idents(seu))

goi <- read.csv(file.path(opt$gene_signature, "26_validated_spatial_genes.csv"))[,2]

seu <- AddModuleScore(seu,
                      features = goi,
                      assay = "SCT",
                      name = "goi")
head(seu@meta.data)
Idents(seu) <- "condition"
my_levels <- c("normal",  "tumor")
Idents(seu) <- factor(Idents(seu), levels= my_levels)
saveRDS(seu, "seu_epi_GC.rds")

VlnPlot(seu, features = c("goi1"))
FeaturePlot(seu, features = c("goi1"))
DimPlot(seu, group.by = "condition")

DefaultAssay(seu) <- "SCT"

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
      if(kruskal_test(module_score ~ ident, data = tmp)$p<0.05) {
        cat("Significant differences between the groups. Running Dunn's test.\n")
        tmp_compare <- tmp  %>% 
          dunn_test(module_score ~ ident,
                    p.adjust.method = "BH",
                    detailed = T) %>%
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


# 04 - Statistical testing ------------------------------------------------

score_pval <- ModScore(module_score = "goi1", seu_ = seu)

score_pval$group1 <- ifelse(score_pval$group1=="normal", "Control", score_pval$group1)
score_pval$group2 <- ifelse(score_pval$group2=="tumor", "Tumor", score_pval$group2)

write.csv(score_pval, file = file.path(opt$output_dir, "T-test TvsN scRNA.csv"), row.names = F)

t_test_tumor_normal_plot <- seu@meta.data %>%
  mutate(condition = ifelse(condition=="normal", "Control", "Tumor")) %>%
  mutate(condition = factor(condition, levels = c("Control", "Tumor"), ordered = T)) %>%
  ggviolin(x = "condition", y = "goi1", fill = "condition") +
  geom_hline(yintercept = 0, lty = 3) +
  stat_pvalue_manual(score_pval, 
                     y.position = 2,
                     step.increase = 0.1,
                     label = "p.signif") +
  labs(subtitle = get_test_label(score_pval)) +
  geom_boxplot(fill = "white", width = 0.3, alpha = 0.7) +
  ggtitle("Tumor vs normal\n(epithelial cells)") +
  scale_fill_manual(values = c("Tumor" = brewer.pal(5, "Reds")[5],
                               "Control" = brewer.pal(5, "Reds")[4])) +
  xlab("") +
  ylab("Module score (26 genes)") +
  labs(fill = "Condition") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

precancer_kruskal_plot <- readRDS(file = "../F01_IW_GIM/temp/Figure_4F.rds")

precancer_kruskal_plot +
  t_test_tumor_normal_plot +
  plot_layout(widths = c(3, 1))

VlnPlot(seu,
        sort = "decreasing",
        assay = "RNA", # This is to use the RNA or SCT assays from the Seurat object
        features = "goi1",
        group.by = "condition",
        pt.size = 0) +
  geom_boxplot(width = 0.3,
               col = "black",
               fill = "white",
               alpha = 0.7) +
  ggtitle("Expression of the 26-gene signature") +
  xlab("Cell type") +
  ylab("Gene module score") +
  scale_fill_brewer(palette = "Reds", 5) +
  # stat_pvalue_manual(score_pval, 
  #                    y.position = 1,
  #                    step.increase = 0.1,
  #                    label = "p.value") +
  theme(legend.position = "none")

write.csv(sig, "modscore_ttest.csv")

