# 00 - Load required R packages -------------------------------------------

library(optparse)
library(tidyverse)
library(Seurat)
library(openxlsx)
library(ComplexHeatmap)
library(harmony)
library(RColorBrewer)
library(rstatix)
library(ggpubr)
library(SingleR)
library(BiocParallel)
library(forcats)

setwd("/mnt/ix1/Projects/M075_201130_GIM_P01/scRNA/F01_IW_GIM/")


# 01 - Set up working directory -------------------------------------------

setwd("./F01_IW_GIM")
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

opt$input_scRNAseq <- "../C01_harmony"
opt$gene_signature <- "../../bulkRNA/D01_20230720_TCGA_validation/temp/26_validated_spatial_genes.csv"
opt$output_dir <- "./temp"
opt$plot_dir <- "./plots"


# 03 - Load scRNAseq data -------------------------------------------------

sc_data <- readRDS(file.path(opt$input_scRNAseq,
                             "aggr.final_celltype.rds"))
table(sc_data@meta.data$celltype,
      sc_data@meta.data$final_celltype)

summary_metadata_scRNAseq <- data.frame(aggregate(sc_data@meta.data$nCount_RNA,
                                                  list(paste(sc_data@meta.data$orig.ident,
                                                             sc_data@meta.data$rev_condition)),
                                                  FUN = mean))

summary_metadata_scRNAseq2 <- data.frame(aggregate(sc_data@meta.data$nFeature_RNA,
                                                   list(paste(sc_data@meta.data$orig.ident,
                                                              sc_data@meta.data$rev_condition)),
                                                   FUN = mean))

summary_metadata_scRNAseq3 <- data.frame(table(paste(sc_data@meta.data$orig.ident,
                                                     sc_data@meta.data$rev_condition)))

if(identical(summary_metadata_scRNAseq$Group.1, summary_metadata_scRNAseq2$Group.1) &
   identical(summary_metadata_scRNAseq$Group.1, as.character(summary_metadata_scRNAseq3$Var1))) {
  summary_metadata_scRNAseq <- data.frame(
    patient = str_split(summary_metadata_scRNAseq$Group.1,
                       " ",
                       simplify = T)[,1],
    condition = str_split(summary_metadata_scRNAseq$Group.1,
                         " ",
                         simplify = T)[,2],
    number_of_cells = summary_metadata_scRNAseq3$Freq,
    mean_number_of_features_per_cell = round(summary_metadata_scRNAseq2$x),
    total_UMI_counts = data.frame(aggregate(sc_data@meta.data$nCount_RNA,
                                        list(paste(sc_data@meta.data$orig.ident,
                                                   sc_data@meta.data$rev_condition)),
                                        FUN = sum))$x,
    mean_UMI_counts_per_cell = round(summary_metadata_scRNAseq$x)
  )
  write.xlsx(summary_metadata_scRNAseq,
             file.path(opt$output_dir,
                       "Summary_metadata.xlsx"))
  rm(summary_metadata_scRNAseq,
     summary_metadata_scRNAseq2,
     summary_metadata_scRNAseq3)
}

sc_metadata <- subset(sc_data@meta.data,
                      select = c(orig.ident, rev_condition))
sc_metadata <- merge(data.frame(table(sc_metadata$orig.ident)),
                     sc_metadata,
                     by.x = "Var1",
                     by.y = "orig.ident",
                     all.x = T)
sc_metadata <- sc_metadata[!duplicated(sc_metadata$Var1),]
sc_metadata$patient <- stringr::str_split(sc_metadata$Var1, "_", simplify = T)[,1]
table(Patient = sc_metadata$patient,
      Condition = sc_metadata$rev_condition)
length(unique(sc_metadata$patient))
# [1] 19

# 04 - SCT ----------------------------------------------------------------

DefaultAssay(sc_data) <- "SCT"
gene_signature <- read.csv(opt$gene_signature)$x
sc_data <- AddModuleScore(sc_data,
                          list(gene_signature),
                          name = "Spatial_genes")
# To remove the score, uncomment the following line
# sc_data@meta.data <- sc_data@meta.data[,str_detect(colnames(sc_data@meta.data), "Spatially", negate = T)]


# 05 - Inspect all lineages -----------------------------------------------


## 5.1 UMAP major lineages ------------------------------------------------

# Format lineage metadata
sc_data@meta.data <- sc_data@meta.data %>%
  mutate(final_celltype = gsub("_", " ", final_celltype)) %>%
  mutate(final_celltype = paste0(toupper(substring(final_celltype, 1, 1)),
                                 substring(final_celltype, 2, nchar(final_celltype)))) %>%
  mutate(final_celltype = gsub("LYZ ", "LYZ-", final_celltype)) %>%
  mutate(final_celltype = gsub("REG3A  ", "REG3A-", final_celltype))

# Truncate axes to make UMAP arrows smaller in R
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

# Generate plot
umap_major_lineages <- DimPlot(sc_data, group.by = 'celltype') +
  scale_color_brewer(palette = "Paired") +
  ggtitle("Major lineages") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme_classic() +
  guides(x = axis,
         y = axis,
         col = guide_legend(ncol = 2,
                            override.aes = list(size = 3))) +
  theme(axis.line = element_line(arrow = arrow(angle = 15,
                                               length = unit(.15,
                                                             "inches"),
                                               type = "closed")),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right")


## 5.2 - UMAP module score ------------------------------------------------

module_major_lineages <- FeaturePlot(sc_data,
                          features =  "Spatial_genes1") +
  scale_color_gradient(low = "skyblue",
                       high = "firebrick") +
  ggtitle("Signature module score across major lineages") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme_classic() +
  guides(x = axis,
         y = axis) +
  theme(axis.line = element_line(arrow = arrow(angle = 15,
                                               length = unit(.15,
                                                             "inches"),
                                               type = "closed")),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right")



## 5.3 - Compare module score between major lineages ----------------------

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

Idents(sc_data) <- "celltype"
score_pval_all_cells <- ModScore(module_score = "Spatial_genes1",
                                 seu_ = sc_data)

score_pval_epi_cells <- subset(score_pval_all_cells,
                               group1%in%"epithelial"|group2%in%"epithelial")

score_pval_epi_cells <- score_pval_epi_cells[c(6,1,5,2,4,3),]

major_lineage_kruskal_plot <- sc_data@meta.data %>%
  ggviolin(x = "celltype",
           y = "Spatial_genes1",
           fill = "celltype") +
  stat_pvalue_manual(score_pval_epi_cells,
                     y.position = 0.8,
                     step.increase = 0.1,
                     label = "p.adj.signif") +
  labs(subtitle = get_test_label(kruskal_test(Spatial_genes1 ~ rev_condition,
                                              data = sc_data@meta.data)),
       caption = "**** Dunn's test FDR-adjusted p ≤ 0.001") +
  ggtitle("Signature module score between major lineages") +
  geom_boxplot(fill = "white",
               width = 0.3,
               alpha = 0.7) +
  geom_hline(yintercept = median(sc_data@meta.data$Spatial_genes1[sc_data@meta.data$celltype=="epithelial"]),
             lty = 3) +
  scale_fill_brewer(type = "seq",
                    palette = "Reds") +
  xlab("") +
  ylab("Module score (26 genes)") +
  labs(fill = "Condition") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.5,
                                   hjust=1))

(umap_major_lineages + module_major_lineages) / major_lineage_kruskal_plot


## 5.4 - Compare module score between minor lineages ----------------------

col <- paletteer::paletteer_d("ggsci::category20_d3")[1:17]
names(col) = unique(sc_data$final_celltype)

col_fun = circlize::colorRamp2(c(-2,
                                 0,
                                 2),
                               c("navy", "white", "firebrick"))


# 06 - Epithelial cell subset ---------------------------------------------


## 6.1 - Subset the dataset -----------------------------------------------

ref_epi <- subset(sc_data,
                  subset = celltype == 'epithelial')


## 6.2 - SCTransform ------------------------------------------------------

ref_epi <- SCTransform(ref_epi) %>% RunPCA()


## 6.3 - Run Harmony for batch correction ---------------------------------

ref_epi <- RunHarmony(ref_epi,
                      group.by.vars = "patientID",
                      assay.use = "SCT")

ref_epi@reductions$harmony


## 6.4 - Run UMAP ---------------------------------------------------------

ref_epi <- RunUMAP(ref_epi,
                   reduction = "harmony",
                   dims = 1:20)

ref_epi <- FindNeighbors(ref_epi,
                         reduction = "harmony",
                         dims = 1:20) %>%
  FindClusters(resolution = 0.6)

FeaturePlot(ref_epi, features =  "Spatial_genes1") +
  scale_color_gradient(low = "skyblue", high = "firebrick") +
  ggtitle("26-gene signature module score:\nEpithelial cell subset") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(angle = 15,
                                               length = unit(.15,
                                                             "inches"),
                                               type = "closed")),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))

Idents(ref_epi) <- ref_epi$rev_condition

score_pval_gim <- ModScore(module_score = "Spatial_genes1",
                           seu_ = ref_epi)

ref_epi@meta.data$rev_condition <- ifelse(ref_epi@meta.data$rev_condition == "normal",
                                          "Normal",
                                          ref_epi@meta.data$rev_condition)

ref_epi@meta.data$rev_condition <- factor(ref_epi@meta.data$rev_condition,
                                          levels= c("Normal",
                                                    "NAG",
                                                    "CAG",
                                                    "GIM",
                                                     "EGC"))

precancer_kruskal_plot <- ref_epi@meta.data %>%
  ggviolin(x = "rev_condition",
           y = "Spatial_genes1",
           fill = "rev_condition") +
  # stat_pvalue_manual(score_pval_gim,
  #                    y.position = 1.3,
  #                    step.increase = 0.1,
  #                    label = "p.adj.signif") +
  labs(subtitle = get_test_label(kruskal_test(Spatial_genes1 ~ rev_condition,
                                              data = ref_epi@meta.data)),
       caption = "All comparisons significant (Dunn's test FDR-adjusted p ≤ 0.001)") +
  ggtitle("Signature module score throughout Correa's cascade (epithelial cells subset)") +
  geom_boxplot(fill = "white",
               width = 0.3,
               alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 3) +
  scale_fill_brewer(type = "seq",
                    palette = "Reds") +
  xlab("") +
  ylab("Module score (26 genes)") +
  labs(fill = "Condition") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.5,
                                   hjust=1))

plot(precancer_kruskal_plot)

write_rds(precancer_kruskal_plot,
          file = file.path(opt$output_dir,"Figure_4F.rds"))
save(precancer_kruskal_plot,
     file = file.path(opt$output_dir,"Figure_4F.rda"))


# 07 - Disaggregate the stem-like cell group ------------------------------


## 7.1 - Import counts from epithelial object -----------------------------

seu <- readRDS("../D01_epithelial/epi_aggr.rds")
seu_counts <- seu@assays$RNA@counts


## 7.2 - Reference mapping ------------------------------------------------

reference <- readRDS("/mnt/ix1/Resources/scRNA_Ref/scRNA_datasets/Busslinger_upperGI_2021/busslinger_gastric_duodenum_epi.rds")
ref_counts <- GetAssayData(object = reference,
                           slot = "data")
# cell_annot <- SingleR(test = seu_counts,
#                       ref = ref_counts, 
#                       labels = reference@meta.data$Celltype,
#                       de.method = "wilcox",
#                       BPPARAM = MulticoreParam(32))

# saveRDS(cell_annot,
#         file.path(opt$output_dir,
#                   "cell_annot_busslinger.rds"),
#         version = 2)

cell_annot <- readRDS("temp/cell_annot_busslinger.rds")

table(Label = cell_annot$labels,
      Lost = is.na(cell_annot$pruned.labels))

if(identical(rownames(cell_annot), rownames(ref_epi@meta.data))) {
  cat("Cell tags from singleR match UMIs from the epithelial cell subset.\n")
  ref_epi@meta.data$specific_cell_type = cell_annot$labels
}

table(ref_epi@meta.data$final_celltype)

ref_epi@meta.data$specific_cell_type <- ifelse(ref_epi@meta.data$final_celltype!="Stem-like",
                                               ref_epi@meta.data$final_celltype,
                                               ref_epi@meta.data$specific_cell_type)


table(ref_epi@meta.data$specific_cell_type)
ref_epi@meta.data$specific_cell_type[ref_epi@meta.data$specific_cell_type=="gastric_Differentiating_cells"] <- "Gastric differentiating cells"
ref_epi@meta.data$specific_cell_type[ref_epi@meta.data$specific_cell_type=="duodenum_Stem_cells"] <- "Duodenum stem cells"
ref_epi@meta.data$specific_cell_type[ref_epi@meta.data$specific_cell_type=="duodenum_Differentiating_stem_cells"] <- "Duodenum differentiating stem cells"
ref_epi@meta.data$specific_cell_type[ref_epi@meta.data$specific_cell_type=="duodenum_Transit_amplyfying_cells"] <- "Duodenum transit amplifying cells"

dim(ref_epi)


## 7.3 - UMAP epithelial cells with granular labels -----------------------

col <- paletteer::paletteer_d("ggsci::category20_d3")[1:17]
names(col) = unique(sc_data$final_celltype)
col_subset <- col[names(col)%in%c("Enterocytes", "Goblet")]
col_subset <- c(brewer.pal(4, "Reds"), col_subset)
cells <- unique(ref_epi@meta.data$specific_cell_type)
names(col_subset) = c("Gastric differentiating cells",
                      "Duodenum transit amplifying cells",
                      "Duodenum differentiating stem cells",
                      "Duodenum stem cells",
                      "Enterocytes",
                      "Goblet")

umap_epi <- DimPlot(ref_epi, group.by = 'specific_cell_type') +
  ggtitle("Reference-mapped epithelial cells\n(45,679 epithelial cells mapped to a normal\nstomach and duodenum reference atlas)") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme_classic() +
  guides(x = axis,
         y = axis) +
  scale_color_manual(values = c(col, col_subset)) +
  theme(axis.line = element_line(arrow = arrow(angle = 15,
                                               length = unit(.15,
                                                             "inches"),
                                               type = "closed")),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right")
  
umap_epi


## 7.4 - Module score (epithelial cells) ----------------------------------

module_epi <- FeaturePlot(ref_epi,
            features =  "Spatial_genes1") +
  scale_color_gradient(low = "skyblue",
                       high = "firebrick") +
  ggtitle("Signature module score in epithelial cells") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme_classic() +
  guides(x = axis,
         y = axis) +
  labs(col = "Module score") +
  theme(axis.line = element_line(arrow = arrow(angle = 15,
                                               length = unit(.15,
                                                             "inches"),
                                               type = "closed")),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5))

umap_epi + module_epi


## 7.5 - Compare module score between epithelial cell types ---------------

Idents(ref_epi) <- "specific_cell_type"

score_pval_celltypes <- ModScore(module_score = "Spatial_genes1",
                                 seu_ = ref_epi)

VlnPlot(object = ref_epi,
        sort = "decreasing",
        features = "Spatial_genes1",
        assay = "RNA",
        group.by = "specific_cell_type",
        pt.size = 0) +
  geom_boxplot(width = 0.3,
               fill = "white",
               alpha = 0.8) +
  geom_hline(yintercept = 0,
             lty = 3) +
  ggtitle("Module score by cell type") +
  labs(subtitle = get_test_label(kruskal_test(Spatial_genes1 ~ specific_cell_type,
                                              data = ref_epi@meta.data)),
       caption = "All comparisons between gastric and intestinal lineages significant\n(Dunn's test FDR-adjusted p ≤ 0.001)") +
  ylab("Module score (26 genes)") +
  xlab("") +
  scale_fill_manual(values = c(col, col_subset)) +
  coord_flip() +
  theme(legend.position = "none")

write.xlsx(list(STable9 = score_pval_celltypes,
                STable10 = score_pval_gim),
           file.path(opt$output_dir, "STables9-10.xlsx"))

## 7.6- Gene signature heatmap --------------------------------------------

dotplot_genes_per_celltype <- DotPlot(object = ref_epi,
                                      features = c(gene_signature, "TFF3"))
prop_df <- dotplot_genes_per_celltype$data
head(prop_df)
exp_df <- prop_df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id,
              values_from = avg.exp.scaled) %>% 
  as.data.frame() 
max(exp_df[,2:ncol(exp_df)]);min(exp_df[,2:ncol(exp_df)])
# exp_df[,2:ncol(exp_df)] <- t(scale(t(exp_df[,2:ncol(exp_df)])))
colnames(exp_df)[1] <- "Gene Symbol"
rownames(exp_df) <- exp_df$`Gene Symbol`
exp_df$`Gene Symbol` <- NULL

cluster_anno <- colnames(exp_df)

column_ha <- HeatmapAnnotation(show_annotation_name = F,
                               `Cell type` = cluster_anno,
                               col = list(`Cell type` = c(col, col_subset)),
                               na_col = "grey")
column_split = ifelse(colnames(exp_df)%in%c("Enterocytes",
                                            "Goblet"),
                      "C",
                      ifelse(colnames(exp_df)%in%c("Duodenum stem cells",
                                                   "Duodenum differentiating stem cells",
                                                   "Duodenum transit amplifying cells"),
                             "B", "A"))

hp <- Heatmap(exp_df[-27,],
              name = "Z-score",
              col = col_fun,
              border = "black",
              show_column_names = F,
              show_row_names = T,
              top_annotation = column_ha,
              clustering_distance_rows = "manhattan",
              clustering_method_rows = "ward.D2",
              clustering_distance_columns = "manhattan",
              clustering_method_columns = "ward.D2",
              cluster_columns = F,
              column_order = c("Parietal",
                               "Chief",
                               "Pit",
                               "LYZ-positive",
                               "Gastric differentiating cells",
                               "REG3A positive",
                               "Neck",
                               "Endocrine",
                               "Isthmus",
                               "Duodenum stem cells",
                               "Duodenum differentiating stem cells",
                               "Duodenum transit amplifying cells",
                               "Enterocytes",
                               "Goblet"),
              column_split = column_split,
              column_title = c("Gastric lineages",
                               "Intestinal\nstem-like",
                               "Differentiated\nintestinal"),
              column_title_gp = gpar(fontsize = 12),
              row_split = 2,
              row_title = c("Expressed by\ngastric enterocytes",
                            "Expressed by\ngastric intestinal\nstem-like cells"),
              row_title_gp = gpar(fontsize = 12),
              row_names_gp = gpar(fontsize = 8),
              show_heatmap_legend = T)

exp_df2 <- exp_df[27,]

hp2 <- Heatmap(exp_df2,
               name = "Z-score",
               col = col_fun,
               border = "black",
               show_column_names = F,
               show_row_names = T,
               clustering_distance_rows = "manhattan",
               clustering_method_rows = "ward.D2",
               clustering_distance_columns = "manhattan",
               clustering_method_columns = "ward.D2",
               cluster_columns = F,
               column_order = c("Parietal",
                                "Chief",
                                "Pit",
                                "LYZ-positive",
                                "Gastric differentiating cells",
                                "REG3A positive",
                                "Neck",
                                "Endocrine",
                                "Isthmus",
                                "Duodenum stem cells",
                                "Duodenum differentiating stem cells",
                                "Duodenum transit amplifying cells",
                                "Enterocytes",
                                "Goblet"),
               column_split = column_split,
               column_title = c("Gastric lineages",
                                "Intestinal\nstem-like",
                               "Differentiated\nintestinal"),
               column_title_gp = gpar(fontsize = 12),
               row_title_gp = gpar(fontsize = 12),
               row_names_gp = gpar(fontsize = 8),
               show_heatmap_legend = FALSE)

# Horizontal concatenation of heatmaps with +
# Vertical concatenation of heatmaps with %v%

ht_list = hp %v% hp2

draw(ht_list,
     padding = unit(c(2,2,6,2),'mm'))
# To add more space before TFF3: ht_gap = unit(c(6, 2), "mm")


## 7.7 - Stacked barplot --------------------------------------------------

cell_by_stage <- data.frame(table(ref_epi$specific_cell_type,
                                  ref_epi$orig.ident))

stage <- data.frame(Var2 = ref_epi$orig.ident,
                    stage = ref_epi$rev_condition)
stage <- stage[!duplicated(stage$Var2),]
cell_by_stage <- merge(cell_by_stage,
                       stage,
                       by = "Var2",
                       al.x = T)

cell_by_stage$stage
order(cell_by_stage$stage)

cell_by_stage <- cell_by_stage[order(cell_by_stage$stage),]

cell_by_stage$Var1 <- as.character(cell_by_stage$Var1)
cell_by_stage$Var1 <- ifelse(cell_by_stage$Var1%in%c("Chief",
                                                     "Endocrine",
                                                     "Gastric differentiating cells",
                                                     "Isthmus",
                                                     "LYZ-positive",
                                                     "Neck",
                                                     "Parietal",
                                                     "Pit",
                                                     "REG3A positive"),
                             "Gastric lineages",
                             cell_by_stage$Var1)

col_subset2 <- col_subset
names(col_subset2)[1] <- "Gastric lineages"
col_subset2[1] <- "lightgoldenrod"
stage <- cell_by_stage$stage

cell_by_stage$Var2 <- factor(cell_by_stage$Var2,
                             levels = c("P6649_21408_normal",
                                        "P7015_21532_normal",
                                        "GSM3954946_NAG1",
                                        "GSM3954947_NAG2",
                                        "GSM3954948_NAG3",
                                        "GSM3954949_CAG1",
                                        "GSM3954950_CAG2",
                                        "GSM3954951_CAG3",
                                        "GSM3954952_IMW1",
                                        "GSM3954953_IMW2",
                                        "GSM3954954_IMS1",
                                        "GSM3954955_IMS2",
                                        "GSM3954956_IMS3",
                                        "GSM3954957_IMS4",
                                        "GSM4546347_Pat25",
                                        "GSM4546351_Pat29",
                                        "P7115_21567_control",
                                        "P6649_21411_GIM",
                                        "P7015_21531_GIM",
                                        "P7016_21533_GIM",
                                        "P7016_21534_GIM",
                                        "GSM3954958_EGC"),
                             ordered = T)
stage <- as.character(cell_by_stage$stage)
stage[cell_by_stage$Var2=="P7016_21534_GIM"] <- "GIM2"

barplot_lab <- unique(paste(str_split(as.character(cell_by_stage$Var2),
                                      "_",
                                      simplify = T)[,1],
                            stage))

cell_by_stage %>%
  mutate(Var1 = factor(Var1,
                       levels = c("Duodenum stem cells",
                                  "Duodenum differentiating stem cells",
                                  "Duodenum transit amplifying cells",
                                  "Enterocytes",
                                  "Goblet",
                                  "Gastric lineages"),
                       ordered = T)) %>%
  ggplot(aes(fill = Var1, y = Freq, x = Var2)) + 
  geom_bar(position="fill",
           stat="identity") +
  ggtitle("Proportion of aberrant epithelial cells\nthroughout Correa's cascade.") +
  scale_fill_manual(values = col_subset2) +
  ylab("Proportion") +
  xlab("") +
  scale_x_discrete(labels = barplot_lab) +
  labs(fill = "Cell type") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust=1))

rm(list = ls())
setwd("../")


