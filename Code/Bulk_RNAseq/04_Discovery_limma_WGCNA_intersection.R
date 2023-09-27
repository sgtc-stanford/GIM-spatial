# 00 - Load required R packages -------------------------------------------

library(optparse)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(VennDiagram)
library(circlize)
library(openxlsx)
library(clusterProfiler)
library(msigdbr)
library(patchwork)


# 01 - Set up working directory -------------------------------------------
setwd("./A06_20230606_limma_WGCNA_intersection/")
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
                     c("-w", "--WGCNA_input_folder"),
                     default=NA,
                     type="character", 
                     help="WGCNA results input folder")
parser <- add_option(parser,
                     c("-l", "--limma_input_folder"),
                     default=NA,
                     type="character",
                     help="limma results input folder")
parser <- add_option(parser,
                     c("-o", "--output_dir"),
                     default=NA,
                     type="character",
                     help="Output directory for analysis results")
parser <- add_option(parser,
                     c("-p", "--plot_dir"),
                     default=NA,
                     type="character",
                     help="Output directory for plots")
opt <- parse_args(parser,
                  args = commandArgs(trailingOnly = TRUE))

opt$limma_input_folder <- "../A04_20230601_Discovery_limma"
opt$WGCNA_input_folder <- "../A05_20230605_Discovery_WGCNA"
opt$output_dir <- "./temp"
opt$plot_dir <- "./plots"


# 03 - Intersect DEA and WGCNA results ------------------------------------
pdf(file.path(opt$plot_dir,
              "DEA_WGCNA_plots.pdf"),
    height = 10.41,
    width = 10.41)

load(file.path(opt$limma_input_folder, "DE_results_discovery.rda"))
DEGs_body <- results_treat$high_vs_low_risk_Bod$SYMBOL[results_treat$high_vs_low_risk_Bod$adj.P.Val<=0.05]
DEGs_antrum <- results_treat$high_vs_low_risk_Ant$SYMBOL[results_treat$high_vs_low_risk_Ant$adj.P.Val<=0.05]
DEGs <- intersect(DEGs_body, DEGs_antrum); rm(DEGs_body, DEGs_antrum)

WGCNA_genes <- read_delim(file.path(opt$WGCNA_input_folder, "./temp/WGCNA_gene_modules.txt"))
WGCNA_genes <- WGCNA_genes[WGCNA_genes$colors%in%c("brown", "purple"),"gene_id"]$gene_id

keep_genes <- intersect(DEGs, WGCNA_genes)


# 04 - Venn diagram -------------------------------------------------------

myCol <- brewer.pal(8, "Paired")[c(2,5)]

venn <- venn.diagram(
  x = list(DEGs, WGCNA_genes),
  category.names = c("Differentially\nexpressed genes\nbody and antrum" ,
                     "WGCNA genes, brown\nand purple modules"),
  filename = NULL,
  # Circles
  lwd = 1,
  # lty = 'blank',
  fill = myCol,
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 25),
  cat.dist = c(0.045, 0.04))

grid.newpage()
grid.draw(venn)


# 05 - Heatmaps -----------------------------------------------------------

metaplasia_markers_complete <- c("LEFTY1", "KLF5", "SALL4", "DMBT1", "LGR5", "TP73", "POU5F1", "SOX2",
                                 "RUNX1", "AXIN2", "VIL1", "BHLHA15", "TNFRSF19", "CCKBR",
                                 "LRIG1", "GAST", "TFF2", "TFF3", "AQP5", "CDX1", "CDX2",
                                 "CFTR", "OLFM4", "IFNG", "LGR4", "LGR6", "SOX9", "CTNNB1",
                                 "MUC6", "KRT7", "KRT17", "PSMA7", "ZFAS1", "ACE2", "EPHB2",
                                 "CD44", "SOX4", "HES1", "EPCAM", "KRT18", "MUC1", "PGA4",
                                 "PGA3", "CHGA", "CHGB", "MKI67", "BRIC", "CEACAM5", "CEACAM6",
                                 "REG2A", "LCN2", "COX7B", "UQCRB", "GKN1", "GKN2", "MUC5A", "PGC",
                                 "MUC2", "ITLN1", "HES6", "CDH17", "COL3A1", "PDGFRB", "REG1A",
                                 "CLDN3", "CKN2A", "RBP4", "FABP1", "TFF1", "SPINK4", "MUC13")

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

expression_data <- read_rds(file.path(opt$limma_input_folder, "Discovery_norm_data.rds"))

side_deg = rowAnnotation(foo = anno_mark(at = which(rownames(expression_data[keep_genes,])%in%metaplasia_markers_complete),
                                         labels = rownames(expression_data[keep_genes,])[rownames(expression_data[keep_genes,])%in%metaplasia_markers_complete]),
                         annotation_name_gp = gpar(fontsize = 4))

ha <- HeatmapAnnotation(`Biopsy site` = factor(ifelse(expression_data$samples$Location=="ANT", "Antrum", "Body")),
                        `OLGIM stage` = expression_data$samples$OLGIM,
                        Lesion = factor(ifelse(expression_data$samples$Severity=="Normal", "Normal",
                                               ifelse(expression_data$samples$Severity=="Mild", "Mild GIM",
                                                      ifelse(expression_data$samples$Severity=="Moderate", "Moderate GIM", "Marked GIM"))),
                                        levels = c("Normal", "Mild GIM", "Moderate GIM", "Marked GIM"),
                                        ordered = T),
                        col = list(`Biopsy site` = c("Antrum" = brewer.pal(8, "Greens")[7],
                                                     "Body" = brewer.pal(8, "Greens")[4]),
                                   `OLGIM stage` = c("OLGIM 0" = brewer.pal(11, "RdYlBu")[7],
                                                     "OLGIM I" = brewer.pal(11, "RdYlBu")[8],
                                                     "OLGIM II" = brewer.pal(11, "RdYlBu")[9],
                                                     "OLGIM III" = brewer.pal(11, "RdYlBu")[10],
                                                     "OLGIM IV" = brewer.pal(11, "RdYlBu")[11]),
                                   Lesion = c("Normal" = brewer.pal(9, "Oranges")[3],
                                              "Mild GIM" = brewer.pal(9, "Oranges")[4],
                                              "Moderate GIM" = brewer.pal(9, "Oranges")[5],
                                              "Marked GIM" = brewer.pal(9, "Oranges")[6])))

draw(Heatmap(t(scale(t(cpm(expression_data$counts[keep_genes,], log = T)))),
             name = "Z-score",
             col = col_fun,
             # row_title = NULL,
             row_split = 5,
             row_title = c("C-1",
                           "C-2",
                           "C-3",
                           "C-4",
                           "C-5"),
             column_split = 3,
             column_title = c("High-risk\nOLGIM",
                              "Low-risk\nbody",
                              "Low-risk\nantrum"),
             show_row_names = F,
             show_column_names = F,
             top_annotation = ha,
             right_annotation = side_deg,
             clustering_distance_columns = "pearson",
             clustering_distance_rows = "pearson",
             clustering_method_columns = "ward.D2",
             clustering_method_rows = "ward.D2",
             show_row_dend = T,
             column_gap = unit(1, "mm"),
             row_gap = unit(0, "mm"),
             border = T))


## 5.1 Extract gene names from each cluster -------------------------------

heat_meta_markers <- Heatmap(t(scale(t(cpm(expression_data$counts[keep_genes,], log = T)))),
                             row_split = 5,
                             clustering_distance_columns = "pearson",
                             clustering_distance_rows = "pearson",
                             clustering_method_columns = "ward.D2",
                             clustering_method_rows = "ward.D2",
                             show_row_dend = F)
heat_genes <- rownames(expression_data$counts[keep_genes,])
clusterlist <- row_order(heat_meta_markers)
cluster_1 <- heat_genes[clusterlist[[1]]]
cluster_2 <- heat_genes[clusterlist[[2]]]
cluster_3 <- heat_genes[clusterlist[[3]]]
cluster_4 <- heat_genes[clusterlist[[4]]]
cluster_5 <- heat_genes[clusterlist[[5]]]

## 5.2 - Plot cluster_5 alone ---------------------------------------------

Heatmap(t(scale(t(cpm(expression_data[cluster_5,], log = T)))),
        name = "Z-score",
        col = col_fun,
        column_split = 2,
        column_title = c("Low-risk body\nand antrum",
                         "High-risk body\nand antrum"),
        row_title = "High-risk OLGIM signature genes",
        show_column_names = F,
        show_row_names = F,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = ha,
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = T)
clusters <- list(cluster_1 = cluster_1,
                 cluster_2 = cluster_2,
                 cluster_3 = cluster_3,
                 cluster_4 = cluster_4,
                 cluster_5 = cluster_5)
write.xlsx(clusters, file.path(opt$output_dir, "Cluster_genes.xlsx"))
save(clusters, file = file.path(opt$output_dir, "Cluster_genes.RData"))


# 06 - Functional enrichment ----------------------------------------------

genelist <- keep_genes


## 6.1 - MSigDB Hallmark gene sets ----------------------------------------

hallmark <- msigdbr(category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

hallmark_results <- subset(enricher(genelist,
                                    universe = intersect(rownames(expression_data),
                                                         hallmark$gene_symbol),
                                    TERM2GENE = hallmark)@result,
                           qvalue <= 0.05)

# Xenobiotic, pancreatic, coagulation, bile acid, kras

hallmark_plot <- hallmark_results %>%
  mutate(logP = -log10(qvalue)) %>%
  mutate(Significant = qvalue<=0.05) %>%
  mutate(ID = gsub("_", " ", ID)) %>%
  head(5) %>%
  ggplot(aes(x = reorder(ID,(logP)), y = logP, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none") +
  scale_fill_gradient(low = "salmon", high = "red4") +
  ggtitle("Hallmark gene sets") +
  xlab("") +
  ylab("-log10 Q-value") +
  theme(axis.text.y = element_text(size = 12))

hallmark_plot


## 6.2 - MSigDB positional gene sets (C1) ---------------------------------

c1 <- msigdbr(category = "C1") %>%
  dplyr::select(gs_name, gene_symbol)

c1_results <- subset(enricher(genelist,
                              universe = intersect(rownames(expression_data),
                                                   c1$gene_symbol),
                              TERM2GENE = c1)@result,
                     qvalue <= 0.05)

# No enriched terms

## 6.3 - MSigDB curated gene sets (C2) ------------------------------------

c2 <- msigdbr(category = "C2") %>%
  dplyr::select(gs_name, gene_symbol)

c2_results <- subset(enricher(genelist,
                              universe = intersect(rownames(expression_data),
                                                   c2$gene_symbol),
                              TERM2GENE = c2)@result,
                     qvalue <= 0.05)

c2_plot <- c2_results %>%
  mutate(logP = -log10(qvalue)) %>%
  mutate(Significant = qvalue<=0.05) %>%
  mutate(ID = gsub("_", " ", ID)) %>%
  head(5) %>%
  ggplot(aes(x = reorder(ID,(logP)), y = logP, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none") +
  scale_fill_gradient(low = "salmon", high = "red4") +
  ggtitle("Curated gene sets") +
  xlab("") +
  ylab("-log10 Q-value") +
  theme(axis.text.y = element_text(size = 12))

c2_plot


## 6.4 - MSigDB regulatory gene sets (C3) ---------------------------------

c3 <- msigdbr(category = "C3") %>%
  dplyr::select(gs_name, gene_symbol)

c3_results <- subset(enricher(genelist,
                              universe = intersect(rownames(expression_data),
                                                   c3$gene_symbol),
                              TERM2GENE = c3)@result,
                     qvalue <= 0.05)

# Enrichment of HNF1 and HNF4 targets

c3_plot <- c3_results %>%
  mutate(logP = -log10(qvalue)) %>%
  mutate(Significant = qvalue<=0.05) %>%
  mutate(ID = gsub("_", " ", ID)) %>%
  head(5) %>%
  ggplot(aes(x = reorder(ID,(logP)), y = logP, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none") +
  scale_fill_gradient(low = "salmon", high = "red4") +
  ggtitle("Regulatory gene sets") +
  xlab("") +
  ylab("-log10 Q-value") +
  theme(axis.text.y = element_text(size = 12))

c3_plot


## 6.5 - MSigDB computational gene sets -----------------------------------

c4 <- msigdbr(category = "C4") %>%
  dplyr::select(gs_name, gene_symbol)

c4_results <- subset(enricher(genelist,
                              universe = intersect(rownames(expression_data),
                                                   c4$gene_symbol),
                              TERM2GENE = c4)@result,
                     qvalue <= 0.05)

# Enrichment of 9 computational gene sets

c4_plot <- c4_results %>%
  mutate(logP = -log10(qvalue)) %>%
  mutate(Significant = qvalue<=0.05) %>%
  mutate(ID = gsub("_", " ", ID)) %>%
  head(5) %>%
  ggplot(aes(x = reorder(ID,(logP)), y = logP, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none") +
  scale_fill_gradient(low = "salmon", high = "red4") +
  ggtitle("Computational gene sets") +
  xlab("") +
  ylab("-log10 Q-value") +
  theme(axis.text.y = element_text(size = 12))

c4_plot


## 6.6 - MSigDB GO (gene ontology) gene sets ------------------------------

c5 <- msigdbr(category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)

c5_results <- subset(enricher(genelist,
                              universe = intersect(rownames(expression_data),
                                                   c5$gene_symbol),
                              TERM2GENE = c5)@result,
                     qvalue <= 0.05)

# Enrichment of terms related to intestinal absorption

c5_plot <- c5_results %>%
  mutate(logP = -log10(qvalue)) %>%
  mutate(Significant = qvalue<=0.05) %>%
  mutate(ID = gsub("_", " ", ID)) %>%
  head(5) %>%
  ggplot(aes(x = reorder(ID,(logP)), y = logP, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none") +
  scale_fill_gradient(low = "salmon", high = "red4") +
  ggtitle("GO gene sets") +
  xlab("") +
  ylab("-log10 Q-value") +
  theme(axis.text.y = element_text(size = 12))

c5_plot


## 6.7 - MSigDB oncogenic signature gene sets -----------------------------

c6 <- msigdbr(category = "C6") %>%
  dplyr::select(gs_name, gene_symbol)

c6_results <- subset(enricher(genelist,
                              universe = intersect(rownames(expression_data),
                                                   c6$gene_symbol),
                              TERM2GENE = c6)@result,
                     qvalue <= 0.05)

# No enriched terms


## 6.8 - MSigDB immunologic gene sets -------------------------------------

c7 <- msigdbr(category = "C7") %>%
  dplyr::select(gs_name, gene_symbol)

c7_results <- subset(enricher(genelist,
                              universe = intersect(rownames(expression_data),
                                                   c7$gene_symbol),
                              TERM2GENE = c7)@result,
                     qvalue <= 0.05)

# EGSE19888_ADENOSINE_A3R_INH_VS_INH_PRETREAT_AND_ACT_WITH_TCELL_MEMBRANES_MAST_CELL_DN
# Baram, 2010, J Immunol., PMID: 20190146

c7_plot <- c7_results %>%
  mutate(logP = -log10(qvalue)) %>%
  mutate(Significant = qvalue<=0.05) %>%
  mutate(ID = gsub("_", " ", ID)) %>%
  head(5) %>%
  ggplot(aes(x = reorder(ID,(logP)), y = logP, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none") +
  scale_fill_gradient(low = "salmon", high = "red4") +
  ggtitle("Immunologic gene sets") +
  xlab("") +
  ylab("-log10 Q-value") +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12))

c7_plot


## 6.9 - MSigDB cell type signature gene sets -----------------------------

c8 <- msigdbr(category = "C8") %>% 
  # subset(str_split(gs_name, "_", simplify = T)[,1]%in%c("BUSSLINGER", "GAO")) %>%
  dplyr::select(gs_name, gene_symbol)
# c8$gs_name[str_detect(c8$gs_name, "ENTEROCYTE")] <- "BUSSLINGER_ENTEROCYTE"

c8_results <- subset(enricher(genelist,
                              universe = intersect(rownames(expression_data),
                                                   c8$gene_symbol),
                              TERM2GENE = c8)@result,
                     qvalue <= 0.05)

c8_plot <- c8_results %>%
  mutate(logP = -log10(qvalue)) %>%
  mutate(Significant = qvalue<=0.05) %>%
  mutate(ID = gsub("_", " ", ID)) %>%
  head(5) %>%
  ggplot(aes(x = reorder(ID,(logP)), y = logP, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none") +
  scale_fill_gradient(low = "salmon", high = "red4") +
  ggtitle("Cell type signature gene sets") +
  xlab("") +
  ylab("-log10 Q-value") +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12))

c8_plot

c3_plot / c5_plot / c8_plot +
  plot_annotation('MSigDB gene sets') 

write_rds(list(cluster_1 = cluster_1,
               cluster_2 = cluster_2,
               cluster_3 = cluster_3,
               cluster_4 = cluster_4,
               cluster_5 = cluster_5),
          file = "Heatmap_clusters.rds")

dev.off()

rm(list = ls())

setwd("../")
