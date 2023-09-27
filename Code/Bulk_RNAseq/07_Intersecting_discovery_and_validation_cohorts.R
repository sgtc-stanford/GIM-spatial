# 00 - Load required R packages -------------------------------------------
library(optparse)
library(tidyverse)
library(edgeR)
library(openxlsx)
library(ggpubr)
library(VennDiagram)
library(ComplexHeatmap)
library(msigdbr)
library(clusterProfiler)
library(RColorBrewer)
library(circlize)
library(patchwork)

# 01 - Set up working directory -------------------------------------------

setwd("./C01_20230608_Intersection_discovery_validation")
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
                     c("-d", "--expression_discovery"),
                     default=NA,
                     type="character", 
                     help="Discovery set normalized expression")
parser <- add_option(parser,
                     c("-v", "--validation_genes"),
                     default=NA,
                     type="character", 
                     help="Validation set expression and results")
parser <- add_option(parser,
                     c("-i", "--intersection_limma_wgcna_discovery"),
                     default=NA,
                     type="character", 
                     help="Intersection between DEA and WGCNA (discovery set)")
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

opt$expression_discovery <- "../A04_20230601_Discovery_limma"
opt$validation_genes <- "../B03_20230607_Validation_limma"
opt$intersection_limma_wgcna_discovery <- "../A06_20230606_limma_WGCNA_intersection"
opt$output_dir <- "./temp"
opt$plot_dir <- "./plots"


# 03 - Intersecting training and testing cohorts --------------------------

pdf(file.path(opt$plot_dir,
              "Validation_plots.pdf"),
    height = 10.41,
    width = 10.41)

discovery_expression <- read_rds(file.path(opt$expression_discovery, "Discovery_norm_data.rds"))
validation_expression <- read_rds(file.path(opt$validation_genes, "Validation_norm_data.rds"))
intersection_genes <- read_rds(file.path(opt$intersection_limma_wgcna_discovery, "Heatmap_clusters.rds"))
load(file.path(opt$expression_discovery, "DE_results_discovery.rda"))
results_discovery <- results_treat
load(file.path(opt$validation_genes, "DE_results_validation.rda"))
results_validation <- results_treat
rm(results_eBayes, results_treat)

discovery_degs <- c(intersection_genes$cluster_1,
                    intersection_genes$cluster_2,
                    intersection_genes$cluster_3,
                    intersection_genes$cluster_4,
                    intersection_genes$cluster_5)

validation_degs <- intersect(subset(results_validation$high_vs_low_risk_Bod,
                                    subset = results_validation$high_vs_low_risk_Bod$adj.P.Val <= 0.05)$SYMBOL,
                             subset(results_validation$high_vs_low_risk_Ant,
                                    subset = results_validation$high_vs_low_risk_Ant$adj.P.Val <= 0.05)$SYMBOL)

validated_genes <- intersect(discovery_degs,
                             validation_degs)

validated_cluster_1 <- intersect(intersection_genes$cluster_1,
                                 validation_degs)
validated_cluster_2 <- intersect(intersection_genes$cluster_2,
                                 validation_degs)
validated_cluster_3 <- intersect(intersection_genes$cluster_3,
                                 validation_degs)
validated_cluster_4 <- intersect(intersection_genes$cluster_4,
                                 validation_degs)
validated_cluster_5 <- intersect(intersection_genes$cluster_5,
                                 validation_degs)
setdiff(intersection_genes$cluster_5, validated_cluster_5)

write.xlsx(list(cluster_1 = validated_cluster_1,
                cluster_2 = validated_cluster_2,
                cluster_3 = validated_cluster_3,
                cluster_4 = validated_cluster_4,
                cluster_5 = validated_cluster_5),
           file.path(opt$output_dir,
                     "validated_genes_by_cluster.xlsx"))


# 04 - Inspect and plot validation by cluster -----------------------------

length(intersection_genes$cluster_1); length(intersect(intersection_genes$cluster_1, validated_genes))
length(intersect(intersection_genes$cluster_1, validated_genes))/length(intersection_genes$cluster_1)
length(intersect(intersection_genes$cluster_1, validated_genes))/length(validated_genes)

length(intersection_genes$cluster_2); length(intersect(intersection_genes$cluster_2, validated_genes))
length(intersect(intersection_genes$cluster_2, validated_genes))/length(intersection_genes$cluster_2)
length(intersect(intersection_genes$cluster_2, validated_genes))/length(validated_genes)

length(intersection_genes$cluster_3); length(intersect(intersection_genes$cluster_3, validated_genes))
length(intersect(intersection_genes$cluster_3, validated_genes))/length(intersection_genes$cluster_3)
length(intersect(intersection_genes$cluster_3, validated_genes))/length(validated_genes)

length(intersection_genes$cluster_4); length(intersect(intersection_genes$cluster_4, validated_genes))
length(intersect(intersection_genes$cluster_4, validated_genes))/length(intersection_genes$cluster_4)
length(intersect(intersection_genes$cluster_4, validated_genes))/length(validated_genes)

length(intersection_genes$cluster_5); length(intersect(intersection_genes$cluster_5, validated_genes))
length(intersect(intersection_genes$cluster_5, validated_genes))/length(intersection_genes$cluster_5)
length(intersect(intersection_genes$cluster_5, validated_genes))/length(validated_genes)

barchart_df <- data.frame(cluster = c("cluster_1",
                                      "cluster_2",
                                      "cluster_3",
                                      "cluster_4",
                                      "cluster_5"),
                          n_validated = c(paste0("Cluster 1: ", length(intersect(intersection_genes$cluster_1, validated_genes)),
                                                 "/",
                                                 length(validated_genes),
                                                 " (",
                                                 round(length(intersect(intersection_genes$cluster_1, validated_genes))/length(validated_genes)*100,2),
                                                 "%)"),
                                          paste0("Cluster 2: ", length(intersect(intersection_genes$cluster_2, validated_genes)),
                                                 "/",
                                                 length(validated_genes),
                                                 " (",
                                                 round(length(intersect(intersection_genes$cluster_2, validated_genes))/length(validated_genes)*100,2),
                                                 "%)"),
                                          paste0("Cluster 3: ", length(intersect(intersection_genes$cluster_3, validated_genes)),
                                                 "/",
                                                 length(validated_genes),
                                                 " (",
                                                 round(length(intersect(intersection_genes$cluster_3, validated_genes))/length(validated_genes)*100,2),
                                                 "%)"),
                                          paste0("Cluster 4: ", length(intersect(intersection_genes$cluster_4, validated_genes)),
                                                 "/",
                                                 length(validated_genes),
                                                 " (",
                                                 round(length(intersect(intersection_genes$cluster_4, validated_genes))/length(validated_genes)*100,2),
                                                 "%)"),
                                          paste0("Cluster 5: ", length(intersect(intersection_genes$cluster_5, validated_genes)),
                                                 "/",
                                                 length(validated_genes),
                                                 " (",
                                                 round(length(intersect(intersection_genes$cluster_5, validated_genes))/length(validated_genes)*100,2),
                                                 "%)")),
                          percentage = as.numeric(c(round(length(intersect(intersection_genes$cluster_1, validated_genes))/length(validated_genes)*100,2),
                                                    round(length(intersect(intersection_genes$cluster_2, validated_genes))/length(validated_genes)*100,2),
                                                    round(length(intersect(intersection_genes$cluster_3, validated_genes))/length(validated_genes)*100,2),
                                                    round(length(intersect(intersection_genes$cluster_4, validated_genes))/length(validated_genes)*100,2),
                                                    round(length(intersect(intersection_genes$cluster_5, validated_genes))/length(validated_genes)*100,2))),
                          status = "Validated genes\nfrom each cluster")

fill.pal <- c(brewer.pal(9, "Blues")[c(8,6)], brewer.pal(9, "Spectral")[c(3,2,1)])

barchart_all <- barchart_df %>%
  ggplot(aes(x = status,
             y = percentage,
             fill = cluster,
             label = n_validated)) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5),
            col = "white",
            fontface = "bold") +
  theme_void() +
  ylab("Percentage") +
  xlab(NULL) +
  scale_fill_manual(values = fill.pal) +
  theme(legend.position = "none")

barchart_df_2 <- data.frame(cluster = rep(c("Cluster 1",
                                            "Cluster 2",
                                            "Cluster 3",
                                            "Cluster 4",
                                            "Cluster 5"), each = 2),
                            n_validated = c(paste0("Validated:\n",
                                                   sum(intersection_genes$cluster_1%in%validated_genes),
                                                   "/",
                                                   length(intersection_genes$cluster_1),
                                                   " (",
                                                   round(sum(intersection_genes$cluster_1%in%validated_genes)/length(intersection_genes$cluster_1)*100,2),
                                                   "%)"),
                                            paste0(#"Not validated:\n",
                                              round((1-sum(intersection_genes$cluster_1%in%validated_genes)/length(intersection_genes$cluster_1))*100,2), "%"),
                                            paste0("Validated:\n",
                                                   sum(intersection_genes$cluster_2%in%validated_genes),
                                                   "/",
                                                   length(intersection_genes$cluster_2),
                                                   " (",
                                                   round(sum(intersection_genes$cluster_2%in%validated_genes)/length(intersection_genes$cluster_2)*100,2),
                                                   "%)"),
                                            paste0(#"Not validated:\n",
                                              round((1-sum(intersection_genes$cluster_2%in%validated_genes)/length(intersection_genes$cluster_2))*100,2), "%"),
                                            paste0("Validated:\n",
                                                   sum(intersection_genes$cluster_3%in%validated_genes),
                                                   "/",
                                                   length(intersection_genes$cluster_3),
                                                   " (",
                                                   round(sum(intersection_genes$cluster_3%in%validated_genes)/length(intersection_genes$cluster_3)*100,2),
                                                   "%)"),
                                            paste0(#"Not validated:\n",
                                              round((1-sum(intersection_genes$cluster_3%in%validated_genes)/length(intersection_genes$cluster_3))*100,2), "%"),
                                            paste0("Validated:\n",
                                                   sum(intersection_genes$cluster_4%in%validated_genes),
                                                   "/",
                                                   length(intersection_genes$cluster_4),
                                                   " (",
                                                   round(sum(intersection_genes$cluster_4%in%validated_genes)/length(intersection_genes$cluster_4)*100,2),
                                                   "%)"),
                                            paste0(#"Not validated:\n",
                                              round((1-sum(intersection_genes$cluster_4%in%validated_genes)/length(intersection_genes$cluster_4))*100,2), "%"),
                                            paste0("Validated:\n",
                                                   sum(intersection_genes$cluster_5%in%validated_genes),
                                                   "/",
                                                   length(intersection_genes$cluster_5),
                                                   " (",
                                                   round(sum(intersection_genes$cluster_5%in%validated_genes)/length(intersection_genes$cluster_5)*100,2),
                                                   "%)"),
                                            paste0(#"Not validated:\n",
                                              round((1-sum(intersection_genes$cluster_5%in%validated_genes)/length(intersection_genes$cluster_5))*100,2), "%")),
                            percentage = c(round(sum(intersection_genes$cluster_1%in%validated_genes)/length(intersection_genes$cluster_1)*100,2),
                                           round((1-sum(intersection_genes$cluster_1%in%validated_genes)/length(intersection_genes$cluster_1))*100,2),
                                           round(sum(intersection_genes$cluster_2%in%validated_genes)/length(intersection_genes$cluster_2)*100,2),
                                           round((1-sum(intersection_genes$cluster_2%in%validated_genes)/length(intersection_genes$cluster_2))*100,2),
                                           round(sum(intersection_genes$cluster_3%in%validated_genes)/length(intersection_genes$cluster_3)*100,2),
                                           round((1-sum(intersection_genes$cluster_3%in%validated_genes)/length(intersection_genes$cluster_3))*100,2),
                                           round(sum(intersection_genes$cluster_4%in%validated_genes)/length(intersection_genes$cluster_4)*100,2),
                                           round((1-sum(intersection_genes$cluster_4%in%validated_genes)/length(intersection_genes$cluster_4))*100,2),
                                           round(sum(intersection_genes$cluster_5%in%validated_genes)/length(intersection_genes$cluster_5)*100,2),
                                           round((1-sum(intersection_genes$cluster_5%in%validated_genes)/length(intersection_genes$cluster_5))*100,2)),
                            status = "validated",
                            col = paste0(c("v",
                                           "nv"),
                                         rep(1:5, each = 2)))
fill.pal2 <- c(rep("grey", 5), fill.pal)

barchart_individual <- barchart_df_2 %>%
  ggplot(aes(x = status,
             y = percentage,
             fill = col,
             label = n_validated)) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5),
            col = "white",
            fontface = "bold") +
  theme_void() +
  facet_wrap(~cluster) +
  scale_fill_manual(values = fill.pal2) +
  theme(legend.position = "none")

print(ggarrange(barchart_all, barchart_individual, ncol = 2))


# 05 - Venn diagram discovery and validation sets -------------------------

myCol2 <- brewer.pal(3, "Paired")[1:3]

grid.newpage()

venn <- venn.diagram(
  x = list(discovery_degs,
           validation_degs,
           intersection_genes$cluster_5),
  category.names = c("Selected genes\ndiscovery set" ,
                     "Selected genes\nvalidation set" ,
                     "High-risk OLGIM signature"),
  filename = NULL,
  # Circles
  lwd = 1,
  lty = 'blank',
  fill = myCol2,
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-17, 17, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

grid.draw(venn)


# 06 - Heatmap discovery and validation sets ------------------------------

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

ha <- HeatmapAnnotation(`Biopsy site` = factor(ifelse(discovery_expression$samples$Location=="ANT", "Antrum", "Body")),
                        `OLGIM stage` = discovery_expression$samples$OLGIM,
                        col = list(`Biopsy site` = c("Antrum" = brewer.pal(8, "Greens")[7],
                                                     "Body" = brewer.pal(8, "Greens")[4]),
                                   `OLGIM stage` = c("OLGIM 0" = brewer.pal(11, "RdYlBu")[7],
                                                     "OLGIM I" = brewer.pal(11, "RdYlBu")[8],
                                                     "OLGIM II" = brewer.pal(11, "RdYlBu")[9],
                                                     "OLGIM III" = brewer.pal(11, "RdYlBu")[10],
                                                     "OLGIM IV" = brewer.pal(11, "RdYlBu")[11])))

side_ha_deg <- rowAnnotation(foo = anno_mark(at = which(validated_cluster_5%in%metaplasia_markers_complete),
                                             labels = rownames(discovery_expression[validated_cluster_5,])[rownames(discovery_expression[validated_cluster_5,])%in%metaplasia_markers_complete]),
                             annotation_name_gp = gpar(fontsize = 4))
dim(cpm(discovery_expression[validated_cluster_5,],
        log = T))

draw(Heatmap(t(scale(t(cpm(discovery_expression[validated_cluster_5,],
                      log = T)))),
        name = "Z-score",
        col = col_fun,
        column_split = 3,
        column_title = c("High-risk body\nand antrum",
                         "Low-risk body",
                         "Low-risk antrum"),
        row_title = "Validated genes from C-5 in the\ndiscovery set (100/105, 95.2%)",
        show_column_names = F,
        show_row_names = F,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = ha,
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = T,
        right_annotation = side_ha_deg))

ha <- HeatmapAnnotation(`Biopsy site` = factor(ifelse(validation_expression$samples$Location=="ANT", "Antrum", "Body")),
                        `OLGIM stage` = validation_expression$samples$OLGIM,
                        col = list(`Biopsy site` = c("Antrum" = brewer.pal(8, "Greens")[7],
                                                     "Body" = brewer.pal(8, "Greens")[4]),
                                   `OLGIM stage` = c("OLGIM 0" = brewer.pal(11, "RdYlBu")[7],
                                                     "OLGIM I" = brewer.pal(11, "RdYlBu")[8],
                                                     "OLGIM II" = brewer.pal(11, "RdYlBu")[9],
                                                     "OLGIM III" = brewer.pal(11, "RdYlBu")[10],
                                                     "OLGIM IV" = brewer.pal(11, "RdYlBu")[11])))

side_ha_deg <- rowAnnotation(foo = anno_mark(at = which(validated_cluster_5%in%metaplasia_markers_complete),
                                             labels = rownames(validation_expression[validated_cluster_5,])[rownames(validation_expression[validated_cluster_5,])%in%metaplasia_markers_complete]),
                             annotation_name_gp = gpar(fontsize = 4))

draw(Heatmap(t(scale(t(cpm(validation_expression[validated_cluster_5,],
                      log = T)))),
        name = "Z-score",
        col = col_fun,
        column_split = 2,
        column_title = c("Low-risk body\nand antrum",
                         "High-risk body\nand antrum"),
        row_title = "Validated genes from C-5 in the\nvalidation set",
        show_column_names = F,
        show_row_names = F,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = ha,
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = T,
        right_annotation = side_ha_deg))


# 07 - MSigDB gene sets ---------------------------------------------------


genelist <- validated_cluster_5


## 7.1 - MSigDB Hallmark gene sets ----------------------------------------

hallmark <- msigdbr(category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

hallmark_results <- subset(enricher(genelist,
                                    TERM2GENE = hallmark)@result,
                           qvalue <= 0.05)

# No enriched terms

## 7.2 - MSigDB positional gene sets (C1) ---------------------------------

c1 <- msigdbr(category = "C1") %>%
  dplyr::select(gs_name, gene_symbol)

c1_results <- subset(enricher(genelist,
                              TERM2GENE = c1)@result,
                     qvalue <= 0.05)

# Only enriched term: chr4q23

## 7.3 - MSigDB curated gene sets (C2) ------------------------------------

c2 <- msigdbr(category = "C2") %>%
  dplyr::select(gs_name, gene_symbol)

c2_results <- subset(enricher(genelist,
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

print(c2_plot)


## 7.4 - MSigDB regulatory gene sets (C3) ---------------------------------

c3 <- msigdbr(category = "C3") %>%
  dplyr::select(gs_name, gene_symbol)

c3_results <- subset(enricher(genelist,
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

print(c3_plot)


## 7.5 - MSigDB computational gene sets -----------------------------------

c4 <- msigdbr(category = "C4") %>%
  dplyr::select(gs_name, gene_symbol)

c4_results <- subset(enricher(genelist,
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

print(c4_plot)


## 6.6 - MSigDB GO (gene ontology) gene sets ------------------------------

c5 <- msigdbr(category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)

c5_results <- subset(enricher(genelist,
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

print(c5_plot)


## 7.7 - MSigDB oncogenic signature gene sets -----------------------------

c6 <- msigdbr(category = "C6") %>%
  dplyr::select(gs_name, gene_symbol)

c6_results <- subset(enricher(genelist,
                              TERM2GENE = c6)@result,
                     qvalue <= 0.05)

# No enriched terms


## 7.8 - MSigDB immunologic gene sets -------------------------------------

c7 <- msigdbr(category = "C7") %>%
  dplyr::select(gs_name, gene_symbol)

c7_results <- subset(enricher(genelist,
                              TERM2GENE = c7)@result,
                     qvalue <= 0.05)

# No enriched terms


## 7.9 - MSigDB cell type signature gene sets -----------------------------

c8 <- msigdbr(category = "C8") %>% 
  # subset(str_split(gs_name, "_", simplify = T)[,1]%in%c("BUSSLINGER", "GAO")) %>%
  dplyr::select(gs_name, gene_symbol)
# c8$gs_name[str_detect(c8$gs_name, "ENTEROCYTE")] <- "BUSSLINGER_ENTEROCYTE"

c8_results <- subset(enricher(genelist,
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

print(c8_plot)

print(c5_plot / c8_plot +
        plot_annotation('MSigDB gene sets') )

print(c5_plot + c8_plot +
        plot_annotation('MSigDB gene sets') )

ora_significant_list <- list(`GO gene sets` = c5_results,
                             `Cell type signature gene sets` = c8_results)

write.xlsx(ora_significant_list,
           file.path(opt$output_dir,
                     "ORA_results.xlsx"))

dev.off()

rm(list = ls())

setwd("../")
