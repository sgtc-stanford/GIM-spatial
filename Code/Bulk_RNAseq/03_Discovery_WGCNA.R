# 00 - Load required R packages -------------------------------------------
library(optparse)
library(tidyverse)
library(edgeR)
library(DESeq2)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 01 - Set up working directory -------------------------------------------

setwd("./A05_20230605_Discovery_WGCNA/")
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
                     c("-i", "--input_DGElist"),
                     default=NA,
                     type="character", 
                     help="DGElist used for WGCNA")
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

opt$input_DGElist <- "../A03_20233105_Discovery_DGEList/"
opt$output_dir <- "./temp"
opt$plot_dir <- "./plots"


# 03 - Import and clean up data for WGCNA ---------------------------------

pdf(file.path(opt$plot_dir,
              "WGCNA_plots.pdf"),
    height = 10.41,
    width = 10.41)

clean_data <- read_rds(file.path(opt$input_DGElist,
                                 "Clean_data_DGEList.rds"))

# Format the "group" column in clean_data$samples
clean_data$samples$group
clean_data$samples$lib.size.1 <- NULL
clean_data$samples$norm.factors.1 <- NULL
dim(clean_data)
# [1] 38643    88

wgcna_exp <- data.frame(rownames(clean_data), clean_data$counts)
wgcna_exp <- as_tibble(wgcna_exp)
names(wgcna_exp)[1] = "SYMBOL"
head(wgcna_exp)

col_sel <- names(wgcna_exp)[-1]     # Get all but first column name

mdata <- wgcna_exp %>%
  tidyr::pivot_longer(
    .,                        # The dot is the the input data
    col = all_of(col_sel)
  )

mdata2 <- data.frame(SampleID = colnames(clean_data),
                     group = clean_data$samples$OLGIM,
                     risk = ifelse(clean_data$samples$OLGIM %in% c("OLGIM IV", "OLGIM III"),
                                   "High-risk",
                                   "Low-risk"))

mdata <- merge(mdata, mdata2, by.x = "name", by.y = "SampleID"); rm(mdata2)

mdata$group = factor(mdata$group,
                     levels = c("OLGIM 0",
                                "OLGIM I",
                                "OLGIM II",
                                "OLGIM III",
                                "OLGIM IV"),
                     ordered = T)

mdata$risk = factor(mdata$risk,
                    levels = c("Low-risk", "High-risk"),
                    ordered = T)

de_input <- clean_data$counts

meta_df <- data.frame(Sample = colnames(de_input),
                      Type = factor(ifelse(clean_data$samples$OLGIM %in% c("OLGIM IV", "OLGIM III"),
                                           "HighRisk",
                                           "LowRisk"),
                                    levels = c("LowRisk", "HighRisk")))

## 3.1 - Preprocess data with DESeq2 --------------------------------------

dds <- DESeqDataSetFromMatrix(round(de_input),
                              meta_df,
                              design = ~Type)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds)
wpn_vsd <- getVarianceStabilizedData(dds)


## 3.2 - Keep 15% most variable genes -------------------------------------

rv_wpn <- rowVars(wpn_vsd)

round(quantile(rv_wpn,
         probs = seq(.1, .9, by = .1)), 3)
#   10%   20%   30%   40%   50%   60%   70%   80%   90% 
# 0.000 0.005 0.032 0.070 0.118 0.182 0.265 0.374 0.639 

q85_wpn <- quantile(rowVars(wpn_vsd), .85)
sum(rv_wpn > q85_wpn)
# [1] 5797

expr_normalized_q85 <- wpn_vsd[rv_wpn > q85_wpn,]
dim(expr_normalized_q85)
# [1] 4809   88
# [1] 5797   88

expr_normalized_q85_df <- data.frame(expr_normalized_q85) %>%
  mutate(
    Gene_id = row.names(expr_normalized_q85)
  ) %>%
  pivot_longer(-Gene_id)


## 3.3 - Soft thresholding ------------------------------------------------

input_mat_85 = t(expr_normalized_q85)

allowWGCNAThreads()

powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

sft_85 = pickSoftThreshold(
  input_mat_85,             # <= Input data
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2))
cex1 = 0.9

plot(sft_85$fitIndices[,1],
     -sign(sft_85$fitIndices[,3]) * sft_85$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste("Scale independence")
)
abline(h = 0.8, col = "black", lty = 2)
text(sft_85$fitIndices[,1],
     -sign(sft_85$fitIndices[,3]) * sft_85$fitIndices[,2],
     labels = powers, cex = cex1, col = "red"
)
sfTop <- -sign(sft_85$fitIndices[,3]) * sft_85$fitIndices[,2]
names(sfTop) <- sft_85$fitIndices[,1]
sfTop
#          1          2          3          4          5          6          7          8          9 
# -0.5318818  0.5563569  0.8209524  0.8132676  0.7932252  0.7984684  0.7940601  0.8121199  0.8048497 
#         10         12         14         16         18         20 
# 0.8057965  0.7965975  0.7945893  0.8108961  0.8115765  0.8236576 
# Power 16 will probably work, at a scale-free topology R^2 >0.8

plot(sft_85$fitIndices[,1],
     sft_85$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
abline(h = 100, col = "black", lty = 2)
text(sft_85$fitIndices[, 1],
     sft_85$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")
meanCon <- sft_85$fitIndices$mean.k.
names(meanCon) <- sft_85$fitIndices[,1]
meanCon
#           1           2           3           4           5           6           7           8 
# 1768.768777  832.846884  474.665756  300.834611  203.966641  144.926492  106.617593   80.583595 
#           9          10          12          14          16          18          20 
#   62.248194   48.961026   31.645801   21.426681   15.048609   10.889445    8.078973
# Power of 16 has a mean connectivity of 15

picked_power = 16


# 04 - WGCNA --------------------------------------------------------------

temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat_85,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor     # Return cor function to original namespace

## 4.1 - WGCNA Plots ------------------------------------------------------

### 4.1.1 - Cluster dendrogram --------------------------------------------

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

dissTOM <- 1-TOMsimilarityFromExpr(input_mat_85, power = 16)

### 4.1.2 - MDS -----------------------------------------------------------

cmd1 <- cmdscale(as.dist(dissTOM),2)
# sizeGrWindow(7,6)
par(mfrow=c(1,1))
plot(cmd1,
     col = mergedColors,
     main = "MDSplot",
     xlab = "ScalingDimension1",
     ylab = "ScalingDimension2")

### 4.1.3 - Downsize gene set to 1500 for plotting ------------------------

restGenes <- (mergedColors != "grey")

# Create vetor to sample 1500 genes (just to downsize the resulting plots)
set.seed(10)

diss1 <- 1-TOMsimilarityFromExpr(input_mat_85[,restGenes], power = 16 )
select = sample(nrow(diss1),
                size = 1500)
# Select colors not from grey module
my_cols <- mergedColors[restGenes]
# Create color palette for heatmaps
myheatcol <- gplots::colorpanel(250,'red',"orange",'lemonchiffon')

### 4.1.4 - TOM plot - similarity matrix ----------------------------------

plotTOM <- diss1[select, select]
hier1 <- hclust(as.dist(plotTOM), method="average" )
diag(diss1) = NA
# sizeGrWindow(9,9)
TOMplot(plotTOM^4,
        hier1,
        as.character(my_cols[select]),
        main = "TOM heatmap plot, module genes",
        col = myheatcol)

### 4.1.5 - Adjacency plot ------------------------------------------------

diss2 <- 1-adjacency(input_mat_85[,restGenes], power = 16)
plotTOM2 <- diss2[select, select]
hier2 <- hclust(as.dist(plotTOM2), method = "average")
diag(diss2) = NA
# sizeGrWindow(9,9)
TOMplot(plotTOM2^4,
        hier2,
        my_cols[select],
        main = "Adjacency heatmap plot, module genes",
        col = myheatcol)

## 4.2 - Write out the results --------------------------------------------

head(module_df[1:5,])
table(module_df$colors)

write_delim(module_df,
            file = file.path(opt$output_dir, "WGCNA_gene_modules.txt"),
            delim = "\t")

## 4.2 Module-trait relationships -----------------------------------------

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat_85, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order <- names(MEs0) %>% gsub("ME","", .)

# Add treatment names using goruped OLGIM
if(identical(rownames(MEs0), colnames(clean_data))) {
  MEs0$treatment <- clean_data$samples$OLGIM  
  # MEs0$loc_sev <- paste(clean_data$samples$Location,
  #                       clean_data$samples$Severity,
  #                       sep = "-")
}

# tidy & plot data
MEs0$treatment <- ifelse(MEs0$treatment%in%c("OLGIM III", "OLGIM IV"), "High-risk", "Low-risk")
MEs0$treatment <- factor(MEs0$treatment, levels = c("Low-risk", "High-risk"))
mME <- MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

subset(mME, subset = name != "grey") %>% ggplot(., aes(x = treatment,
                                                       y = name,
                                                       fill = scale(value))) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(min(scale(mME$value)), max(scale(mME$value)))) +
  # theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships",
       x = "Group",
       y = "Modules",
       fill = "Scaled\nEigen gene")

# brown, greenyellow, magenta, red, purple, tan and pink modules look interesting


# 05 - Inspect hierarchical clustering with gene modules ------------------


# Heatmaps exported as pdf. Dimensions (height X width): 10.41 X 6.25 

with(module_df, table(colors[colors!="grey"]))
# black        blue       brown       green greenyellow     magenta        pink      purple 
#   136        1219         826         460          49          87         127          69 
# red      salmon         tan   turquoise      yellow 
# 207          40          48        1621         505 

## 5.1 - brown module -----------------------------------------------------

# N = 826
modules_of_interest = "brown"

submod <- module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) <- module_df$gene_id

subexpr <- expr_normalized_q85[submod$gene_id,]

norm_data <- read_rds("../A04_20230601_Discovery_limma/Discovery_norm_data.rds")
keep_wgcna <- rownames(norm_data)%in%submod$gene_id

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

side_ha = rowAnnotation(foo = anno_mark(at = which(rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete),
                                        labels = rownames(norm_data[keep_wgcna,])[rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete]),
                        annotation_name_gp = gpar(fontsize = 4))

ha <- HeatmapAnnotation(`Biopsy site` = factor(ifelse(norm_data$samples$Location=="ANT", "Antrum", "Body")),
                        `OLGIM stage` = norm_data$samples$OLGIM,
                        Lesion = factor(ifelse(norm_data$samples$Severity=="Normal", "Normal",
                                               ifelse(norm_data$samples$Severity=="Mild", "Mild GIM",
                                                      ifelse(norm_data$samples$Severity=="Moderate", "Moderate GIM", "Marked GIM"))),
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

Heatmap(t(scale(t(cpm(norm_data[keep_wgcna,], log = T)))),
        name = "Z-score",
        col = col_fun,
        column_split = 3,
        column_title = c("Low-risk\nbody",
                         "Low-risk\nantrum",
                         "High-risk\nOLGIM"),
        row_title = "Genes from brown module",
        show_row_names = F,
        show_column_names = F,
        top_annotation = ha,
        right_annotation = side_ha,
        row_names_gp = gpar(fontsize = 2),
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = F)

# brown module captures most of the metaplasia markers and separates high- from low-risk OLGIM stages

## 5.2 - greenyellow module -----------------------------------------------

# N = 49

modules_of_interest = "greenyellow"

submod <- module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) <- module_df$gene_id

subexpr <- expr_normalized_q85[submod$gene_id,]

keep_wgcna <- rownames(norm_data)%in%submod$gene_id

side_ha = rowAnnotation(foo = anno_mark(at = which(rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete),
                                        labels = rownames(norm_data[keep_wgcna,])[rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete]),
                        annotation_name_gp = gpar(fontsize = 4))

Heatmap(t(scale(t(cpm(norm_data[keep_wgcna,], log = T)))),
        name = "Z-score",
        col = col_fun,
        row_title = "Genes from greenyellow module",
        show_row_names = F,
        show_column_names = F,
        top_annotation = ha,
        right_annotation = side_ha,
        row_names_gp = gpar(fontsize = 2),
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = F)
# greenyellow module is not informative on GIM

## 5.3 - magenta module ---------------------------------------------------

# N = 87

modules_of_interest = "magenta"

submod <- module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) <- module_df$gene_id

subexpr <- expr_normalized_q85[submod$gene_id,]

keep_wgcna <- rownames(norm_data)%in%submod$gene_id

side_ha = rowAnnotation(foo = anno_mark(at = which(rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete),
                                        labels = rownames(norm_data[keep_wgcna,])[rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete]),
                        annotation_name_gp = gpar(fontsize = 4))

Heatmap(t(scale(t(cpm(norm_data[keep_wgcna,], log = T)))),
        name = "Z-score",
        col = col_fun,
        row_title = "Genes from magenta module",
        show_row_names = F,
        show_column_names = F,
        top_annotation = ha,
        right_annotation = side_ha,
        row_names_gp = gpar(fontsize = 2),
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = F)
# magenta module is not informative on GIM

## 5.4 - red module -----------------------------------------------------

# N = 207

modules_of_interest = "red"

submod <- module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) <- module_df$gene_id

subexpr <- expr_normalized_q85[submod$gene_id,]

keep_wgcna <- rownames(norm_data)%in%submod$gene_id

side_ha = rowAnnotation(foo = anno_mark(at = which(rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete),
                                        labels = rownames(norm_data[keep_wgcna,])[rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete]),
                        annotation_name_gp = gpar(fontsize = 4))

Heatmap(t(scale(t(cpm(norm_data[keep_wgcna,], log = T)))),
        name = "Z-score",
        col = col_fun,
        row_title = "Genes from red module",
        show_row_names = F,
        show_column_names = F,
        top_annotation = ha,
        right_annotation = side_ha,
        row_names_gp = gpar(fontsize = 2),
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = F)
# red module is mostly antrum vs corpus

## 5.5 - purple module ----------------------------------------------------

# N = 69

modules_of_interest = "purple"

submod <- module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) <- module_df$gene_id

subexpr <- expr_normalized_q85[submod$gene_id,]

keep_wgcna <- rownames(norm_data)%in%submod$gene_id

side_ha = rowAnnotation(foo = anno_mark(at = which(rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete),
                                        labels = rownames(norm_data[keep_wgcna,])[rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete]),
                        annotation_name_gp = gpar(fontsize = 4))

Heatmap(t(scale(t(cpm(norm_data[keep_wgcna,], log = T)))),
        name = "Z-score",
        col = col_fun,
        column_split = 3,
        column_title = c("Low-risk\nantrum",
                         "High-risk\nOLGIM",
                         "Low-risk\nbody"),
        row_title = "Genes from purple module",
        show_row_names = F,
        show_column_names = F,
        top_annotation = ha,
        right_annotation = side_ha,
        row_names_gp = gpar(fontsize = 2),
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = F)
# purple module represents genes that are down-regulated in advanced OLGIM

## 5.6 - tan module -------------------------------------------------------

# N = 48

modules_of_interest = "tan"

submod <- module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) <- module_df$gene_id

subexpr <- expr_normalized_q85[submod$gene_id,]

keep_wgcna <- rownames(norm_data)%in%submod$gene_id

side_ha = rowAnnotation(foo = anno_mark(at = which(rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete),
                                        labels = rownames(norm_data[keep_wgcna,])[rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete]),
                        annotation_name_gp = gpar(fontsize = 4))

Heatmap(t(scale(t(cpm(norm_data[keep_wgcna,], log = T)))),
        name = "Z-score",
        col = col_fun,
        row_title = "Genes from tan module",
        show_row_names = F,
        show_column_names = F,
        top_annotation = ha,
        right_annotation = side_ha,
        row_names_gp = gpar(fontsize = 2),
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = F)
# tan module is mostly noisy

## 5.7 - pink module ------------------------------------------------------

# N = 127

modules_of_interest = "pink"

submod <- module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) <- module_df$gene_id

subexpr <- expr_normalized_q85[submod$gene_id,]

keep_wgcna <- rownames(norm_data)%in%submod$gene_id

side_ha = rowAnnotation(foo = anno_mark(at = which(rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete),
                                        labels = rownames(norm_data[keep_wgcna,])[rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete]),
                        annotation_name_gp = gpar(fontsize = 4))

Heatmap(t(scale(t(cpm(norm_data[keep_wgcna,], log = T)))),
        name = "Z-score",
        col = col_fun,
        row_title = "Genes from pink module",
        show_row_names = F,
        show_column_names = F,
        top_annotation = ha,
        right_annotation = side_ha,
        row_names_gp = gpar(fontsize = 2),
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = F)
# pink module is mostly noise

## 5.8 - Plot brown and purple module heatmap -----------------------------

modules_of_interest = c("brown", "purple")

# Pull out list of genes in that module
submod <- module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) <- module_df$gene_id

subexpr <- expr_normalized_q85[submod$gene_id,]

keep_wgcna <- rownames(norm_data)%in%submod$gene_id

side_ha = rowAnnotation(foo = anno_mark(at = which(rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete),
                                        labels = rownames(norm_data[keep_wgcna,])[rownames(norm_data[keep_wgcna,])%in%metaplasia_markers_complete]),
                        annotation_name_gp = gpar(fontsize = 4))

Heatmap(t(scale(t(cpm(norm_data[keep_wgcna,], log = T)))),
        name = "Z-score",
        col = col_fun,
        row_title = "Genes from brown and\npurple modules",
        column_split = 2,
        column_title = c("Low-risk\nOLGIM",
                         "High-risk\nOLGIM"),
        show_row_names = F,
        show_column_names = F,
        top_annotation = ha,
        right_annotation = side_ha,
        row_names_gp = gpar(fontsize = 2),
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = F)

write_rds(keep_wgcna, "WGCNA_genes.rds")

dev.off()

rm(list = ls())

setwd("../")
