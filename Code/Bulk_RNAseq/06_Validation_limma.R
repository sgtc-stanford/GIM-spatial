# 00 - Load required R packages -------------------------------------------
library(optparse)
library(edgeR)
library(ComplexHeatmap)
library(RColorBrewer)
library(factoextra)
library(tidyverse)
library(tidyquant)
library(ggpubr)
library(ggdensity)
library(circlize)
library(openxlsx)
library(ggrepel)
library(patchwork)

# 01 - Set up working directory -------------------------------------------

setwd("./B03_20230607_Validation_limma")
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
                     help="DGElist used for differential expression analysis with limma-voom")
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

opt$input_DGElist <- "../B02_20230607_Validation_DGEList"
opt$output_dir <- "./temp"
opt$plot_dir <- "./plots"


# 03 - Preprocess data ----------------------------------------------------
pdf(file.path(opt$plot_dir,
              "Validation_limma.pdf"),
    height = 10.41,
    width = 10.41)
raw_data <- read_rds(file.path(opt$input_DGElist, "Raw_data_DGEList.rds"))

# Look at the mean and median library size
c(mean(raw_data$samples$lib.size) * 1e-6, median(raw_data$samples$lib.size) * 1e-6)
# [1] 49.87928 47.73916
# See how many genes are not expressed in any of the samples

table(rowSums(raw_data$counts==0)==ncol(raw_data))
# FALSE
# 32994

par(mfrow = c(1,2))
plotDensities(log1p(raw_data$counts),
              legend = F,
              main = "A. Raw data")

# Filter lowly expressed genes
keep.exprs <- filterByExpr(raw_data, group=raw_data$samples$group)
table(keep.exprs)
# keep.exprs
# FALSE  TRUE 
# 13922 19072 
filtered_data <- raw_data[keep.exprs,, keep.lib.sizes=FALSE]
dim(filtered_data)
# [1] 19072   215
plotDensities(log1p(filtered_data$counts),
              legend = F,
              main = "B. Filtered data")

# Normalizing the data
norm_data <- calcNormFactors(filtered_data, method = "TMM")
par(mfrow=c(1,2))
lcpm <- cpm(filtered_data, log = T)
boxplot(lcpm, las = 2, main = "A. Unnormalized data", ylab="Log-cpm")
lcpm <- cpm(norm_data, log = T)
boxplot(lcpm, las = 2, main = "B. Normalized data", ylab="Log-cpm")

par(mfrow = c(1,1))

cor_samples <- cor(cpm(norm_data$counts, log = T))

plot(hclust(as.dist(cor_samples)))

ha <- HeatmapAnnotation(`Biopsy site` = factor(ifelse(norm_data$samples$Location=="ANT", "Antrum", "Body")),
                        `OLGIM stage` = norm_data$samples$OLGIM,
                        Batch = factor(norm_data$samples$Batch),
                        col = list(`Biopsy site` = c("Antrum" = brewer.pal(8, "Greens")[7],
                                                     "Body" = brewer.pal(8, "Greens")[4]),
                                   `OLGIM stage` = c("OLGIM 0" = brewer.pal(11, "RdYlBu")[7],
                                                     "OLGIM I" = brewer.pal(11, "RdYlBu")[8],
                                                     "OLGIM II" = brewer.pal(11, "RdYlBu")[9],
                                                     "OLGIM III" = brewer.pal(11, "RdYlBu")[10],
                                                     "OLGIM IV" = brewer.pal(11, "RdYlBu")[11]),
                                   Batch = c("B2" = brewer.pal(3, "Dark2")[1],
                                             "B3" = brewer.pal(3, "Dark2")[2],
                                             "B4" = brewer.pal(3, "Dark2")[3])))

draw(Heatmap(cor_samples,
        name = "Pearson's correlation",
        top_annotation = ha,
        show_row_names = F,
        show_column_names = F,
        show_row_dend = F,
        show_column_dend = F,
        column_dend_reorder = T))


# 04 - PCA ----------------------------------------------------------------

group <- paste(gsub(" ", "", norm_data$samples$OLGIM), norm_data$samples$Location, sep = "_")
pc <- prcomp(t(norm_data$counts), scale. = T)
exp_variance <- round(100*(pc$sdev^2/sum(pc$sdev^2)), 2)
print(fviz_screeplot(pc, addlabels=T))
pc_df <- as.data.frame(pc$x)
location <- factor(ifelse(norm_data$samples$Location=="BOD", "Body", "Antrum"))
batch <- factor(norm_data$samples$Batch)
pal <- c(brewer.pal(9, "YlOrRd")[c(1,3,5,7,9)], brewer.pal(11, "RdYlBu")[7:11])
pal <- pal[c(1,6,2,7,3,8,4,9,5,10)]
group
pc_df$group <- group
pc_df$group2 <- ifelse(str_detect(pc_df$group, "OLGIMIII")|str_detect(pc_df$group, "OLGIMIV"), "OLGIM III-IV", "OLGIM 0-I-II")
pc_df$location <- location

pc_validation_batch <- ggplot(pc_df, aes(PC1, PC2, fill = batch)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  # geom_hdr(probs = c(0.5), alpha = 0.5) +
  geom_point(size = 2, shape = 21, colour = "black", alpha = 0.8) +
  xlab(paste0("PC1: ", exp_variance[1], "%")) +
  ylab(paste0("PC2: ", exp_variance[2], "%")) +
  scale_fill_viridis_d() +
  theme_tq() +
  theme(legend.position = "bottom")

pc_validation_all <- ggplot(pc_df, aes(PC1, PC2, fill = group)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  # geom_hdr(probs = c(0.5), alpha = 0.5) +
  geom_point(size = 2, shape = 21, colour = "black", alpha = 0.8) +
  xlab(paste0("PC1: ", exp_variance[1], "%")) +
  ylab(paste0("PC2: ", exp_variance[2], "%")) +
  scale_fill_manual(labels = c('OLGIM 0 - Antrum',
                               'OLGIM 0 - Body',
                               'OLGIM I - Antrum',
                               'OLGIM I - Body',
                               'OLGIM II - Antrum',
                               'OLGIM II - Body',
                               'OLGIM III - Antrum',
                               'OLGIM III - Body',
                               'OLGIM IV - Antrum',
                               'OLGIM IV - Body'),
                    values = pal) +
  theme_tq() +
  theme(legend.position = "bottom")

pc_validation_split <- ggplot(pc_df, aes(PC1, PC2, fill = group)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  # geom_hdr(probs = c(0.5), alpha = 0.5) +
  geom_point(size = 2, shape = 21, colour = "black", alpha = 0.8) +
  facet_grid(rows = vars(group2), cols = vars(location)) +
  xlab(paste0("PC1: ", exp_variance[1], "%")) +
  ylab(paste0("PC2: ", exp_variance[2], "%")) +
  labs(fill = "") +
  scale_fill_manual(labels = c('OLGIM 0 - Antrum',
                               'OLGIM 0 - Body',
                               'OLGIM I - Antrum',
                               'OLGIM I - Body',
                               'OLGIM II - Antrum',
                               'OLGIM II - Body',
                               'OLGIM III - Antrum',
                               'OLGIM III - Body',
                               'OLGIM IV - Antrum',
                               'OLGIM IV - Body'),
                    values = pal) +
  theme_tq() +
  theme(legend.position = "none")

print((pc_validation_batch + pc_validation_all) / pc_validation_split)

# 05 - Differential expression analysis -----------------------------------

cond <- paste(norm_data$samples$Location,
              ifelse(norm_data$samples$OLGIM%in%c("OLGIM 0", "OLGIM I", "OLGIM II"), "low_risk", "high_risk"),
              sep = "_")
cond <- factor(cond,
               levels = c("BOD_low_risk", "BOD_high_risk",
                          "ANT_low_risk", "ANT_high_risk"),
               ordered = T)
design <- model.matrix(~0+cond)
colnames(design) <- levels(cond)
contr.matrix <- makeContrasts(high_vs_low_risk_Bod = BOD_high_risk - BOD_low_risk,
                              high_vs_low_risk_Ant = ANT_high_risk - ANT_low_risk,
                              Interaction_BodAnt = (BOD_high_risk - BOD_low_risk) - (ANT_high_risk - ANT_low_risk),
                              levels = design)
v <- voom(norm_data, design)
corfit <- duplicateCorrelation(v,
                               design,
                               block = str_split(norm_data$samples$id, "_", simplify = T)[,1])
corfit$consensus
# [1] 0.3087526

par(mfrow=c(1,2))
v <- voom(norm_data, design, block = norm_data$samples$patient, correlation =
            corfit$consensus, plot = T)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
par(mfrow=c(1,1))

summary(decideTests(efit))
#        high_vs_low_risk_Bod high_vs_low_risk_Ant Interaction_BodAnt
# Down                   1337                 1573                 39
# NotSig                16161                15469              19021
# Up                     1574                 2030                 12

results_eBayes <- list(high_vs_low_risk_Bod = topTable(efit, 1, n = Inf, confint = T),
                       high_vs_low_risk_Ant = topTable(efit, 2, n = Inf, confint = T),
                       Interaction_BodAnt = topTable(efit, 3, n = Inf, confint = T))

tfit <- treat(vfit, lfc = log2(1.25))
summary(decideTests(tfit))
#        high_vs_low_risk_Bod high_vs_low_risk_Ant Interaction_BodAnt
# Down                    130                  223                 12
# NotSig                18449                18154              19060
# Up                      493                  695                  0

results_treat <- list(high_vs_low_risk_Bod = topTreat(tfit, 1, n = Inf, confint = T),
                      high_vs_low_risk_Ant = topTreat(tfit, 2, n = Inf, confint = T),
                      Interaction_BodAnt = topTreat(tfit, 3, n = Inf, confint = T))

write.xlsx(results_eBayes, file.path(opt$output_dir, "eBayes_grouped_OLGIM_multi_level_design.xlsx"))
write.xlsx(results_treat, file.path(opt$output_dir, "Treat_1.25_grouped_OLGIM_multi_level_design.xlsx"))

# 06 - Plots --------------------------------------------------------------

differential_bod <- rownames(subset(results_treat$high_vs_low_risk_Bod, adj.P.Val< 0.05))
differential_ant <- rownames(subset(results_treat$high_vs_low_risk_Ant, adj.P.Val< 0.05))
differential_bod_ant <- rownames(subset(results_eBayes$Interaction_BodAnt, adj.P.Val< 0.05))

degs <- differential_ant[differential_ant%in%differential_bod]

# Remove genes with significant interaction term
degs <- setdiff(degs, differential_bod_ant) 

vennDiagram(vennCounts(decideTests(efit)))
vennDiagram(vennCounts(decideTests(tfit)))

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
interest_genes <- metaplasia_markers_complete[metaplasia_markers_complete%in%degs]
all(interest_genes%in%norm_data$genes$SYMBOL)

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))


## 6.1 - Body -------------------------------------------------------------

results_treat$high_vs_low_risk_Bod$label <- results_treat$high_vs_low_risk_Bod$SYMBOL%in%interest_genes
results_treat$high_vs_low_risk_Bod$label <- ifelse(results_treat$high_vs_low_risk_Bod$label, results_treat$high_vs_low_risk_Bod$SYMBOL, "")


### 6.1.1 - MA plot -------------------------------------------------------

print(results_treat$high_vs_low_risk_Bod %>%
  mutate(Sign = ifelse(logFC>0 & adj.P.Val<=0.05, "Up", ifelse(logFC<0 & adj.P.Val<=0.05, "Down", "Non-DEG"))) %>%
  mutate(trans = ifelse(Sign =="Non-DEG", 0.2, 0.8)) %>%
  ggplot(aes(AveExpr, logFC, col = Sign, label = label, alpha = trans)) +
  ggtitle("Body") +
  geom_point(size = 2) +
  geom_text_repel(col = "black", max.overlaps = Inf) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Average expression (lcpm)") +
  ylab(expression(log[2]~"FC")) +
  scale_color_manual(values = c("Up" = "red", "Down" = "navy", "Non-DEG" = "grey"),
                     name = "Status",
                     labels = c("Up-regulated", "Down-regulated", "Non-differentially expressed")) +  scale_alpha(guide = "none") +
  theme_classic() +
  theme(legend.position = "bottom"))


### 6.1.2 - Volcano plot --------------------------------------------------

print(results_treat$high_vs_low_risk_Bod %>%
  mutate(Sign = ifelse(logFC>0 & adj.P.Val<=0.05, "Up", ifelse(logFC<0 & adj.P.Val<=0.05, "Down", "Non-DEG"))) %>%
  mutate(trans = ifelse(Sign =="Non-DEG", 0.2, 0.8)) %>%
  mutate(logp = -log10(adj.P.Val)) %>%
  ggplot(aes(logFC, logp, col = Sign, label = label, alpha = trans)) +
  ggtitle("Body") +
  geom_point(size = 2) +
  geom_text_repel(col = "black", max.overlaps = Inf) +
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  xlab(expression(log[2]~"FC")) +
  ylab(expression("-log"[10]~"adjusted P-val")) +
  scale_color_manual(values = c("Up" = "red", "Down" = "navy", "Non-DEG" = "grey"),
                     name = "Status",
                     labels = c("Up-regulated", "Down-regulated", "Non-differentially expressed")) +
  scale_alpha(guide = "none") +
  theme_classic() +
  theme(legend.position = "bottom"))


## 6.2 - Antrum -----------------------------------------------------------

results_treat$high_vs_low_risk_Ant$label <- results_treat$high_vs_low_risk_Ant$SYMBOL%in%interest_genes
results_treat$high_vs_low_risk_Ant$label <- ifelse(results_treat$high_vs_low_risk_Ant$label, results_treat$high_vs_low_risk_Ant$SYMBOL, "")


### 6.2.1 MA plot ---------------------------------------------------------

print(results_treat$high_vs_low_risk_Ant %>%
  mutate(Sign = ifelse(logFC>0 & adj.P.Val<=0.05, "Up", ifelse(logFC<0 & adj.P.Val<=0.05, "Down", "Non-DEG"))) %>%
  mutate(trans = ifelse(Sign =="Non-DEG", 0.2, 0.8)) %>%
  ggplot(aes(AveExpr, logFC, col = Sign, label = label, alpha = trans)) +
  ggtitle("Antrum") +
  geom_point(size = 2) +
  geom_text_repel(col = "black", max.overlaps = Inf) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Average expression (lcpm)") +
  ylab(expression(log[2]~"FC")) +
  scale_color_manual(values = c("Up" = "red", "Down" = "navy", "Non-DEG" = "grey"),
                     name = "Status",
                     labels = c("Up-regulated", "Down-regulated", "Non-differentially expressed")) +  scale_alpha(guide = "none") +
  scale_alpha(guide = "none") +
  theme_classic() +
  theme(legend.position = "bottom"))


### 6.2.2 Volcano plot ----------------------------------------------------

print(results_treat$high_vs_low_risk_Ant %>%
  mutate(Sign = ifelse(logFC>0 & adj.P.Val<=0.05, "Up", ifelse(logFC<0 & adj.P.Val<=0.05, "Down", "Non-DEG"))) %>%
  mutate(trans = ifelse(Sign =="Non-DEG", 0.2, 0.8)) %>%
  mutate(logp = -log10(adj.P.Val)) %>%
  ggplot(aes(logFC, logp, col = Sign, label = label, alpha = trans)) +
  ggtitle("Antrum") +
  geom_point(size = 2) +
  geom_text_repel(col = "black", max.overlaps = Inf) +
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  xlab(expression(log[2]~"FC")) +
  ylab(expression("-log"[10]~"adjusted P-val")) +
  scale_color_manual(values = c("Up" = "red", "Down" = "navy", "Non-DEG" = "grey"),
                     name = "Status",
                     labels = c("Up-regulated", "Down-regulated", "Non-differentially expressed")) +
  scale_alpha(guide = "none") +
  theme_classic() +
  theme(legend.position = "bottom"))


## 6.3 Venn Diagram -------------------------------------------------------

myCol <- brewer.pal(3, "Pastel2")

venn <- VennDiagram::venn.diagram(
  x = list(differential_bod, differential_ant, differential_bod_ant),
  category.names = c("Differentially expressed\ngenes body (FC threshold 1.25)" ,
                     "Differentially expressed\ngenes antrum (FC threshold 1.25)" ,
                     "Interaction (no FC threshold)"),
  filename = NULL,
  # Circles
  lwd = 1,
  lty = 'blank',
  fill = myCol,
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
grid.newpage()
grid.draw(venn)


## 6.4 Heatmap ------------------------------------------------------------

# Heatmap exported as pdf. Dimensions (height X width): 10.41 X 6.25 
side_ha_deg = rowAnnotation(foo = anno_mark(at = which(rownames(norm_data[degs,])%in%interest_genes),
                                            labels = rownames(norm_data[degs,])[rownames(norm_data[degs,])%in%interest_genes]),
                            annotation_name_gp = gpar(fontsize = 4))

draw(Heatmap(t(scale(t(cpm(norm_data[degs,], log = T)))),
        name = "Z-score",
        col = col_fun,
        column_split = 3,
        column_title = c("High-risk antrum\nand body",
                         "Low-risk\nbody",
                         "Low-risk\nantrum"),
        row_title = "Common differentially expressed\ngenes in body and antrum",
        show_column_names = F,
        show_row_names = F,
        top_annotation = ha,
        right_annotation = side_ha_deg,
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        show_row_dend = F))

save(discovery_ebayes = results_eBayes,
     descovery_treat = results_treat,
     file = "DE_results_validation.rda")

saveRDS(norm_data, "Validation_norm_data.rds")

dev.off()

rm(list = ls())

setwd("../")
