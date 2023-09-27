# 00 - Set up the environment ####

library(optparse)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
library(edgeR)
library(openxlsx)
library(factoextra)
library(ggsci)
library(tidyquant)
library(patchwork)
library(ggrepel)

# 01 - Set up working directory -------------------------------------------

setwd("./D01_20230720_TCGA_validation/")
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
                     c("-t", "--TCGA_STAD"),
                     default=NA,
                     type="character", 
                     help="Expression data for TCGA STAD cohort")
parser <- add_option(parser,
                     c("-s", "--validated_genes"),
                     default=NA,
                     type="character", 
                     help="List of 36 validated genes from the spatial assay")
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

opt$TCGA_STAD <- "/mnt/ix1/Projects_lite/20230308_IW_TCGAbiolinks/A01_RNAseq/STAD_RNAseq.rds"
opt$validated_genes <- "../../Visium/C03_pseudobulk_DE/35_overexpressed_signature_genes.csv"
opt$output_dir <- "./temp"
opt$plot_dir <- "./plots"


# 03 - Import TCGA STAD data ----------------------------------------------

tcga_stad <- readRDS(opt$TCGA_STAD)

## 3.1 - Extract expression data from assay -------------------------------

tcga_exp <- assay(tcga_stad)

# Rename gene ids
if(identical(rowData(tcga_stad)$gene_id, rownames(tcga_exp))) {
  cat("Switching ENSEMBL ids for HUGO symbols.\n")
  rownames(tcga_exp) <- rowData(tcga_stad)$gene_name
}

## 3.2 - Extract sample metadata from assay -------------------------------

tcga_samples <- as.data.frame(colData(tcga_stad))

if (sum(table(tcga_samples$paper_Molecular.Subtype)) <
    sum(table(TCGA_MolecularSubtype(tcga_samples$patient)$subtype$subtype))) {
  cat("Molecular subtypes from the original publication are outdated.\n")
  temp <- subset(merge(tcga_samples,
                       TCGA_MolecularSubtype(tcga_samples$patient)$subtype,
                       by.x = "patient",
                       by.y = "patients",
                       all = T), !is.na(patient))
  rownames(temp) <- temp$barcode
  temp <- temp[rownames(tcga_samples),]
  if(identical(rownames(temp), rownames(tcga_samples))) {
    cat(" Adding molecular subtypes from PanCancer study to new column \"PanCancer_Molecular.Subtype\".")
    temp = temp %>%
      mutate(subtype = gsub("GI.", "", subtype)) %>%
      dplyr::rename(PanCancer_Molecular.Subtype = subtype) %>%
      mutate(samples = NULL) %>%
      mutate(color = NULL)
  }
  tcga_samples = temp
  rm(temp)
}

tcga_samples$PanCancer_Molecular.Subtype[is.na(tcga_samples$PanCancer_Molecular.Subtype)] <- "Unknown"
tcga_samples$paper_Lauren.Class <- as.character(tcga_samples$paper_Lauren.Class)
tcga_samples$paper_Lauren.Class[is.na(tcga_samples$paper_Lauren.Class)] <- "Unknown"
tcga_samples$paper_Lauren.Class <- ifelse(tcga_samples$paper_Lauren.Class=="NA", "Unknown", tcga_samples$paper_Lauren.Class)

intestinal_tumor <- tcga_samples$patient[tcga_samples$paper_Lauren.Class=="Intestinal"&tcga_samples$shortLetterCode=="TP"]

for (i in 1:nrow(tcga_samples)) {
  if (tcga_samples$shortLetterCode[i]=="NT"&tcga_samples$patient[i]%in%intestinal_tumor) {
    tcga_samples$paper_Lauren.Class[i] = "Intestinal"
  }
  rm(i)
}

tcga_samples$PanCancer_Molecular.Subtype <- ifelse(tcga_samples$PanCancer_Molecular.Subtype=="HM-SNV", "MSI", tcga_samples$PanCancer_Molecular.Subtype)
tcga_samples$PanCancer_Molecular.Subtype <- as.character(tcga_samples$PanCancer_Molecular.Subtype)
tcga_samples$PanCancer_Molecular.Subtype[is.na(tcga_samples$PanCancer_Molecular.Subtype)] <- "Unkown"

for (i in 1:nrow(tcga_samples)) {
  if (tcga_samples$PanCancer_Molecular.Subtype[i]=="EBV"&tcga_samples$patient[i]%in%intestinal_tumor) {
    tcga_samples$PanCancer_Molecular.Subtype[i] = "EBV"
  } else if (tcga_samples$PanCancer_Molecular.Subtype[i]=="CIN"&tcga_samples$patient[i]%in%intestinal_tumor) {
    tcga_samples$PanCancer_Molecular.Subtype[i] = "CIN"
  } else if (tcga_samples$PanCancer_Molecular.Subtype[i]=="MSI"&tcga_samples$patient[i]%in%intestinal_tumor) {
    tcga_samples$PanCancer_Molecular.Subtype[i] = "MSI"
  } else if (tcga_samples$PanCancer_Molecular.Subtype[i]=="GS"&tcga_samples$patient[i]%in%intestinal_tumor) {
    tcga_samples$PanCancer_Molecular.Subtype[i] = "GS"
  }
  rm(i)
}
table(tcga_samples$paper_Lauren.Class)
table(tcga_samples$PanCancer_Molecular.Subtype)

rm(intestinal_tumor)

tcga_samples <- tcga_samples %>%
  rename(Lauren = paper_Lauren.Class) %>%
  rename(Molecular_subtype = PanCancer_Molecular.Subtype)


## 3.3 - Generate the DGElist for differential expression -----------------

tcga_stad <- DGEList(counts = tcga_exp,
                     samples = tcga_samples,
                     group = factor(tcga_samples$shortLetterCode,
                                    levels = c("NT", "TP"),
                                    ordered = T),
                     genes = rowData(tcga_stad))

tcga_stad <- tcga_stad[!duplicated(rownames(tcga_stad)),]


### 3.3.1 - Subset to intestinal-type GC ----------------------------------

tcga_stad <- tcga_stad[,tcga_stad$samples$Lauren == "Intestinal"]
table(group = tcga_stad$samples$group,
      Mol_subtype = tcga_stad$samples$Molecular_subtype)
#      Mol_subtype
# group CIN EBV GS MSI Unknown
#   NT  11   0  0   6       1
#   TP  94  13  8  47      18


### 3.3.2 - Save metadata for the paper -----------------------------------

TCGA_metadata <- data.frame(TCGA_barcode = tcga_stad$samples$barcode,
                            patient_id = tcga_stad$samples$patient,
                            sample_type = tcga_stad$samples$group,
                            mapped_reads = tcga_stad$samples$lib.size,
                            Lauren_class = tcga_stad$samples$Lauren,
                            location = tcga_stad$samples$paper_Anatomic.Region,
                            molecular_subtype = tcga_stad$samples$Molecular_subtype)

write.xlsx(list(TCGA_metadata = TCGA_metadata),
           file.path(opt$output_dir,
                     "STable4_TCGA_cohort_metadata.xlsx"))

## 3.4 - Preprocessing the data -------------------------------------------

par(mfrow = c(1,2))
plotDensities(cpm(tcga_stad$counts, log = T),
              legend = F,
              main = "A. Raw data (TCGA)")


### 3.4.1 - Filter lowly expressed genes ----------------------------------

keep.exprs <- filterByExpr(tcga_stad, group = tcga_stad$samples$group)
table(keep.exprs)
# keep.exprs
# FALSE  TRUE 
# 25916 33511 
tcga_stad_filtered <- tcga_stad[keep.exprs,, keep.lib.sizes=FALSE]
dim(tcga_stad_filtered)
# [1] 33511   198
plotDensities(cpm(tcga_stad_filtered$counts, log = T),
              legend = F,
              main = "B. Filtered data (TCGA)")

### 3.4.2 - Normalizing the data ------------------------------------------

tcga_stad_norm <- calcNormFactors(tcga_stad_filtered, method = "TMM")
par(mfrow=c(1,2))

tcga_stad_lcpm <- cpm(tcga_stad_filtered, log = T)
boxplot(tcga_stad_lcpm, las = 2, main = "A. Unnormalized data (TCGA)",
        ylab="Log-cpm")

tcga_stad_lcpm <- cpm(tcga_stad_norm, log = T)
boxplot(tcga_stad_lcpm, las = 2, main = "B. Normalized data (TCGA)",
        ylab="Log-cpm")

par(mfrow = c(1,1))


## 3.5 - PCA --------------------------------------------------------------

tcga_pc <- prcomp(t(tcga_stad_norm$counts), scale. = T)
exp_variance_tcga <- round(100*(tcga_pc$sdev^2/sum(tcga_pc$sdev^2)), 2)
fviz_screeplot(tcga_pc, addlabels=T)
tcga_pc_df <- as.data.frame(tcga_pc$x)
tcga_pc_df$group <- tcga_stad_norm$samples$group


tcga_group <- ggplot(tcga_pc_df, aes(PC1, PC2, fill = group)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(size = 3, shape = 21, colour = "black", alpha = 0.8) +
  xlab(paste0("PC1: ", exp_variance_tcga[1], "%")) +
  ylab(paste0("PC2: ", exp_variance_tcga[2], "%")) +
  labs(fill = "Sample type:") +
  scale_fill_aaas() +
  theme_tq() +
  theme(legend.position = "bottom")

tcga_subtype <- ggplot(tcga_pc_df, aes(PC1, PC2, fill = tcga_stad_norm$samples$Molecular_subtype)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(size = 3, shape = 21, colour = "black", alpha = 0.8) +
  xlab(paste0("PC1: ", exp_variance_tcga[1], "%")) +
  ylab(paste0("PC2: ", exp_variance_tcga[2], "%")) +
  labs(fill = "Molecular subtype:") +
  scale_fill_brewer(palette = "Dark2") +
  theme_tq() +
  theme(legend.position = "bottom")

tcga_group + tcga_subtype


# 04 - Differential expression analysis -----------------------------------

group <- tcga_stad_norm$samples$group
patient <- tcga_stad_norm$samples$patient

design <- model.matrix(~0+group+patient)
colnames(design) <- gsub("group", "", gsub("patient",  "",  colnames(design)))
colnames(design) <- make.names(colnames(design))
rownames(design) <- colnames(tcga_stad_norm)

cont <- makeContrasts(TvsN = TP - NT,
                      levels = design)

par(mfrow = c(1,2))

v <- voom(tcga_stad_norm, design, plot = T)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = cont)
efit <- eBayes(vfit)
plotSA(efit, main = "Final model: Mean-variance trend")
tfit <- treat(vfit, lfc = log2(1.25))
plotSA(tfit, main = "Final model: Mean-variance trend")

res <- list(TvsN_eBayes = topTable(efit, 1, confint = T, number = Inf)[,-c(1:4,8:10)],
            TvsN_treat = topTreat(tfit, 1, confint = T, number = Inf))

write.xlsx(res, file.path(opt$output_dir,
                          "TCGA_STAD_resuts_intestinal_subset_TvsN.xlsx"))

## 4.1 - GIM-specific genes -----------------------------------------------

spatially_resolved_genes_base <- read.xlsx("/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/C03_pseudobulk_DE/metaplasia_overexpressed_genes.xlsx")
spatially_resolved_genes_pit <- read.xlsx("/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/C03_pseudobulk_DE/metaplasia_overexpressed_genes.xlsx",
                                          sheet = 2)
identical(spatially_resolved_genes_base$symbol,
          spatially_resolved_genes_pit$symbol)
spatially_resolved_genes_average <- data.frame(symbol = spatially_resolved_genes_base$symbol,
                                               fc_GIMvsbase = spatially_resolved_genes_base$logFC,
                                               P.adj_base = spatially_resolved_genes_base$adj.P.Val,
                                               fc_GIMvspit = spatially_resolved_genes_pit$logFC,
                                               P.adj_pit = spatially_resolved_genes_pit$adj.P.Val)
spatially_resolved_genes_average$aveFC <- (spatially_resolved_genes_average$fc_GIMvsbase+spatially_resolved_genes_average$fc_GIMvspit)/2
rownames(spatially_resolved_genes_average) <- spatially_resolved_genes_average$symbol
write.xlsx(spatially_resolved_genes_average, file.path(opt$output_dir, "Visium_results_all_gim_genes.xlsx"))

spatially_resolved_genes_signature <- c(read.xlsx("/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/C03_pseudobulk_DE/35_genes.xlsx")$symbol, "ADGRG7")
validated_spatially_resolved_genes <- spatially_resolved_genes_average[spatially_resolved_genes_signature,]
write.xlsx(validated_spatially_resolved_genes, file.path(opt$output_dir, "Visium_results_36_genes.xlsx"))

tcga_spatially_resolved <- res$TvsN_eBayes[spatially_resolved_genes_signature,]
tcga_spatially_resolved$Status <- ifelse(tcga_spatially_resolved$logFC>=0, "Included", "Excluded")
write.xlsx(tcga_spatially_resolved, file.path(opt$output_dir, "TCGA_results_spatial_signature_36_genes.xlsx"))

tcga_spatially_resolved %>%
  mutate(DE_status = ifelse(adj.P.Val<=0.05&logFC>0, "Up, DEG",
                            ifelse(adj.P.Val<=0.05&logFC<0, "Down, DEG",
                                   ifelse(logFC>0, "Up, non-DEG", "Down, non-DEG")))) %>%
  mutate(spatial_fc = validated_spatially_resolved_genes$aveFC) %>%
  ggplot(aes(spatial_fc, logFC, fill = DE_status, label = gene_name)) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf) +
  xlab(expression(log[2]~"FC (spatially-resolved metaplasia vs normal glands)")) +
  ylab(expression(log[2]~"FC (TCGA tumor vs tumor-adjacent controls)")) +
  scale_fill_manual(
    values = c("Down, DEG" = "navy",
               "Up, DEG" = "firebrick",
               "Down, non-DEG" = "blue",
               "Up, non-DEG" = "red")) +
  labs(col = "Differential expression status:") +
  coord_cartesian(xlim = c(0, 7), ylim = c(-4, 6), expand = TRUE) +
  labs(fill = "DE status:") +
  theme_classic() +
  theme(legend.position = "bottom")

write.csv(tcga_spatially_resolved$gene_name[tcga_spatially_resolved$Status=="Included"],
          file.path(opt$output_dir, "26_validated_spatial_genes.csv"))

rm(list = ls())

setwd("../")
