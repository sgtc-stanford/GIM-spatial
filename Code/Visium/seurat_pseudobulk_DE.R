#################################
# Pseudobulk DE gene expression #
#################################
# Based on: https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

library(ComplexHeatmap)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(optparse)
library(RColorBrewer)
library(Seurat)
library(tidyverse)


# Setup environment -------------------------------------------------------

set.seed(42)
sessionInfo()

# Parse command line arguments  -------------------------------------------

parser <- OptionParser()
parser <- add_option(parser, c('-c', '--counts'), default=NA, type='character', 
                     help='Path to counts matrix (genes x samples)')
parser <- add_option(parser, c('-s', '--samples'), default=NA, type='character',
                     help='Path to sample metadata table(.csv)')
parser <- add_option(parser, c('-o', '--output_dir'), default=NA, type='character',
                     help='Output directory')
opt <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

opt$counts <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/C01_pseudobulk/pseudobulk_A08_manual_spot_annotation_raw.txt'
opt$samples <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/samples.csv'
opt$output_dir <- '/mnt/ix1/Projects/M075_201130_GIM_P01/Visium/C03_pseudobulk_DE'


# Load data ---------------------------------------------------------------

counts <- read_tsv(opt$counts)
samples <- read_csv(opt$samples)

signature_genes <- read_csv("/mnt/ix1/Projects/M075_201130_GIM_P01/bulkRNA/A04_221222_DEA_WGCNA_validation/validated_genes_c5_for_spatial_mapping.csv")
signature_genes <- signature_genes %>% pull(x)

fn <- tools::file_path_sans_ext(basename(opt$counts)) 
prefix <- paste(opt$output_dir, fn, sep='/')

# metadata <- tibble(
#   files = list.files(path=data.dir, pattern = '*.rds'),
#   sample = tools::file_path_sans_ext(files),
# )
# 
# metadata <- left_join(metadata, samples, by='sample')

# Create DEG list ---------------------------------------------------------

counts_df <- counts %>% 
  select(-(Gene)) %>% 
  as.data.frame()
rownames(counts_df) <- counts %>% pull(Gene)

pseudobulk_metadata <- tibble(
  pseudobulk = colnames(counts_df)
)

pseudobulk_metadata <- pseudobulk_metadata %>% 
  separate(pseudobulk, c('array', 'LIMS_dissection', 'group'), sep='_', remove=FALSE) %>%
  mutate(LIMS_dissection = as.double(LIMS_dissection)) %>%
  left_join(samples, by='LIMS_dissection')

group <- factor(pseudobulk_metadata$group,
                levels=c('Base', 'Pit', 'Metaplasia'),
                ordered=TRUE)
y <- DGEList(counts=counts_df, group = group, samples = pseudobulk_metadata %>% 
               mutate(dissection = factor(LIMS_dissection)) %>% as.data.frame())
y$samples$array <- factor(pseudobulk_metadata$array.x)
y$samples$dissection <- factor(pseudobulk_metadata$LIMS_dissection)
y$samples$region <- factor(pseudobulk_metadata$region)
y$samples$metaplasia <- factor(pseudobulk_metadata$metaplasia)
y$samples$OLGIM <- factor(pseudobulk_metadata$OLGIM)
y$samples


# Vis library size --------------------------------------------------------

col <- as.vector(recode_factor(y$samples$group, 
                               Base = "#4DAF4A", 
                               Pit = "#377EB8", 
                               Metaplasia = "#FF7F00"))


ggplot(data = y$samples, mapping = aes(x=group, y=lib.size, fill=dissection)) +
  geom_bar(stat='identity') +
  labs(x = element_blank(),
       y = 'Library size',
       fill = 'Patient') +
  theme(text = element_text(size=20))
ggsave(paste(prefix, 'library_size_barplot.jpeg', sep='_'), dpi = 300)


# Convert raw counts to cpm -----------------------------------------------

counts_cpm <- cpm(y)
counts_lcpm <- cpm(y, log=TRUE)

L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
c(L, M)

boxplot(y$samples$lib.size)
abline(h=mean(y$samples$lib.size), lty=2)

with(y$samples, boxplot(lib.size~sample+group), las =2 )

boxplot(log1p(y$counts))
boxplot(counts_lcpm)

summary(counts_lcpm)

# Filter out low expressed genes ------------------------------------------

# Count number of genes that are absent in all samples
table(rowSums(y$counts == 0)==nrow(y$samples))

keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)

keep.genes <- keep[keep == T]
length(keep.genes)

# Plot
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(y)

jpeg(paste(prefix, 'filter_low_expression_genes.jpeg', sep='_'), width = 1600, height = 800, res=150)
par(mfrow=c(1,2))
plot(density(counts_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main='A. Raw data',
     xlab = 'Log-cpm')
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(counts_lcpm[,i])
  print(col[i])
  lines(den$x, den$y, col=col[i], lwd=2)  
}
legend('topright', colnames(y), text.col = col, bty='n')
counts_lcpm <- cpm(y, log=TRUE)
plot(density(counts_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(counts_lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend('topright', colnames(y), text.col = col, bty='n')
dev.off()

# Normalise gene expression distribution ----------------------------------

# Normalise by library size
y <- calcNormFactors(y, method = 'TMM')
y$samples

jpeg(paste(prefix, 'TMM_normalisation.jpeg', sep='_'), width = 1600, height = 800, res=150)
par(mfrow=c(1,2))
boxplot(counts_lcpm, las=2, col=col, main='')
title(main='A. Unnormalised data', ylab='Log-cpm')
counts_lcpm <- cpm(y, log=TRUE)
boxplot(counts_lcpm, las=2, col=col, main='')
title(main='B. Normalised data', ylab='Log-cpm')
dev.off()

# Unsupervised clustering of samples --------------------------------------
col.dissection <- y$samples$dissection
levels(col.dissection) <- brewer.pal(nlevels(col.dissection), 'Set2')
col.dissection <- as.character(col.dissection)

jpeg(paste(prefix, 'MDS_plot_groups_and_patients.jpeg', sep='_'), width = 1600, height = 800, res=150)
par(mfrow=c(1,2))
plotMDS(counts_lcpm, labels=group, col=col, dim=c(1,2))
abline(h=0, v=0, lty=2)
title(main='A. Sample groups')
plotMDS(counts_lcpm, labels=y$samples$dissection, col=col.dissection, dim=c(1,2))
abline(h=0, v=0, lty=2)
title(main='B. Patient')
dev.off()

jpeg(paste(prefix, 'MDS_plot_region.jpeg', sep='_'), width = 1600, height = 800, res=150)
par(mfrow=c(1,2))
col.region <- y$samples$region 
levels(col.region) <- brewer.pal(nlevels(col.region), 'Set1')
col.region <- as.character(col.region)
plotMDS(counts_lcpm, labels=y$samples$region, col=col.region, dim=c(1,2))
plotMDS(counts_lcpm, labels=y$samples$region, col=col.region, dim=c(2,3))
dev.off()

jpeg(paste(prefix, 'MDS_plot_metaplasia_severity.jpeg', sep='_'), width = 1600, height = 800, res=150)
par(mfrow=c(1,2))
col.metaplasia <- y$samples$metaplasia
levels(col.metaplasia) <- brewer.pal(nlevels(col.metaplasia), 'Set1')
col.metaplasia <- as.character(col.metaplasia)
plotMDS(counts_lcpm, labels=y$samples$metaplasia, col=col.metaplasia, dim=c(1,2))
plotMDS(counts_lcpm, labels=y$samples$metaplasia, col=col.metaplasia, dim=c(2,3))
dev.off()

pc <- prcomp(t(counts_lcpm))
pc$sdev
exp.var <- round((pc$sdev^2 / sum(pc$sdev^2)) * 100, 2)
par(mfrow=c(1,1))
screeplot(pc)
pc.df <- as.data.frame(pc$x)

ggplot(data=pc.df, mapping=aes(x=PC1, y=PC2, col=group)) +
  geom_point(size=3) + 
  geom_vline(xintercept=0, lty=2) + 
  geom_hline(yintercept=0, lty=2) +
  ggtitle('Conditions') +
  scale_color_manual(values=c('Base' = col[1], 'Pit' = col[3], 'Metaplasia' = col[2])) +
  theme(legend.position = 'bottom') +
  theme_classic()
ggsave(paste(prefix, 'PCA_groups.jpeg', sep='_'), dpi = 300)

# Select covariates -------------------------------------------------------

counts_lcpm_corrected <- removeBatchEffect(counts_lcpm, batch = y$samples$LIMS_dissection)

jpeg(paste(prefix, 'MDS_plot_batch_correction_by_patient.jpeg', sep='_'), width = 1600, height = 800, res=150)
par(mfrow=c(1,2))
plotMDS(counts_lcpm, labels=y$samples$dissection, col=col.dissection, dim=c(1,2))
abline(h=0, v=0, lty=2)
title(main='A. Patients')
plotMDS(counts_lcpm_corrected, labels=y$samples$dissection, col=col.dissection, dim=c(1,2))
abline(h=0, v=0, lty=2)
title(main='B. Batch corrected')
dev.off()

counts_lcpm_corrected <- removeBatchEffect(counts_lcpm, batch = y$samples$region)

jpeg(paste(prefix, 'MDS_plot_batch_correction_by_region.jpeg', sep='_'), width = 1600, height = 800, res=150)
par(mfrow=c(1,2))
plotMDS(counts_lcpm, labels=y$samples$region, col=col.region, dim=c(1,2))
abline(h=0, v=0, lty=2)
title(main='A. Regions')
plotMDS(counts_lcpm_corrected, labels=y$samples$region, col=col.region, dim=c(1,2))
abline(h=0, v=0, lty=2)
title(main='B. Batch corrected')
dev.off()

# Create design matrix and contrasts --------------------------------------

dissection <- y$samples$dissection
region <- y$samples$region
design <- model.matrix(~0+group + dissection + region)
colnames(design) <- gsub('group', '', colnames(design))
design

contr.matrix <- makeContrasts(
  MetaplasiavsBase = Metaplasia - Base,
  MetaplasiavsPit = Metaplasia - Pit,
  BasevsPit = Base - Pit,
  levels = colnames(design)
)
contr.matrix

# Removing heteroscedascity from count data -------------------------------

jpeg(paste(prefix, 'remove_heteroscedascity.jpeg', sep='_'), width = 1600, height = 800, res=150)
par(mfrow=c(1,2))
v <- voom(y, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main='Final model: Mean-variance trend')
dev.off()



# Examining the no. of DE genes -------------------------------------------
summary(decideTests(efit))
tfit <- treat(vfit, lfc=0.0)
dt <- decideTests(tfit)
summary(dt)

saveRDS(vfit, paste(prefix, 'vfit.rds', sep='_'))
saveRDS(efit, paste(prefix, 'efit.rds', sep='_'))
saveRDS(tfit, paste(prefix, 'tfit.rds', sep='_'))
saveRDS(dt, paste(prefix, 'dt.rds', sep='_'))

jpeg(paste(prefix, 'venn_diagram_metaplasia_DE_genes.jpeg', sep='_'), width = 1200, height = 1200, res=150)
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
vennDiagram(dt[, 1:2], circle.col = c('turquoise', 'salmon'), cex=0.8)
dev.off()

jpeg(paste(prefix, 'venn_diagram_metaplasia_DE_overexpressed_genes.jpeg', sep='_'), width = 1200, height = 1200, res=150)
dt.upregulated <- dt[which(dt[,1] > 0 & dt[,2]>0), ]
vennDiagram(dt.upregulated[, 1:2], circle.col = c('turquoise', 'salmon'), cex=1)
dev.off()

plotMD(tfit, column=1, status=dt[, 1], main=colnames(tfit)[1],
       xlim=c(-8,13))

# Export DE genes ---------------------------------------------------------
metaplasia.vs.base <- topTreat(tfit,confint = T, coef=1, n=Inf)
metaplasia.vs.pit <- topTreat(tfit, confint = T, coef=2, n=Inf)
base.vs.pit <- topTreat(tfit, confint = T, coef=3, n=Inf)

head(metaplasia.vs.base, 20)
head(metaplasia.vs.pit, 20)
head(base.vs.pit, 20)

write.csv(as_tibble(metaplasia.vs.base, rownames='Gene'), 
          file = paste(prefix, 'metaplasia_vs_base.csv', sep='_'), row.names=FALSE)
write.csv(as_tibble(metaplasia.vs.pit, rownames='Gene'), 
          file= paste(prefix, 'metaplasia_vs_pit.csv', sep='_'), row.names=FALSE)
write.csv(as_tibble(base.vs.pit, rownames='Gene'), 
          file= paste(prefix, 'base_vs_pit.csv', sep ='_'), row.names=FALSE)


# Find genes significantly overexpressed in metaplasia ------------------------
common.degs <- intersect(rownames(metaplasia.vs.base)[metaplasia.vs.base$logFC > 0 &
                                                        metaplasia.vs.base$adj.P.Val <= 0.05],
                         rownames(metaplasia.vs.pit)[metaplasia.vs.pit$logFC > 0 &
                                                       metaplasia.vs.pit$adj.P.Val <= 0.05])


length(common.degs)
common.degs

# Export
metaplasia_degs <- list(metaplasia.vs.base = data.frame(symbol = common.degs, metaplasia.vs.base[common.degs,]),
                        metaplasia.vs.pit = data.frame(symbol = common.degs, metaplasia.vs.pit[common.degs,]),
                        base.vs.pit = data.frame(symbol=common.degs, base.vs.pit[common.degs,]))
write.xlsx(metaplasia_degs, file.path(opt$output_dir, 'metaplasia_overexpressed_genes.xlsx'))

# Heatmap of all sig. overexpressed genes in metaplasia -------------------

# Heatmap annotations with sample metadata
ha_samples = y$samples %>%
  as_tibble() %>%
  mutate(dissection = as.double(as.character(dissection)),
         Lesion = metaplasia,
         `Biopsy site`=region,
         `OLGIM stage`=OLGIM)


ha_samples$Lesion <- recode_factor(ha_samples$Lesion, 'moderate'='Moderate GIM', 'severe'='Marked GIM')
ha_samples$`OLGIM stage` <- addNA(ha_samples$`OLGIM stage`)
levels(ha_samples$`OLGIM stage`) <- c(levels(ha_samples$OLGIM), 'Tumor')
ha_samples$`OLGIM stage` <- recode_factor(ha_samples$`OLGIM stage`, `2`='OLGIM II', `3`='OLGIM III')


column_ha <- HeatmapAnnotation(df = as.data.frame(select(ha_samples, Lesion, atrophy, gastritis, `OLGIM stage`, `Biopsy site`)),
                               col = list(Lesion = c('Moderate GIM'= brewer.pal(9, "Oranges")[5], 
                                                     'Marked GIM'=brewer.pal(9, "Oranges")[6]),
                                          `Biopsy site` = c('antrum' = brewer.pal(8, "Greens")[7], 
                                                     'incisura' = brewer.pal(8, "Greens")[6]),
                                          atrophy = c('neg'= 'green', 
                                                      'mild'=brewer.pal(9, "Oranges")[4], 
                                                      'moderate' = brewer.pal(9, "Oranges")[5]),
                                          gastritis = c('mild_chronic' = brewer.pal(9, "Oranges")[4], 
                                                        'moderate_chronic'=brewer.pal(9, "Oranges")[5]), 
                                          `OLGIM stage` = c('OLGIM II' = brewer.pal(11, "RdYlBu")[9],
                                                            'OLGIM III' = brewer.pal(11, "RdYlBu")[10],
                                                            'Tumor' = brewer.pal(9, 'Reds')[4])
                               ),
                               simple_anno_size = unit(3, 'mm'),
                               annotation_name_gp = gpar(fontsize = 10)
)

counts_lcpm_z <- t(scale(t(counts_lcpm)))

# features <- c('CPS1', 'CDX1', 'HOXA10', 'CLDN3', 'CDH17', 'ANXA13', 'HKDC1', 
#               'DMBT1', 'PRAP1', 'GPD1', 'SLC39A5', 'OLFM4', 'ANPEP', 'SLC13A2',
#               'ONECUT2')

RNAscope_genes <- c('TFF3', 'ANXA13', 'HKDC1', 'DMBT1', 'OLFM4', 'CPS1', 'ANPEP', 'SLC39A5',
                    'ONECUT2', 'CLDN3', 'CDH17', 'CDX1')

setdiff(RNAscope_genes, common.degs)

labels <- RNAscope_genes
labels_indices <- which(rownames(counts_lcpm_z[common.degs,]) %in% labels)

ha_labels = rowAnnotation(foo = anno_mark(at=labels_indices, labels = labels))

hmap_common <- Heatmap(counts_lcpm_z[common.degs,], name='z-score', 
        row_title = 'Genes', row_names_gp = gpar(fontsize=9), show_row_names = F,
        column_title = 'Pseudobulk Samples', column_title_side = 'top', column_names_gp = gpar(fontsize=11, col= col),
        clustering_method_rows = 'ward.D2', clustering_distance_rows = 'manhattan', row_dend_width = unit(1.5, 'cm'),
        clustering_method_columns = 'ward.D2', clustering_distance_columns = 'manhattan',
        cluster_columns= T,
        column_split = 2,
        top_annotation = column_ha)

draw(hmap_common)

# pdf(file = file.path(opt$output_dir, 'Heatmap_metaplasia_overexpressed_genes.pdf'), height=6.25, width=10.14)
# draw(hmap_common)
# dev.off()


# Find and export signature genes significantly overexpressed in metaplasia ---
length(signature_genes)
genes_35 <- intersect(signature_genes, common.degs)
genes_35_lst <- list(metaplasia.vs.base = data.frame(symbol = genes_35,  metaplasia.vs.base[genes_35,]),
                 metaplasia.vs.pit = data.frame(symbol = genes_35, metaplasia.vs.pit[genes_35,]),
                 base.vs.pit = data.frame(symbol = genes_35, base.vs.pit[genes_35,]))
openxlsx::write.xlsx(genes_35_lst, file.path(opt$output_dir, '35_genes.xlsx'))


# Heatmap of 35 genes -----------------------------------------------------

hmap_35 <- Heatmap(counts_lcpm_z[genes_35,], name='z-score', 
                       row_title = 'Genes', row_names_gp = gpar(fontsize=8), show_row_names = T,
                       column_title = 'Pseudobulk Samples', column_title_side = 'top', column_names_gp = gpar(fontsize=11, col= col),
                       clustering_method_rows = 'ward.D2', clustering_distance_rows = 'manhattan', row_dend_width = unit(1.5, 'cm'),
                       clustering_method_columns = 'ward.D2', clustering_distance_columns = 'manhattan',
                       cluster_columns= T,
                       column_split = 2,
                       top_annotation = column_ha)

draw(hmap_35)

# pdf(file = file.path(opt$output_dir, 'Heatmap_35_overexpressed_signature_genes.pdf'), height=6.25, width=10.14)
# draw(hmap_35)
# dev.off()


# Export gene lists -------------------------------------------------------

common.degs
genes_35

write.csv(common.degs, file = file.path(opt$output_dir, 'metaplasia_overexpressed_genes.csv'), row.names = F)
write.csv(genes_35, file = file.path(opt$output_dir, '35_overexpressed_signature_genes.csv'), row.names = F)

# library(gplots)
metaplasia.vs.base.topgenes <- rownames(metaplasia.vs.base)[1:100]
de_gene_pairs <- list(metaplasia.vs.base, metaplasia.vs.pit, base.vs.pit)
# sig_genes <- c()
#for (pair in de_gene_pairs) {
#  print(head(pair))
#  print(pair)
#  print(pair[pair$adj.P.val < 0.05, ])
#  append(sig_genes, rownames(pair[pair$adj.P.val < 0.05]))
#}


# Heatmap -----------------------------------------------------------------



features <- c('CPS1', 'CDX1', 'HOXA10', 'CLDN3', 'CDH17', 'ANXA13', 'HKDC1', 
                           'DMBT1', 'PRAP1', 'GPD1', 'SLC39A5', 'OLFM4', 'ANPEP', 'SLC13A2',
                           'ONECUT2')

# Heatmap of signature genes from bulk

# i <- which(rownames(counts_lcpm) %in% cluster_5)
i <- which(rownames(counts_lcpm %in% signature_genes))

mycol <- colorpanel(1000, 'blue', 'white', 'red')

signature_genes
rownames(counts_lcpm)




pdf(file = file.path(opt$output_dir, 'Heatmap_39_genes.pdf'), height=10.14, width=10.14)
Heatmap(counts_lcpm_z[intersect(signature_genes, common.degs),], name='z-score', 
        row_title = 'Genes', row_names_gp = gpar(fontsize=9), show_row_names = T,
        column_title = 'Pseudobulk Samples', column_title_side = 'top', column_names_gp = gpar(fontsize=11, col= col),
        clustering_method_rows = 'ward.D2', row_dend_width = unit(1.5, 'cm'),
        clustering_method_columns = 'ward.D2',
        cluster_columns= T,
        column_split = 2,
        top_annotation = column_ha)
dev.off()

length(intersect(intersect(signature_genes, rownames(counts_lcpm)), features))
setdiff(features, intersect(signature_genes, rownames(counts_lcpm)))

setdiff(features, common.degs)

sort(common.degs)
# Heatmap of significant genes 



intersect(signature_genes, common.degs)

length(intersect(common.degs, signature_genes))
setdiff(features, intersect(signature_genes, common.degs))

metaplasia.vs.base.sig <- metaplasia.vs.base[metaplasia.vs.base$adj.P.Val < 0.05, ]
metaplasia.vs.pit.sig <- metaplasia.vs.pit[metaplasia.vs.pit$adj.P.Val < 0.05, ]
base.vs.pit.sig <- base.vs.pit[base.vs.pit$adj.P.Val < 0.05, ]
sig.genes <-  c(rownames(metaplasia.vs.base.sig),
                rownames(metaplasia.vs.pit.sig),
                rownames(base.vs.pit.sig)) %>%
  unique()

i <- which(rownames(v) %in% sig.genes)


fc_thresh = 0

metaplasia.vs.base.up <- metaplasia.vs.base %>%
  as_tibble(rownames = 'gene') %>%
  filter(logFC>fc_thresh, adj.P.Val < 0.05)
metaplasia.vs.pit.up <- metaplasia.vs.pit %>%
  as_tibble(rownames = 'gene') %>%
  filter(logFC>fc_thresh, adj.P.Val < 0.05)

# sig.genes <-  c(pull(metaplasia.vs.base.up, gene),
#                 pull(metaplasia.vs.pit.up, gene)) %>%
#   unique()

sig.genes <-  intersect(pull(metaplasia.vs.base.up, gene),
                pull(metaplasia.vs.pit.up, gene)) %>%
  unique()

sig.genes %>% 
  enframe(name = NULL, value = 'metaplasia_overexpressed') %>%
  write_csv(paste(prefix, 'metaplasia_overexpressed_genes.csv', sep='_'))

cluster_5_up <- intersect(sig.genes, cluster_5)

cluster_5_up %>%
  enframe(name = NULL, value = 'metaplasia_overexpressed_C5') %>%
  write_csv(paste(prefix, 'metaplasia_overexpressed_C5_genes.csv', sep='_'))

metaplasia.vs.base.up %>%
  filter(gene %in% cluster_5) %>%
ggplot(mapping = aes(x=AveExpr, y=logFC)) +
  geom_point()

i <- which(rownames(counts_lcpm) %in% sig.genes)

mycol <- colorpanel(1000, 'blue', 'white', 'red')

counts_lcpm_z <- t(scale(t(counts_lcpm)))

labels <- intersect(sig.genes, cluster_5)
labels_indices <- which(rownames(counts_lcpm_z[i,]) %in% labels)

ha_labels = rowAnnotation(foo = anno_mark(at=labels_indices, labels = labels))

Heatmap(counts_lcpm_z[i,], name='z-score', 
        row_title = 'Genes', row_names_gp = gpar(fontsize=9), show_row_names = F,
        column_title = 'Pseudobulk Samples', column_title_side = 'top', column_names_gp = gpar(fontsize=11, col= col),
        clustering_method_rows = 'ward.D2', row_dend_width = unit(1.5, 'cm'),
        clustering_method_columns = 'ward.D2',
        cluster_columns= T,
        column_split = 2,
        top_annotation = column_ha,
        right_annotation = ha_labels)


i <- which(rownames(counts_lcpm) %in% cluster_5_up)

Heatmap(counts_lcpm_z[i,], name='z-score', 
        row_title = 'Genes', row_names_gp = gpar(fontsize=9), show_row_names = T,
        column_title = 'Pseudobulk Samples', column_title_side = 'top', column_names_gp = gpar(fontsize=11, col= col),
        clustering_method_rows = 'ward.D2', row_dend_width = unit(1.5, 'cm'),
        clustering_method_columns = 'ward.D2',
        cluster_columns= T,
        column_split = 2,
        top_annotation = column_ha)

ht <- Heatmap(counts_lcpm_z[i,], name='z-score', 
#        heatmap_width = unit(15, 'cm'), heatmap_height = unit(9, 'cm'),
        row_title = 'Genes', row_names_gp = gpar(fontsize=9),
        column_title = '%s', column_title_side = 'bottom', 
        show_column_names = F, column_names_gp = gpar(fontsize=11, col= col),
        clustering_method_rows = 'ward.D2', row_dend_width = unit(1.5, 'cm'),
        clustering_method_columns = 'ward.D2',
        cluster_columns= T,
        column_split = y$samples$group, column_gap = unit(0.8, 'mm'),
        row_split = 3, row_gap = unit(0.8, 'mm'),
        top_annotation = column_ha)
draw(ht, column_title = 'Pseudobulk samples')


write.fit(tfit, dt, file="pseudobulk_de_limma_lfc1.txt")




