# 00 - Load required R packages -------------------------------------------
library(optparse)
library(openxlsx)
library(edgeR)
library(tidyverse)

# 01 - Set up working directory -------------------------------------------

setwd("./B02_20230607_Validation_DGEList")
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
                     c("-c", "--counts"),
                     default=NA,
                     type="character", 
                     help="Path to count files (*.counts$)")
parser <- add_option(parser,
                     c("-e", "--external_files"),
                     default=NA,
                     type="character",
                     help="gtf file used for gene annotations")
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

opt$counts <- "../B01_221201_htSeq_tpm"
opt$external_files <- "../A00_External_files/"
opt$output_dir <- "./temp"
opt$plot_dir <- "./plots"

# 03 - Import raw count files ---------------------------------------------

pdf(file.path(opt$plot_dir,
              "Validation_DGEList_plots.pdf"),
    height = 10.41,
    width = 10.41)

files <- list.files(opt$counts, # List all files from this folder
                    pattern = "*.counts.txt$", # Listed must en with ".counts"
                    full.names = T) # Files will be listed with relative path from this folder
raw_data <- readDGE(files, # Read all files listed in object files into a DGEList type object
                    h = F) # The first line is not a header
rm(files)

# 04 - Format sample metadata ---------------------------------------------
# Extract sample names from filenames, store them in name column and add them as column names
# We do this with str_split function from stringr package
raw_data$samples$name <- str_split(str_split(raw_data$samples$files,
                                             "/", simplify = T)[,3],
                                   "\\.", simplify = T)[,1]
colnames(raw_data) <- raw_data$samples$name

# Extract LIMS barcode from sample names and store in barcode column
raw_data$samples$barcode <- str_split(raw_data$samples$name,
                                      "_",
                                      simplify = T)[,2]
# Import external metadata, that is matched to this raw data by the LIMS barcode
batch2_metadata <- read.xlsx(file.path(opt$external_files, "20220830_GAPS_B2_info.xlsx"),
                             detectDates = T)
batch2_metadata[,8:10] <- NULL
batch2_metadata$Batch <- "B2"

batch3_metadata <- read.xlsx(file.path(opt$external_files, "20220901_GAPS_B3_info.xlsx"),
                             detectDates = T)
batch3_metadata[,8:10] <- NULL
batch3_metadata$Batch <- "B3"

batch4_metadata <- read.xlsx(file.path(opt$external_files, "20220922_GAPS_B4_info.xlsx"),
                             detectDates = T)
batch4_metadata[,8:10] <- NULL
batch4_metadata$Batch <- "B4"

sample_metadata <- rbind.data.frame(batch2_metadata,
                                    batch3_metadata,
                                    batch4_metadata)

rm(batch2_metadata, batch3_metadata, batch4_metadata)

# All LIMS barcodes have an .R01 suffix that will be removed to match with the barcode extracted from the count files
if (all(str_split(sample_metadata$Barcode,
                  "\\.",
                  simplify = T)[,2]=="R01")) {
  sample_metadata$Barcode = str_split(sample_metadata$Barcode,
                                      "\\.",
                                      simplify = T)[,1]
}

aux <- data.frame(Sample = raw_data$samples$name,
                  Barcode = raw_data$samples$barcode)

sample_metadata <- merge(aux,
                         sample_metadata,
                         by = "Barcode"); rm(aux)
table(is.na(sample_metadata$OLGIM))
# 
# FALSE  TRUE 
#   229    59

rownames(sample_metadata) <- sample_metadata$Sample

if(all(rownames(sample_metadata)%in%colnames(raw_data))&all(colnames(raw_data)%in%rownames(sample_metadata))) {
  sample_metadata <- sample_metadata[colnames(raw_data),]
}

if(identical(colnames(raw_data), rownames(sample_metadata))){
  raw_data$samples <- merge(raw_data$samples,
                            sample_metadata,
                            by.x = "name",
                            by.y = "Sample",
                            all = T)
}

identical(raw_data$samples$name, colnames(raw_data))
dim(raw_data)
# [1] 54633   288

raw_data <- raw_data[,!is.na(raw_data$samples$OLGIM)]

dim(raw_data)
# [1] 54633   229

raw_data$samples$group <- factor(ifelse(raw_data$samples$OLGIM%in%c(0,1,2),
                                        "low_risk",
                                        "high_risk"),
                                 levels = c("low_risk",
                                            "high_risk"),
                                 ordered = T)

raw_data$samples$OLGIM <- factor(ifelse(raw_data$samples$OLGIM == 0,
                                        "OLGIM 0", ifelse(
                                          raw_data$samples$OLGIM == 1,
                                          "OLGIM I", ifelse(
                                            raw_data$samples$OLGIM == 2,
                                            "OLGIM II", ifelse(
                                              raw_data$samples$OLGIM == 3,
                                              "OLGIM III", "OLGIM IV")
                                          ))),
                                 levels = c("OLGIM 0",
                                            "OLGIM I",
                                            "OLGIM II",
                                            "OLGIM III",
                                            "OLGIM IV"),
                                 ordered = T)

raw_data$samples$OLGIM

raw_data$samples <- data.frame(id = str_split(raw_data$samples$name, "_", simplify = T)[,1],
                               raw_data$samples)

raw_data$samples$id <- paste(raw_data$samples$id,
                             tolower(raw_data$samples$Location),
                             sep = "_")

# Inspect the type of lesions
table(raw_data$samples$Type)
#   Adenocarcinoma Dysplasia/Cancer              GIM           Normal 
#                1                3               62              163 
# For differential expression analysis, remove patients with dysplasia or cancer from the analysis

cancer_ids <- str_split(raw_data$samples$name,
                        "_",
                        simplify = T)[raw_data$samples$Type%in%c("Adenocarcinoma",
                                                                 "Dysplasia/Cancer"),
                                      1]
table(!str_split(colnames(raw_data), "_", simplify = T)[,1]%in%cancer_ids)
# FALSE  TRUE 
#    11   218

dim(raw_data)
# [1] 54633   229

all(substring(colnames(raw_data), 1, 6)==substring(raw_data$samples$id, 1, 6))
# [1] TRUE

raw_data <- raw_data[,!str_split(colnames(raw_data),
                                 "_",
                                 simplify = T)[,1]%in%cancer_ids]
dim(raw_data)
# [1] 54633   218

dups <- raw_data$samples$id[duplicated(raw_data$samples$id)]
dups
# [1] "P07115_ant" "P07115_bod" "P09417_bod"
# These patients have more than 1 biopsy per anatomic site.
# We will keep the worse diagnosis for OLGIM cases

table(raw_data$samples$id[raw_data$samples$id%in%dups],
      raw_data$samples$Type[raw_data$samples$id%in%dups])

raw_data <- raw_data[,!duplicated(raw_data$samples$id)]
dim(raw_data)
# [1] 54633   215

# Plot mosaic plot to explore composition of cases (high-risk) and controls (low-risk) between the batches
with(raw_data$samples,
     mosaicplot(table(Batch,
                      ifelse(OLGIM%in%c("OLGIM III", "OLGIM IV"),
                             "High-risk", "Low-risk")),
                color = T,
                main = "Mosaic plot:\nStrata by batch",
                las = 2))

with(raw_data$samples,
     mosaicplot(table(Batch,
                      OLGIM),
                color = T,
                main = "Mosaic plot:\nStrata by batch",
                las = 2))

# Examine distribution of cases (high-risk) between batches
table(raw_data$samples$Batch, ifelse(raw_data$samples$OLGIM%in%c(0,1,2), "low_risk", "high-risk"))
chisq.test(table(raw_data$samples$Batch, raw_data$samples$OLGIM%in%c(0,1,2)))
# 	Chi-squared test for given probabilities
# data:  table(raw_data$samples$Batch, raw_data$samples$OLGIM %in% c(0,     1, 2))
# X-squared = 4.1674, df = 2, p-value = 0.1245

# There's no significant association between the number of cases and sequencing batch

# 05 - Update gene annotations --------------------------------------------
# Load the original annotation file as a data.frame
gtf <- data.frame(rtracklayer::import(file.path(opt$external_files,
                                                "hg38.ncbiRefSeq.gtf")))
# Filter duplicated gene IDs, we are not interested in the genomic coordinates used for annotations of different transcripts
colnames(gtf)
gtf <- subset(gtf,
              subset = !duplicated(gtf$gene_id),
              select = c(gene_id,
                         seqnames))
# We'll only keep sequences mapped to chromosomes 1:23, X, Y and mitochondrial genes
chr <- paste0("chr", c(1:23, "X", "Y", "M"))
gtf <- gtf[gtf$seqnames%in%chr,]
# Since we have no stable identifiers, we will create unique IDs for each gene by concatenating the symbol or alias (we don't know yet) using the chromosome info
gtf$common_id <- with(gtf,
                      paste(gene_id,
                            seqnames,
                            sep = "_"))

# We'll next upload a list of gene ids to update the annotations in genenames.org
# Upload gene list to multi-symbol checker at https://www.genenames.org/tools/multi-symbol-checker/
genes <- gtf$gene_id
write.table(genes,
            file.path(opt$output_dir,
                      "Original_gtf_genes_with_Chr.txt"),
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F)

# Load hgnc multi-symbol checker output file
genes <- read.csv(file.path(opt$external_files, "hgnc-symbol-check.csv"),
                  stringsAsFactors = F,
                  skip = 1)
# Filter any unmatched or withdrawn genes
table(genes$Match.type)
#    Alias symbol Approved symbol Entry withdrawn Previous symbol       Unmatched 
#             551           38572              10             292           11970 
genes <- genes[!genes$Match.type%in%c("Entry withdrawn", "Unmatched"),]
sum(genes$Match.type!="Approved symbol")
# [1] 843
# Create common id to match with the gtf file using the original gene symbol
genes$chr <- str_split(genes$Location,
                       "\\.",
                       simplify = T)[,1]
genes$chr <- gsub("p",
                  "q",
                  genes$chr)
genes$chr <- str_split(genes$chr,
                       "q",
                       simplify = T)[,1]
# Check that all chromosomes are correctly annotated
table(genes$chr)
genes[genes$chr=="",]
# LEFTY3P does not have an associated chromosome from genenames, but matched LEFTY3P from chr1 in our annotation file
genes$chr[genes$chr==""] <- "1"
# Check for duplicates
# Identify chromosome location
genes$chr <- gsub("14 un",
                  "14",
                  genes$chr)
genes$chr <- gsub("mitochondria",
                  "M",
                  genes$chr)
genes$chr <- paste0("chr",
                    genes$chr)
table(genes$chr)
genes$common_id <- paste(genes$Input,
                         genes$chr,
                         sep = "_")

table(genes$common_id%in%gtf$common_id)
# FALSE  TRUE 
#   658 38757 
table(unique(gtf$common_id)%in%genes$common_id)
# FALSE  TRUE 
# 12025 38652 

genes2 <- merge(genes,
                gtf,
                by = "common_id")
length(setdiff(genes$common_id, genes2$common_id))
# 655 IDs were discarded because they matched gene symbols located at different chromosomes than those from the original annotation file

genes <- genes2; rm(genes2); rm(gtf)

# Check for duplicates
table(duplicated(genes$Input))
# FALSE  TRUE 
# 38652   105 

updated <- genes[genes$Match.type!="Approved symbol",]
genes <- genes[genes$Match.type=="Approved symbol",]
updated <- updated[!updated$Input%in%genes$Input,]

if (all(!duplicated(updated$Input))&all(!duplicated(genes$Input))) {
  genes <- rbind.data.frame(genes, updated)
  rm(updated)
}

raw_data <- raw_data[genes$Input,]
identical(rownames(raw_data), genes$Input)
colnames(genes)
genes[,c("common_id", "HGNC.ID", "gene_id", "seqnames")] <- NULL

# Identify any gene that was duplicated after updating symbol
table(duplicated(genes$Approved.symbol))
genes$Approved.symbol[duplicated(genes$Approved.symbol)]
# [1] "IDSP1"       "SEPTIN14P20" "SEPTIN14P6"  "SEPTIN14P23" "SUCLA2"      "MALAT1"     
# [7] "NABP1"       "RPL36AL"     "MYO18A"

# Summarize the expression levels of duplicates by adding them up under their approved symbol
raw_data$counts <- aggregate(raw_data$counts,
                             list(genes$Approved.symbol),
                             FUN = "sum")
dim(raw_data$counts)
# [1] 38643   216

# Create the final gene annotation sheet
summarized_from <- aggregate(genes$Input,
                             list(genes$Approved.symbol),
                             FUN = "paste", sep = ", ")
summarized_from$x <- as.character(summarized_from$x)
genes <- subset(genes,
                subset = !duplicated(Approved.symbol),
                select = c(Approved.symbol,
                           Approved.name,
                           Location,
                           chr))
genes <- merge(genes,
               summarized_from,
               by.x = "Approved.symbol",
               by.y = "Group.1")
colnames(genes) <- c("SYMBOL",
                     "GENENAME",
                     "LOCATION",
                     "CHR",
                     "PREVIOUS_SYMBOL")

# Remove all "c()" from the PREVIOUS_SYMBOL column
genes$PREVIOUS_SYMBOL <- gsub("[^A-Z , 0-9]", "", genes$PREVIOUS_SYMBOL)

if(all(genes$SYMBOL %in% raw_data$counts$Group.1)&all(raw_data$counts$Group.1 %in% genes$SYMBOL)) {
  cat("Updating gene annotations in DGEList object.\n")
  rownames(raw_data$counts) <- raw_data$counts$Group.1
  raw_data$counts$Group.1 <- NULL
  raw_data$counts <- data.matrix(raw_data$counts)
  raw_data$counts <- raw_data$counts[genes$SYMBOL,]
}

# 06 - Verify that DGEList object is correctly formatted ------------------

# Verify that the sample metadata and gene expression matrices match
if(identical(colnames(raw_data), raw_data$samples$id)&
   identical(genes$SYMBOL, rownames(raw_data))) {
  cat("DGEList is correctly formatted.\n")
}

clean_data <- DGEList(counts = raw_data$counts,
                      samples = raw_data$samples,
                      genes = genes,
                      group = as.numeric(raw_data$samples$OLGIM),
                      remove.zeros = F)

raw_data <- DGEList(counts = raw_data$counts,
                    samples = raw_data$samples,
                    genes = genes,
                    group = as.numeric(raw_data$samples$OLGIM),
                    remove.zeros = T)

raw_data$samples$group <- raw_data$samples$OLGIM
raw_data$samples$lib.size.1 <- NULL
raw_data$samples$norm.factors.1 <- NULL
raw_data$samples$group.1 <- NULL
dim(raw_data)
# [1] 32994   215
genes <- raw_data$genes
write_delim(genes, file.path(opt$output_dir, "Gene_annotations.tsv"))
raw_data$genes <- subset(raw_data$genes,
                         select = c(SYMBOL, CHR))

saveRDS(raw_data, file = "Raw_data_DGEList.rds")
saveRDS(clean_data, file = "Clean_data_DGEList.rds")
rm(list = ls())

dev.off()

setwd("../")
