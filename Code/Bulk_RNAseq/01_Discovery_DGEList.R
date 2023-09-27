# 00 - Load required R packages -------------------------------------------
library(optparse)
library(openxlsx)
library(edgeR)
library(tidyverse)

# 01 - Set up working directory -------------------------------------------

setwd("./A03_20233105_Discovery_DGEList/")
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

opt$counts <- "../A02_220428_htSeq_counts/"
opt$external_files <- "../A00_External_files/"
opt$output_dir <- "./temp"
opt$plot_dir <- "./plots"

# 03 - Import raw count files ---------------------------------------------

pdf(file.path(opt$plot_dir, "DGEList_plots.pdf"),
    onefile = F,
    height = 10.41,
    width = 10.41)

files <- list.files(opt$counts, # List all files from this folder
                    pattern = "*.counts$", # Listed must en with ".counts"
                    full.names = T) # Files will be listed with relative path from this folder
raw_data <- readDGE(files, # Read all files listed in object files into a DGEList type object
                    h = F) # The first line is not a header
rm(files)

# 04 - Format sample metadata ---------------------------------------------
# Extract sample names from filenames, store them in name column and add them as column names
# We do this with str_split function from stringr package
raw_data$samples$name <- str_split(str_split(raw_data$samples$files,
                                             "//", simplify = T)[,2],
                                   "\\.", simplify = T)[,1]
colnames(raw_data) <- str_split(str_split(colnames(raw_data),
                                          "//",
                                          simplify = T)[,2],
                                "\\.", simplify = T)[,1]

# Extract LIMS barcode from sample names and store in barcode column
raw_data$samples$barcode <- str_split(raw_data$samples$name,
                                      "_",
                                      simplify = T)[,2]
# Import external metadata, that is matched to this raw data by the LIMS barcode
sample_metadata <- read.xlsx(file.path(opt$external_files, "sample_metadata.xlsx"),
                             detectDates = T)
sample_metadata2 <- subset(read.xlsx(file.path(opt$external_files, "sample_metadata_original.xlsx"),
                                     detectDates = T),
                           select = c(Barcode, OLGA, OLGIM))
sample_metadata <- merge(sample_metadata,
                         sample_metadata2,
                         by = "Barcode")
rm(sample_metadata2)

colnames(sample_metadata)[colnames(sample_metadata)=="Barcode"] <- "barcode"
# Column `37` has no actual data on it, it is just a remainder from the excel spreadsheet.
# Column "GIM" is also a remainder from the excel spreadsheet and will be removed
# We'll remove them.
sample_metadata[,c("37","GIM")] <- NULL

# Merge external metadata with samples using the LIMS barcode
dim(raw_data$samples)
# [1] 95  6

raw_data$samples <- merge(raw_data$samples,
                          sample_metadata,
                          by = "barcode",
                          all = T)
rm(sample_metadata)
# View(raw_data$samples)
# Sample barcoded 22654 was not sequenced, and will be removed from samples dataframe
raw_data$samples <- raw_data$samples[!is.na(raw_data$samples$files),]
raw_data$samples <- data.frame(id = raw_data$samples$name,
                               raw_data$samples)
raw_data$samples$name <- NULL
raw_data$samples$patient <- str_split(raw_data$samples$id,
                                      "_",
                                      simplify = T)[,1]
type <- raw_data$samples$Type
table(type)

# Update the annotation of normal samples on a different column, depending on whether:
# - they are normal controls
# - they are, in fact, tumor adjacent mucosa (TAM)
# - they raise suspicion on their class as they have history of GIM (GIM surveillance)
# - they tested positive for H. pylori (Hp)

raw_data$samples$Type[raw_data$samples$Comments=="Antrum has signet ring carcinoma"] <- "TAM"
raw_data$samples$Type[raw_data$samples$Type=="GIM"&raw_data$samples$Comments=="Antrum has GIM and signet ring carcinoma"] <- "TAM"
raw_data$samples$Type[raw_data$samples$Comments=="(No) H pylori seen?"] <- "Hp"
# Update the annotation of the tumor samples, which are either recurrent or diffuse (signet ring)
raw_data$samples$Type[raw_data$samples$Comments=="Recurrent tumor"] <- "ReccurrentGC"
raw_data$samples$Comments[raw_data$samples$Type=="Adenocarcinoma"]
# [1] "Antrum has GIM and signet ring carcinoma"
raw_data$samples$Type[raw_data$samples$Type=="Adenocarcinoma"] <- "DiffuseGC"
raw_data$samples$Type[raw_data$samples$Type=="Signet ring"] <- "DiffuseGC"
table(raw_data$samples$Type, type)
table(raw_data$samples$Type)
with(raw_data$samples,
     table(OLGIM, Type))

# Match column order with that of sample metadata
table(raw_data$samples$id%in%colnames(raw_data))
raw_data$counts <- raw_data$counts[,raw_data$samples$id]
identical(raw_data$samples$id,colnames(raw_data))
with(raw_data$samples, table(OLGIM, Type))
# Let's overwrite the filter to keep only normal and GIM samples
filt <- raw_data$samples$Type%in%c("Normal", "GIM")
table(filt)
# FALSE  TRUE 
#     6    89 
raw_data <- raw_data[,filt]
dim(raw_data)
# Check the remaining columns from the metadata
with(raw_data$samples, table(Location))
# ANT DRY will be re-assigned as ANT, and L1 will be removed from the dataset
raw_data$samples$Location[raw_data$samples$Location=="ANT DRY"] <- "ANT"
filt <- raw_data$samples$Location!="L1"
table(filt)
raw_data <- raw_data[,filt]
table(raw_data$samples$Location)
# Check for paired samples
paired <- as.data.frame(table(raw_data$samples$patient))
paired <- paired$Var1[paired$Freq==2]
raw_data$samples$paired <- FALSE
raw_data$samples$paired[raw_data$samples$patient%in%paired] <- TRUE
with(raw_data$samples, table(paired))
with(raw_data$samples, table(paired, Location))
with(raw_data$samples, table(Location, Severity))
raw_data$samples$Severity[is.na(raw_data$samples$Severity)] <- "Normal"

# Format OLGA
raw_data$samples$OLGA <- ifelse(raw_data$samples$OLGA==0, "OLGA 0",
                                 ifelse(raw_data$samples$OLGA==1, "OLGA I",
                                        ifelse(raw_data$samples$OLGA==2, "OLGA II",
                                               ifelse(raw_data$samples$OLGA==3, "OLGA III", "OLGA IV"))))

raw_data$samples$OLGA <- factor(raw_data$samples$OLGA,
                                 levels = c("OLGA 0", "OLGA I", "OLGA II", "OLGA III", "OLGA IV"),
                                 ordered = T)

# Format OLGIM
raw_data$samples$OLGIM <- ifelse(raw_data$samples$OLGIM==0, "OLGIM 0",
                                 ifelse(raw_data$samples$OLGIM==1, "OLGIM I",
                                        ifelse(raw_data$samples$OLGIM==2, "OLGIM II",
                                               ifelse(raw_data$samples$OLGIM==3, "OLGIM III", "OLGIM IV"))))

raw_data$samples$OLGIM <- factor(raw_data$samples$OLGIM,
                                 levels = c("OLGIM 0", "OLGIM I", "OLGIM II", "OLGIM III", "OLGIM IV"),
                                 ordered = T)

# Clean the workspace
rm(filt, paired, type)

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
# [1] 38643    89

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
# [1] 32055    88
genes <- raw_data$genes
write_delim(genes, file.path(opt$output_dir, "Gene_annotations.tsv"))
raw_data$genes <- subset(raw_data$genes,
                         select = c(SYMBOL, CHR))

saveRDS(raw_data, file = "Raw_data_DGEList.rds")
saveRDS(clean_data, file = "Clean_data_DGEList.rds")

dev.off()

rm(list = ls())



setwd("../")
