# Set environment variables to force save and saveRDS to use format readable by R<3.5
Sys.setenv(R_DEFAULT_SAVE_VERSION = 2)
Sys.setenv(R_DEFAULT_SERIALIZE_VERSION = 2)

R_SAVE_VERSION=2  #Not sure if this is needed?

# Set cellranger reference version# (from environment variable), default to 3.0.0 if not provided
# Use ref version 1.2.1 for cellranger 1 or 2; ref version 3.0.0 for cellranger 3
cr_ref_version = Sys.getenv('CR_VERSION')
if ( cr_ref_version < '1' ) { cr_ref_version = '3.0.0' }

cr_genome = Sys.getenv('GENOME')
if ( cr_genome < 'a' ) { cr_genome = 'GRCh38' }

##General methods for working with scRNA Seurat objects
GENE_XREF = paste0('/mnt/ix1/Resources/10X_Resources/cellranger_rna/refdata-cellranger-', cr_ref_version, '/', cr_genome, '/genes/genes_xref.txt')
GENE_POS_XREF = paste0('/mnt/ix1/Resources/10X_Resources/cellranger_rna/refdata-cellranger-', cr_ref_version, '/', cr_genome, '/genes/gene_pos_xref.tsv')

MIN_CPG = 3
MIN_GPC = 200

filesep = .Platform$file.sep

# Create biomaRt gene reference dataset
biomartGeneXref <- function(genes, gene_fld='ensembl_gene_id', organism='hsapiens', build=38, version=NULL) {
  library(biomaRt)
  ensembl = "ENSEMBL_MART_ENSEMBL"
  org_genes = paste0(organism,'_gene_ensembl')
    
  #Host names (in order of preference) below  
  build37_host = c('feb2014.archive', 'grch37')
  build38_host = c('apr2018.archive', 'dec2017.archive', 'dec2016.archive', 'grch38')
  
  #Build list of available versions and associated host names
  archives <- as.data.frame(listEnsemblArchives())
  
  if (!is.null(version)) {
    url = archives[archives$version==paste('Ensembl',version),'url']
    gene_mart = useMart(ensembl, dataset=org_genes, host=substr(url,8,100))
    
  } else if (build==38) {
    gene_mart = useMart(ensembl, dataset=org_genes, host=paste0(build38_host[1],'.ensembl.org'))
  
  } else if (build==37) {   
    gene_mart = useMart(ensembl, dataset=org_genes, host=paste0(build37_host[1],'.ensembl.org'))

  } else {
    print(paste("Build",build,"not supported. Only supporting build 37,38."))
    return()
  }
  ds = useDataset(org_genes, mart=gene_mart)
  res <- getBM(attributes = c("ensembl_gene_id", "entrezgene", "hgnc_id", "hgnc_symbol", "chromosome_name", 
                              "start_position", "end_position", "description"),
               filter=gene_fld, values=genes, mart=ds)
  return(res) 
}

##Determine format of gene name rows
# GeneIDs in the seurat object will be defined in one of four ways:
#   - Ensembl gene name only
#   - Hugo gene name only
#   - Ensembl + Hugo names in format ENSxxxxxxx::Hugo (full Ensembl and Hugo gene names)
#   - Hugo name, with last 3 digits of Ensembl gene in format:  Hugo::nnn
chkGeneNameFormat <- function(genesV) {
  rows_with_delim = grep("::", genesV)
  gene1_first3 = substr(genesV[1],1,3)

  if (length(rows_with_delim) > 0) {
    fmt <- ifelse(gene1_first3 == 'ENS', 'both', 'hybrid') 
  } else {
    fmt <- ifelse(gene1_first3 == 'ENS', 'ensembl', 'hugo')
  }
  return (fmt)
}

# gene_xref file is used for all gene mapping functions 
#   col1 = Ensembl ID, col2 = Hugo ID, col3 = Hugo unique ID, col4 = HybridID 
#   hybridIDs will be hugo gene, last 3 digits of Ensembl ID, with '::' separator

##Get a unique Hugo ID for each Ensembl ID
mapEnsembl2Hugo <- function(ensemblIDs, gene_xref=GENE_XREF){
  gmap = read.table(gene_xref, sep="\t", header=F, check.names = F, stringsAsFactors = F)
  gidx = match(ensemblIDs, gmap$V1);  
  return (gmap$V3[gidx])
}

mapHugo2Ensembl <- function(hugoIDs, gene_xref=GENE_XREF) {
  gmap = read.table(gene_xref, sep="\t", header=F, check.names = F, stringsAsFactors = F)
  gidx = match(hugoIDs,gmap$V3);  
  return (gmap$V1[gidx])
}

mapEntrez2Hugo <- function(entrezIDs, gene_xref=GENE_XREF) {
  gmap <- read.table(gene_xref, sep="\t", header=F, check.names=F, stringsAsFactors=F)
  gidx <- match(entrezIDs, gmap$V5)
  return (gmap$V3[gidx])
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Read and write 10X sparse matrix trios
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
read10Xtrio <- function(dirname='.', prefix10x, geneID='hugo') {
  prefix10x <- ifelse(prefix10x=="", prefix10x, paste0(prefix10x,"."))
  fn_prefix = paste0(dirname, filesep, prefix10x)
  fn_suffix = ifelse(file.exists(paste0(fn_prefix, 'matrix.mtx.gz')), '.gz', '')

  barcode.path <- paste0(fn_prefix, "barcodes.tsv", fn_suffix)
  matrix.path <- paste0(fn_prefix, "matrix.mtx", fn_suffix)
  expr_matrix <- readMM(matrix.path)
  
  features_tsv <- paste0(fn_prefix, "features.tsv", fn_suffix)
  genes_tsv    <- paste0(fn_prefix, "genes.tsv", fn_suffix)
  features.path <- ifelse(file.exists(features_tsv), features_tsv, genes_tsv)
  feature.names = read.delim(features.path, sep="\t", header=F, stringsAsFactors=F)
  rownames(expr_matrix)= parseGeneNames(feature.names)
  
  barcode.names = read.delim(barcode.path, sep="\t", header=F, stringsAsFactors=F) 
  colnames(expr_matrix) = barcode.names$V1
  
  return(expr_matrix)
}

parseGeneNames <- function(feature.names, geneID='hugo') {
  if (geneID == 'ensembl') {
    gene_names = feature.names$V1
  } else if (geneID == 'both') {
    gene_names = paste(feature.names$V1, feature.names$V2, sep="::")
  } else {
    gene_names = mapEnsembl2Hugo(feature.names$V1, GENE_XREF)
	#gene_names = make.unique(as.vector(genes$V2))
  }
  return(gene_names)
} 

modifyBarcodeSuffix <- function(barcodes, cell.suffix) {
  bc_has_suffix = (length(grep("-", barcodes)[1]) > 0)
  
  if (cell.suffix == 'N' && bc_has_suffix) {
    new_barcodes = gsub("-\\d", '', barcodes)  #Delete barcode suffix
	
  } else if (cell.suffix %in% c('K','N')) {
    new_barcodes = barcodes
	
  } else if (bc_has_suffix) {
    new_bc = paste0('-',cell.suffix)
    new_barcodes = gsub("-\\d", new_bc, barcodes)  #Replace barcode suffix
	
  } else {
    new_barcodes = paste(barcodes,c('-',suffix), sep="")
  }
  return(new_barcodes)
}  
  
addSampleId <- function(seu_) {
  barcode_split <- sapply(colnames(seu_), strsplit,"-")
  barcode_ids   <- as.matrix(sapply(barcode_split, function(x) x[2]))
  dimnames(barcode_ids) <- list(colnames(seu_), c('sample.id'))
  seu_$sample.id <- barcode_ids
  return(seu_)
}

write10Xtrio <- function(seu.data, prefix, overwrite=T){
  if (!overwrite && file.exists(paste0(prefix,".matrix.mtx"))) {
    print(paste0(prefix, "_matrix.mtx already exists and overwrite set to FALSE."))
    return()
  }
  
  if (length(grep("::", rownames(seu.data)[1]) == 0)) {
    genes = rownames(seu.data)
  } else {
    genes = sapply(rownames(seu.data), strsplit,"::")
    genes = t(as.data.frame(genes));
  }
  
  writeMM(as(seu.data,"sparseMatrix"), file=paste0(prefix,".matrix.mtx"), col.names=F )
  write.table(genes, file=paste0(prefix,".genes.tsv"), quote=F, row.names=F, col.names=F, sep="\t" )
  write.table(colnames(seu.data), file=paste0(prefix, ".barcodes.tsv"), quote=F, row.names=F, col.names=F )
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Create or return Seurat objects
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
##Create Seurat object
create10X_Seurat <- function(expr_mtx, min.cells=MIN_CPG, min.genes=MIN_GPC, max.mito=0.1, organism='Homo sapiens', qc=TRUE, projectName) {
  seu_     <- CreateSeuratObject(counts=expr_mtx, min.cells=min.cells, min.features=min.genes, project=projectName)
  seu_     <- addSampleId(seu_)
  if (qc) {
    seu_   <- annotateMito(seu_, organism=organism)
	  seu_   <- annotateCellCycle(seu_, organism=organism)
    }
  return(seu_)
}

##Read Rdata/Robj file and return contents (allows object to be given new name)
loadRObj <- function(robj) {
  load(robj)
  get(ls()[ls() != "robj"])
}

##Find name of seurat object in current environment, and return contents
getSeuratObj <- function(obj_list) {
  print(obj_list)
  seu_ = NULL
  
  if ( 'seurat_obj' %in% obj_list ) {
    print('Have Seurat object: seurat_obj')
    seu_ = seurat_obj
	
  } else if ( 'seurat_tsne' %in% obj_list ) {
    print('Have Seurat object: seurat_tsne')
    seu_ = seurat_tsne
	
  } else {
    for ( obj in obj_list ) {
      temp_obj = get(obj)
      if ( class(temp_obj)[1] == 'seurat' ) {
        print(paste0('Have Seurat object:',obj))
        seu_ <- temp_obj
	break
      }
    }
  }
  return(seu_)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Cell annotation
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
annotateCellCycle <- function(seu_, organism='Homo sapiens') {
  # Human genes available within Seurat as cc.genes, override below for mouse
  if (organism == 'Mus musculus') {
    cc.genes <- readRDS('/mnt/ix1/Resources/scRNA_Ref/CellCycle/mouse_cell_cycle_genes.rds')
  }	
  seu_ <- CellCycleScoring(object=seu_, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=F)
  seu_$CC.Difference <- seu_$S.Score - seu_$G2M.Score
  return(seu_)
}

##Annotate with percent mitochondrial gene/UMI counts
annotateMito <- function(seu_, subsetMito=F, max.mito=0.1, geneID='hugo', organism='Homo sapiens') {
  if (organism == 'Mus musculus') {
    mito_prefix = 'mt-'
  } else {
    mito_prefix = 'MT-'
  }	
  if (geneID == 'both') {
    mito_regex = paste0("^ENS.{12}::", mito_prefix)
  } else {
    mito_regex = paste0("^", mito_prefix)
  }
  	
  seu_[["percent.mito"]] <- PercentageFeatureSet(seu_, pattern=mito_regex)
  if (subsetMito) {
    seu_ <- subset(seu_, subset = percent.mito < max.mito)
  }
  return(seu_)
}  

##True/False annotation using a file with list of barcodes
annotateBarcodes <- function(seu_, barcodes_tsv, meta_name, numeric=F) {
  barcodes = read.table(barcodes_tsv)
  seu_ <- annotateBarcodeList(seu_, barcodes$V1, meta_name, numeric=numeric)
  return(seu_)
}

annotateBarcodeList <- function(seu_, barcodes, meta_name, numeric=F) {
  meta_df = colnames(seu_) %in% barcodes  #Values will be TRUE/FALSE
  if (numeric) {
    meta_df = sapply(meta_df, as.numeric) #Translate TRUE/FALSE to 1/0
  }
  meta_matrix = as.matrix(meta_df)
  dimnames(meta_matrix) <- list(colnames(seu_), c(meta_name))
  seu_[[meta_name]] <- meta_matrix
  return(seu_)
}

##Annotate with cellity QC score (0=low quality, 1=high quality)
annotateCellity <- function(seu_, hq_barcodes_tsv) {
  return (annotateBarcodes(seu_, hq_barcodes_tsv, 'cellity.qc', numeric=T))
}

##Annotate from inferCNV type file, with barcodes and CNV group number
annotateCNV <- function(seu_, barcodes_tsv, meta_name) {
  annot_barcodes = read.table(barcodes_tsv, sep=' ')
  seu_barcodes = colnames(seu_)
  
  df <- seu_barcodes %in% annot_barcodes$V1
  df <- sapply(df, as.numeric)
  idx = match(annot_barcodes$V1, seu_barcodes)
  df[idx] = annot_barcodes$V2
  
  mtx = as.matrix(df)
  dimnames(mtx) = list(seu_barcodes, meta_name)
  seu_[[meta_name]] <- mtx
  return (seu_)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Seurat processing
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
##Normalization and scaling
transform_Seurat <- function(seu_, scale.factor=10000, vars.to.regress=c("nCount_RNA", "percent.mito")) {							  
  seu_ <- NormalizeData(object=seu_, normalization.method="LogNormalize", scale.factor=scale.factor)
  seu_ <- FindVariableFeatures(object=seu_, selection.method = "vst", nfeatures = 2000)
  
  all.genes <- rownames(seu_)
  seu_ <- ScaleData(object=seu_, vars.to.regress=vars.to.regress, features=all.genes) 
  
  return(seu_)
}

##Alternate/new transformation method: sctransform, use directly from Seurat:
#seu_ <- SCTransform(seu_, vars.to.regress=c("nCount_RNA", "percent.mito"), verbose = FALSE)
  
  
##Dimensionality reduction using PCA
pca_Seurat <- function(seu_, pca.pct=90) {
  seu_ <- RunPCA(seu_, features=VariableFeatures(object=seu_))
  
  print(seu_[["pca"]], dims=1:5, nfeatures=5)
  k_dims <- pcaPctVar(seu_, pca.pct=pca.pct)
  print(paste(k_dims, "PCs explain", pca.pct, "% of variance"))
  
  return(seu_)
}

pcaPctVar <- function(seu_, pca.pct=90) {
  sdev = seu_@reductions$pca@stdev
  sdev = 100*sdev/sum(sdev)
  k=1
  while(sum(sdev[1:k]) < pca.pct) {
    k=k+1
  }
  return(k)
}  

##Seurat clustering and tSNE
tsne_Seurat <- function(seu_, do.fast=TRUE, nr_pca=20, minpct_varexplained=90, check_duplicates=F, 
                        perplexity=30, clustering=T, cluster_res=0.6) {

  if (minpct_varexplained & minpct_varexplained > 0) {
    k_dims <- pcaPctVar(seu_, pca.pct=minpct_varexplained)
    print(paste(k_dims, "PCs explain", minpct_varexplained, "% of variance"))
  } 
  if (nr_pca & nr_pca > 0) { nr_dims = nr_pca } else { nr_dims = k_dims }  
  
  if (clustering) {
    print(paste("Finding clusters using", nr_dims, "dimensions..."))
    #defaults for FindClusters: k.param=30, k.scale=25 (lower k (nearest neighbor distance) => more clusters)
    seu_ <- FindClusters(object=seu_, reduction.type="pca", dims.use=1:nr_dims, k.param=30, 
	          resolution=cluster_res, print.output=0, save.SNN=TRUE)
  }

  seu_ <- RunTSNE(object=seu_, dims.use=1:nr_dims, do.fast=do.fast, check_duplicates=check_duplicates, perplexity=perplexity) 
  return(seu_)
}

##Seurat clustering and UMAP
umap_Seurat <- function(seu_, dims.reduce, cluster.res, do.reorder=T) {
  seu_ <- FindNeighbors(seu_, dims=1:dims.reduce)
  seu_ <- FindClusters(seu_, resolution=cluster.res)
  if (do.reorder) {
    seu_ <- BuildClusterTree(seu_, reorder=T, reorder.numeric=T, dims=1:dims.reduce)
  }	
  seu_ <- RunUMAP(seu_, dims=1:dims.reduce, do.label=T)
  return(seu_)
}  

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Cluster markers and candidate genes for ligand/receptor identification
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
getClusterMarkers <- function(seu_, cluster_info, dir='.', only.pos=T, min.pct=0.25, csv=T) {
  markers <- FindMarkers(seu_, ident.1=cluster_info$cident, only.pos=only.pos, min.pct=min.pct)
  markers$gene = rownames(markers)
  rownames(markers) = c()
  if (csv) { write.csv(markers, file=paste0(dir,'/markers_',cluster_info$cname,'.csv'), row.names=F) }
  markers$cluster = cluster_info$cname
  markers$cident  = cluster_info$cident
  return(markers)
}

matchCandidateGenes <- function(cluster_markers, candidate_genes) {
  markers_found <- cluster_markers[rownames(cluster_markers) %in% candidate_genes,]
  return(markers_found)
}

matchLigandReceptors <- function(ligand_markers, receptor_markers, ref_pairs, dir='.') {
  ligand_names = unique(ligand_markers$gene)
  receptor_names = unique(receptor_markers$gene)

  ligand_receptors = ref_pairs[ref_pairs$Ligand.ApprovedSymbol %in% ligand_names & ref_pairs$Receptor.ApprovedSymbol %in% receptor_names,]
  write.csv(ligand_receptors, file=paste0(dir, '/', "ligand_", as.character(ligand_markers$cluster[1]),"_receptor_", as.character(receptor_markers$cluster[1]), ".csv"))  
  return(ligand_receptors)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# inferCNV drilldown 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
addGenePositions <- function(df, gene_pos_xref=GENE_POS_XREF, geneID='hugo') {
  gene_coords = read.table(gene_pos_xref, sep='\t', stringsAsFactors=F)
  if (geneID == 'ensembl') {
    idx = match(rownames(df), gene_coords$V1)
  } else {
    idx = match(rownames(df), gene_coords$V5)
  }
  df$chr = gene_coords$V2[idx]
  df$start_pos = gene_coords$V3[idx]
  df$end_pos   = gene_coords$V4[idx]
  return(df)
}

dfGeneCoords <- function(gene_list) {
  df = data.frame(matrix(nrow=length(gene_list), ncol=3))
  rownames(df) = gene_list
  colnames(df) = c('chr','start_pos','end_pos')
  return(addGenePositions(df))
}

genesFromCoords <- function(chr, start_pos, end_pos, gene_coords) {
  target_genes = gene_coords[gene_coords$V2 == chr & gene_coords$V3 > start_pos & gene_coords$V4 < end_pos,]
  return(target_genes)
  }

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Summary counts
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
cellCountSummary <- function(seu_, col2='sample.id', nr.decimal=4) {
  tbl = table(Idents(seu_), seu_@meta.data[, col2])
  tbl = cbind(tbl, prop.table(x=tbl, margin=2))
  #first half of columns are counts; second half are proportions.  Round proportions cols.
  c1 = ncol(tbl)/2+1
  tbl[,c1:ncol(tbl)] <- mapply(round, tbl[,c1:ncol(tbl)], nr.decimal)
  return (tbl)
}

aggrBySample <- function(seu_) {
  df_cells = Idents(seu_)
  df_cells = cbind(df_cells, seu_$sample.id)
  aggr_ct  = aggregate(df_cells, by=list(df_cells[,2], df_cells[,1]), FUN=length)
  colnames(aggr_ct) = c('sample.id', 'cluster', 'nr.cells')
  aggr_ct[,2] = as.numeric(aggr_ct[,2])-1
  aggr_ct   <- aggr_ct[order(aggr_ct$cluster),1:3]
  rownames(aggr_ct) <- NULL
  return (aggr_ct)
}
