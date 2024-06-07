setwd("/data/R01/liym223/projects/skin/lib1-5/")

library(Signac)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(scater)
library(celda)
library(DropletUtils)
library(Matrix)
library(DoubletFinder)
set.seed(1234)

library(dplyr)
library(GenomeInfoDb)
library(hdf5r)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(utils)
library(tidyr)

#load("~/projects/skin/Merged_wound_scar/1010wp.RData")

######################################################################################

#....................................................................................
# CREATE SEURAT OBJECT
#....................................................................................
######################################################################################
#### lib5
# file path 
counts <- Read10X_h5("/data/R01/liym223/projects/skin/raw_lib5/new_raw_data/skin5_new/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/R01/liym223/projects/skin/raw_lib5/new_raw_data/skin5_new/outs/atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
GenomeInfoDb::seqlevelsStyle(annotation) <- "UCSC"


# create a Seurat object containing the RNA adata
lib5 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)
# create ATAC assay and add it to the object
lib5[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

lib5

######################################################################################
# FEATURE FILTERING
######################################################################################
##Quality control
setwd(dir = "/data/R01/liym223/projects/skin/lib1-5/lib5/")

DefaultAssay(lib5) <- "ATAC"

# calculate the score of NS and TSS
lib5 <- NucleosomeSignal(lib5)
lib5 <- TSSEnrichment(lib5)

#calculate percent.mt
per_barcode_metrics <- read.csv("/data/R01/liym223/projects/skin/raw_lib5/new_raw_data/skin5_new/outs/per_barcode_metrics.csv")
per_barcode_metrics_cell <- per_barcode_metrics[which(per_barcode_metrics$is_cell==1),]
per_barcode_metrics_cell$percent.mt <- per_barcode_metrics_cell$atac_mitochondrial_reads*100/per_barcode_metrics_cell$atac_raw_reads

lib5[["percent.mt.atac.reads"]] <- as.data.frame(as.data.frame(per_barcode_metrics_cell[,-c(1:31)]) %>% `rownames<-`(per_barcode_metrics_cell$barcode))
nrow(lib5[["percent.mt.atac.reads"]])
nrow(lib5@meta.data)
View(lib5[["percent.mt.atac.reads"]])



DefaultAssay(lib5) <- "RNA"
lib5[["percent.mt"]] <- PercentageFeatureSet(lib5, pattern = "^MT-")
View(lib5[["percent.mt"]])

VlnPlot(
  object = lib5,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal","nFeature_RNA","percent.mt", "percent.mt.atac.reads"),
  ncol = 4,
  pt.size = 0
)

densityRNA<-plot(density(lib5@meta.data$nCount_RNA),xlim=c(0,5000))
densityRNA_feature<-plot(density(lib5@meta.data$nFeature_RNA))
densityATAC<-plot(density(lib5@meta.data$nCount_ATAC),xlim=c(0,10000))

View(lib5@meta.data)

#raw data summary
raw_data_lib5 <- as.data.frame(lib5@meta.data)
raw_data_lib5 <- raw_data_lib5[,-1]
rmax <- vector()
rmin <- vector()
rmean <- vector()
rmed <- vector()
rvar <- vector()
rsd <- vector()
for (i in 1:ncol(raw_data_lib5))
{
  rmax = c(rmax,max(raw_data_lib5[,i]))
  rmin = c(rmin,min(raw_data_lib5[,i]))
  rmean = c(rmean,mean(raw_data_lib5[,i]))
  rmed = c(rmed,median(raw_data_lib5[,i]))
}
raw_data_lib5_summary <- as.data.frame(rbind(rmin,rmax,rmed))
colnames(raw_data_lib5_summary) <- colnames(raw_data_lib5)

write.csv(raw_data_lib5_summary, "./QC/raw_data_lib5_summary.csv")

#filter out low quality cells
#reference: signac

lib5_filtered <- subset(
  x = lib5,
  subset = nCount_ATAC < 100000 & nCount_ATAC > 400 &
    nCount_RNA < 25000 & nCount_RNA > 400 &
    nucleosome_signal < 2 &
    TSS.enrichment > 3 &
    percent.mt < 6
)

lib5_filtered

VlnPlot(
  object = lib5_filtered,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal","nFeature_RNA","percent.mt", "percent.mt.atac.reads"),
  ncol = 4,
  pt.size = 0
)


#QC data summary
QC_data_lib5 <- as.data.frame(lib5_filtered@meta.data)
QC_data_lib5 <- QC_data_lib5[,-1]
QC_data_lib5 <- QC_data_lib5[,c(1:10)]
rmax <- vector()
rmin <- vector()
rmed <- vector()

for (i in 1:ncol(QC_data_lib5))
{
  rmin = c(rmin,min(QC_data_lib5[,i]))
  rmax = c(rmax,max(QC_data_lib5[,i]))
  rmed = c(rmed,median(QC_data_lib5[,i]))
}

QC_data_lib5_summary <- as.data.frame(rbind(rmin,rmax,rmed))
colnames(QC_data_lib5_summary) <- colnames(raw_data_lib5)
write.csv(QC_data_lib5_summary, './QC/QC_data_lib5_summary.csv')

##
medianFeature<- median(lib5$nFeature_RNA)
madFeature <- mad(lib5$nFeature_RNA)
lowerBound<- min(lib5$nFeature_RNA)
upperBound <- medianFeature+madFeature*3

# PLOT: Cells per gene
cells_per_gene_lib5=ggplot(lib5_filtered@meta.data, aes(nFeature_RNA)) + geom_histogram(binwidth = 50) +
  geom_vline(xintercept = c(medianFeature), color='black') + coord_cartesian(xlim = c(0,4000)) +
  #geom_vline(xintercept = c(lowerBound, upperBound), color='red')+
  ggtitle(paste0(lib5_filtered@project.name,": median=",medianFeature,
                 #", lowerBound=",lowerBound,", upperBound=",upperBound,
                 ", Cells pre-filter=",nrow(raw_data_lib5),
                 ", Cells post-filter=",nrow(QC_data_lib5)))
cells_per_gene_lib5



######################################################################################
# DOUBLETFINDER
######################################################################################
library("DoubletFinder")
## Pre-process Seurat object -------------------------------------------------------------------------------------------------
DefaultAssay(lib5_filtered) <- "RNA"
lib5_filtered <- NormalizeData(lib5_filtered)
lib5_filtered <- ScaleData(lib5_filtered)
lib5_filtered <- SCTransform(lib5_filtered)
lib5_filtered <- RunPCA(lib5_filtered)
lib5_filtered <- RunTSNE(lib5_filtered, dims = 1:50, verbose = F)
#lib5_filtered <- RunPCA(lib5_filtered, pc.genes = lib5_filtered@var.genes, pcs.print = 0)
#lib5_filtered <- RunTSNE(lib5_filtered, dims.use = 1:10, verbose=TRUE)
lib5_filtered <- FindNeighbors(lib5_filtered, dims = 1:50, verbose = F)
lib5_filtered <- FindClusters(lib5_filtered, verbose = T)
lib5_filtered <- RunUMAP(lib5_filtered, dims = 1:50, verbose = F)
lib5_filtered <- RunUMAP(lib5_filtered, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# PLOT: umap/tsne
DimPlot(lib5_filtered, reduction = "umap", label = TRUE)
DimPlot(lib5_filtered, reduction = "tsne", label = TRUE)

# Compute expected doublet rate
cellN=nrow(lib5_filtered@meta.data)
expDoubletRate = (cellN*0.0008)*0.01#按每增加1000个细胞，双细胞比率增加千分之8来计算

## pK Identification ---------------------------------------------------------------------------------------------------------
sweep.res.list_lib5 <- paramSweep_v3(lib5_filtered, PCs = 1:30, sct = T)
sweep.stats_lib5 <- summarizeSweep(sweep.res.list_lib5, GT = FALSE)
bcmvn_lib5 <- find.pK(sweep.stats_lib5)
bcmvn_lib5$pK <- as.numeric(as.character(bcmvn_lib5$pK))
pK1=bcmvn_lib5$pK[bcmvn_lib5$BCmetric==max(bcmvn_lib5$BCmetric)]
print(head(pK1))#0.25
#NULL. the github author replied: "NULL" is just an output of the code as it prints the bcMVN x pK plot into the global environment. Does not influence the results.

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- lib5_filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(expDoubletRate*nrow(lib5_filtered@meta.data))  ##Default: Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
# use nExp_poi(934)
lib5_filtered <- doubletFinder_v3(lib5_filtered, PCs = lib5_filtered@commands$RunUMAP.SCT.pca$dims, pN = 0.25, pK = pK1, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
# use nExp_poi.adj(836)
lib5_filtered <- doubletFinder_v3(lib5_filtered, PCs = lib5_filtered@commands$RunUMAP.SCT.pca$dims, pN = 0.25, pK = pK1, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)

#lib5_filtered <- doubletFinder_v3(lib5_filtered, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_11062338", sct = T)

## Plot results --------------------------------------------------------------------------------------------------------------
lib5_filtered@meta.data$DF.classifications_0.25_0.25_934
lib5_filtered@meta.data$DF.classifications_0.25_0.25_836

table(lib5_filtered@meta.data$DF.classifications_0.25_0.25_934)
table(lib5_filtered@meta.data$DF.classifications_0.25_0.25_836)

lib5_filtered@meta.data[,"DF_hi.lo"] <- lib5_filtered@meta.data$DF.classifications_0.25_0.25_934
lib5_filtered@meta.data$DF_hi.lo[which(lib5_filtered@meta.data$DF_hi.lo == "Doublet" & lib5_filtered@meta.data$DF.classifications_0.25_0.25_836 == "Singlet")] <- "Doublet_lo"
lib5_filtered@meta.data$DF_hi.lo[which(lib5_filtered@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
table(lib5_filtered@meta.data$DF_hi.lo)
# Doublet_hi Doublet_lo    Singlet 
# 836         98       9872 

DimPlot(lib5_filtered, reduction = "umap", group.by ="DF_hi.lo",cols =c("black","red","gold"))
DimPlot(lib5_filtered, reduction = "wnn.umap", group.by ="DF_hi.lo",cols =c("black","red","gold"))

#Deleting the doublets
lib5_filtered_DoubletRM <- subset(
  x = lib5_filtered,
  subset = DF_hi.lo=="Singlet"
)

#################################################################################
##### Peak Calling 
#################################################################################
# call peaks using MACS2
peaks <- Signac::CallPeaks(lib5_filtered_DoubletRM,assay = "ATAC", macs2.path = "/data/R01/liym223/anaconda3/bin/macs2",outdir = "/data/R01/liym223/projects/skin/lib1-5/lib5/")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
DefaultAssay(lib5_filtered_DoubletRM) <- "ATAC"
macs2_counts <- FeatureMatrix(
  fragments = Fragments(lib5_filtered_DoubletRM),
  features = peaks,
  cells = colnames(lib5_filtered_DoubletRM))

# create a new assay using the MACS2 peak set and add it to the Seurat object
lib5_filtered_DoubletRM[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

saveRDS(lib5_filtered_DoubletRM, file = 'lib5_filtered_DoubletRM.RDS')


