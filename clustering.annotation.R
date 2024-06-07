######################## Clustering and annotation  ###################
#### (1) marker list ####
Skin_cellType_markers <- read.csv("/data/R01/liym223/projects/skin/merged_lib1-4/cell_type_annotation/skin_celltype_marker_new.csv", header = T)

########################################################
#### (2) clustering ####
# build a joint neighbor graph using both assays
#We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering
DefaultAssay(lib12345) <- "integratedRNA"
lib12345 <- RunPCA(lib12345, npcs = 50, verbose = FALSE)

# Test FindMultiModalNeighbors dimension setting
# Selection criteria
# 1.Cumulative contribution of principal components is greater than 90%
# 2.PC itself contributes less than 5% to the variance
# 3.Less than 0.1% difference between two consecutive PCs
pct <- lib12345[["pca"]]@stdev / sum( lib12345[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC 
cumu <- cumsum(pct) 
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5 
co1 <- which(cumu > 90 & pct < 5)[1] 
co1 #41
# Determine the difference between variation of PC and subsequent PC 
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%. 
co2 #16
# Minimum of the two calculation 
pcs <- min(co1, co2) 
pcs #16
# Create a dataframe with values 
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct)) 
# Elbow plot to visualize
pdf("./Routput/cluster/plotdf_pca.pdf") 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") + 
  theme_bw()
dev.off()


# Dimensional calculations are also performed for lsi
lst <- lib12345 [["lsi"]]@stdev / sum( lib12345[["lsi"]]@stdev) * 100
# Calculate cumulative percents for each PC 
cumu2 <- cumsum(lst) 
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5 
co1 <- which(cumu2 > 90 & lst < 5)[1] 
co1 #39
# Determine the difference between variation of PC and subsequent PC 
co2 <- sort(which((lst[1:length(lst) - 1] - lst[2:length(lst)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%. 
co2 #13
# Minimum of the two calculation 
pcs2 <- min(co1, co2) 
pcs2 #13
# Create a dataframe with values 
plot_df2 <- data.frame(pct = lst, cumu = cumu2, rank = 1:length(lst)) 
# Elbow plot to visualize
pdf("./Routput/cluster/plotdf_lsi.pdf") 
ggplot(plot_df2, aes(cumu2, lst, label = rank, color = rank > pcs2)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(lst[lst > 5]), color = "grey") + 
  theme_bw()
dev.off()

lib12345 <- FindMultiModalNeighbors(
  object = lib12345,
  reduction.list = list("pca", "integratedLSI"), 
  dims.list = list(1:16, 2:16), # The dim can usually be adjusted to get better clustering results
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

lib12345 <- FindNeighbors(object = lib12345, reduction = 'pca', dims = 1:40)

#### (3) clustree ####
# Use clustree to see which resolution is better
library(clustree)
lib12345 <- FindClusters(object = lib12345, verbose = FALSE,algorithm = 3, resolution = c(seq(.1,1.5,.1)))
#colnames(lib12345@meta.data)
pdf("./Routput/cluster/clustree_lib12345_unsupervised.pdf",width=14, height=14)
clustree(lib12345@meta.data, prefix = "integratedRNA_snn_res.")
dev.off()

lib12345 <- FindClusters(object = lib12345, reduction.type="umap", verbose = FALSE, 
                            resolution =0.5, algorithm = 3)

# Build a joint UMAP visualization
lib12345 <- RunUMAP(
  object = lib12345,
  nn.name = "weighted.nn",
  assay = "integratedRNA",
  verbose = TRUE
)

######### (4)orig.ident DimPlot #########
Idents(lib12345)<-lib12345$orig.ident
lib12345$orig.ident<-factor(lib12345$orig.ident,levels=c("scar1","scar2","wound1","lib3","lib4"))

pdf("./Routput/cluster/lib12345_cluster_splitbysample.pdf",width=15, height=8)
DimPlot(lib12345, label = TRUE, repel = TRUE, reduction = "umap",split.by = "orig.ident")
dev.off()
pdf("./Routput/cluster/lib12345_cluster_sample_consistence.pdf")
DimPlot(lib12345, label = TRUE, repel = TRUE, reduction = "umap",group.by = "orig.ident")
DimPlot(lib12345, label = TRUE, repel = TRUE, reduction = "umap",group.by = "seurat_clusters")
dev.off()

# Save the integrated lib12345 object
saveRDS(lib12345,"./Routput/lib12345_integrated.RDS")


#### (5) compute wnn ####
# build a joint neighbor graph using both assays
#We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

lib12345 <- RunUMAP(lib12345, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
lib12345 <- FindClusters(lib12345, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

lib12345 <- FindMultiModalNeighbors(lib12345, reduction.list = list("pca", "lsi"), dims.list = list(1:42, 2:40))

#use wsnn to cluster
#lib12345$wsnn_res.0.8
pdf("./Routput/cluster/groupbySNN_wnnumap_res0.8.pdf")
DimPlot(lib12345, label = TRUE, repel = TRUE, reduction = "umap",group.by = "integratedRNA_snn_res.0.8")
dev.off()
#30 clusters
#get the bubble plot
pdf("./Routput/cluster/bubbleplot_res0.8.pdf", width=16, height=14)
DotPlot(lib12345,features = markers_update,assay = 'SCT',group.by = "wsnn_res.0.8")+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.85,hjust = 0.75))
dev.off()
#res=0.8, ?cluster
#NA markers: IGJ, Cd79A, GLDN5, VMF, CD31, CD45, Cytokeratin-20, NCAM, HLA-DRB, CD56

#### (6) Rename seurat_clusters####
newClusterIDs <- c(
  '1' = 'Endothelial cell','7' = 'Endothelial cell','19' = 'Endothelial cell','18_0' = 'Endothelial cell','18_6' = 'Endothelial cell',
  '0' = 'Fibroblast','4' = 'Fibroblast','8' = 'Fibroblast','9' = 'Fibroblast','14' = 'Fibroblast','21' = 'Fibroblast',
  '2' = 'Keratinocyte','10' = 'Keratinocyte','11' = 'Keratinocyte','18_1' = 'Keratinocyte','18_3' = 'Keratinocyte', '18_4' = 'Keratinocyte','18_5' = 'Keratinocyte',
  '13' = 'Smooth muscle cell','5' = 'Smooth muscle cell','17' = 'Smooth muscle cell',
  '12' = 'Immune cell', '15' = 'Immune cell','18_2' = 'Immune cell',
  '23' = 'Melanocyte',
  '24' = 'Schwann cell','16' = 'Schwann cell',
  '3' = 'Sweat and sebaceous gland cell',
  '6' = 'Hair follicle cell', '20' = 'Hair follicle cell', '22' = 'Hair follicle cell'
  )

lib12345 <- RenameIdents(lib12345, newClusterIDs)
lib12345$celltype <- lib12345@active.ident


