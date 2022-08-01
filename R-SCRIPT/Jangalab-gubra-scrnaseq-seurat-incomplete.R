#load required packages & set path
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

library(devtools)
install_github('immunogenomics/presto')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
install.packages(Polychrome)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dittoSeq")

library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(patchwork)
library(presto)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(Polychrome)
library(cowplot)
library(dittoSeq)

Polychrome::swatch(viridis(20))

#load data
lean_1_chow<- read.table("GSM4795784_Lean-chow1.txt", header = TRUE, row.names = 1)
lean_2_chow<-read.table("GSM4795785_Lean-chow2.txt", header = TRUE, row.names = 1)
lean_3_chow<-read.table("GSM4795786_Lean-chow3.txt", header = TRUE, row.names = 1)
lean_4_chow<-read.table("GSM4795787_Lean-chow4.txt", header = TRUE, row.names = 1)
nash_1_dio<-read.table("GSM4795788_DIO-NASH1.txt", header = TRUE, row.names = 1)
nash_2_dio<-read.table("GSM4795789_DIO-NASH2.txt", header = TRUE, row.names = 1)

#CREATE SEURAT OBJECT LEAN-CHOW & DIO-NASH

CreateSeuratObject(
  counts=lean_1_chow, 
  project="jangalabproject", 
  min.cells = 3, 
  min.features=200
) -> lean_1_chow

lean_1_chow

CreateSeuratObject(
  counts=lean_2_chow, 
  project="jangalabproject", 
  min.cells = 3, 
  min.features=200
) -> lean_2_chow

lean_2_chow

CreateSeuratObject(
  counts=lean_3_chow, 
  project="jangalabproject", 
  min.cells = 3, 
  min.features=200
) -> lean_3_chow

lean_3_chow

CreateSeuratObject(
  counts=lean_4_chow, 
  project="jangalabproject", 
  min.cells = 3, 
  min.features=200
) -> lean_4_chow

lean_4_chow

CreateSeuratObject(
  counts= nash_1_dio, 
  project="jangalabproject", 
  min.cells = 3, 
  min.features=200
) -> nash_1_dio

nash_1_dio

CreateSeuratObject(
  counts=nash_2_dio, 
  project="jangalabproject", 
  min.cells = 3, 
  min.features=200
) -> nash_2_dio

nash_2_dio


#merge lean seurat object
merged <- merge(lean_1_chow,  y=c(lean_2_chow, lean_3_chow, lean_4_chow, nash_1_dio, nash_2_dio), add.cell.ids = c("lean_chow", "lean_chow", "lean_chow", "lean_chow", "nash_dio", "nash_dio"), project = "mergelean_nash")
merged 
rm(lean_1_chow, lean_2_chow, lean_3_chow, lean_4_chow, nash_1_dio, nash_2_dio)

#qc&filtering-------------

#create sample column
merged$sample<- rownames(merged@meta.data)
#split sample column
merged@meta.data<-separate(merged@meta.data, col = "sample", into = c("mouse", "type", "barcode"), sep = "_")

#observe merged seurat object
head(colnames(merged))
tail(colnames(merged))
unique(sapply(X = strsplit(colnames(merged), split = "_"), FUN = "[", 1))
table(merged$orig.ident)

#QC1-calculate mitochondran%
merged$MTpercent<- PercentageFeatureSet(merged, pattern = "^MT-")

#QC2- % of Largest Gene
apply(
  merged@assays$RNA@counts,
  2,
  function(x)(100*max(x))/sum(x)
) -> merged$Percent.Largest.Gene

head(merged$Percent.Largest.Gene)


#plot QC metrics
VlnPlot(merged, features=c("nCount_RNA","MTpercent"))
VlnPlot(merged, features=c("nCount_RNA","MTpercent")) + scale_y_log10()
FeatureScatter(merged,feature1 = "nCount_RNA", feature2 = "Percent.Largest.Gene")

#plot QC metrices in ggplot2 & base R
as_tibble(
  merged[[c("nCount_RNA","nFeature_RNA","Percent.Largest.Gene")]],
  rownames="Cell.Barcode"
) -> qc.metrics

qc.metrics

qc.metrics %>%
  arrange(Percent.Largest.Gene) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=Percent.Largest.Gene)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Example of plotting QC metrics") +
  geom_hline(yintercept = 500) +
  geom_hline(yintercept = 800) 

qc.metrics %>%
  ggplot(aes(Percent.Largest.Gene)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Largest gene") +
  geom_vline(xintercept = 10)

qc.metrics %>%
  ggplot(aes(Percent.Largest.Gene)) + 
  geom_histogram(binwidth = 0.7, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Largest Gene") +
  geom_vline(xintercept = 20)

#Filtering data to remove unusual QC metrices
merged_filtered <- subset(merged, subset = nCount_RNA > 800 &
                            nFeature_RNA > 300 &
                            MTpercent < 10)

merged_filtered
merged
FeatureScatter(merged_filtered,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_filtered, split.by = 'mouse')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
merged.integrated <- IntegrateData(anchorset = anchors)


# perform integration analysis--------
merged.integrated <- NormalizeData(object = merged.integrated)
merged.integrated <- FindVariableFeatures(object = merged.integrated)
merged.integrated <- ScaleData(object = merged.integrated)
merged.integrated <- RunPCA(object = merged.integrated)
ElbowPlot(merged.integrated)
merged.integrated <- FindNeighbors(object = merged.integrated, dims = 1:23)
merged.integrated <- FindClusters(object = merged.integrated)
merged.integrated <- RunUMAP(object = merged.integrated, dims = 1:23)

# plot
p1 <- DimPlot(merged.integrated, reduction = 'umap', group.by = 'mouse')
p2 <- DimPlot(merged.integrated, reduction = 'umap', label = TRUE, repel = TRUE)
p1 + p2

DimPlot(merged.integrated, reduction = "umap", split.by = "mouse")


#Finding Markers for each Cluster
#Single Prediction
view(FindMarkers(merged.integrated,ident.1 = 0, min.pct = 0.25))
VlnPlot(merged.integrated,features="Fcmr")
#Multiple Prediction
# This loop just runs the FindMarkers function on all of the clusters
lapply(
  levels(merged.integrated[["seurat_clusters"]][[1]]),
  function(x)FindMarkers(merged.integrated,ident.1 = x,min.pct = 0.25)
) -> cluster.markers

# This simply adds the cluster number to the results of FindMarkers
sapply(0:(length(cluster.markers)-1),function(x) {
  cluster.markers[[x+1]]$gene <<- rownames(cluster.markers[[x+1]])
  cluster.markers[[x+1]]$cluster <<- x
})

#Finally we collapse the list of hits down to a single table and sort it by FDR to put the most significant ones first
as_tibble(do.call(rbind,cluster.markers)) %>% arrange(p_val_adj) -> cluster.markers
cluster.markers

cluster.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1) %>%
  pull(gene) -> best.wilcox.gene.per.cluster

best.wilcox.gene.per.cluster
VlnPlot(merged.integrated,features=best.wilcox.gene.per.cluster)

#find cell type
cluster0_markers<-FindMarkers(merged.integrated,ident.1 = 0 , only.pos = TRUE, min.pct = 0.25)
cluster1_markers<-FindMarkers(merged.integrated,ident.1 = 1 , only.pos = TRUE, min.pct = 0.25)
cluster2_markers<-FindMarkers(merged.integrated,ident.1 = 2 , only.pos = TRUE, min.pct = 0.25)
cluster3_markers<-FindMarkers(merged.integrated,ident.1 = 3 , only.pos = TRUE, min.pct = 0.25)
cluster4_markers<-FindMarkers(merged.integrated,ident.1 = 4 , only.pos = TRUE, min.pct = 0.25)
cluster5_markers<-FindMarkers(merged.integrated,ident.1 = 5 , only.pos = TRUE, min.pct = 0.25)
cluster6_markers<-FindMarkers(merged.integrated,ident.1 = 6 , only.pos = TRUE, min.pct = 0.25)
cluster7_markers<-FindMarkers(merged.integrated,ident.1 = 7 , only.pos = TRUE, min.pct = 0.25)
cluster8_markers<-FindMarkers(merged.integrated,ident.1 = 8 , only.pos = TRUE, min.pct = 0.25)
cluster9_markers<-FindMarkers(merged.integrated,ident.1 = 9 , only.pos = TRUE, min.pct = 0.25)
cluster10_markers<-FindMarkers(merged.integrated,ident.1 = 10 , only.pos = TRUE, min.pct = 0.25)
cluster11_markers<-FindMarkers(merged.integrated,ident.1 = 11 , only.pos = TRUE, min.pct = 0.25)
cluster12_markers<-FindMarkers(merged.integrated,ident.1 = 12 , only.pos = TRUE, min.pct = 0.25)
cluster13_markers<-FindMarkers(merged.integrated,ident.1 = 13 , only.pos = TRUE, min.pct = 0.25)
cluster14_markers<-FindMarkers(merged.integrated,ident.1 = 14 , only.pos = TRUE, min.pct = 0.25)
cluster15_markers<-FindMarkers(merged.integrated,ident.1 = 15 , only.pos = TRUE, min.pct = 0.25)
cluster16_markers<-FindMarkers(merged.integrated,ident.1 = 16 , only.pos = TRUE, min.pct = 0.25)
cluster17_markers<-FindMarkers(merged.integrated,ident.1 = 17 , only.pos = TRUE, min.pct = 0.25)
cluster18_markers<-FindMarkers(merged.integrated,ident.1 = 18 , only.pos = TRUE, min.pct = 0.25)
cluster19_markers<-FindMarkers(merged.integrated,ident.1 = 19 , only.pos = TRUE, min.pct = 0.25)
cluster20_markers<-FindMarkers(merged.integrated,ident.1 = 20 , only.pos = TRUE, min.pct = 0.25)


head(cluster11_markers, n=10)
#save results
write.table(cluster0_markers, file='cluster0.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster1_markers, file='cluster1.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster2_markers, file='cluster2.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster3_markers, file='cluster3.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster4_markers, file='cluster4.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster5_markers, file='cluster5.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster6_markers, file='cluster6.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster7_markers, file='cluster7.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster8_markers, file='cluster8.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster9_markers, file='cluster9.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster10_markers, file='cluster10.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster11_markers, file='cluster11.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster12_markers, file='cluster12.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster13_markers, file='cluster13.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster14_markers, file='cluster14.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster15_markers, file='cluster15.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster16_markers, file='cluster16.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster17_markers, file='cluster17.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster18_markers, file='cluster18.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster19_markers, file='cluster19.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(cluster20_markers, file='cluster20.txt', quote=FALSE, sep='\t', col.names = NA)

#Assigning cell type identity to clusters
new.cluster.ids <- c('B cells', 'Endothelial cells','Macrophages', 'T cells', 'Macrophages', 'Macrophages','T cells', 'Epithelial cells', 'Dendritic cells', 'Macrophages', 'NK','Endothelial cells', 'Endothelial cells','Macrophages', 'Endothelial cells', 'Dividing cells', 'Dendritic cells', 'Hepatocytes','Dendritic cells', 'HSC', 'Granulocytes' )
names(new.cluster.ids) <- levels(merged.integrated)
merged.integrated <- RenameIdents(merged.integrated, new.cluster.ids)
DimPlot(merged.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) 
merged.integrated$seurat_clusters
#
all_markers
your_markers <- c("Cd79a","Ly6d","Cd79b" , "Aqp1" ,  "Kdr"  ,    "Gpihbp1" , "Ccl5"  ,   "Nkg7" ,    "Klrk1" ,   "Tm4sf1"  ,
                  "Fabp4"  ,  "Cavin2","Rspo3",    "Ptprb" ,   "Ear2" ,    "Napsa"  ,  "Pglyrp1" , "Plpp3"  ,  "Dnase1l3", "Stmn1"  ,  "Birc5"  ,  "Pclaf"  ,  "Siglech" , "Cox6a2",
                  "Fabp1"  ,  "Wfdc21" ,  "Car3"  ,   "Gadd45b" , "Fscn1",    "Relb"   ,  "Dcn"   ,   "Col3a1" ,  "Col14a1" , "Ms4a6c"  , "Alox5ap" , "S100a4",  
 "Ifitm1" ,  "Ccl9"  ,   "Hdc"    ,  "Ms4a4b" ,  "Cd3d"    , "Cd3g"   ,  "C1qb"   ,  "C1qc"  ,   "C1qa"    , "Spp1" ,    "Clu"   ,   "Tm4sf4"  ,
 "Plbd1"  ,  "H2-DMb1" )

dittoDimPlot(merged.integrated, "mouse")
#active.ident
dittoBarPlot(merged.integrated, "cell_type", group.by = "mouse")
dittoBarPlot(merged.integrated.sce, "cell_type", group.by = "mouse",
             scale = "count")
dittoBoxPlot(merged.integrated, "Cd74", group.by = "mouse")
#
top2 = cluster.markers %>% dplyr::group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
plot = DotPlot(object = merged.integrated, features = unique(top2$gene), assay="RNA")
plot + coord_flip()

DotPlot(merged.integrated, features = unique(top2$gene), cols = c("blue", "red"), dot.scale = 8, split.by = "mouse") +
  RotatedAxis()

genes<-unique(top2$gene)


dittoHeatmap(merged.integrated, genes,
             colors = 1:11,
             color.panel = c("#F8766D", "#DB8E00", "#AEA200", "#64B200", "#00BD5C", "#00C1A7", "#00BADE", "#00A6FF", "#B385FF", "#EF67EB", "#FF63B6"),
             annot.by = c("cell_type", "mouse"),
             scaled.to.max = TRUE,
             max.color = "#0A2F51", min.color = "#DEEDCF",
             )

#retrieved tables
cluster.averages <- AverageExpression(merged.integrated)
write.table(merged.integrated@active.ident, file='Convert_UMI_Label3.tsv',quote=TRUE, sep='\t', col.names = NA)
write.table(merged.integrated@assays[["RNA"]]@counts, file='Gene_Count_per_Cell3.tsv', quote=FALSE, sep='\t', col.names = NA)
merged.integrated[["RNA"]]@counts

#create reference file for seurat

celltype<- as.data.frame(merged.integrated@active.ident, colnames= c(' barcodes, cell_types'))
head(celltype)

merged.integrated <- AddMetaData(
  object = merged.integrated,
  metadata = merged.integrated@active.ident,
  col.name = 'cell_type')

merged.integrated@meta.data
mouseSigMatrixGUBRA<- data.frame(merged.integrated[["RNA"]]@counts)

celltype_labels<-t(merged.integrated@meta.data$cell_type)
colnames(mouseSigMatrixGUBRA)<-celltype_labels
head(mouseSigMatrixGUBRA)
write.table(mouseSigMatrixGUBRA, file = 'mousesignaturematixGUBRA.txt', quote = FALSE,sep='\t', col.names = NA)


#Pull number of cells in a seurat object

library(data.table)
library(magrittr)
## extract meta data
md <- merged.integrated@meta.data %>% as.data.table
# the resulting md object has one "row" per cell

## count the number of cells per unique combinations of "mouse" and "cell_type"
viewmd<-md[, .N, by = c("mouse", "cell_type")]

view(viewmd)

#Alternative
sum(merged.integrated$cell_type == "Macrophages")
sum(merged.integrated$cell_type == "Hepatocytes")

