#
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#Load libraries
library(reshape)
library(ggplot2)
library(readxl)
library(pheatmap)
library(wesanderson)
library(tidyverse)
library(data.table)
library(dittoSeq)
library(magrittr)
###

##Xiongx & gubra -sigcell data

merged.integrated
###
dittoDimPlot(merged.integrated, "mouse")
dittoBarPlot(merged.integrated, "cell_type", group.by = "mouse")
dittoBarPlot(merged.integrated, "cell_type", group.by = "mouse",
             scale = "count")

#Pull number of cells in a seurat object
## extract meta data
md <- merged.integrated@meta.data %>% as.data.table
# the resulting md object has one "row" per cell

## count the number of cells per unique combinations of "mouse" and "cell_type"
viewmd<-md[, .N, by = c("mouse", "cell_type")]
view(viewmd)

#R-Scrip for visualizations

#####plot correlation matrix for single cell data

#Gubra-data
GUBRAsigmatix<- read.table('mousesignaturematixGUBRA_inferred_phenoclasses.txt', header = TRUE, sep = "\t",  row.names =1)
head(GUBRAsigmatix)
#correlation
cor_GUBRA<- cor(GUBRAsigmatix)
#visualize-2
cor_GUBRA2 <- melt(cor_GUBRA)
cor_GUBRA2$value<- round(cor_GUBRA2$value, digits = 2)
head(cor_GUBRA2)
#
heatmap2<- ggplot(cor_GUBRA2, aes(x = X1, y = X2, fill = value)) +
  geom_tile(color = "white", lwd= 1.5, linetype= 1) +
  scale_fill_gradient2(low = "#2C7BB6",
                       mid = "white",
                       high = "#D7191C", breaks=seq(-1,1,0.2), limits= c(-1,1)) +
  coord_fixed() + 
  theme_minimal(base_family="Helvetica")+
  guides(fill = guide_colourbar(barwidth = 0.5,barheight = 20))
 
heatmap2 + 
  geom_text(aes(X2, X1, label = value), color = "black", size = 4) +
  theme(axis.ticks=element_blank())+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

## XIONGX-DATA
Xiongsigmatix<- read.table('mousesignaturematixX_inferred_phenoclasses.txt', header = TRUE, sep = "\t",  row.names =1)
head(Xiongsigmatix)
#correlation
cor_Xiong<- cor(Xiongsigmatix)
#visualize-2
cor_Xiong2 <- melt(cor_Xiong)
cor_Xiong2$value<- round(cor_Xiong2$value, digits = 2)
head(cor_Xiong2)
#
heatmap2<- ggplot(cor_Xiong2, aes(x = X1, y = X2, fill = value)) +
  geom_tile(color = "white", lwd= 1.5, linetype= 1) +
  scale_fill_gradient2(low = "#2C7BB6",
                       mid = "white",
                       high = "#D7191C", breaks=seq(-1,1,0.2), limits= c(-1,1)) +
  coord_fixed() + 
  theme_minimal(base_family="Helvetica")+
  guides(fill = guide_colourbar(barwidth = 0.5,barheight = 20))

heatmap2 + 
  geom_text(aes(X2, X1, label = value), color = "black", size = 4) +
  theme(axis.ticks=element_blank())+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))


#P-VALUE
#STATS
## to identify the differences between two independent  conditions.  
#However, the continuous response variable of interest is not normally distributed.

wilcox.test(gubraresult$Macrophages~gubraresult$Type)
wilcox.test(gubraresult$Hepatocytes~gubraresult$Type)
wilcox.test(gubraresult$`Endothelial cells`~gubraresult$Type)
wilcox.test(gubraresult$`Epithelial cells`~gubraresult$Type)
wilcox.test(gubraresult$HSC~gubraresult$Type)



######DECONVOLUTION BAR PLOTS
#barplot celltype with gubra
gubraresult<- read_excel("Gubradeconvoluted.xlsx", sheet = "Sheet1")
gubraresult<- as.data.frame(gubraresult)
#reshape
gubraresult2 <- melt(gubraresult)
gubraresult2<- as.data.frame(gubraresult2)
gubraresult2['new']<-c('DIO-NASH','DIO-NASH','DIO-NASH','DIO-NASH', 'CHOW','CHOW', 'CHOW', 'CHOW', 'CHOW' )

gubraplot<- ggplot(gubraresult2, aes(fill=variable, y= value, x= Mixture))+ theme(axis.ticks=element_blank())+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  geom_bar(position="stack", stat="identity")+facet_grid(. ~ new, drop=TRUE,scale="free",space="free_x")+xlab("Samples")+ ylab("Cell type Proportions")
 
#plot with standard error bars
gubraresult
gubraresult$Type<- c('DIO-NASH','DIO-NASH','DIO-NASH','DIO-NASH','Lean','Lean','Lean','Lean','Lean' )
GG <- aggregate(cbind(Macrophages,Hepatocytes,`Endothelial cells`, `Epithelial cells`, HSC)~Type , data=gubraresult , mean)
rownames(GG) <- GG[,1]
GG <- as.matrix(GG[,-1])
Gg<-as.data.frame(GG)
#Plot boundaries
lim <- 1.2*max(GG)
#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
#Then I calculate the standard deviation:
stdev <- aggregate(cbind(Macrophages,Hepatocytes,`Endothelial cells`, `Epithelial cells`, HSC)~Type , data=gubraresult , sd)
rownames(stdev) <- stdev[,1]
stdev <- as.matrix(stdev[,-1]) * 1.96 / 10
#I am ready to add the error bar on the plot using my "error bar" function !
ze_barplot <- barplot(GG , beside=T , legend.text=T,col=c("blue" , "skyblue") , ylim=c(0,lim) , ylab="height")
error.bar(ze_barplot,GG, stdev)

#P-VALUE
#STATS
## to identify the differences between two independent  conditions.  
#However, the continuous response variable of interest is not normally distributed.

wilcox.test(gubraresult$Macrophages~gubraresult$Type)
wilcox.test(gubraresult$Hepatocytes~gubraresult$Type)
wilcox.test(gubraresult$`Endothelial cells`~gubraresult$Type)
wilcox.test(gubraresult$`Epithelial cells`~gubraresult$Type)
wilcox.test(gubraresult$HSC~gubraresult$Type)

####

#barplot celltype with JANGALABDATA
Jangaresult<- read_excel("CIBERSORTx_Job178_Results.xlsx", sheet = "Sheet1")
Jangaresult<- as.data.frame(Jangaresult)
#reshape
Jangaresult2 <- melt(Jangaresult)
Jangaresult2<- as.data.frame(Jangaresult2)
Jangaresult2['new']<-c('Healthy','Healthy','NASH','NASH' )

Jangaplot<- ggplot(Jangaresult2, aes(fill=variable, y= value, x= Mixture))+ theme(axis.ticks=element_blank())+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  geom_bar(position="stack", stat="identity")+facet_grid(. ~ new, drop=TRUE,scale="free",space="free_x")+xlab("Samples")+ ylab("Cell type Proportions")

#plot with standard error bars
Jangaresult
Jangaresult$Type<- c('Healthy','Healthy','NASH','NASH' )

GG <- aggregate(cbind(Macrophages,Hepatocytes,`Endothelial cells`, Cholangiocytes, HSC)~Type , data=Jangaresult , mean)
rownames(GG) <- GG[,1]
GG <- as.matrix(GG[,-1])
Gg<-as.data.frame(GG)

#Plot boundaries
lim <- 1.2*max(GG)

#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#Then I calculate the standard deviation for each specie and condition :
stdev <- aggregate(cbind(Macrophages,Hepatocytes,`Endothelial cells`, Cholangiocytes, HSC)~Type , data=Jangaresult , sd)
rownames(stdev) <- stdev[,1]
stdev <- as.matrix(stdev[,-1]) * 1.96 / 10

#I am ready to add the error bar on the plot using my "error bar" function !
ze_barplot <- barplot(GG , beside=T , legend.text=T,col=c("blue" , "skyblue") , ylim=c(0,lim) , ylab="height")
error.bar(ze_barplot,GG, stdev)

#STATS
## to identify the differences between two independent  conditions.  
#However, the continuous response variable of interest is not normally distributed.

wilcox.test(Jangaresult$Macrophages~Jangaresult$Type)
wilcox.test(Jangaresult$Hepatocytes~Jangaresult$Type)
wilcox.test(Jangaresult$`Endothelial cells`~Jangaresult$Type)
wilcox.test(Jangaresult$Cholangiocytes~Jangaresult$Type)
wilcox.test(Jangaresult$HSC~Jangaresult$Type)
####
#barplot celltype with furutaDATA
#furuta et al
furutaresult<- read_excel("/Users/kiteomoru/Downloads/JangaJOB/CIBERSORTxResults/DeconvolutionResult/CIBERSORTx_Job181_Results.xlsx", sheet = "Sheet1")
furutaresult<- as.data.frame(furutaresult)

#reshape
furutaresult<- as.data.frame(furutaresult)
furutaresult2['new']<-c('CHOW','CHOW', 'CHOW','FFC-NASH','FFC-NASH','FFC-NASH')
furutaresult2 <- melt(furutaresult)

furutaplot<- ggplot(furutaresult2, aes(fill=variable, y= value, x= Mixture))+ theme(axis.ticks=element_blank())+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  geom_bar(position="stack", stat="identity")+facet_grid(. ~ new, drop=TRUE,scale="free",space="free_x")+xlab("Samples")+ ylab("Cell type Proportions")

#plot with standard error bars
furutaresult2
furutaresult2$Type<- c('FFC-NASH','FFC-NASH','FFC-NASH', 'CHOW','CHOW', 'CHOW' )

GG <- aggregate(cbind(Macrophages,Hepatocytes,`Endothelial cells`, Cholangiocytes, HSC)~new , data=furutaresult , mean)
rownames(GG) <- GG[,1]
GG <- as.matrix(GG[,-1])
Gg<-as.data.frame(GG)

#Plot boundaries
lim <- 1.2*max(GG)

#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#Then I calculate the standard deviation for each specie and condition :
stdev <- aggregate(cbind(Macrophages,Hepatocytes,`Endothelial cells`, Cholangiocytes, HSC)~new , data=furutaresult , sd)
rownames(stdev) <- stdev[,1]
stdev <- as.matrix(stdev[,-1]) * 1.96 / 10

#I am ready to add the error bar on the plot using my "error bar" function !
ze_barplot <- barplot(GG , beside=T , legend.text=T,col=c("blue" , "skyblue") , ylim=c(0,lim) , ylab="height")
error.bar(ze_barplot,GG, stdev)

#
#STATS
## to identify the differences between two independent  conditions.  
#However, the continuous response variable of interest is not normally distributed.

wilcox.test(furutaresult$Macrophages~furutaresult$new)
wilcox.test(furutaresult$Hepatocytes~furutaresult$new)
wilcox.test(furutaresult$`Endothelial cells`~furutaresult$new)
wilcox.test(furutaresult$Cholangiocytes~furutaresult$new)
wilcox.test(furutaresult$HSC~furutaresult$new)
###
####


#Sgnature ganes

Hsc<-  c('Dcn', 'Gsn', 'Ecm1', 'Cxcl12', 'Bgn', 'Colec11', 'Col14a1', 'Col3a1', 'Pam', 'Ifitm1', 'Rgs5', 'Prelp', 'Sod3', 'Reln', 'Lrat', 'Igfbp5', 'Steap4', 'Col1a2', 'Angptl6', 'Lum', 'Slpi', 'Col1a1', 'Cygb', 'Igfbp3', 'Abcc9', 'Dpt', 'Tagln', 'Prss23', 'Igfbp6', 'Colec10', 'Pth1r', 'C1s1', 'Rspo3', 'Mustn1', 'Vipr1', 'Mgp', 'Gpx3', 'Pcolce', 'Htra3', 'Col6a1', 'Smoc2', 'Col6a3', 'Sparcl1', 'Aebp1', 'Clec3b', 'Nrxn1', 'Pdgfrb', 'Tcf21', 'Tnxb', 'Abi3bp')
Hepatocytes<-c( 'Alb','Fabp1', 'Apoa2', 'Mup20', 'Apoc3', 'Mup3', 'Apoa1', 'Aldob', 'Gsta3', 'Car3', 'Serpina1b', 'Bhmt', 'Scd1', 'Cdo1', 'Serpina1a', 'Acaa1b', 'Ass1', 'Akr1c6', 'Phyh', 'Apoc4', 'Serpina1e', 'Wfdc21', 'Ttc36', 'Gnmt', 'Rgn', 'Hp', 'Serpina3k', 'Uox', 'Hpx', 'Hpd', 'Hamp', 'Rida', 'Fbp1', 'Orm1', 'Sord', 'Angptl3', 'Fgg', 'Serpina1d', 'Fgb', 'Mat1a', 'Fabp2', 'H2.Q10', 'Tdo2', 'Arg1', 'Acox1', 'Otc', 'Apoh', 'Pck1', 'Cyp3a11', 'Serpina12', 'Gamt')
Cholangiocytes<- c('Spp1', 'Clu', 'Tm4sf4', 'Cp', 'Tstd1', 'Epcam', 'Atp1b1', 'Ankrd1', 'Sorbs2', 'Tspan8', 'Ddit4l', 'Onecut2', 'Sox9', 'Cldn3', 'Fxyd3', 'Cd24a', 'Dsg2', 'Pdzk1ip1', 'Anxa4', 'Cldn7', 'Kifc3', 'Scara3', 'Cdh1', 'Fgfr3', 'Dcdc2a', 'Smim22', 'Prox1', 'Gsta3', 'Plet1', 'Krt19', 'Pkhd1', 'Sftpd', 'Wwc1', 'Aldob', 'Clmn', 'Ptprf', 'Slc5a1', 'Msmo1', 'Cxadr', 'Slc2a2', 'Bche', 'Wfdc2', 'Fdft1', 'Tesc', 'Enc1', 'Kcne3', 'Vegfa', 'X1700011H14Rik', 'Bex4', 'Wfdc15b')
Kc <- c('Lyz2', 'Wfdc17','Tyrobp', 'Fcer1g', 'Ctss', 'C1qb', 'C1qc', 'Cd5l', 'C1qa', 'Cd52', 'Laptm5', 'Clec4f', 'Lst1', 'Cybb', 'Vsig4', 'Aif1', 'Csf1r', 'Mpeg1', 'Hmox1', 'Ly86', 'Lcp1','Coro1a', 'Mafb', 'Ccl6', 'S100a4', 'Fyb', 'Ms4a6c', 'Cd44', 'Ptprc', 'Il1b', 'Cd53', 'Spi1','Cfp', 'Alox5ap', 'H2.DMa', 'Itgal', 'Rgs1', 'Fcgr4', 'Ear2', 'Fam49b', 'Ifi30', 'Ccl4', 'Bcl2a1b', 'H2.DMb1', 'AF251705', 'Cd68', 'Acp5', 'Rac2', 'Tnfaip2', 'Clec4a3')
Endo<-  c('Aqp1', 'Dnase1l3', 'Fcgr2b', 'Kdr', 'Fabp4', 'Clec4g', 'Gpihbp1', 'Gpr182', 'Stab2',  'Cyp26b1', 'Cldn5', 'Stab1', 'Cd55', 'Flt1', 'Oit3', 'Adgrl4', 'Ntn4', 'Fam167b', 'Adam23', 'Tmem2', 'Plxnc1', 'Pde2a', 'Cd300lg', 'Elk3', 'Mmrn2', 'Pecam1', 'Flt4', 'Cyyr1', 'Cyp4b1', 'Tpbgl', 'Rasgrp3', 'Bgn', 'Gimap4', 'Acer2', 'Tek', 'Dpp4', 'Nrp2', 'Sox18', 'Ctla2a', 'F8', 'Plk2', 'Ap1b1', 'Myct1', 'Gimap6', 'Fam124a', 'Agpat5', 'Jam2', 'Arhgap31', 'Tie1', 'X1810011O10Rik')


#function for scaling
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

#create metadata jangala&furuta
rep<- c(1,2,3,4,5 ,1,2,3,4,5)
protocol<-c('Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy','NASH', 'NASH', 'NASH', 'NASH','NASH')
tissue<- c('liver','liver','liver','liver','liver','liver','liver','liver','liver', 'liver')
organism<-c('Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse', 'Mouse' )

metadata_df<-data_frame(rep, protocol, tissue, organism)
 rw<-c('CHOW1', 'CHOW2', 'CHOW3', 'Liver1','Liver2', 'FFC1', 'FFC2', 'FFC3', 'NASH-R1', 'HASH-R2')
 metadata_df<-as.data.frame(metadata_df) 
 rownames(metadata_df)<-rw 

##
#Preparing the data for object
rownames(metadata_df)
#check if the rownames and column names fo the two data sets match
rownames(metadata3_df) == colnames(jangaFurutaMIXTURE)


#create metadata gubradata
#create metadata
rep<- c(1,2,3,4,1,2,3,4,5)
protocol<-c('NASH', 'NASH', 'NASH', 'NASH', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
tissue<- c('liver','liver','liver','liver','liver','liver','liver','liver','liver')
organism<-c('Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse' )

metadatag_df<-data_frame(rep, protocol, tissue, organism)
rwg<-c('DIONASH1', 'DIONASH2', 'DIONASH3', 'DIONASH4','Lean1', 'Lean2', 'Lean3', 'Lean4', 'Lean5')
metadatag_df<-as.data.frame(metadatag_df) 

rownames(metadatag_df)<-rwg
colnames(metadatag_df)

metadatag_df<- metadatag_df[,-5]
##
#Preparing the data for object
rownames(metadatag_df)
colnames(GUBRAMIXTURE)
#check if the rownames and column names fo the two data sets match
rownames(metadatag_df) == colnames(GUBRAMIXTURE)


#heatmaps for signature genes expression

#MIXTUREDATA-jangalab_furuta
jangaFurutaMIXTURE<-read.delim('ALLRAW_COUNTS_JANGA_FURUTA copy.txt', header = TRUE , sep = "\t")
jangaFurutaMIXTURE= as.data.frame(jangaFurutaMIXTURE)
jangaFurutaMIXTURE= jangaFurutaMIXTURE[!duplicated(jangaFurutaMIXTURE$Gene),]
row.names(jangaFurutaMIXTURE)<- jangaFurutaMIXTURE$Gene
head(jangaFurutaMIXTURE)
str(jangaFurutaMIXTURE)
jangaFurutaMIXTURE= jangaFurutaMIXTURE[,2:11]

# transform to Z-scale
jangaFurutaMIXTURE <- t(scale(t(data.matrix(jangaFurutaMIXTURE))))
jangaFurutaMIXTURE= as.matrix(jangaFurutaMIXTURE)

#REMOVE NA VALUES
jangaFurutaMIXTURE<-jangaFurutaMIXTURE %>%
  na.omit()
#PLOT HEATMAP
pheatmap(jangaFurutaMIXTURE [rownames(jangaFurutaMIXTURE) %in% Hepatocytes ,],labels_row =  Hepatocytes ,  
         scale="row", main = "Hepatocytes", annotation = dplyr::select(metadata3_df, protocol), color = pal)

pheatmap(jangaFurutaMIXTURE [rownames(jangaFurutaMIXTURE) %in% Hsc ,],labels_row =  Hsc ,  
         scale="row", main = "Hsc", annotation = dplyr::select(metadata3_df, protocol), color = pal)

pheatmap(jangaFurutaMIXTURE [rownames(jangaFurutaMIXTURE) %in% Cholangiocytes ,],labels_row =  Cholangiocytes ,  
         scale="row", main = "Cholangiocytes", annotation = dplyr::select(metadata3_df, protocol), color = pal)

pheatmap(jangaFurutaMIXTURE [rownames(jangaFurutaMIXTURE) %in% Endo ,],labels_row =  Endo ,  
         scale="row", main = "Endo", annotation = dplyr::select(metadata3_df, protocol), color = pal)

pheatmap(jangaFurutaMIXTURE [rownames(jangaFurutaMIXTURE) %in% Kc ,],labels_row =  Kc ,  
         scale="row", main = "Kc", annotation = dplyr::select(metadata3_df, protocol), color = pal)

#GUBRAMIXTURE
#MIXTUREDATA-jangalab_furuta
GUBRAMIXTURE<-read.delim('ALLRAW_COUNTS_JANGA_FURUTA copy.txt', header = TRUE , sep = "\t")
GUBRAMIXTURE= as.data.frame(GUBRAMIXTURE)
GUBRAMIXTURE= GUBRAMIXTURE[!duplicated(GUBRAMIXTURE$Gene),]
row.names(GUBRAMIXTURE)<- GUBRAMIXTURE$Gene
head(GUBRAMIXTURE)
str(GUBRAMIXTURE)
GUBRAMIXTURE= GUBRAMIXTURE[,2:11]

# transform to Z-scale
GUBRAMIXTURE <- t(scale(t(data.matrix(GUBRAMIXTURE))))
GUBRAMIXTURE= as.matrix(GUBRAMIXTURE)

#REMOVE NA VALUES
GUBRAMIXTURE<-GUBRAMIXTURE %>%
  na.omit()
#PLOT HEATMAP
pheatmap(GUBRAMIXTURE [rownames(GUBRAMIXTURE) %in% Hepatocytes ,],labels_row =  Hepatocytes ,  
         scale="row", main = "Hepatocytes", annotation = dplyr::select(metadatag_df, protocol), color = pal)

pheatmap(GUBRAMIXTURE [rownames(GUBRAMIXTURE) %in% Hsc ,],labels_row =  Hsc ,  
         scale="row", main = "Hsc", annotation = dplyr::select(metadatag_df, protocol), color = pal)

pheatmap(GUBRAMIXTURE [rownames(GUBRAMIXTURE) %in% Cholangiocytes ,],labels_row =  Cholangiocytes ,  
         scale="row", main = "Cholangiocytes", annotation = dplyr::select(metadatag_df, protocol), color = pal)

pheatmap(GUBRAMIXTURE [rownames(GUBRAMIXTURE) %in% Endo ,],labels_row =  Endo ,  
         scale="row", main = "Endo", annotation = dplyr::select(metadatag_df, protocol), color = pal)

pheatmap(GUBRAMIXTURE [rownames(jangaFurutaMIXTURE) %in% Kc ,],labels_row =  Kc ,  
         scale="row", main = "Kc", annotation = dplyr::select(metadatag_df, protocol), color = pal)



#NASH Associations heatmap

associations<- read_excel("NASHAssociations.xlsx", sheet = "Sheet2")
associations%>%
  as.data.frame()->df_associations

rownames(df_associations) = df_associations$Gene
df_associations %>%
  dplyr::select(3:21)-> df_ass

df_ass_norm <- t(apply(df_ass, 1, cal_z_score))

df_associations %>%
  dplyr::select(2)-> nashasso_ano

##/Users/kiteomoru/Desktop/HTP NOTES/Copy of TPMSALLLIVER.xlsx
associations2<- read_excel("/Users/kiteomoru/Desktop/HTP NOTES/Copy of TPMSALLLIVER.xlsx", sheet = "Sheet1")
associations2%>%
  as.data.frame()->associations2

rownames(associations2) = associations2$GeneSymbol
associations2 %>%
  dplyr::select(3:21)-> df_ass

df_ass_norm <- t(apply(df_ass, 1, cal_z_score))

associations2 %>%
  dplyr::select(2)-> nashasso_ano

##### create metadata
colnames(df_ass_norm)
rep<- c(1,2,3,4,5,6,7,8,9 ,1,2,3,4,5,6,7,8,9,10)
protocol<-c('NASH', 'NASH', 'NASH', 'NASH','NASH', 'NASH','NASH','NASH','NASH','Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy','Healthy','Healthy','Healthy','Healthy','Healthy')
tissue<- c('liver','liver','liver','liver','liver','liver','liver','liver','liver', 'liver','liver','liver','liver','liver','liver','liver','liver','liver','liver')
organism<-c('Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse', 'Mouse' , 'Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse','Mouse')

metadataAll<-data_frame(rep, protocol, tissue, organism)
rw<-c('srr_3_ffc3', 'srr_2_ffc2', 'srr_1_ffc1',  'DIONASH1'  ,   'DIONASH2'  ,   'DIONASH3'   ,  'DIONASH4' ,'NASH-R2','NASH-R1',  'liver1_count', 'liver2_count','srr_5_chow2', 'srr_4_chow1', 'srr_6_chow3','Lean1'    , 'Lean2' ,      
       'Lean3'   ,     'Lean4'    ,    'Lean5' )
metadataAll<-as.data.frame(metadataAll) 
##########
rownames(metadataAll)<-rw 

rownames(metadataAll) == colnames(df_ass_norm)


metadataAll
nashcolanno_df = data.frame("Protocol" = metadataAll$protocol)
rownames(nashcolanno_df) = rownames(metadataAll)

nashrowanno_df = data.frame(associations = factor(nashasso_ano$Association, ordered = TRUE))
rownames(nashrowanno_df) = rownames(nashasso_ano)

pheatmap(df_ass_norm , annotation_row = nashrowanno_df , annotation_col = nashcolanno_df,cluster_rows = F, scale = "row", color = pal2, cutree_rows = 5 )


#####3Functional analysis plots
DEGs<- read_excel("", sheet = "")
DEGs_plot<- ggplot(DEGs, aes(reorder(x=Term, -log10(`PValue`)), y=-log10(`PValue`), fill=  Expression))+
  geom_bar(stat="identity")+
  scale_fill_brewer(palette= 'Dark2')+
  theme_minimal()+
  coord_flip()+
  facet_grid( Category~., scale="free")



