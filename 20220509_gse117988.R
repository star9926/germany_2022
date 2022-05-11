#20220509gse117988
rm(list = ls()) 
options(warn=-1) 
suppressMessages(library(Seurat))
# 读取表达矩阵
start_time <- Sys.time()
raw_dataPBMC <- read.csv('./GSE117988_raw.expMatrix_PBMC.csv.gz', 
                         header = TRUE, 
                         row.names = 1)
end_time <- Sys.time()
end_time - start_time
#Time difference of 1.486653 mins

dim(raw_dataPBMC) 

#normalization
dataPBMC <- log2(1 + sweep(raw_dataPBMC, 2, 
                           median(colSums(raw_dataPBMC))/colSums(raw_dataPBMC), '*'))
head(colnames(dataPBMC))
tail(colnames(dataPBMC))

timepoints<-sapply(colnames(dataPBMC), 
                   function(x) unlist(strsplit(x,'\\.'))[2])
table(timepoints)

timepoints<-ifelse(timepoints=='1','PBMC_Pre',
                   ifelse(timepoints=='2','PBMC_EarlyD27',
                          ifelse(timepoints=='3','PBMC_RespD376', 'PBMC_ARD614')))
table(timepoints)

#quality control
fivenum(apply(dataPBMC,1,function(x) sum(x>0)))
boxplot(apply(dataPBMC,1,function(x) sum(x>0)))
#75% of genes only express in 207 cells 
fivenum(apply(dataPBMC,2,function(x) sum(x>0)))
hist(apply(dataPBMC,2,function(x) sum(x>0)))
#75% of cells express less than 481 genes in total

#create seurat object
PBMC<- CreateSeuratObject(dataPBMC,
                          min.cells = 1,
                          min.features = 0,
                          project = '10x_PBMC')
PBMC


#add metadata
PBMC<-AddMetaData(PBMC,
                  metadata = apply(raw_dataPBMC,2,sum),
                  col.name = 'nUMI_raw')
PBMC<-AddMetaData(PBMC,
                  metadata = timepoints,
                  col.name = 'TimePoint')

#scale data to normalize
PBMC<-ScaleData(PBMC,
                vars.to.regress = c('nUMI_raw'),
                model.use = 'linear',
                use.umi = F)

#find high variable genes
#PBMC <- FindVariableGenes(object = PBMC, 
                          mean.function = ExpMean, 
                          dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, 
                          x.high.cutoff = 3,
                          y.cutoff = 0.5)

#more proper for seurat v3 or higher version
PBMC <- FindVariableFeatures(object = PBMC, 
                             mean.function = ExpMean, 
                             dispersion.function = LogVMR, 
                             mean.cutoff = c(0.0125,3), 
                             dispersion.cutoff = c(0.5,Inf))

#use high variable genes for dimension reduction
PBMC <- RunPCA(object = PBMC, pc.genes = PBMC@var.genes)

#find clusters based on the results of pca
PBMC<-RunUMAP(object = PBMC,
              reduction = 'pca',
              dims = 1:50)

#PBMC <- FindClusters(object = PBMC, 
                     reduction.type = "pca", 
                     dims.use = 1:10, 
                     resolution = 1, 
                     k.param = 35, 
                     save.SNN = TRUE) # 13 clusters

#more proper for seurat 3 or higher version
PBMC <- FindNeighbors(PBMC,
                      reduction = "pca",
                      dims = 1:10,
                      k.param = 35)

PBMC <- FindClusters(PBMC, 
                     resolution = 0.9, 
                     verbose=F) 

#use tsne to reduce dimension
PBMC <- RunTSNE(object = PBMC, 
                dims.use = 1:10)
#visualization
TSNEPlot(PBMC, 
         cols = c('green4', 'pink', '#FF7F00', 'orchid', 
                        '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 
                        'grey60', 'grey', 'red', '#FB9A99', 'black'))

#more proper for seurat 3 or higher version
colp <- c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 
        'dodgerblue2', 'grey30', 'yellow', 'grey60', 
        'grey', 'red', '#FB9A99', 'black') 

DimPlot(PBMC, 
        cols = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 
                       'dodgerblue2', 'grey30', 'yellow', 'grey60', 
                       'grey', 'red', '#FB9A99', 'black'))


save(PBMC,file = 'patient1.PBMC.output.Rdata')

table(PBMC@meta.data$res.1)
class(PBMC@meta.data$res.1)
table(as.numeric(PBMC@meta.data$res.1))

TSNEPlot(PBMC, group.by = "TimePoints")

#table(PBMC@meta.data$TimePoint,PBMC@ident)
table(PBMC@meta.data$TimePoint,Idents(PBMC))
PBMC_pre@active.ident

allGenes = row.names(PBMC@raw.data)
markGenes<-c('CD3D','CD3E','TRAC','IL7R','GZMA','FCGR3A',
             'CD14','MS4A1','FCER1A')

FeaturePlot(PBMC,
            features=markGenes,
            cols=c('grey','blue'),
            reduction='tsne')

a<-read.table('celltype-patient1-pbmc.txt')
new.cluster.ids <- as.character(a[,2])
names(new.cluster.ids)<-levels(PBMC)
PBMC<-RenameIdents(PBMC,new.cluster.ids)

DimPlot(PBMC,
        reduction = 'tsne',
        label = T,
        pt.size = 0.5,
        cols = colp)+NoLegend()


timepoints<-PBMC@meta.data$TimePoint
table(timepoints)

PBMC_pre<-subset(PBMC,TimePoint=='PBMC_RespD376')
table(Idents(PBMC_pre))


PBMC_pre_deg<-subset(PBMC_pre,
                     PBMC_pre@active.ident %in% c('CD8+ effector T cells ',
                                             'CD8+ cytotoxic T cells'))
#None of the requested variables were found:
#table get 201 Trues
BiocManager::install('monocle')
library(monocle)
#droped for now

start_time<-Sys.time()
raw_datatumor<-read.csv('GSE117988_raw.expMatrix_Tumor.csv.gz',
                        header = T,row.names = 1)
end_time<-Sys.time()
end_time-start_time
dim(raw_datatumor)

#
datatumor<-log2(1+sweep(raw_datatumor,2,
                        median(colSums(raw_datatumor))/colSums(raw_datatumor),
                        '*'))
head(colnames(datatumor))

#
celltype<-sapply(colnames(datatumor),function(x) unlist(strsplit(x,'\\.'))[2])
#for seurat v3
#cellTypes <- sapply(colnames(datatumor), function(x) ExtractField(x, 2, '[.]'))
#for seurat v2

celltype<-ifelse(celltype=='1','tumor_brfore','tumor_acquiredResistance')
table(celltype)

#
fivenum(apply(datatumor,1,function(x) sum(x>0)))
#75% genes are expressed on 566 cells 
fivenum(apply(datatumor,2,function(x) sum(x>0)))
#75% cells express 1971 genes

tumor<-CreateSeuratObject(datatumor,
                          min.cells = 1,
                          min.features = 0,
                          project = '10x_tumor')
tumor








