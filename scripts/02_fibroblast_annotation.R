table(pbmc$celltype)
Idents(pbmc) <- 'celltype'
bcell.filt=subset(pbmc,idents = 'Fibroblast')
table(bcell.filt$group)
library(harmony)
set.seed(15764)
table(bcell.filt$orig.ident)
if(T){
  source('./scRNA_scripts/harmony.R')
  sce.all_int = run_harmony(bcell.filt)
}
print(dim(bcell.filt))
print(dim(sce.all_int))
sce.all_int <- readRDS('sce.all_int.rds')
sce.all_int@meta.data$seurat_clusters <- sce.all_int@meta.data$RNA_snn_res.0.2
Idents(sce.all_int) <- "seurat_clusters"
DimPlot(sce.all_int, reduction = "umap")
set.seed(15865)
library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check
DotPlot(sce.all_int, features = unique(genes_to_check),
        assay='RNA' )  + coord_flip()
set.seed(1566)
#ggsave('check_last_markers.pdf',height = 11,width = 11)
mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)
umap =DimPlot(sce.all_int, reduction = "umap",cols = mycolors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.2",label = T,label.box = T)
umap
ggsave('group_umap_endo_sub.pdf',width = 8,height = 7)
tsne =DimPlot(sce.all_int, reduction = "tsne",cols = mycolors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.1",label = T,label.box = T)
tsne
ggsave('group_tsne_endo_sub.pdf',width = 8,height = 7)
set.seed(1524567)
#####细胞生物学命名
celltype=data.frame(ClusterID=0:11,
                    celltype= 0:11) 
celltype[celltype$ClusterID %in% c(0,1),2]='CAF_MYOC'
celltype[celltype$ClusterID %in% c(2),2]='CAF_MMP11'
celltype[celltype$ClusterID %in% c(3,8),2]='Fibroblast_CCL14'
celltype[celltype$ClusterID %in% c(4),2]='CAF_PCPE2'
celltype[celltype$ClusterID %in% c(5),2]='CAF_IL1B'
celltype[celltype$ClusterID %in% c(6),2]='CAF_CD24'
celltype[celltype$ClusterID %in% c(7),2]='Fibroblast_IGHG2'
celltype[celltype$ClusterID %in% c(9),2]='CAF_UBE2C'
celltype[celltype$ClusterID %in% c(10),2]='Fibroblast_DES'
celltype[celltype$ClusterID %in% c(11),2]='Fibroblast_NCAM'
table(sce.all_int@meta.data$RNA_snn_res.0.2)
table(sce.all_int$celltype)
sce.all_int@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all_int@meta.data[which(sce.all_int@meta.data$RNA_snn_res.0.2 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all_int@meta.data$celltype)
th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 

Idents(sce.all_int)<-sce.all_int$celltype
DimPlot(sce.all_int, reduction = "umap", group.by = "celltype",label = TRUE, pt.size = 0.5) + NoLegend()
library(patchwork)
set.seed(1568)
celltype_umap =DimPlot(sce.all_int, reduction = "umap",cols = mycolors,pt.size = 0.3,
                       group.by = "celltype",label = F)
celltype_umap
DimPlot(sce.all_int, reduction = "tsne",cols = mycolors,pt.size = 0.3,
        group.by = "celltype",label = T)
