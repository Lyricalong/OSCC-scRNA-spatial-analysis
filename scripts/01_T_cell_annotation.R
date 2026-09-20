
table(pbmc$celltype)
Idents(pbmc) <- 'celltype'
bcell.filt=subset(pbmc,idents = 'T')
library(harmony)
set.seed(15264)
table(bcell.filt$orig.ident)
if(T){
  source('./scRNA_scripts/harmony.R')
  sce.all_int = run_harmony(bcell.filt)
}
print(dim(bcell.filt))
print(dim(sce.all_int))
sce.all_int <- readRDS('sce.all_int.rds')
sce.all_int@meta.data$seurat_clusters <- sce.all_int@meta.data$RNA_snn_res.0.5
Idents(sce.all_int) <- "seurat_clusters"
DimPlot(sce.all_int, reduction = "umap")
ggsave('tsne_by_sce.all_int.pdf',width = 8,height = 7)
set.seed(15265)
genes_to_check = c(
  "CD4", #T_check_CD4
  "CCR7" , #Naive_CD4_T
  "RORA", #Memory_CD4
  "FOXP3","CTLA4", #Treg_CD4
  "CXCL13", #Th1_CD4
  "BCL6", #Tfh_CD4
  'CD8A','CD8B', #CD8
  'GZMK', 'MKI67', #CD8
  'EOMES','TOX'
) 
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
              group.by = "RNA_snn_res.0.5",label = T,label.box = T)
umap
ggsave('group_umap_endo_sub.pdf',width = 8,height = 7)
tsne =DimPlot(sce.all_int, reduction = "tsne",cols = mycolors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.1",label = T,label.box = T)
tsne
ggsave('group_tsne_endo_sub.pdf',width = 8,height = 7)
set.seed(1524567)
celltype=data.frame(ClusterID=0:12,
                    celltype= 0:12) 
celltype[celltype$ClusterID %in% c(0),2]='CD8_Tex_GZMK'
celltype[celltype$ClusterID %in% c(1),2]='CD8_Tn_CCR7'
celltype[celltype$ClusterID %in% c(2),2]='CD8_CXCL13'
celltype[celltype$ClusterID %in% c(3,4),2]='CD4_Treg_FOXP3'
celltype[celltype$ClusterID %in% c(5),2]='CD4_Tm_RORA'
celltype[celltype$ClusterID %in% c(8),2]='CD4_Tfh_BCL6'
celltype[celltype$ClusterID %in% c(7,9),2]='CD8_MKI67'
celltype[celltype$ClusterID %in% c(6,10),2]='CD8_Tm_GZMK'
celltype[celltype$ClusterID %in% c(11,12),2]='CD4_MKI67'
table(sce.all_int@meta.data$RNA_snn_res.0.5)
table(sce.all_int$celltype)
sce.all_int@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all_int@meta.data[which(sce.all_int@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all_int@meta.data$celltype)
th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 

Idents(sce.all_int)<-sce.all_int$celltype
DimPlot(sce.all_int, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
library(patchwork)

set.seed(1568)
celltype_umap =DimPlot(sce.all_int, reduction = "umap",cols = mycolors,pt.size = 0.3,
                       group.by = "celltype",label = T)
celltype_umap
ggsave('umap_by_endo_celltype.pdf',width = 8,height = 7)
DimPlot(sce.all_int, reduction = "tsne",cols = mycolors,pt.size = 0.3,
        group.by = "celltype",label = T)
ggsave('tsne_by_endo_celltype.pdf',width = 8,height = 7)
set.seed(14562)
source('./scRNA_scripts/mycolors.R')
celltype_umap =DimPlot(pbmc, reduction = "umap",cols = mycolors,pt.size = 0.3,
                       group.by = "celltype",label = T)
celltype_umap
table(sce.all_int$celltype)
table(pbmc$celltype)
set.seed(12543570)
Idents(pbmc) = pbmc$celltype
pbmc$celltype = as.character(Idents(pbmc))
pbmc$celltype = ifelse(pbmc$celltype=="T",
                       sce.all_int$celltype[match(colnames(pbmc),colnames(sce.all_int))],
                       pbmc$celltype)
DimPlot(pbmc, reduction = "umap",cols = mycolors,pt.size = 0.3,
        group.by = "celltype",label = T)
DimPlot(pbmc, reduction = "tsne",cols = mycolors,pt.size = 0.3,
        group.by = "celltype",label = T)
table(pbmc$celltype)
table(pbmc$group)