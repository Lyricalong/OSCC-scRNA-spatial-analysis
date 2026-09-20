#devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
mycolors <-c('#e96153','#ef8a46' ,'#faa95a'  ,'#fed06e',
             '#fee7b5'  ,'#abdce0'  ,'#72bcd5','#528fac',
             '#386795','#1f466f'
)
table(scRNA$patient)
Idents(scRNA) <- scRNA$cluster
cellRatio<-cellRatioPlot(object = scRNA,                         
                         sample.name = "patient",                         
                         celltype.name = "cluster",                         
                         flow.curve = 0.5,fill.col = mycolors)+  
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
cellRatio
#width = 8,height = 6