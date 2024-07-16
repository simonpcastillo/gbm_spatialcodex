path_data= '/Users/spcastillo/Downloads/gbm_fromseadragon_20231109/output_test'
list_pred = paste(path_data, list.files(path_data, pattern = "feat*"), sep = '/')
links_all = data.frame()
for (j in list_pred) {
  df0 = read.csv(j)[, c('cell.id', 'phenotype2')]
  roiid = str_remove(str_split(j, '/')[[1]][length(str_split(j, '/')[[1]])], 'features_')
  links0 = read.csv(paste(path_data, paste0('links_30um_', roiid), sep = '/'))
  merged1=merge(data.frame(cell.id=links0[,1]), df0, by='cell.id'); names(merged1) = c('id1', 'ph1')
  merged2=merge(data.frame(cell.id=links0[,2]), df0, by='cell.id'); names(merged2) = c('id2', 'ph2')
  phlinks = cbind(merged1, merged2)
  phlinks$links = paste(phlinks$ph1, phlinks$ph2, sep='-')
  phlinks$loop = phlinks$id1 == phlinks$id2
  phlinks = phlinks[phlinks$loop==FALSE,]
  phlinks$time= str_extract(roiid, "[P,R]")
  print(str_extract(roiid, "[P,R]"))
  phlinks$roi= str_remove(roiid, ".csv")
  links_all=rbind(links_all, phlinks)
}

links_p = links_all[links_all$time=='P',]

lt = nrow(links_p[links_p$links %in% c('Lymphoid-Tumor', 'Tumor-Lymphoid'),])/ nrow(links_p)
mt = nrow(links_p[links_p$links %in% c('Myeloid-Tumor', 'Tumor-Myeloid'),])/ nrow(links_p)
ot= nrow(links_p[links_p$links  %in% c('Other-Tumor','Tumor-Other'),])/ nrow(links_p)
ml = nrow(links_p[links_p$links  %in% c('Myeloid-Lymphoid','Lymphoid-Myeloid'),])/ nrow(links_p)
ol = nrow(links_p[links_p$links  %in% c('Other-Lymphoid','Lymphoid-Other'), ])/ nrow(links_p)
om = nrow(links_p[links_p$links  %in% c('Other-Myeloid','Myeloid-Other'),])/ nrow(links_p)
ll = nrow(links_p[links_p$links  %in% c('Lymphoid-Lymphoid'),])/ nrow(links_p)
mm = nrow(links_p[links_p$links  %in% c('Myeloid-Myeloid'),])/ nrow(links_p)
tt = nrow(links_p[links_p$links  %in% c('Tumor-Tumor'),])/ nrow(links_p)
oo = nrow(links_p[links_p$links  %in% c('Other-Other'),])/ nrow(links_p)

nodes = data.frame(id=1:4, label = c('Tumor', 'Lymphoid', 'Myeloid', 'Other'))
edgesP = data.frame(from=c(1,1,1,1,2,2,2,3,3,4), to=c(1,2,3,4,2,3,4,3,4,4), weight = c(tt,lt, mt, ot,ll, ml, ol,mm, om, oo)); edgesP$sign = 1


links_r = links_all[links_all$time=='R',]

lt = nrow(links_r[links_r$links %in% c('Lymphoid-Tumor', 'Tumor-Lymphoid'),])/ nrow(links_r)
mt = nrow(links_r[links_r$links == 'Myeloid-Tumor' | links_r$links =='Tumor-Myeloid',])/ nrow(links_r)
ot= nrow(links_r[links_r$links == 'Other-Tumor' |links_r$links =='Tumor-Other',])/ nrow(links_r)
ml = nrow(links_r[links_r$links == 'Myeloid-Lymphoid' | links_r$links =='Lymphoid-Myeloid',])/ nrow(links_r)
ol = nrow(links_r[links_r$links == 'Other-Lymphoid' | links_r$links =='Lymphoid-Other', ])/ nrow(links_r)
om = nrow(links_r[links_r$links == 'Other-Myeloid' | links_r$links =='Myeloid-Other',])/ nrow(links_r)
ll = nrow(links_r[links_r$links  %in% c('Lymphoid-Lymphoid'),])/ nrow(links_r)
mm = nrow(links_r[links_r$links  %in% c('Myeloid-Myeloid'),])/ nrow(links_r)
tt = nrow(links_r[links_r$links  %in% c('Tumor-Tumor'),])/ nrow(links_r)
oo = nrow(links_r[links_r$links  %in% c('Other-Other'),])/ nrow(links_r)

nodes = data.frame(id=1:4, label = c('Tumor', 'Lymphoid', 'Myeloid', 'Other'))
edgesR = data.frame(from=c(1,1,1,1,2,2,2,3,3,4), to=c(1,2,3,4,2,3,4,3,4,4), weight = c(tt,lt, mt, ot,ll, ml, ol,mm, om, oo)); edgesR$sign = -1

edges = rbind(edgesP, edgesR)



lt = nrow(links_p[links_p$links %in% c('Lymphoid-Tumor', 'Tumor-Lymphoid'),])//nrow(links_r[links_r$links %in% c('Lymphoid-Tumor', 'Tumor-Lymphoid'),])
mt = nrow(links_r[links_r$links == 'Myeloid-Tumor' | links_r$links =='Tumor-Myeloid',])
ot= nrow(links_r[links_r$links == 'Other-Tumor' |links_r$links =='Tumor-Other',])/ nrow(links_r)
ml = nrow(links_r[links_r$links == 'Myeloid-Lymphoid' | links_r$links =='Lymphoid-Myeloid',])/ nrow(links_r)
ol = nrow(links_r[links_r$links == 'Other-Lymphoid' | links_r$links =='Lymphoid-Other', ])/ nrow(links_r)
om = nrow(links_r[links_r$links == 'Other-Myeloid' | links_r$links =='Myeloid-Other',])/ nrow(links_r)
ll = nrow(links_r[links_r$links  %in% c('Lymphoid-Lymphoid'),])/ nrow(links_r)
mm = nrow(links_r[links_r$links  %in% c('Myeloid-Myeloid'),])/ nrow(links_r)
tt = nrow(links_r[links_r$links  %in% c('Tumor-Tumor'),])/ nrow(links_r)
oo = nrow(links_r[links_r$links  %in% c('Other-Other'),])/ nrow(links_r)

nodes = data.frame(id=1:4, label = c('Tumor', 'Lymphoid', 'Myeloid', 'Other'))
edgesRatio = data.frame(from=c(1,1,1,1,2,2,2,3,3,4), to=c(1,2,3,4,2,3,4,3,4,4), weight = c(tt,lt, mt, ot,ll, ml, ol,mm, om, oo)); edgesR$sign = -1



edgesRatio = data.frame(edgesP, weightR=edgesR$weight)
edgesRatio$weightratio = edgesRatio$weight /edgesRatio$weightR

library(tidygraph)
net.tidy.p <- tbl_graph(
  nodes = nodes, edges = edgesRatio, directed = T
)

require(ggraph)
svglite('plots/interaction_freq_legend.svg', width = 3.5, height = 1.5)
print(ggraph(net.tidy.p, layout = "linear") + 
  geom_edge_arc(aes(width = weightratio, color=weightratio<1), alpha = 0.8, strength=-1) + 
  #scale_edge_width(limits = c(0, 10),range = c(0,3)) +
  geom_edge_loop(aes(width = weightratio,color=weightratio<1, direction=0, span=30),alpha = 0.8, force_flip = F)+ 
  scale_edge_color_manual(values = c('#3DC2C1', '#C23D3E'))+
  #scale_edge_color_gradient(low="gray", high="#3DC2C1",limits=c(0, 10)) +
  geom_node_text(aes(label = label), repel = F) +
  labs(edge_width = "Fold change", edge_color ='pGBM < rGBM') +
  theme_graph()
)
dev.off()
