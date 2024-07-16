## Author: Simon Castillo spcastillo@mdanderson.org
## Last modified (ddmmyyyy): 13092024

pacman::p_load(igraph, dplyr, stringr)
if(!dir.exists('step2_output-network')) dir.create('step2_output-network')


f0 = Sys.glob('step1_output-spatial/*.csv')

data_pcf = data.frame()
for(i in 1:length(f0)){
  df0 = read.csv(f0[i])
  data_pcf = rbind(data_pcf, df0)
}

data_pcf %>%
  mutate(tempid = paste(sample, region, path_region, pair, sep='---')) %>%
  group_by(tempid) %>%
  summarise( mGmean30= mean(Gmean30),
            mGmean100= mean(Gmean100) )%>%
  ungroup() -> data_pcf2


data_pcf2$sample= lapply(str_split(data_pcf2$tempid,pattern= '---'), function(x) x[1])
data_pcf2$region= lapply(str_split(data_pcf2$tempid,pattern= '---'), function(x) x[2])
data_pcf2$path_region= lapply(str_split(data_pcf2$tempid,pattern= '---'), function(x) x[3])
data_pcf2$pair= lapply(str_split(data_pcf2$tempid,pattern= '---'), function(x) x[4])
data_pcf2$cell1= lapply(str_split(data_pcf2$pair,pattern= '-vs-'), function(x) x[1])
data_pcf2$cell2= lapply(str_split(data_pcf2$pair,pattern= '-vs-'), function(x) x[2])

data_pcf2 %>%
  dplyr::select(sample, region, path_region, pair, cell1, cell2,mGmean30, mGmean100) -> data_pcf2

pathreg= 'CT'
metrics_net = data.frame()
metrics_node = data.frame()

samples = unlist(unique(data_pcf2$sample))
for(sam in samples){
  for(reg in unique(data_pcf2[data_pcf2$sample==sam,]$region)){
    for(d in c(30, 100)){
      print(reg)
      targetg = paste0('mGmean',d)
      reg_df = 
        data_pcf2 %>%
        filter(sample== sam, region== reg)
      reg_mat = matrix(data=NA, nrow = 7, ncol = 7)
      nods = c('DC', 'Microglia',  'Macrophage', 'Tumor', 'TCD4', 'TCD8', 'Monocyte')
      colnames(reg_mat) <- rownames(reg_mat)<- nods
      
      for(i in 1:nrow(reg_df)){
        target_cell1 = reg_df[i,]$cell1[[1]]
        target_cell2 = reg_df[i,]$cell2[[1]]
        mGmedval = unlist(reg_df[i, targetg])
        reg_mat[target_cell1,target_cell2]<-mGmedval
        reg_mat[target_cell2,target_cell1]<-mGmedval
      }
      
      reg_mat[reg_mat<1]<- 0
      reg_mat[is.na(reg_mat)]<- 0
      
      net <- graph_from_adjacency_matrix(reg_mat, weighted = T, diag = F,mode = 'upper' )
      clos = closeness(net, normalized = TRUE, mode = 'all')
      str = strength(net)
      bet = betweenness(net, directed = FALSE, normalized = TRUE)
      deg = degree(net)
      transiti = transitivity(net, type = 'weighted', isolates = 'zero')
      df00 = data.frame(names(V(net)), deg, clos, str, bet, transiti)
      df00$sample = sam
      df00$reg = reg
      df00$path_reg = pathreg
      df00$time = str_extract(sam, 'recurrent|primary')
      df00$d = d
      assor_degree = assortativity_degree(net)
      wtc <- cluster_walktrap(net)
      mod = modularity(wtc)
      trans = mean(transitivity(net, type = 'barrat'), na.rm=T)
      trans_uw = mean(transitivity(net, type = 'global'), na.rm=T)
      dens= edge_density(net)
      df00_net = data.frame(sam, reg, pathreg, d,assor_degree, mod, trans,trans_uw, dens)
      
      metrics_node = rbind(metrics_node, df00)
      metrics_net = rbind(metrics_net, df00_net)
      
    }
  }
}

metrics_net$time = str_extract(metrics_net$sam, 'recurrent|primary')
metrics_net$pat = str_extract(metrics_net$sam, '\\d{2}')


write.csv(metrics_net, 'step2_output-network/metrics_networks.csv')
