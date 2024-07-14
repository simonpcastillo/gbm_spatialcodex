require(dplyr); require(raster); require(spatstat)

pacman::p_load(foreach,doParallel,parallel, dplyr, reshape2,stringr,ClusterR)
## Author: Simon Castillo spcastillo@mdanderson.org
## Last modified (ddmmyyyy): 14092024


#BiocManager::install("BiocNeighbors")
require(BiocNeighbors)
library(data.table)

load("/data/dataafterQC_CTIT_phenotyped_20240311.rdata")

markers<- c('Olig2','Nestin','CD45', 'CD3e', 'CD68')
normalize<- TRUE


links_all1 = data.frame()
links_all10 = data.frame()
links_all30 = data_frame()
features_all = data.frame()

ss_all_gat = phenotyped_data %>% filter(acquisition_id %in% c("s338_c003_v002_r001_reg008", 's338_c005_v001_r001_reg005','s338_c006_v001_r001_reg001'))

for(sam in unique(ss_all_gat$sample)){
  print(sam)
  ss_all_gat %>%
    filter(sample==sam) -> d0
  
  for(reg in unique(d0$region)){
    d0 %>% filter(region == reg, path_region == 'Cellular tumor') -> d1
    d2= d1 %>% dplyr::select(study_id, acquisition_id, x, y, size)
    
    for (j in markers) {
      m0 = d1[, j]
      if(normalize == FALSE){
        dat = data.frame(m0)
      }
      
      if(normalize == TRUE){
        dat = data.frame((m0- mean(m0))/sd(m0))
      }
      colnames(dat)<- j
      d2 = cbind(d2, dat)
    }
    d2$phenotype = d1$phenotype
    
    others_df = d2 %>% filter(phenotype=='Other')
    others_df = others_df %>% sample_frac(size = 0.25, replace = F)
    others_df$phenotype2 = others_df$phenotype
    
    d2 %>% filter(phenotype!='Other') %>% 
      mutate(phenotype2 = case_when(phenotype %in% c('TCD4', 'TCD8')~'Lymphoid',
                                    phenotype %in% c('DC', 'Macrophage', 'Microglia', 'Monocyte') ~ 'Myeloid',
                                    phenotype == 'Tumor' ~ 'Tumor')
      )->nonothers_df
    
    if(nrow(nonothers_df)==0) next
    
    allreg = rbind(nonothers_df, others_df)
    allreg$cell_id <- 1:nrow(allreg)
    allreg$cell.id <- paste0(allreg$acquisition_id, '.', allreg$cell_id)
    
    k.out <- findNeighbors(allreg[, c('x', 'y')], threshold=26.497 * 3) # threshold= dis * 0.3774 #from pixels to um. 26.497 px ~ 10um
    id = allreg$cell.id
    names(k.out$index)<- id
    mm=sapply(k.out$index, function(x) unlist(str_split(x, pattern = ',')))
    nf=sapply(mm, function(x) length(x))
    linksdf30=data.frame(focal=rep(names(nf), nf),target=paste0(str_split(names(nf)[1], "[.]")[[1]][1],"." ,unname(unlist(mm))))
    
    linksdf30<- data.frame(unique(as.data.table(linksdf30)[, c("focal", "target") := list(pmin(focal, target),
                                                                                          pmax(focal, target))], by = c("focal", "target")))
    
    k.out <- findNeighbors(allreg[, c('x', 'y')], threshold=26.497) 
    id = allreg$cell.id
    names(k.out$index)<- id
    mm=sapply(k.out$index, function(x) unlist(str_split(x, pattern = ',')))
    nf=sapply(mm, function(x) length(x))
    linksdf10=data.frame(focal=rep(names(nf), nf),target=paste0(str_split(names(nf)[1], "[.]")[[1]][1],"." ,unname(unlist(mm))))
    
    linksdf10<- data.frame(unique(as.data.table(linksdf10)[, c("focal", "target") := list(pmin(focal, target),
                                                                                          pmax(focal, target))], by = c("focal", "target")))
    
    
    k.out <- findNeighbors(allreg[, c('x', 'y')], threshold=26.497*0.1)
    id = allreg$cell.id
    names(k.out$index)<- id
    mm=sapply(k.out$index, function(x) unlist(str_split(x, pattern = ',')))
    nf=sapply(mm, function(x) length(x))
    linksdf1=data.frame(focal=rep(names(nf), nf),target=paste0(str_split(names(nf)[1], "[.]")[[1]][1],"." ,unname(unlist(mm))))
    
    linksdf1<- data.frame(unique(as.data.table(linksdf1)[, c("focal", "target") := list(pmin(focal, target),
                                                                                        pmax(focal, target))], by = c("focal", "target")))
    
    
    
    features = allreg %>% dplyr::select(cell.id, size, all_of(markers), phenotype,phenotype2)
    
    
    write.csv(linksdf1, paste0('output_test/links_1um_',sam,"_",reg, ".csv"), row.names = F)
    write.csv(linksdf10, paste0('output_test/links_10um_',sam,"_",reg, ".csv"), row.names = F)
    write.csv(linksdf30, paste0('output_test/links_30um_',sam,"_",reg, ".csv"), row.names = F)
    write.csv(features, paste0('output_test/features_',sam,"_",reg, ".csv"), row.names = F)

  }

}

