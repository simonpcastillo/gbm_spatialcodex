## Author: Simon Castillo spcastillo@mdanderson.org
## Last modified (ddmmyyyy): 13092024

pacman::p_load(dplyr, raster, spatstat)

load('step0_output-qc/dataafterQC_CTIT_phenotyped.rdata')
if(!dir.exists('step1_output-spatial')) dir.create('step1_output-spatial')

alldata_afterQC$path_region %>% table()
CT_pcf= data.frame()

pacman::p_load(foreach,doParallel,parallel, dplyr, reshape2,stringr,ClusterR)



for(sam in unique(alldata_afterQC$sample)){
  alldata_afterQC %>%
    filter(sample==sam) -> d0
  
  for(reg in unique(d0$region)){
    iw_masks <- raster(Sys.glob(paste0('data/masks/', sam,'/*' , reg ,'-labels.png')))
    m0=as.matrix(iw_masks)
    m0 = m0 ==3 #only cellular tumour
    m0 <- apply(m0, 2, rev)
    
    if(sum(m0)==0) next
    d0 %>%
      filter(region== reg, path_region == 'Cellular tumor') -> d1
    
    pp1=ppp(d1$x, d1$y, window=owin(mask = m0), marks = as.factor(d1$phenotype))
    comlevels = combn(unique(d1$phenotype)[!unique(d1$phenotype) %in% 'Other'], 2)
    print(paste('Processing ', sam, reg, sep='...'))
    
    cl <- makePSOCKcluster(3)
    registerDoParallel(cl) 
    clusterEvalQ(cl, .libPaths("/home/spcastillo/R/ubuntu/4.3.1"))
    
    reg0 = foreach (j = 1:(length(comlevels)/2), .combine= rbind) %dopar%{
      require(dplyr); require(raster); require(spatstat)
      levelstocomp = comlevels[, j]
      lvl1 = levelstocomp[1]
      lvl2 = levelstocomp[2]
      pcfcross1 <- pcfcross(pp1, lvl1, lvl2, stoyan=0.1,correction = 'translation', r = 0:300)
      
      pcfcross1 -> p2
      
      p2$radius= p2$r * 0.3774 #from pixels to um
      max(pcfcross1[!pcfcross1$trans %in% c(Inf, NaN, NA),]$trans) -> Gmax
      min(pcfcross1[!pcfcross1$trans %in% c(Inf, NaN, NA),]$trans) -> Gmin
      mean(p2[!p2$trans %in% c(Inf, NaN, NA) & p2$radius> 5 & p2$radius < 10,]$trans) -> Gmean10
      mean(p2[!p2$trans %in% c(Inf, NaN, NA) & p2$radius> 5 & p2$radius < 50,]$trans) -> Gmean50
      mean(p2[!p2$trans %in% c(Inf, NaN, NA) & p2$radius> 5 & p2$radius < 100,]$trans) -> Gmean100
      
      
      p2 = p2[!p2$trans%in% c(Inf, NaN, NA),]
      p2[p2$trans == Gmax,]$radius -> Radmax
      p2[p2$trans == Gmin,]$radius -> Radmin
      
      df_out = data.frame(sample = sam, region = reg, path_region = 'CT',
                          pair= paste0(levelstocomp, collapse = '-vs-'),
                          Gmax= Gmax, Dmax = Radmax, Gmin = Gmin, Dmin = Radmin,
                          Gmean10=Gmean10, Gmean50=Gmean50, Gmean100 = Gmean100)
      df_out
      write.csv(df_out, paste0('step1_output-spatial/',paste(sam, reg, paste0(levelstocomp, collapse = '-vs-'), sep='_'), '.csv'))
      
    }
    stopCluster(cl) 
    
    CT_pcf = rbind(CT_pcf, reg0)
    
  }
  
  
}
