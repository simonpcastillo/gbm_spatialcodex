## Author: Simon Castillo spcastillo@mdanderson.org
## Last modified (ddmmyyyy): 13092024

if(!require("pacman")) install.packages("pacman")
source('/misc/misc_functions.R')
pacman::p_load(foreach,doParallel,parallel, dplyr, reshape2,stringr,ClusterR, SpatialMap, raster)

if(!dir.exists('step0_output-qc')) dir.create('step0_output-qc')

#begins parallel computing
cl <- makePSOCKcluster(8)
registerDoParallel(cl) 
files = Sys.glob('data/processed_rdata_enablemedicine/*_unfiltered.rda')
clusterEvalQ(cl, .libPaths("/home/spcastillo/R/ubuntu/4.3.1"))

# Expression and coordinates for each cell passing QC
expression<- foreach (k=1:length(files), .combine = rbind) %dopar% {
  print(paste0('Processing ', k, ' of ', length(files)))
  pacman::p_load(foreach,doParallel,parallel, dplyr, reshape2, stringr, ClusterR)
  load(files[k])
  
  sample_id = str_extract(files[k], 'IW[P,R]\\d{2}')
  sample_data = data.frame()
  
  for (i in names(rawdata_sm@'regions')) {
    upper_z = 1.5
    lower_z = -1.5
    
    a<-rawdata_sm@'regions'[[i]]@'coordinates'
    df0 = rawdata_sm@'regions'[[i]]@'Data'
    df0bg = rawdata_sm@'regions'[[i]]@'bgData'
    
    # single-cell QC: signal intensity
    signalsum = colSums(df0)
    signalsum_modz = (signalsum- median(signalsum))/mad(signalsum) # modified z score
    filter_ssum = c(which(c(signalsum_modz > upper_z) %in% 'TRUE'),
                    which(c(signalsum_modz < lower_z) %in% 'TRUE'))
    filter_ssum_id = colnames(df0)[filter_ssum]
    
    # single-cell QC: cell size
    size = a$size
    size_modz = (size- median(size))/mad(size)
    filter_size = c(which(c(size_modz > upper_z) %in% 'TRUE'),
                    which(c(size_modz < lower_z) %in% 'TRUE'))
    filter_size_id = rownames(a)[filter_size]
    
    # single-cell QC: coefficient of variation in signal intensity
    cv=sapply(data.frame(df0), function(x){sd(as.numeric(x), na.rm=TRUE)/mean(as.numeric(x), na.rm=TRUE)})
    cv_modz = (cv- median(cv))/mad(cv)
    filter_cv = c(which(c(cv_modz < lower_z) %in% 'TRUE'))
    filter_cv_id = colnames(df0)[filter_cv]
    
    
    # single-cell QC: DAPI
    df0bg_dapi = df0bg[grep('DAPI',rownames(df0bg)),]
    dna = sapply(data.frame(df0bg_dapi), function(x) sum(x))
    dna_modz = (dna- median(dna))/mad(dna)
    filter_dna = c(which(c(dna_modz > upper_z) %in% 'TRUE'),
                   which(c(dna_modz < lower_z) %in% 'TRUE'))
    
    filter_dna_id = colnames(df0bg_dapi)[filter_dna]
    
    # single-cell QC: all cells ids that are considered outliers
    filter_out_cells = unique(c(filter_ssum_id, filter_dna_id, filter_size_id, filter_cv_id))
    
    df0_post_filter = df0[,!colnames(df0) %in% filter_out_cells ] #data without outliers
    df0_post_filter = data.frame(t(df0_post_filter))
    colids = colnames(df0_post_filter)
    rowids = rownames(df0_post_filter)
    df0_post_filter = as.data.frame(df0_post_filter)
    colnames(df0_post_filter)= colids
    rownames(df0_post_filter) = rowids
    
    df0_post_filter$region = i
    df0_post_filter$sample = sample_id
    sample_data = rbind(sample_data,df0_post_filter)
  }
  
  sample_data
}

coordinates<- foreach (k=1:length(files), .combine = rbind) %dopar% {
  print(paste0('Processing ', k, ' of ', length(files)))
  pacman::p_load(foreach,doParallel,parallel, dplyr, reshape2,
                 stringr)
  load(files[k])
  
  sample_id = str_extract(files[k], 'IW[P,R]\\d{2}')
  
  
  sample_data = data.frame()
  coordinates_data = data.frame()
  
  for (i in names(rawdata_sm@'regions')) {
    upper_z = 1.5
    lower_z = -1.5
    
    a<-rawdata_sm@'regions'[[i]]@'coordinates'
    df0 = rawdata_sm@'regions'[[i]]@'Data'
    df0bg = rawdata_sm@'regions'[[i]]@'bgData'
    
    # single-cell QC: signal intensity
    signalsum = colSums(df0)
    signalsum_modz = (signalsum- median(signalsum))/mad(signalsum)
    filter_ssum = c(which(c(signalsum_modz > upper_z) %in% 'TRUE'),
                    which(c(signalsum_modz < lower_z) %in% 'TRUE'))    
    filter_ssum_id = colnames(df0)[filter_ssum]
    
    # single-cell QC: cell size
    size = a$size
    size_modz = (size- median(size))/mad(size)
    filter_size = c(which(c(size_modz > upper_z) %in% 'TRUE'),
                    which(c(size_modz < lower_z) %in% 'TRUE'))
    filter_size_id = rownames(a)[filter_size]
    
    # single-cell QC: coefficient of variation in signal intensity
    cv=sapply(data.frame(df0), function(x){sd(as.numeric(x), na.rm=TRUE)/mean(as.numeric(x), na.rm=TRUE)})
    cv_modz = (cv- median(cv))/mad(cv)
    filter_cv = c(which(c(cv_modz < lower_z) %in% 'TRUE'))
    filter_cv_id = colnames(df0)[filter_cv]
    
    
    # single-cell QC: DAPI
    df0bg_dapi = df0bg[grep('DAPI',rownames(df0bg)),]
    dna = sapply(data.frame(df0bg_dapi), function(x) sum(x))
    dna_modz = (dna- median(dna))/mad(dna)
    filter_dna = c(which(c(dna_modz > upper_z) %in% 'TRUE'),
                   which(c(dna_modz < lower_z) %in% 'TRUE'))
    filter_dna_id = colnames(df0bg_dapi)[filter_dna]
    
    
    filter_out_cells = unique(c(filter_ssum_id, filter_dna_id, filter_size_id, filter_cv_id)) 
    coords_post_filter = a[!rownames(a) %in% filter_out_cells,]
    coordinates_data = rbind(coordinates_data, coords_post_filter)
  }
  
  coordinates_data
}

stopCluster(cl) # ends parallel

merged_afterQC = cbind(coordinates, expression); merged_afterQC=merged_afterQC[,-57]
merged_afterQC= merged_afterQC %>% mutate(region = paste0('reg',str_pad(region, 3, pad='0'))) 

merged_path = data.frame()
correctxy = read.csv('data/correctioncellpos.csv') # cell coordinate correction necessary to match path mask and cell coordinates
num_ids = do.call(paste0, expand.grid( c('P','R'), str_pad(5:11, 2, pad='0')))


for(num_id in num_ids){
  iw_masks <- Sys.glob(paste0('data/masks/','IW',num_id, '/*.png'))
  
  for(i in 1:length(iw_masks)){
    t0 = raster(iw_masks[i])
    reg = str_extract(iw_masks[i], 'reg\\d{3}')
    sam0 = paste0('IW',num_id)
    print(paste(sam0, reg, sep='-'))
    
    dfx<- merged_afterQC[merged_afterQC$sample == sam0 & merged_afterQC$region == reg,]
    
    corrdf = correctxy[correctxy$sample==sam0,]
    offsety  =  as.numeric(corrdf[corrdf$region==reg,]$offsety)
    offsetx = as.numeric(corrdf[corrdf$region==reg,]$offsetx)
    
    if(sam0 %in% c('IWP05', 'IWP09','IWP10', 'IWP11', 'IWP08', 'IWP07')){
      dfx %>% mutate(y =offsety - y) -> dfx
    }
    if(sam0 %in% c('IWR05', 'IWR07', 'IWR10')){
      dfx %>% mutate(y =offsety + y) -> dfx
    }
    if(sam0 %in% c('IWP06', 'IWR08')){
      dfx %>% mutate(x = offsetx - x, y = offsety - y) -> dfx
    }
    
    if(sam0 %in% c('IWR06','IWR09', 'IWR11')){
      dfx %>% mutate(x2 = x, y2 = y, x = offsetx - y2, y = x2) %>% dplyr::select(-c(y2, x2)) -> dfx
    }
    
    dfx %>% mutate(y = dim(t0)[1]-y) -> dfx
    dfx$in_region = extract(t0, dfx[,c('x', 'y')]) #extract pathology region for each cell and assign as new cell feature
    
    dfx=dfx %>%
      mutate(path_region = case_when( in_region == '1' ~ 'Normal',
                                      in_region == '2' ~ 'Infiltrating tumor',
                                      in_region == '3' ~ 'Cellular tumor',
                                      in_region == '4' ~ 'Necrosis',
                                      in_region == '5' ~ 'Necrosis_pp',
                                      in_region == '6' ~ 'Necrosis_treatment',
                                      in_region == '7' ~ 'Depopulated tumor',
                                      in_region == '255' ~ 'background')) 
    merged_path = rbind(merged_path, dfx)
  }
  
}

targetmarkers = colnames(merged_path)[15:56]
markerspos_before = list()
markersll = list()
markerspos_after2 = data.frame()

for (aqid in unique(merged_path$acquisition_id)) {
  markerspos_after = list()
  
  for (j in targetmarkers) {
    print(j)
    m0 = merged_path[merged_path$acquisition_id == aqid, j]
    m0 = m0[!is.na(m0)]
    dat = data.frame((m0- mean(m0))/sd(m0)) # zscore normalization
    gmm = GMM(dat, 2, dist_mode = "eucl_dist", "random_subset", 100, 100)
    pr = predict(gmm, newdata = dat)
    dx = data.frame(numcluster = pr) %>% 
      mutate(cluster = case_when(numcluster == match(max(gmm$centroids), gmm$centroids) ~ TRUE, # marker positive: TRUE
                                 numcluster != match(max(gmm$centroids), gmm$centroids) ~ FALSE )  # marker positive: FALSE
      )
    lldf = data.frame(gmm$Log_likelihood)
    names(dx)[2] <- paste0(j, 'pos')
    markerspos_after[[paste0(j, 'pos')]] = dx[,2]
    
  }
  markerspos_after = data.frame(markerspos_after)
  markerspos_after2 = rbind(markerspos_after2, markerspos_after)
}

markerspos_after = data.frame(markerspos_after2)

markerspos_after%>%
  mutate(
    TCD4 = case_when(CD45pos == TRUE & CD3epos == TRUE & CD4pos == TRUE ~ TRUE),
    TCD8 = case_when(CD45pos == TRUE & CD3epos == TRUE & CD8pos == TRUE ~ TRUE),
    Monocyte = case_when(CD3epos == FALSE & CD14pos == TRUE & CD16pos == TRUE & CD45pos == TRUE ~ TRUE),
    Macrophage =  case_when(CD3epos == FALSE& CD45pos == TRUE & CD68pos == TRUE ~ TRUE), 
    DC = case_when(CD3epos == FALSE & CD45pos == TRUE & CD11cpos == TRUE & CD370pos == TRUE & CD1cpos == TRUE ~ TRUE), 
    Microglia = case_when(CD3epos == FALSE & CD45pos == FALSE & CX3CR1pos == TRUE ~ TRUE),
    Tumor = case_when(CD45pos == FALSE & Olig2pos == TRUE | CD45pos == FALSE & Nestinpos == TRUE ~ TRUE)) -> cellphenotypes_after


cellphenotypes_after[is.na(cellphenotypes_after)]<- FALSE

cells_after=apply(cellphenotypes_after[, 43:49], 1, function(x) names(x)[which(x %in% TRUE)]) 
lng_a=sapply(cells_after, function(x) length(x))
multiph_a=which(lng_a %in% c(2,3))
phenotyped_data = cbind(merged_path, cellphenotypes_after)

phenotyped_data$multiph <- NA
phenotyped_data[multiph_a,]$multiph <-'Yes'
phenotyped_data[is.na(phenotyped_data$multiph),]$multiph <-'No'
cells_after = apply(cellphenotypes_after[, 43:49], 1, function(x) names(x)[match(TRUE, x)]) 
phenotyped_data$phenotype = cells_after
phenotyped_data[is.na(phenotyped_data$phenotype),]$phenotype <-'Other'
phenotyped_data[phenotyped_data$TCD4 == TRUE & alldata_afterQC$TCD8 == TRUE,]$multiph = FALSE


save(phenotyped_data,file=paste0('step0_output-qc/dataafterQC_CTIT_phenotyped.rdata')) 
