# Demo data
exam_scores <- data.frame(
  row.names = c("Student.1", "Student.2", "Student.3"),
  Biology = c(7.9, 3.9, 9.4),
  Physics = c(10, 20, 0),
  Maths = c(3.7, 11.5, 2.5),
  Sport = c(8.7, 20, 4),
  English = c(7.9, 7.2, 12.4),
  Geography = c(6.4, 10.5, 6.5),
  Art = c(2.4, 0.2, 9.8),
  Programming = c(0, 0, 20),
  Music = c(20, 20, 20)
)
exam_scores

library(fmsb)
# Define the variable ranges: maximum and minimum
max_min <- data.frame(
  Biology = c(20, 0), Physics = c(20, 0), Maths = c(20, 0),
  Sport = c(20, 0), English = c(20, 0), Geography = c(20, 0),
  Art = c(20, 0), Programming = c(20, 0), Music = c(20, 0)
)
rownames(max_min) <- c("Max", "Min")

# Bind the variable ranges to the data
df <- rbind(max_min, exam_scores)
df


create_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

BEAUTIFUL RADAR CHART IN R USING FMSB AND GGPLOT PACKAGES
HOMEDATA VISUALIZATIONBEAUTIFUL RADAR CHART IN R USING FMSB AND GGPLOT PACKAGES
Search

radar-chart-in-r-customized-fmstb-radar-chart-1.png



12 Dec
Beautiful Radar Chart in R using FMSB and GGPlot Packages
Alboukadel |  Circular Plot |  Data Visualization |  2
A radar chart, also known as a spider plot is used to visualize the values or scores assigned to an individual over multiple quantitative variables, where each variable corresponds to a specific axis.

This article describes how to create a radar chart in R using two different packages: the fmsb or the ggradar R packages.

Note that, the fmsb radar chart is an R base plot. The ggradar package builds a ggplot spider plot.

You will learn:
  
  how to create a beautiful fmsb radar chart
how to create ggplot radar chart
alternatives to radar charts



Contents:
  
  Demo data
fmsb radar chart
Prerequisites
Data preparation
Basic radar plot
Customize the radar charts
Create radar charts for multiple individuals
Compare every profile to an average profile
ggplot radar chart using the ggradar R package
Prerequisites
Key function and arguments
Data preparation
Basic radar plot
Customize radar charts
Radar chart with multiple individuals or groups
Alternatives to radar charts
Case when all quantitative variables have the same scale
Case when you have a lot of individuals to plot or if your variables have different scales
Conclusion
Demo data
We’ll use a demo data containing exam scores for 3 students on 9 topics (Biology, Physics, etc). The scores range from 0 to 20. Columns are quantitative variables and rows are individuals.

# Demo data
exam_scores <- data.frame(
  row.names = c("Student.1", "Student.2", "Student.3"),
  Biology = c(7.9, 3.9, 9.4),
  Physics = c(10, 20, 0),
  Maths = c(3.7, 11.5, 2.5),
  Sport = c(8.7, 20, 4),
  English = c(7.9, 7.2, 12.4),
  Geography = c(6.4, 10.5, 6.5),
  Art = c(2.4, 0.2, 9.8),
  Programming = c(0, 0, 20),
  Music = c(20, 20, 20)
)
exam_scores
##           Biology Physics Maths Sport English Geography Art Programming Music
## Student.1     7.9      10   3.7   8.7     7.9       6.4 2.4           0    20
## Student.2     3.9      20  11.5  20.0     7.2      10.5 0.2           0    20
## Student.3     9.4       0   2.5   4.0    12.4       6.5 9.8          20    20
fmsb radar chart
Prerequisites
Install the fmsb R package:
  
  install.packages("fmsb")
Load the package:
  
  library(fmsb)
Data preparation
The data should be organized as follow:
  
  The row 1 must contain the maximum values for each variable
The row 2 must contain the minimum values for each variable
Data for cases or individuals should be given starting from row 3
The number of columns or variables must be more than 2.
# Define the variable ranges: maximum and minimum
max_min <- data.frame(
  Biology = c(20, 0), Physics = c(20, 0), Maths = c(20, 0),
  Sport = c(20, 0), English = c(20, 0), Geography = c(20, 0),
  Art = c(20, 0), Programming = c(20, 0), Music = c(20, 0)
)
rownames(max_min) <- c("Max", "Min")

# Bind the variable ranges to the data
df <- rbind(max_min, exam_scores)

ggradar(
  df[3:5,], 
  values.radar = c("0", "10", "20"),
  grid.min = 0, grid.mid = 10, grid.max = 20,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#00AFBB", "#E7B800", "#FC4E07"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
)


for (samp in unique(alldata_afterQC$sample)){
  samp = 'IWR06'
df_radar = alldata_afterQC %>%
  filter(sample == samp, !phenotype %in% c('Olig2prolif', 'Olig2nonprolif', 'Other')) %>%
  group_by(path_region) %>%
  summarise(cd4t=length(phenotype[phenotype=='TCD4']),
            cd8t=length(phenotype[phenotype=='TCD8']),
            cd4cd8t = length(phenotype[phenotype=='TCD4CD8']),
            monocyte=length(phenotype[phenotype=='Monocyte']),
            microglia = length(phenotype[phenotype=='Microglia']),
            Mac=length(phenotype[phenotype=='Macrophage']),
            cd4Mac=length(phenotype[phenotype=='CD4Macrophage']), 
            DC = length(phenotype[phenotype=='DC']), 
            # tumprol= length(phenotype[phenotype=='Olig2prolif']),
            # tumnonprol= length(phenotype[phenotype=='Olig2nonprolif']),
            total = length(phenotype)) %>%
  mutate(across(cd4t:DC, ~ ./total)) %>%
  select(-c(total)) %>%
  ungroup() %>%
  mutate(path_region = case_when(path_region == 'Cellular tumor' ~ 'CT',
                                 path_region == 'Infiltrating tumor' ~ 'IT'))


svglite(file=paste0('plots/radar/immune_',samp, '.svg'),width = 3, height = 3)
print(ggradar(
  df_radar, 
  values.radar = c(""),
  grid.min = 0, grid.mid = 0.5, grid.max = 0.75,
  # Polygons
  group.line.width = 1, 
  group.point.size = 2,legend.position = 'bottom'
)
)
dev.off()
}



require(vegan)

?diversity

data(BCI)
diversity(df_radar[,-1])

divdf=data.frame()
for (samp in unique(alldata_afterQC$sample)){
  df_radar = alldata_afterQC %>%
    filter(sample == samp, !phenotype %in% c('Olig2prolif', 'Olig2nonprolif', 'Other')) %>%
    group_by(path_region) %>%
    summarise(cd4t=length(phenotype[phenotype=='TCD4']),
              cd8t=length(phenotype[phenotype=='TCD8']),
              cd4cd8t = length(phenotype[phenotype=='TCD4CD8']),
              monocyte=length(phenotype[phenotype=='Monocyte']),
              microglia = length(phenotype[phenotype=='Microglia']),
              Mac=length(phenotype[phenotype=='Macrophage']),
              cd4Mac=length(phenotype[phenotype=='CD4Macrophage']), 
              DC = length(phenotype[phenotype=='DC']), 
              # tumprol= length(phenotype[phenotype=='Olig2prolif']),
              # tumnonprol= length(phenotype[phenotype=='Olig2nonprolif']),
              total = length(phenotype)) %>%
    mutate(across(cd4t:DC, ~ ./total)) %>%
    select(-c(total)) %>%
    ungroup() %>%
    mutate(path_region = case_when(path_region == 'Cellular tumor' ~ 'CT',
                                   path_region == 'Infiltrating tumor' ~ 'IT'))
  
  
  dv0=data.frame(shannon=diversity(df_radar[,-1], 'simpson'))
  dv0$path_region = df_radar$path_region
  dv0$sample = samp
divdf = rbind(divdf, dv0)
}


svglite(file=paste0('plots/immunediversity.svg'),width = 4, height = 2.3)
print(
divdf %>%
  mutate(time = str_extract(sample, '[R,P]'))%>%
ggplot(aes(path_region, shannon)) +
  geom_boxplot(aes(fill= path_region))+
  labs(x='', y='Simpson\'s diversity')+
  scale_fill_manual(values = c("#FC4E07", "#E7B800"))+
  scale_y_continuous(limits = c(0,1))+
  geom_point() +
  facet_wrap(.~time)+
  stat_anova_test() +
  theme_minimal()+
  guides(fill = 'none')
)
dev.off()



selectdf %>%
  ggplot(aes(path_region, nMicroglia))+
  geom_boxplot(aes(fill=path_region))+  
  scale_fill_manual(values = c("#FC4E07", "#E7B800"))+
  geom_point() +
  facet_wrap(.~time)+
  stat_anova_test() +
  theme_minimal()


svglite(file=paste0('plots/Microglia.svg'),width = 3, height = 1.8)
print(
  selectdf %>%
    mutate(path_region = case_when(path_region == 'Cellular tumor' ~ 'CT',
                                   path_region == 'Infiltrating tumor' ~ 'IT'))%>%
    ggplot(aes(time, nMicroglia))+
    geom_boxplot(aes(fill=path_region))+  
    scale_fill_manual(values = c("#FC4E07", "#E7B800"))+
    geom_quasirandom(size=0.25) +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    labs(y= 'Microglia density', x='')+
    facet_wrap(.~path_region)+
    stat_anova_test() +
    theme_minimal()+
    theme(text = element_text(size=8))+
    guides(fill = 'none')
)
  dev.off()
  
  
  
  require(Hmisc); library(corrplot)
  immune<- selectdf %>%
    filter(path_region == 'Infiltrating tumor', time == 'R') %>%
    select(nTCD4:nOlig2Ki67neg)
  
  res = cor(immune)
  res2 <- rcorr(as.matrix(immune), )
  
  
  corrplot(res2$r, type = "upper",method = 'ellipse',
           tl.col = "black", tl.srt = 45, p.mat = res2$P, insig = 'blank', diag = FALSE, sig.level = 0.01)
  
  
  library(network)
  library(sna)
  library(ggplot2)
  
  
  pvalmat = as.matrix(res2$P <0.01)
  pvalmat[pvalmat==TRUE] <-1
  diag(pvalmat)<- 0
  
  netplot = res2$r * pvalmat

 
  library(igraph)
  net <- graph_from_adjacency_matrix(netplot, weighted = T, diag = F,mode = 'upper' )
  E(net)$width = abs(E(net)$weight)*5
  E(net)$arrow.size = 0
  E(net)$color <- ifelse(E(net)$weight>0, 'purple', 'blue')
  V(net)$color <- DiscretePalette(ncol(netplot), palette = NULL, shuffle = FALSE)
  graph_attr(net, "layout") <- layout.fruchterman.reingold
  
  svglite('plots/networks/IT_Rec.svg',width = 4, height = 4)
  plot(net, edge.curved=0.2, vertex.size=1)
  dev.off()
    
  
  
  alldata_afterQC %>%
    filter(phenotype %in% c('Other')) %>%
    nrow()
  