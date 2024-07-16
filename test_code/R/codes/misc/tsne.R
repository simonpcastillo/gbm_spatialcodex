require(ggplot2)


alldata_afterQC_path %>%
  filter(sample== 'IWP05', region =='reg001') %>%
  ggplot(aes(x,y,colour=path_region)) +
  geom_point()
  

require(ggplot2)

alldata_afterQC_path %>%
  filter(sample == 'IWP05', region== 'reg001')%>%
  ggplot(aes(x,y, colour=path_region)) +
  geom_point()


alldata_afterQC%>%
  group_by(acquisition_id) %>%
  summarise(sample = unique(sample),
            region = unique(region),
            totalcells= n(),
            nSALL1pos = length(SALL1pos[SALL1pos=='High']),
            nICAM1pos = length(ICAM1pos[ICAM1pos=='High']),
            nCD163pos = length(CD163pos[CD163pos=='High']),
            nCD16pos = length(CD16pos[CD16pos=='High']),
            nCD66bpos = length(CD66bpos[CD66bpos=='High']),
            nPD1pos = length(PD1pos[PD1pos=='High']),
            nAHRpos = length(AHRpos[AHRpos=='High']),
            nCD8pos = length(CD8pos[CD8pos=='High']),
            nTMEM119pos = length(TMEM119pos[TMEM119pos=='High']),
            nAXLpos = length(AXLpos[AXLpos=='High']),
            nCD107apos = length(CD107apos[CD107apos=='High']),
            nVISTApos = length(VISTApos[VISTApos=='High']),
            nPerforinpos = length(Perforinpos[Perforinpos=='High']),
            nLysozymepos = length(Lysozymepos[Lysozymepos=='High']),
            nPDL1pos = length(PDL1pos[PDL1pos=='High']),
            nGal9pos = length(Gal9pos[Gal9pos=='High']),
            nCD49dpos = length(CD49dpos[CD49dpos=='High']),
            nTIM3pos = length(TIM3pos[TIM3pos=='High']),
            nGranulysinpos = length(Granulysinpos[Granulysinpos=='High']),
            nCD56pos = length(CD56pos[CD56pos=='High']),
            nCD11cpos = length(CD11cpos[CD11cpos=='High']),
            nHLA.DRpos = length(HLA.DRpos[HLA.DRpos=='High']),
            nCD14pos = length(CD14pos[CD14pos=='High']),
            nFoxP3pos = length(FoxP3pos[FoxP3pos=='High']),
            nLYVE1pos = length(LYVE1pos[LYVE1pos=='High']),
            nGranzymeBpos = length(GranzymeBpos[GranzymeBpos=='High']),
            nCD45pos = length(CD45pos[CD45pos=='High']),
            nCD4pos = length(CD4pos[CD4pos=='High']),
            nOlig2pos = length(Olig2pos[Olig2pos=='High']),
            nNestinpos = length(Nestinpos[Nestinpos=='High']),
            nFoxP2pos = length(FoxP2pos[FoxP2pos=='High']),
            nCD370pos = length(CD370pos[CD370pos=='High']),
            nB2Mpos = length(B2Mpos[B2Mpos=='High']),
            nKi67pos = length(Ki67pos[Ki67pos=='High']),
            nCCR2pos = length(CCR2pos[CCR2pos=='High']),
            nCX3CR1pos = length(CX3CR1pos[CX3CR1pos=='High']),
            nCD68pos = length(CD68pos[CD68pos=='High']),
            nCD1cpos = length(CD1cpos[CD1cpos=='High'])) ->summarycounts


dfratio <- summarycounts %>%
  mutate(across(nSALL1pos:nCD1cpos, ~ ./totalcells))

write.csv(dfratio, 'sumamry_relativeabundance_allmarkers.csv')

library(tsne); require(plotly); require(dplyr)

dfratio_only <- alldata_afterQC_path %>%
  filter(in_region != 255) %>%
  filter(path_region == 'Cellular tumor') %>%
  group_by(acquisition_id) %>%
  sample_frac(size= 0.001)

dfratio <- dfratio_only[,15:56]

set.seed(0)
tsnedf <- tsne(dfratio, initial_dims = 3)
tsnedf<- data.frame(tsnedf)
tsnedf$sample = dfratio_only$sample
tsnedf$region = dfratio_only$region
tsnedf$acquisition_id = dfratio_only$acquisition_id

fig <-  plot_ly(data = tsnedf ,x =  ~X1, y = ~X2, type = 'scatter', mode = 'markers', split = ~sample)

fig <- fig %>%
  layout(
    plot_bgcolor = "#e5ecf6"
  )

fig


fig <-  plot_ly(data = tsnedf ,x =  ~X1, y = ~X2, z = ~X3, color = ~sample, colors = c('#636EFA','#EF553B','#00CC96') ) %>% 
  add_markers(size = 8) %>%
  layout( 
    xaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'), 
    yaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    scene =list(bgcolor = "#e5ecf6"))
fig


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("BEclear")
library(BEclear)


t_dfratio = t(dfratio)

ex.samples = dfratio_only %>% select(cell.id, acquisition_id)

colnames(ex.samples) <- c('sample_id', 'batch_id')

colnames(t_dfratio) = ex.samples$sample_id

batchEffect <- calcBatchEffects(
  data = t_dfratio, samples = ex.samples,
  adjusted = TRUE, method = "fdr"
)

mdifs <- batchEffect$med
pvals <- batchEffect$pval
summary <- calcSummary(medians = mdifs, pvalues = pvals)


library(kBET)
batch.estimate <- kBET(dfratio[, 'CD8'], dfratio_only$acquisition_id)
