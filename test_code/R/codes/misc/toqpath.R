size_px = c(8064,8064)

toqpath <- name0[, c('x', 'y')]

toqpath %>%
  mutate(y = size_px[2]-y) -> toqpath

write.table(toqpath, file='/Volumes/proj5/GBM_MDACC/validation/cellpos/s338_c004_v001_r001_reg001.tsv', quote=FALSE, sep='\t', row.names =F)
