library(MODISTools)

rm(list=ls())

site_info <- read.csv('AmeriFlux_FLUXNET_GEDI_forest_sites_GEE.csv')

band <- "Fpar_500m"   # Lai_500m, Fpar_500m, FparLai_QC, LaiStdDev_500m
products <- mt_products()
bands <- mt_bands(product = "MOD15A2H")
subsets <- mt_batch_subset(df = site_info,
                           product = "MOD15A2H",
                           band = band,   
                           internal = TRUE,
                           km_lr = 1,
                           km_ab = 1,
                           out_dir = 'Data/',
                           start = "2000-01-01",
                           end = "2023-12-31")
write.csv(subsets, paste0('Data/MOD15A2H_',band,'_ICOS_FLUXNET.csv'), row.names = FALSE)
