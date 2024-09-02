library(MODISTools)

rm(list=ls())

setwd("C://Users/ppushpendra/OneDrive - The University of Alabama/Summer_Visit_to_OSU/CLM_Input_Data_Preparation/Prepare_LAI_SAI_Data/")
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
write.csv(subsets, paste0('Data/MOD15A2H_',band,'_AmeriFlux_FLUXNET_GEDI_Forest.csv'), row.names = FALSE)
