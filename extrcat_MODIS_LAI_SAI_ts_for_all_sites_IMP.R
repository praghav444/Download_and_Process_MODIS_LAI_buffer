library(zoo)
library(lubridate)
library(dplyr)
library(tidyr)
#clear R environment
rm(list=ls(all=TRUE))

# Read MODIS data that were downloaded using MODIStool 
df_lai <- read.csv('Data/MOD15A2H_Lai_500m_ICOS_FLUXNET.csv', header=TRUE)
df_qc <- read.csv('Data/MOD15A2H_FparLai_QC_ICOS_FLUXNET.csv', header=TRUE)
df_sd <- read.csv('Data/MOD15A2H_LaiStdDev_500m_ICOS_FLUXNET.csv', header=TRUE)
site_codes <- unique(df_lai$site)

"
MODISTools R package only lets you to extract 1km around site coordinates.
This is greater than the average site footprint.
Code returns data for 5 x 5 pixels, manually verified than the pixel numbers
provided (1-25) go by rows, i.e.

 [ 1 ] [ 2 ]  [ 3 ] [ 4 ]  [ 5 ]
 [ 6 ] [ 7*]  [ 8*] [ 9*]  [10 ]
 [11 ] [12*]  [ce*] [14*]  [15 ]
 [16 ] [17*]  [18*] [19*]  [20 ]
 [21 ] [22 ]  [23 ] [24 ]  [25 ]

 Extracting pixels in the centre and immediately around it (*)
 These correspond to a radius of 500m around site coordinates
"

pixel_no <- c(7, 8, 9, 12, 13, 14, 17, 18, 19)  # Pixels just around the centre pixel

#Use only good quality data
#Random bit integer format, ask Martin if need to work these out again...

qc_flags <- c(0, 2, 24 ,26, 32, 34, 56, 58)

#Loop through sites
for (s in 1:length(site_codes)) {
  
  #Exception for US-ORv, wetland site with no MODIS LAI available
  if (site_codes[s] == "US-ORv") next
  print(paste0('----Processing for site ',site_codes[s],'----'))
  #Read data
  lai <- subset(df_lai, site==site_codes[s])
  qc  <- subset(df_qc, site==site_codes[s])
  sd  <- subset(df_sd, site==site_codes[s])
  
  
  #Number of time steps
  
  #If data files have been acquired on different days from MODIS
  #server, they might be different lengths so take minimum. adjust further down
  
  no_tsteps <- min(nrow(lai), nrow(sd), nrow(qc)) / max(lai$pixel)
  
  
  #Extract 3 x3 pixels
  lai_pixel <- matrix(nrow=no_tsteps, ncol=length(pixel_no))
  sd_pixel  <- matrix(nrow=no_tsteps, ncol=length(pixel_no))
  qc_pixel  <- matrix(nrow=no_tsteps, ncol=length(pixel_no))
  
  #Save time stamps
  lai_time <- as.Date(lai$calendar_date[which(lai$pixel == pixel_no[1])])
  
  #Loop through pixels
  for (p in 1:length(pixel_no)) {
    
    #Get time series for pixel and scale using scale factor (and adjust for different lengths with min_dim)
    lai_pixel[,p] <- lai$value[which(lai$pixel == pixel_no[p])][1:no_tsteps] * lai$scale[1]
    sd_pixel[,p]  <- sd$value[which(sd$pixel == pixel_no[p])][1:no_tsteps] * sd$scale[1]
    qc_pixel[,p]  <- qc$value[which(sd$pixel == pixel_no[p])][1:no_tsteps] 
    
  }
  
  ##################################
  ### Mask out poor quality data ###
  ##################################
  
  #Mask out where QC flag 
  lai_pixel <- replace(lai_pixel, !(qc_pixel %in% qc_flags), NA)
  sd_pixel  <- replace(sd_pixel, !(qc_pixel %in% qc_flags), NA)
  
  #Also mask out where sd really low (likely cloud effects)
  sd_pixel  <- replace(sd_pixel, sd_pixel < 0.1, NA)
  lai_pixel <- replace(lai_pixel, is.na(sd_pixel), NA)
  
  #Set fill values to missing
  lai_pixel <- replace(lai_pixel, lai_pixel > 10, NA)
  
  
  ### Average each time step ###
  
  #Initialise lai time series
  lai_ts <- vector(length=no_tsteps)
  #lai_ts_mean <- vector(length=no_tsteps)
  
  #Grid maximum LAI
  grid_max_lai_ts <- vector(length=no_tsteps)
  #Loop through time steps
  for (t in 1:no_tsteps) {
    
    #If no values available
    if (all(is.na(lai_pixel[t,]))) {
      
      lai_ts[t] <- NA
      grid_max_lai_ts[t] <- NA
      
      #If values available
    } else {
      
      #Weight grid cell estimates by their standard deviation
      #Following Martin's method (https://github.com/mdekauwe/get_MODIS_LAI_australia/blob/master/build_modis_climatology.py),
      #but normalising by sum of standard deviations
      
      sd_vals <- sd_pixel[t,]
      
      weights = (1/sd_vals**2) / sum(1/sd_vals**2, na.rm=TRUE)
      
      #Check that weights sum up to 1 (because of a precision issue presumably,
      #rounding to 5 decimals, otherwise might not equal 1 even when correct)
      if (round(sum(weights, na.rm=TRUE), 5) != 1) stop("Weighting not correct")
      
      #Calculate weighted average
      lai_ts[t] <- weighted.mean(lai_pixel[t,], w=weights, na.rm=TRUE)
      grid_max_lai_ts[t] <- max(lai_pixel[t,], na.rm=T)
      #lai_ts_mean[t] <- mean(lai_pixel[t,], na.rm=T)
    }
  }
  
  
  ######################################
  ### Gapfill and smooth with spline ###
  ######################################
  
  lai_ts_df <- data.frame(DateTime=lai_time, LAI=lai_ts)
  lai_ts_df$DateTime <- as.POSIXct(lai_ts_df$DateTime, format="%Y-%m-%d")
  
  hr <- zoo(lai_ts_df$LAI, lai_ts_df$DateTime)
  rng <- range(time(hr))
  tt <- seq(rng[1], rng[2], by = "hour")
  
  #Gapfill with spline (and cap negative values)
  temp_file <- list.files(path = 'C://Users/ppushpendra/OneDrive - The University of Alabama/Summer_Visit_to_OSU/CLM_Input_Data_Preparation/FLUXNET_DATA/AmeriFlux/RAW_FILES/',
                          pattern = paste0('AMF_',site_codes[s],'_FLUXNET_FULLSET_*'), recursive = TRUE, full.names = TRUE)
  temp_df <- read.csv(temp_file)
  temp_df$DateTime <- strptime(temp_df$TIMESTAMP_START, format="%Y%m%d%H%M")
  cols_to_delete <- c("Date", "Hour", "LAI_smooth", "LAI")
  temp_df <- temp_df %>%
    select(-one_of(cols_to_delete))
  
  #Define spline function
  func1 = splinefun(x=lai_ts_df$DateTime, y=lai_ts_df$LAI, method="fmm",  ties = mean)
  #Gapfill with spline (and cap negative values)
  lai_spline <- func1(tt)
  lai_spline[lai_spline < 0] <- 0
  
  #Smooth with spline (and cap negative values)
  smooth_lai_ts = smooth.spline(tt, lai_spline)$y
  smooth_lai_ts[smooth_lai_ts < 0] <- 0
  
  
  smooth_lai_ts <- data.frame(DateTime=tt, LAI_smooth=smooth_lai_ts, LAI=lai_spline)
  smooth_lai_ts$DateTime <- strptime(smooth_lai_ts$DateTime, format="%Y-%m-%d %H:%M")
  smooth_lai_ts$Date <- date(smooth_lai_ts$DateTime)
  smooth_lai_ts$Hour <- hour(smooth_lai_ts$DateTime)
  
  temp_df$Date <- date(temp_df$DateTime)
  temp_df$Hour <- hour(temp_df$DateTime)
  
  temp_df <- temp_df[,!(names(temp_df) %in% c('DateTime'))]
  smooth_lai_ts <- smooth_lai_ts[,!(names(smooth_lai_ts) %in% c('DateTime'))]
  temp_df <- left_join(temp_df, smooth_lai_ts, by=c('Date','Hour'))
  temp_df <- temp_df[!(duplicated(temp_df$TIMESTAMP_START)),]
  
  temp_df <- temp_df %>%
    fill(LAI, .direction = "downup") %>%
    fill(LAI_smooth, .direction = "downup")
  
  # Delete columns starting with "RECO" and "NEE" and ending with QC
  temp_df <- temp_df %>%
    select(-matches("^RECO"), -matches("^NEE"), -matches("QC$"))
  # #Test
  plot(lai_spline, type='l')
  lines(smooth_lai_ts$LAI, col='blue')
  
  #----Add Stem Area Index (SAI)-------------Zeng et. al. 2002, Journal of Climate, p1835-------
  # Get SAI from LAI (https://github.com/ESCOMP/CTSM/blob/master/src/biogeochem/CNVegStructUpdateMod.F90)
  dtsmonth <- 2592000 # Seconds on 1 month (30 days)
  dt <- 60 * 60       # seconds in 1 hour
  tsai_min <- 1       # or 0.5 
  tsai_alpha <- 1.0 - 0.5 * dt / dtsmonth
  LAI <- temp_df$LAI_smooth
  SAI <- numeric(length(LAI))
  SAI[1] <- max(LAI[0], tsai_min)
  for (t in 2:length(LAI)) {
    LAI_old <- LAI[t - 1]
    LAI_current <- LAI[t]
    SAI_old <- SAI[t - 1]
    SAI[t] <- max(tsai_alpha * SAI_old + max(LAI_old - LAI_current, 0), tsai_min)
  }
  plot(SAI[1:50000])
  temp_df$SAI_smooth <- SAI
  #Save
  write.csv(temp_df, temp_file, row.names = FALSE)
}
