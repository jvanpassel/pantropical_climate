library(sp)
library(raster)
library(sf)
library(geosphere)
library(dplyr)
library(tsbox)
library(tidyr)
library(BIOMASS)
library(scales)
library(zoo)
library(RColorBrewer)
library(terra)
library(ggplot2)
library(nlme)
library(MuMIn)
library(lubridate)
library(ggpubr)
library(SPEI)
library(MASS)
library(reticulate)
library(maps)
library(tidyverse)
library(readr)
library(stringr)
library(xts)
library(ggtern)
library(cowplot)
library(scales)
library(patchwork)
library(rnaturalearth)
library(scales)
library(reshape2)

# 1) Calculate monthly precipitation and temperature anomalies ####
prec_3imerg<-stack('df.3IMERG.Tropics.climate.rspd.0024.prec.mask50.tif')
prec_chirps<-stack('df.chirps.Tropics.climate.rspd.0024.prec.mask50.tif')
prec_cru<-stack('df.CRU.Tropics.climate.rspd.0024.prec.mask50.tif')
prec_era5<-stack('df.ERA5.Tropics.climate.rspd.0024.prec.mask50.tif')
prec_gldas<-stack('df.GLDAS.Tropics.climate.rspd.0024.prec.mask50.tif')
prec_gpcc<-stack('df.GPCC.Tropics.climate.rspd.0024.prec.mask50.tif')
prec_merra2<-stack('df.MERRA2.Tropics.climate.rspd.0024.prec.mask50.tif')
prec_mswep<-stack('df.MSWEP.Tropics.climate.rspd.0024.prec.mask50.tif')
prec_ncep<-stack('df.NCEP.Tropics.climate.rspd.0024.prec.mask50.tif')
prec_weighted_av<-rast('df.Tropics.climate.rspd.0024.prec.mean.weightedR2.tif')

tav_cru<-stack('df.CRU.Tropics.climate.rspd.0024.tav.mask50.tif')
tav_era5<-stack('df.ERA5.Tropics.climate.rspd.0024.tav.mask50.tif')
tav_gldas<-stack('df.GLDAS.Tropics.climate.rspd.0024.tav.mask50.tif')
tav_merra2<-stack('df.MERRA2.Tropics.climate.rspd.0024.tav.mask50.tif')
tav_berk<-stack('df.Berk.Tropics.climate.rspd.0024.tav.mask50.tif')
tav_cams<-stack('df.CAMS.Tropics.climate.rspd.0024.tav.mask50.tif')
tav_ncep<-stack('df.NCEP.Tropics.climate.rspd.0024.tav.mask50.tif')
tav_weighted_av<-rast('df.Tropics.climate.rspd.0024.tav.mean.weightedR2.tif')

#run for each product
prec_mask<-prec_weighted_av
#don't include 2015-16-23-24 for calculating monthly mean and SD
prec_mask_removeEN<-prec_mask[[c(1:180,205:276)]]
months_yr_total<-rep(seq(1:12),25)
months_yr_removeEN<-rep(seq(1:12),21)
#monthly mean
prec_monthlymean<-stackApply(prec_mask_removeEN,fun = mean,months_yr_removeEN)
prec_monthlymean_full<-prec_monthlymean[[months_yr_total]]
#monthly SD
prec_monthlysd<-stackApply(prec_mask_removeEN,fun = sd,months_yr_removeEN)
prec_monthlysd_full<-prec_monthlysd[[months_yr_total]]
#standardised anomalies
prec_st_anomalies<-(prec_mask-prec_monthlymean_full)/prec_monthlysd_full
prec_st_anomalies[is.na(prec_st_anomalies)]<-0
writeRaster(prec_st_anomalies,'df.Tropics.climate.rspd.0024.prec.mean.weightedR2.st.anomalies.removeEN.tif',overwrite=T)

# 2) Use 1SD thresholds to define spatial and temporal coverage of EN ####
#Dry events ####
prec_an<-prec_st_anomalies
prec_an_df<-as.data.frame(prec_an,xy=T)
prec_an_df_nona<-prec_an_df[!is.na(prec_an_df[,3]),]
nonarows<-which(!is.na(prec_an_df[,3]))
neganomaly<-(-1) #change to how big anomaly should be ~ SD
find_dry_periods<-function(x){
  if(!all(is.na(x))){
    negindex<-which(x<0) # Find positions of all negative anomalies
    signnegindex<-which(x<neganomaly) # Find positions of negative anomalies
    if(length(signnegindex)!=0){
      droughtstart<-split(signnegindex, cumsum(c(1, diff(signnegindex) != 1))) # Make list of all separate significant droughts
      droughtlist<-list()
      for(i in 1:length(droughtstart)){
        if(length(droughtstart[[i]])>1){ # Drought should start with at least two months of sign neg anomaly
          start<-which(negindex==droughtstart[[i]][[1]]) # Find position of drought start
          negindex_2<-negindex[-c(1:(start-1))] # Delete everything before this drought
          dry<-cumsum(c(1, diff(negindex_2) != 1)) # See how many months were dry after drought start
          if(sum(dry==1)>1){ # Should be at least 1 month afterwards (sum at least 2)
            dur<-sum(dry==1)
            droughtlist[[i]]<-c(negindex[start]:(negindex[start]+dur-1)) # Get positions of total droughts
          }
        }
      }
      if(length(droughtlist)!=0){
        if(length(which(sapply(droughtlist, is.null)))>0){
          droughtlist<-droughtlist[-which(sapply(droughtlist, is.null))] # Remove empty list elements (due to droughtstart = 1 or sum = 1)
        }
        if(length(droughtlist)!=0){
          droughts<-unlist(droughtlist)
          y<-x
          y[-droughts]<-NA # Only droughts are not NA
        } 
      } else {
        y<-rep(NA,length(x))
      }
    } else { 
      y<-rep(NA,length(x))
    }
  } else {
    y<-rep(NA,length(x))
  }
  return(y)
}
dry_periods_df<-as.data.frame(t(apply(prec_an_df[3:302],MARGIN = 1,FUN = find_dry_periods)))
#convert csv to tif of 2023-2024
prec_3imerg_an<-stack('df.3IMERG.Tropics.climate.rspd.0024.prec.st.anomalies.removeEN.tif')
prec_an_df<-as.data.frame(prec_3imerg_an,xy=T)
prec_coord<-prec_an_df[,1:3]
#only look at 2023-2024
prec_dry1_2324<-dry_periods_df[,277:300]
#find area of dry period
prec_dry1_2324_nonapixels<-which(rowSums(!is.na(prec_dry1_2324))>0)
prec_dry_coord<-prec_coord
prec_dry_coord[!is.na(prec_dry_coord[,3]),3]<-0
prec_dry_coord[prec_dry1_2324_nonapixels,3]<-1
length(prec_dry1_2324_nonapixels)/sum(!is.na(prec_dry_coord[,3]))*100
#get percentages per continent
prec_dry_coord_am<-prec_dry_coord[prec_dry_coord$x<(-25),]
prec_dry_coord_af<-prec_dry_coord[prec_dry_coord$x>(-25) & prec_dry_coord$x<50,]
prec_dry_coord_as<-prec_dry_coord[prec_dry_coord$x>50,]
sum(prec_dry_coord_am[,3]>0,na.rm=T)/sum(!is.na(prec_dry_coord_am[,3]))*100
sum(prec_dry_coord_af[,3]>0,na.rm=T)/sum(!is.na(prec_dry_coord_af[,3]))*100
sum(prec_dry_coord_as[,3]>0,na.rm=T)/sum(!is.na(prec_dry_coord_as[,3]))*100
colnames(prec_dry_coord)<-c('x','y','dry')
#convert to 0-1 raster to show where dry periods occurred during EN
prec_dry <- rasterFromXYZ(prec_dry_coord[, c('x','y','dry')])
crs(prec_dry)<-'+proj=longlat +datum=WGS84 +no_defs'
writeRaster(prec_dry,'df.Tropics.climate.rspd.0024.prec.mean.weightedR2.dry.removeEN.tif',overwrite=T)

#Hot events ####
tav_an<-tav_st_anomalies
tav_an_df<-as.data.frame(tav_an,xy=T)
posanomaly<-1 #change to how big anomaly should be ~ SD
find_hot_periods<-function(x){
  if(!all(is.na(x))){
    posindex<-which(x>0) # Find positions of all positive anomalies
    signposindex<-which(x>posanomaly) # Find positions of positive anomalies
    if(length(signposindex)!=0){
      hotstart<-split(signposindex, cumsum(c(1, diff(signposindex) != 1))) # Make list of all separate significant droughts
      hotlist<-list()
      for(i in 1:length(hotstart)){
        if(length(hotstart[[i]])>1){ # Drought should start with at least two months of sign neg anomaly
          start<-which(posindex==hotstart[[i]][[1]]) # Find position of drought start
          posindex_2<-posindex[-c(1:(start-1))] # Delete everything before this drought
          hot<-cumsum(c(1, diff(posindex_2) != 1)) # See how many months were dry after drought start
          if(sum(hot==1)>1){ # Should be at least 1 month afterwards (sum at least 2)
            dur<-sum(hot==1)
            hotlist[[i]]<-c(posindex[start]:(posindex[start]+dur-1)) # Get positions of total droughts
          }
        }
      }
      if(length(hotlist)!=0){
        if(length(which(sapply(hotlist, is.null)))>0){
          hotlist<-hotlist[-which(sapply(hotlist, is.null))] # Remove empty list elements (due to droughtstart = 1 or sum = 1)
        }
        if(length(hotlist)!=0){
          hots<-unlist(hotlist)
          y<-x
          y[-hots]<-NA # Only droughts are not NA
        } 
      } else {
        y<-rep(NA,length(x))
      }
    } else { 
      y<-rep(NA,length(x))
    }
  } else {
    y<-rep(NA,length(x))
  }
  return(y)
}
hot_periods_df<-as.data.frame(t(apply(tav_an_df[3:302],MARGIN = 1,FUN = find_hot_periods)))
#convert csv to tif of 2023-2024
prec_3imerg_an<-stack('./RDSfiles/df.3IMERG.Tropics.climate.rspd.0024.prec.st.anomalies.tif')
prec_an_df<-as.data.frame(prec_3imerg_an,xy=T)
prec_coord<-prec_an_df[,1:3]
#only look at 2023-2024
tav_hot1_2324<-hot_periods_df[,277:300]
#find area of dry period
tav_hot1_2324_nonapixels<-which(rowSums(!is.na(tav_hot1_2324))>0)
tav_hot_coord<-prec_coord
tav_hot_coord[!is.na(tav_hot_coord[,3]),3]<-0
tav_hot_coord[tav_hot1_2324_nonapixels,3]<-1
length(tav_hot1_2324_nonapixels)/sum(!is.na(tav_hot_coord[,3]))*100
#get percentages per continent
tav_hot_coord_am<-tav_hot_coord[tav_hot_coord$x<(-25),]
tav_hot_coord_af<-tav_hot_coord[tav_hot_coord$x>(-25) & tav_hot_coord$x<50,]
tav_hot_coord_as<-tav_hot_coord[tav_hot_coord$x>50,]
sum(tav_hot_coord_am[,3]>0,na.rm=T)/sum(!is.na(tav_hot_coord_am[,3]))*100
sum(tav_hot_coord_af[,3]>0,na.rm=T)/sum(!is.na(tav_hot_coord_af[,3]))*100
sum(tav_hot_coord_as[,3]>0,na.rm=T)/sum(!is.na(tav_hot_coord_as[,3]))*100
colnames(tav_hot_coord)<-c('x','y','hot')
#convert to 0-1 raster to show where hot periods occurred during EN
tav_hot <- rasterFromXYZ(tav_hot_coord[, c('x','y','hot')])
crs(tav_hot)<-'+proj=longlat +datum=WGS84 +no_defs'
writeRaster(tav_hot,'./RDSfiles/df.Tropics.climate.rspd.0024.tav.mean.weightedR2.hot.removeEN.tif',overwrite=T)

#Hot and dry events ####
#first combine them into new df per product with sum if both events co-occuring and NA if not
dry1_df<-dry_periods_df
hot1_df<-hot_periods_df
hotdry_periods_df<-dry1_df+hot1_df
write.csv(hotdry_periods_df,'./RDSfiles/df.Tropics.climate.rspd.0024.mean.weightedR2.hotdry.1SD.removeEN.csv')
#only look at 2023-2024
hotdry1_2324<-hotdry_periods_df[,277:300]
hot1_2324<-hot1_df[,277:300]
dry1_2324<-dry1_df[,277:300]
#find area of hot and dry period
prec_3imerg_an<-stack('./RDSfiles/df.3IMERG.Tropics.climate.rspd.0024.prec.st.anomalies.tif')
prec_an_df<-as.data.frame(prec_3imerg_an,xy=T)
prec_coord<-prec_an_df[,1:5]
hotdry1_2324_nonapixels<-which(rowSums(!is.na(hotdry1_2324))>0)
hot1_2324_nonapixels<-which(rowSums(!is.na(hot1_2324))>0)
dry1_2324_nonapixels<-which(rowSums(!is.na(dry1_2324))>0)
hotdry_coord<-prec_coord
hotdry_coord[!is.na(hotdry_coord[,3]),3]<-0
hotdry_coord[!is.na(hotdry_coord[,3]),4]<-0
hotdry_coord[!is.na(hotdry_coord[,3]),5]<-0
hotdry_coord[hot1_2324_nonapixels,3]<-1
hotdry_coord[dry1_2324_nonapixels,4]<-1
hotdry_coord[hotdry1_2324_nonapixels,5]<-1
#get percentages per continent
sum(hotdry_coord[,5]>0,na.rm=T)/sum(!is.na(hotdry_coord[,5]))*100
hotdry_coord_am<-hotdry_coord[hotdry_coord$x<(-25),]
hotdry_coord_af<-hotdry_coord[hotdry_coord$x>(-25) & hotdry_coord$x<50,]
hotdry_coord_as<-hotdry_coord[hotdry_coord$x>50,]
sum(hotdry_coord_am[,5]>0,na.rm=T)/sum(!is.na(hotdry_coord_am[,5]))*100
sum(hotdry_coord_af[,5]>0,na.rm=T)/sum(!is.na(hotdry_coord_af[,5]))*100
sum(hotdry_coord_as[,5]>0,na.rm=T)/sum(!is.na(hotdry_coord_as[,5]))*100
#convert to 0-1 raster to show where hot periods occurred during EN
colnames(hotdry_coord)<-c('x','y','hot','dry','dry_hot')
hotdry <- rasterFromXYZ(hotdry_coord[, c('x','y','dry_hot')])
crs(hotdry)<-'+proj=longlat +datum=WGS84 +no_defs'
writeRaster(hotdry,'df.Tropics.climate.rspd.0024.mean.weightedR2.hotdry.removeEN.tif',overwrite=T)


# 3) Look at temporal variation in coverage between products #####
prec_3imerg_an<-stack('df.3IMERG.Tropics.climate.rspd.0024.prec.st.anomalies.removeEN.tif')
prec_an_df<-as.data.frame(prec_3imerg_an,xy=T)
nonarows<-which(!is.na(prec_an_df[,3]))
prec_coord<-prec_an_df[nonarows,1:2]

prec_coord_am<-prec_coord
prec_coord_am<-prec_coord_am[prec_coord_am$x<(-25)&prec_coord_am$x>(-120),]
am_rows<-which(prec_coord$x<(-25)&prec_coord$x>(-120))
prec_coord_af<-prec_coord
prec_coord_af<-prec_coord_af[prec_coord_af$x<(50)&prec_coord_af$x>(-25),]
af_rows<-which(prec_coord$x<(50)&prec_coord$x>(-25))
prec_coord_as<-prec_coord
prec_coord_as<-prec_coord_as[prec_coord_as$x<(180)&prec_coord_as$x>(50),]
as_rows<-which(prec_coord$x<(180)&prec_coord$x>(50))

#take continental % of locations with dry/hot/hot and dry events
prec_prod1<-dry_periods_df[nonarows,] #or hot_periods_df or hotdry_periods_df
prec_prod1_am<-prec_prod1[am_rows,]
prec_prod1_af<-prec_prod1[af_rows,]
prec_prod1_as<-prec_prod1[as_rows,]
prec_prod1_perc<-as.numeric(apply(prec_prod1,MARGIN = 2,FUN = function(x) sum(!is.na(x))/length(x)*100))
plot(prec_prod1_perc,type='l',ylab='% of area experiencing dry/hot event',xlab='Months since Jan 2000')
prec_prod1_am_perc<-as.numeric(apply(prec_prod1_am,MARGIN = 2,FUN = function(x) sum(!is.na(x))/length(x)*100))
plot(prec_prod1_am_perc,type='l',ylab='% of area experiencing dry/hot event',xlab='Months since Jan 2000')
prec_prod1_af_perc<-as.numeric(apply(prec_prod1_af,MARGIN = 2,FUN = function(x) sum(!is.na(x))/length(x)*100))
plot(prec_prod1_af_perc,type='l',ylab='% of area experiencing dry/hot event',xlab='Months since Jan 2000')
prec_prod1_as_perc<-as.numeric(apply(prec_prod1_as,MARGIN = 2,FUN = function(x) sum(!is.na(x))/length(x)*100))
plot(prec_prod1_as_perc,type='l',ylab='% of area experiencing dry/hot event',xlab='Months since Jan 2000')
prec_prod1_cont<-rbind(prec_prod1_am_perc,prec_prod1_af_perc,prec_prod1_as_perc,prec_prod1_perc)
write.csv(prec_prod1_cont,'./RDSfiles/df.Tropics.climate.rspd.0024.prec.mean.weightedR2.dry.1SD.removeEN.perc.ts.csv')



# 4) Plot variation in temporal coverage over 2000-2024 ####
#Precipitation anomalies ####
imerg_1<-read.csv('df.3IMERG.Tropics.climate.rspd.0024.dry.1SD.removeEN.perc.ts.csv')[,-1]
chirps_1<-read.csv('df.chirps.Tropics.climate.rspd.0024.dry.1SD.removeEN.perc.ts.csv')[,-1]
cru_1<-read.csv('df.CRU.Tropics.climate.rspd.0024.dry.1SD.removeEN.perc.ts.csv')[,-1]
era5_1<-read.csv('df.ERA5.Tropics.climate.rspd.0024.dry.1SD.removeEN.perc.ts.csv')[,-1]
gldas_1<-read.csv('df.GLDAS.Tropics.climate.rspd.0024.dry.1SD.removeEN.perc.ts.csv')[,-1]
gpcc_1<-read.csv('df.GPCC.Tropics.climate.rspd.0024.dry.1SD.removeEN.perc.ts.csv')[,-1]
merra2_1<-read.csv('df.MERRA2.Tropics.climate.rspd.0024.dry.1SD.removeEN.perc.ts.csv')[,-1]
mswep_1<-read.csv('df.MSWEP.Tropics.climate.rspd.0024.dry.1SD.removeEN.perc.ts.csv')[,-1]
ncep_1<-read.csv('df.NCEP.Tropics.climate.rspd.0024.dry.1SD.removeEN.perc.ts.csv')[,-1]
av_1<-read.csv('df.Tropics.climate.rspd.0024.prec.mean.weightedR2.dry.1SD.removeEN.perc.ts.csv')[,-1]

EN_timing<-c(32:39,55:57,59:60,62:63,80:85,115:116,118:124,185:197,224:226,230:231,239,283:291)
dates<-seq(as.Date('2000-01-01'),as.Date('2024-12-01'),by='month')
dry_1_am<-rbind(imerg_1[1,],chirps_1[1,],cru_1[1,],era5_1[1,],gldas_1[1,],gpcc_1[1,],
                mswep_1[1,],ncep_1[1,])
dry_1_af<-rbind(imerg_1[2,],chirps_1[2,],cru_1[2,],era5_1[2,],gldas_1[2,],gpcc_1[2,],
                mswep_1[2,],ncep_1[2,])
dry_1_as<-rbind(imerg_1[3,],chirps_1[3,],cru_1[3,],era5_1[3,],gldas_1[3,],gpcc_1[3,],
                mswep_1[3,],ncep_1[2,])
dry_1_pan<-rbind(imerg_1[4,],chirps_1[4,],cru_1[4,],era5_1[4,],gldas_1[4,],gpcc_1[4,],
                 mswep_1[4,],ncep_1[4,])
dry_av_am<-as.numeric(av_1[1,])
dry_av_af<-as.numeric(av_1[2,])
dry_av_as<-as.numeric(av_1[3,])
dry_av_pan<-as.numeric(av_1[4,])

dry_1_am_mean<-colMeans(dry_1_am)
dry_1_am_min<-apply(dry_1_am,MARGIN = 2,FUN = function(x) min(x))
dry_1_am_max<-apply(dry_1_am,MARGIN = 2,FUN = function(x) max(x))
max(dry_av_am)
which(dry_av_am==max(dry_av_am))
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = dry_1_am_min, ymax = dry_1_am_max),fill='darkseagreen')+
  geom_line(aes(x = dates,y = dry_1_am_mean),col='darkgreen') +
  geom_line(aes(x = dates,y = dry_av_am),col='darkgreen',linetype='dashed') +
  xlab('Time')+
  ylab('% of continent experiencing dry event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dry_1_af_mean<-colMeans(dry_1_af)
dry_1_af_min<-apply(dry_1_af,MARGIN = 2,FUN = function(x) min(x))
dry_1_af_max<-apply(dry_1_af,MARGIN = 2,FUN = function(x) max(x))
max(dry_av_af)
which(dry_av_af==max(dry_av_af))
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = dry_1_af_min, ymax = dry_1_af_max),fill='darkseagreen')+
  geom_line(aes(x = dates,y = dry_1_af_mean),col='darkgreen') +
  geom_line(aes(x = dates,y = dry_av_af),col='darkgreen',linetype='dashed') +
  xlab('Time')+
  ylab('% of continent experiencing dry event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dry_1_as_mean<-colMeans(dry_1_as)
dry_1_as_min<-apply(dry_1_as,MARGIN = 2,FUN = function(x) min(x))
dry_1_as_max<-apply(dry_1_as,MARGIN = 2,FUN = function(x) max(x))
max(dry_av_as)
which(dry_av_as==max(dry_av_as))
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = dry_1_as_min, ymax = dry_1_as_max),fill='darkseagreen')+
  geom_line(aes(x = dates,y = dry_1_as_mean),col='darkgreen') +
  geom_line(aes(x = dates,y = dry_av_as),col='darkgreen',linetype='dashed') +
  xlab('Time')+ 
  ylab('% of continent experiencing dry event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Temperature anomalies ####
cru_1<-read.csv('df.CRU.Tropics.climate.rspd.0024.hot.1SD.removeEN.perc.ts.csv')[,-1]
era5_1<-read.csv('df.ERA5.Tropics.climate.rspd.0024.hot.1SD.removeEN.perc.ts.csv')[,-1]
gldas_1<-read.csv('df.GLDAS.Tropics.climate.rspd.0024.hot.1SD.removeEN.perc.ts.csv')[,-1]
merra2_1<-read.csv('df.MERRA2.Tropics.climate.rspd.0024.hot.1SD.removeEN.perc.ts.csv')[,-1]
berk_1<-read.csv('df.Berk.Tropics.climate.rspd.0024.hot.1SD.removeEN.perc.ts.csv')[,-1]
cams_1<-read.csv('df.CAMS.Tropics.climate.rspd.0024.hot.1SD.removeEN.perc.ts.csv')[,-1]
ncep_1<-read.csv('df.NCEP.Tropics.climate.rspd.0024.hot.1SD.removeEN.perc.ts.csv')[,-1]
av_1<-read.csv('df.Tropics.climate.rspd.0024.tav.mean.weightedR2.hot.1SD.removeEN.perc.ts.csv')[,-1]

EN_timing<-c(32:39,55:57,59:60,62:63,80:85,115:116,118:124,185:197,224:226,230:231,239,283:291)
dates<-seq(as.Date('2000-01-01'),as.Date('2024-12-01'),by='month')
hot_1_am<-rbind(cru_1[1,],era5_1[1,],gldas_1[1,],berk_1[1,],cams_1[1,],ncep_1[1,])
hot_1_af<-rbind(cru_1[2,],era5_1[2,],gldas_1[2,],berk_1[2,],cams_1[2,],ncep_1[2,])
hot_1_as<-rbind(cru_1[3,],era5_1[3,],gldas_1[3,],berk_1[3,],cams_1[3,],ncep_1[3,])
hot_1_pan<-rbind(cru_1[4,],era5_1[4,],gldas_1[4,],berk_1[4,],cams_1[4,],ncep_1[4,])
hot_av_am<-as.numeric(av_1[1,])
hot_av_af<-as.numeric(av_1[2,])
hot_av_as<-as.numeric(av_1[3,])
hot_av_pan<-as.numeric(av_1[4,])

hot_1_am_mean<-colMeans(hot_1_am)
hot_1_am_min<-apply(hot_1_am,MARGIN = 2,FUN = function(x) min(x))
hot_1_am_max<-apply(hot_1_am,MARGIN = 2,FUN = function(x) max(x))
max(hot_1_am_max-hot_1_am_min)
which(c(hot_1_am_max-hot_1_am_min)==max(hot_1_am_max-hot_1_am_min))
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = hot_1_am_min, ymax = hot_1_am_max),fill='mistyrose')+
  geom_line(aes(x = dates,y = hot_1_am_mean),col='darkred') +
  geom_line(aes(x = dates,y = hot_av_am),col='darkred',linetype='dashed') +
  xlab('Time')+
  ylab('% of continent experiencing hot event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hot_1_af_mean<-colMeans(hot_1_af)
hot_1_af_min<-apply(hot_1_af,MARGIN = 2,FUN = function(x) min(x))
hot_1_af_max<-apply(hot_1_af,MARGIN = 2,FUN = function(x) max(x))
max(hot_1_af_max-hot_1_af_min)
which(c(hot_1_af_max-hot_1_af_min)==max(hot_1_af_max-hot_1_af_min))
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = hot_1_af_min, ymax = hot_1_af_max),fill='mistyrose')+
  geom_line(aes(x = dates,y = hot_1_af_mean),col='darkred') +
  geom_line(aes(x = dates,y = hot_av_af),col='darkred',linetype='dashed') +
  xlab('Time')+
  ylab('% of continent experiencing hot event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hot_1_as_mean<-colMeans(hot_1_as)
hot_1_as_min<-apply(hot_1_as,MARGIN = 2,FUN = function(x) min(x))
hot_1_as_max<-apply(hot_1_as,MARGIN = 2,FUN = function(x) max(x))
max(hot_1_as_max-hot_1_as_min)
which(c(hot_1_as_max-hot_1_as_min)==max(hot_1_as_max-hot_1_as_min))
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = hot_1_as_min, ymax = hot_1_as_max),fill='mistyrose')+
  geom_line(aes(x = dates,y = hot_1_as_mean),col='darkred') +
  geom_line(aes(x = dates,y = hot_av_as),col='darkred',linetype='dashed') +
  xlab('Time')+ 
  ylab('% of continent experiencing dry event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Hot-dry anomalies ####
#load all
imerg_1<-read.csv('df.3IMERG.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.perc.ts.csv')[,-1]
chirps_1<-read.csv('df.chirps.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.perc.ts.csv')[,-1]
cru_1<-read.csv('df.CRU.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.perc.ts.csv')[,-1]
era5_1<-read.csv('df.ERA5.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.perc.ts.csv')[,-1]
gldas_1<-read.csv('df.GLDAS.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.perc.ts.csv')[,-1]
gpcc_1<-read.csv('df.GPCC.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.perc.ts.csv')[,-1]
mswep_1<-read.csv('df.MSWEP.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.perc.ts.csv')[,-1]
berk_1<-read.csv('df.Berk.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.perc.ts.csv')[,-1]
cams_1<-read.csv('df.CAMS.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.perc.ts.csv')[,-1]
ncep_1<-read.csv('df.NCEP.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.perc.ts.csv')[,-1]
av_1<-read.csv('df.Tropics.climate.rspd.0024.mean.weightedR2.hotdry.1SD.removeEN.perc.ts.csv')[,-1]

EN_timing<-c(32:39,55:57,59:60,62:63,80:85,115:116,118:124,185:197,224:226,230:231,239,283:291)
dates<-seq(as.Date('2000-01-01'),as.Date('2024-12-01'),by='month')
hotdry_1_am<-rbind(imerg_1[1,],chirps_1[1,],cru_1[1,],era5_1[1,],gldas_1[1,],gpcc_1[1,],mswep_1[1,],
                   berk_1[1,],cams_1[1,],ncep_1[1,])
hotdry_1_af<-rbind(imerg_1[2,],chirps_1[2,],cru_1[2,],era5_1[2,],gldas_1[2,],gpcc_1[2,],mswep_1[2,],
                   berk_1[2,],cams_1[2,],ncep_1[2,])
hotdry_1_as<-rbind(imerg_1[3,],chirps_1[3,],cru_1[3,],era5_1[3,],gldas_1[3,],gpcc_1[3,],mswep_1[3,],
                   berk_1[3,],cams_1[3,],ncep_1[3,])
hotdry_1_pan<-rbind(imerg_1[4,],chirps_1[4,],cru_1[4,],era5_1[4,],gldas_1[4,],gpcc_1[4,],mswep_1[4,],
                    berk_1[4,],cams_1[4,],ncep_1[4,])
hotdry_av_am<-as.numeric(av_1[1,])
hotdry_av_af<-as.numeric(av_1[2,])
hotdry_av_as<-as.numeric(av_1[3,])
hotdry_av_pan<-as.numeric(av_1[4,])

hotdry_1_am_mean<-colMeans(hotdry_1_am)
hotdry_1_am_min<-apply(hotdry_1_am,MARGIN = 2,FUN = function(x) min(x))
hotdry_1_am_max<-apply(hotdry_1_am,MARGIN = 2,FUN = function(x) max(x))
max(hotdry_1_am_max-hotdry_1_am_min)
which(c(hotdry_1_am_max-hotdry_1_am_min)==max(hotdry_1_am_max-hotdry_1_am_min))
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = hotdry_1_am_min, ymax = hotdry_1_am_max),fill='peachpuff')+
  geom_line(aes(x = dates,y = hotdry_1_am_mean),col='darkorange2') +
  geom_line(aes(x = dates,y = hotdry_av_am),col='darkorange2',linetype='dashed') +
  xlab('Time')+
  ylab('% of continent experiencing hot event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hotdry_1_af_mean<-colMeans(hotdry_1_af)
hotdry_1_af_min<-apply(hotdry_1_af,MARGIN = 2,FUN = function(x) min(x))
hotdry_1_af_max<-apply(hotdry_1_af,MARGIN = 2,FUN = function(x) max(x))
max(hotdry_1_af_max-hotdry_1_af_min)
which(c(hotdry_1_af_max-hotdry_1_af_min)==max(hotdry_1_af_max-hotdry_1_af_min))
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = hotdry_1_af_min, ymax = hotdry_1_af_max),fill='peachpuff')+
  geom_line(aes(x = dates,y = hotdry_1_af_mean),col='darkorange2') +
  geom_line(aes(x = dates,y = hotdry_av_af),col='darkorange2',linetype='dashed') +
  xlab('Time')+
  ylab('% of continent experiencing hot event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hotdry_1_as_mean<-colMeans(hotdry_1_as)
hotdry_1_as_min<-apply(hotdry_1_as,MARGIN = 2,FUN = function(x) min(x))
hotdry_1_as_max<-apply(hotdry_1_as,MARGIN = 2,FUN = function(x) max(x))
max(hotdry_1_as_max-hotdry_1_as_min)
which(c(hotdry_1_as_max-hotdry_1_as_min)==max(hotdry_1_as_max-hotdry_1_as_min))
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = hotdry_1_as_min, ymax = hotdry_1_as_max),fill='peachpuff')+
  geom_line(aes(x = dates,y = hotdry_1_as_mean),col='darkorange2') +
  geom_line(aes(x = dates,y = hotdry_av_as),col='darkorange2',linetype='dashed') +
  xlab('Time')+ 
  ylab('% of continent experiencing dry event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Combine dry/hot/hot-dry in one plot (FIGURE 1) ####
#American forests
year<-rep(2000:2024,each=12)
dry_1_am_mean<-colMeans(dry_1_am)
dry_1_am_min<-apply(dry_1_am,MARGIN = 2,FUN = function(x) min(x))
dry_1_am_max<-apply(dry_1_am,MARGIN = 2,FUN = function(x) max(x))
max(dry_av_am) #65%
year[which(dry_av_am==max(dry_av_am))] #298-2024
hot_1_am_mean<-colMeans(hot_1_am)
hot_1_am_min<-apply(hot_1_am,MARGIN = 2,FUN = function(x) min(x))
hot_1_am_max<-apply(hot_1_am,MARGIN = 2,FUN = function(x) max(x))
max(hot_av_am) #99%
year[which(hot_av_am==max(hot_av_am))] #285-2023
hotdry_1_am_mean<-colMeans(hotdry_1_am)
hotdry_1_am_min<-apply(hotdry_1_am,MARGIN = 2,FUN = function(x) min(x))
hotdry_1_am_max<-apply(hotdry_1_am,MARGIN = 2,FUN = function(x) max(x))
max(hotdry_av_am) #65%
year[which(hotdry_av_am==max(hotdry_av_am))] #298-2024
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = dry_1_am_min, ymax = dry_1_am_max),fill='darkseagreen',alpha=0.5)+
  geom_ribbon(aes(x = dates,ymin = hot_1_am_min, ymax = hot_1_am_max),fill='mistyrose',alpha=0.5)+
  geom_ribbon(aes(x = dates,ymin = hotdry_1_am_min, ymax = hotdry_1_am_max),fill='papayawhip',alpha=0.5)+
  geom_line(aes(x = dates,y = dry_av_am),col='darkgreen') +
  geom_line(aes(x = dates,y = hot_av_am),col='darkred') +
  geom_line(aes(x = dates,y = hotdry_av_am),col='darkorange2') +
  xlab('Time')+
  ylab('% of continent experiencing extreme event') +
  ylim(0,100) +
  #for yearly ticks but only labels per 10 years
  scale_x_date(breaks=seq(as.Date('1999-01-01'),as.Date('2025-12-01'),by='1 year'),
               labels = function(x) {
                 yrs <- as.numeric(format(x, "%Y"))
                 ifelse(yrs %% 10 == 0, yrs, "")},expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_line(aes(x = dates,y = dry_1_am_mean),col='darkgreen') +
  geom_line(aes(x = dates,y = hot_1_am_mean),col='darkred') +
  geom_line(aes(x = dates,y = hotdry_1_am_mean),col='darkorange2') +
  geom_line(aes(x = dates,y = dry_av_am),col='darkgreen',linetype='dashed') +
  geom_line(aes(x = dates,y = hot_av_am),col='darkred',linetype='dashed') +
  geom_line(aes(x = dates,y = hotdry_av_am),col='darkorange2',linetype='dashed') +
  xlab('Time')+
  ylab('% of continent experiencing extreme event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#African forests
dry_1_af_mean<-colMeans(dry_1_af)
dry_1_af_min<-apply(dry_1_af,MARGIN = 2,FUN = function(x) min(x))
dry_1_af_max<-apply(dry_1_af,MARGIN = 2,FUN = function(x) max(x))
max(dry_av_af) #80%
year[which(dry_av_af==max(dry_av_af))] #time 299-2024
hot_1_af_mean<-colMeans(hot_1_af)
hot_1_af_min<-apply(hot_1_af,MARGIN = 2,FUN = function(x) min(x))
hot_1_af_max<-apply(hot_1_af,MARGIN = 2,FUN = function(x) max(x))
max(hot_av_af) #99%
year[which(hot_av_af==max(hot_av_af))] #time 232-2019
hotdry_1_af_mean<-colMeans(hotdry_1_af)
hotdry_1_af_min<-apply(hotdry_1_af,MARGIN = 2,FUN = function(x) min(x))
hotdry_1_af_max<-apply(hotdry_1_af,MARGIN = 2,FUN = function(x) max(x))
max(hotdry_av_af) #54%
year[which(hotdry_av_af==max(hotdry_av_af))] #time 299-2024
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = dry_1_af_min, ymax = dry_1_af_max),fill='darkseagreen',alpha=0.5)+
  geom_ribbon(aes(x = dates,ymin = hot_1_af_min, ymax = hot_1_af_max),fill='mistyrose',alpha=0.5)+
  geom_ribbon(aes(x = dates,ymin = hotdry_1_af_min, ymax = hotdry_1_af_max),fill='papayawhip',alpha=0.5)+
  geom_line(aes(x = dates,y = dry_av_af),col='darkgreen') +
  geom_line(aes(x = dates,y = hot_av_af),col='darkred') +
  geom_line(aes(x = dates,y = hotdry_av_af),col='darkorange2') +
  xlab('Time')+
  ylab('% of continent experiencing extreme event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_line(aes(x = dates,y = dry_1_af_mean),col='darkgreen') +
  geom_line(aes(x = dates,y = hot_1_af_mean),col='darkred') +
  geom_line(aes(x = dates,y = hotdry_1_af_mean),col='darkorange2') +
  geom_line(aes(x = dates,y = dry_av_af),col='darkgreen',linetype='dashed') +
  geom_line(aes(x = dates,y = hot_av_af),col='darkred',linetype='dashed') +
  geom_line(aes(x = dates,y = hotdry_av_af),col='darkorange2',linetype='dashed') +
  xlab('Time')+
  ylab('% of continent experiencing extreme event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#Asian forests
dry_1_as_mean<-colMeans(dry_1_as)
dry_1_as_min<-apply(dry_1_as,MARGIN = 2,FUN = function(x) min(x))
dry_1_as_max<-apply(dry_1_as,MARGIN = 2,FUN = function(x) max(x))
max(dry_av_as) #51%
year[which(dry_av_as==max(dry_av_as))] #time 190-2015
hot_1_as_mean<-colMeans(hot_1_as)
hot_1_as_min<-apply(hot_1_as,MARGIN = 2,FUN = function(x) min(x))
hot_1_as_max<-apply(hot_1_as,MARGIN = 2,FUN = function(x) max(x))
max(hot_av_as) #95%
year[which(hot_av_as==max(hot_av_as))] #time 299-2024
hotdry_1_as_mean<-colMeans(hotdry_1_as)
hotdry_1_as_min<-apply(hotdry_1_as,MARGIN = 2,FUN = function(x) min(x))
hotdry_1_as_max<-apply(hotdry_1_as,MARGIN = 2,FUN = function(x) max(x))
max(hotdry_av_as) #33%
year[which(hotdry_av_as==max(hotdry_av_as))] #time 189-2015
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = dry_1_as_min, ymax = dry_1_as_max),fill='darkseagreen',alpha=0.5)+
  geom_ribbon(aes(x = dates,ymin = hot_1_as_min, ymax = hot_1_as_max),fill='mistyrose',alpha=0.5)+
  geom_ribbon(aes(x = dates,ymin = hotdry_1_as_min, ymax = hotdry_1_as_max),fill='papayawhip',alpha=0.5)+
  geom_line(aes(x = dates,y = dry_av_as),col='darkgreen') +
  geom_line(aes(x = dates,y = hot_av_as),col='darkred') +
  geom_line(aes(x = dates,y = hotdry_av_as),col='darkorange2') +
  xlab('Time')+
  ylab('% of continent experiencing extreme event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_line(aes(x = dates,y = dry_1_as_mean),col='darkgreen') +
  geom_line(aes(x = dates,y = hot_1_as_mean),col='darkred') +
  geom_line(aes(x = dates,y = hotdry_1_as_mean),col='darkorange2') +
  geom_line(aes(x = dates,y = dry_av_as),col='darkgreen',linetype='dashed') +
  geom_line(aes(x = dates,y = hot_av_as),col='darkred',linetype='dashed') +
  geom_line(aes(x = dates,y = hotdry_av_as),col='darkorange2',linetype='dashed') +
  xlab('Time')+
  ylab('% of continent experiencing extreme event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#pantropical forests
dry_1_pan_mean<-colMeans(dry_1_pan)
dry_1_pan_min<-apply(dry_1_pan,MARGIN = 2,FUN = function(x) min(x))
dry_1_pan_max<-apply(dry_1_pan,MARGIN = 2,FUN = function(x) max(x))
dry_1_pan_diff<-dry_1_pan_max-dry_1_pan_min
max(dry_av_pan) #51%
year[which(dry_av_pan==max(dry_av_pan))] #time 299-2024
hot_1_pan_mean<-colMeans(hot_1_pan)
hot_1_pan_min<-apply(hot_1_pan,MARGIN = 2,FUN = function(x) min(x))
hot_1_pan_max<-apply(hot_1_pan,MARGIN = 2,FUN = function(x) max(x))
hot_1_pan_diff<-hot_1_pan_max-hot_1_pan_min
max(hot_av_pan) #97%
year[which(hot_av_pan==max(hot_av_pan))] #time 297-2024
hotdry_1_pan_mean<-colMeans(hotdry_1_pan)
hotdry_1_pan_min<-apply(hotdry_1_pan,MARGIN = 2,FUN = function(x) min(x))
hotdry_1_pan_max<-apply(hotdry_1_pan,MARGIN = 2,FUN = function(x) max(x))
hotdry_1_pan_diff<-hotdry_1_pan_max-hotdry_1_pan_min
max(hotdry_av_pan) #44%
year[which(hotdry_av_pan==max(hotdry_av_pan))] #time 190-2015
ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_ribbon(aes(x = dates,ymin = dry_1_pan_min, ymax = dry_1_pan_max),fill='darkseagreen',alpha=0.5)+
  geom_ribbon(aes(x = dates,ymin = hot_1_pan_min, ymax = hot_1_pan_max),fill='mistyrose',alpha=0.5)+
  geom_ribbon(aes(x = dates,ymin = hotdry_1_pan_min, ymax = hotdry_1_pan_max),fill='papayawhip',alpha=0.5)+
  geom_line(aes(x = dates,y = dry_av_pan),col='darkgreen') +
  geom_line(aes(x = dates,y = hot_av_pan),col='darkred') +
  geom_line(aes(x = dates,y = hotdry_av_pan),col='darkorange2') +
  xlab('Time')+
  ylab('% of continent experiencing extreme event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot()+
  geom_rect(aes(xmin=dates[32], xmax=dates[39], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[55], xmax=dates[57], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[59], xmax=dates[60], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[62], xmax=dates[63], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[80], xmax=dates[85], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[115], xmax=dates[116], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[118], xmax=dates[124], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[185], xmax=dates[197], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[224], xmax=dates[226], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[230], xmax=dates[231], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[239], xmax=dates[239], ymin=0, ymax=100),fill='lightgrey')+
  geom_rect(aes(xmin=dates[283], xmax=dates[291], ymin=0, ymax=100),fill='lightgrey')+
  geom_line(aes(x = dates,y = dry_1_pan_mean),col='darkgreen') +
  geom_line(aes(x = dates,y = hot_1_pan_mean),col='darkred') +
  geom_line(aes(x = dates,y = hotdry_1_pan_mean),col='darkorange2') +
  geom_line(aes(x = dates,y = dry_av_pan),col='darkgreen',linetype='dashed') +
  geom_line(aes(x = dates,y = hot_av_pan),col='darkred',linetype='dashed') +
  geom_line(aes(x = dates,y = hotdry_av_pan),col='darkorange2',linetype='dashed') +
  xlab('Time')+
  ylab('% of continent experiencing extreme event') +
  ylim(0,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



# 5) Make map of event occurrence during 2023-24 for all products (FIGURE 2) ####
prec_3imerg_an<-stack('df.3IMERG.Tropics.climate.rspd.0024.prec.st.anomalies.removeEN.tif')
prec_an_df<-as.data.frame(prec_3imerg_an,xy=T)
prec_coord<-prec_an_df[,1:13]
#Dry events ####
prec_3imerg_dry1_df<-read.csv('df.3IMERG.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_chirps_dry1_df<-read.csv('df.chirps.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_cru_dry1_df<-read.csv('df.CRU.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_era5_dry1_df<-read.csv('df.ERA5.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_gldas_dry1_df<-read.csv('df.GLDAS.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_gpcc_dry1_df<-read.csv('df.GPCC.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_merra2_dry1_df<-read.csv('df.MERRA2.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_mswep_dry1_df<-read.csv('df.MSWEP.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_ncep_dry1_df<-read.csv('df.MSWEP.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_av_dry1_df<-read.csv('df.Tropics.climate.rspd.0024.prec.mean.weightedR2.dry.1SD.removeEN.csv')[,-1]

#only look at 2023-2024
prec_3imerg_dry1_df<-prec_3imerg_dry1_df[,277:300]
prec_chirps_dry1_df<-prec_chirps_dry1_df[,277:300]
prec_cru_dry1_df<-prec_cru_dry1_df[,277:300]
prec_era5_dry1_df<-prec_era5_dry1_df[,277:300]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,277:300]
prec_gpcc_dry1_df<-prec_gpcc_dry1_df[,277:300]
prec_mswep_dry1_df<-prec_mswep_dry1_df[,277:300]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,277:300]
prec_av_dry1_df<-prec_av_dry1_df[,277:300]
#OR only look at 2015-2016
prec_3imerg_dry1_df<-prec_3imerg_dry1_df[,169:192]
prec_chirps_dry1_df<-prec_chirps_dry1_df[,169:192]
prec_cru_dry1_df<-prec_cru_dry1_df[,169:192]
prec_era5_dry1_df<-prec_era5_dry1_df[,169:192]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,169:192]
prec_gpcc_dry1_df<-prec_gpcc_dry1_df[,169:192]
prec_mswep_dry1_df<-prec_mswep_dry1_df[,169:192]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,169:192]
prec_av_dry1_df<-prec_av_dry1_df[,169:192]
#find area
prec_3imerg_nonapixels<-which(rowSums(!is.na(prec_3imerg_dry1_df))>0)
prec_chirps_nonapixels<-which(rowSums(!is.na(prec_chirps_dry1_df))>0)
prec_cru_nonapixels<-which(rowSums(!is.na(prec_cru_dry1_df))>0)
prec_era5_nonapixels<-which(rowSums(!is.na(prec_era5_dry1_df))>0)
prec_gldas_nonapixels<-which(rowSums(!is.na(prec_gldas_dry1_df))>0)
prec_gpcc_nonapixels<-which(rowSums(!is.na(prec_gpcc_dry1_df))>0)
prec_mswep_nonapixels<-which(rowSums(!is.na(prec_mswep_dry1_df))>0)
prec_ncep_nonapixels<-which(rowSums(!is.na(prec_ncep_dry1_df))>0)
prec_av_nonapixels<-which(rowSums(!is.na(prec_av_dry1_df))>0)
prec_dry_coord<-prec_coord
prec_dry_coord[!is.na(prec_dry_coord[,3]),3:11]<-0
prec_dry_coord[prec_3imerg_nonapixels,3]<-1
prec_dry_coord[prec_chirps_nonapixels,4]<-1
prec_dry_coord[prec_cru_nonapixels,5]<-1
prec_dry_coord[prec_era5_nonapixels,6]<-1
prec_dry_coord[prec_gldas_nonapixels,7]<-1
prec_dry_coord[prec_gpcc_nonapixels,8]<-1
prec_dry_coord[prec_mswep_nonapixels,9]<-1
prec_dry_coord[prec_ncep_nonapixels,10]<-1
prec_dry_coord[prec_av_nonapixels,11]<-1
colnames(prec_dry_coord)[1:11]<-c('x','y',c(paste0('dry_',1:9)))
#convert to 0-1 raster to show where dry periods occurred during EN
prec_dry <- rasterFromXYZ(prec_dry_coord[, c('x','y',c(paste0('dry_',1:9)))])
crs(prec_dry)<-'+proj=longlat +datum=WGS84 +no_defs'
prec_dry<-rast(prec_dry)
count_ones <- app(prec_dry[[1:8]], function(x) sum(x == 1, na.rm = F))
plot(count_ones)
writeRaster(count_ones,'dry_events_1516_countones.tif',overwrite=T)
plot(prec_dry[[9]])
writeRaster(prec_dry[[9]],'dry_events_1516_occ_weightedav.tif')
poly <- as.polygons(count_ones)
poly_sf <- st_as_sf(poly)
names(poly_sf)[1] <- "count"
av_poly <- as.polygons(prec_dry[[9]])
av_poly_sf <- st_as_sf(av_poly)
stipple_sf <- av_poly_sf %>% filter(dry_9 == 1)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
poly_sf$fill_color <- alpha(scales::col_numeric(c("honeydew2", "darkgreen"), 
                                                domain = c(0,8))(poly_sf$count),
                            poly_sf$count / 8)
ggplot() +
  geom_sf(data = poly_sf, aes(fill = fill_color), color = NA) +
  scale_fill_identity() +
  geom_sf(data = stipple_sf, linewidth = 0.3, color = "black",fill='NA') +
  geom_sf(data = world,fill = NA, color = "grey40", linewidth = 0.1) +
  theme_minimal()
plot(count_ones,col=colorRampPalette(c("honeydew2", "darkgreen"))(8))

#Hot events ####
prec_cru_dry1_df<-read.csv('df.CRU.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_era5_dry1_df<-read.csv('df.ERA5.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_gldas_dry1_df<-read.csv('df.GLDAS.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_merra2_dry1_df<-read.csv('df.MERRA2.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_berk_dry1_df<-read.csv('df.Berk.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_cams_dry1_df<-read.csv('df.CAMS.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_ncep_dry1_df<-read.csv('df.NCEP.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_av_dry1_df<-read.csv('df.Tropics.climate.rspd.0024.tav.mean.weightedR2.hot.1SD.removeEN.csv')[,-1]
#only look at 2023-2024
prec_cru_dry1_df<-prec_cru_dry1_df[,277:300]
prec_era5_dry1_df<-prec_era5_dry1_df[,277:300]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,277:300]
prec_berk_dry1_df<-prec_berk_dry1_df[,277:300]
prec_cams_dry1_df<-prec_cams_dry1_df[,277:300]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,277:300]
prec_av_dry1_df<-prec_av_dry1_df[,277:300]
#OR only look at 2015-2016
prec_cru_dry1_df<-prec_cru_dry1_df[,169:192]
prec_era5_dry1_df<-prec_era5_dry1_df[,169:192]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,169:192]
prec_berk_dry1_df<-prec_berk_dry1_df[,169:192]
prec_cams_dry1_df<-prec_cams_dry1_df[,169:192]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,169:192]
prec_av_dry1_df<-prec_av_dry1_df[,169:192]

prec_cru_nonapixels<-which(rowSums(!is.na(prec_cru_dry1_df))>0)
prec_era5_nonapixels<-which(rowSums(!is.na(prec_era5_dry1_df))>0)
prec_gldas_nonapixels<-which(rowSums(!is.na(prec_gldas_dry1_df))>0)
prec_berk_nonapixels<-which(rowSums(!is.na(prec_berk_dry1_df))>0)
prec_cams_nonapixels<-which(rowSums(!is.na(prec_cams_dry1_df))>0)
prec_ncep_nonapixels<-which(rowSums(!is.na(prec_ncep_dry1_df))>0)
prec_av_nonapixels<-which(rowSums(!is.na(prec_av_dry1_df))>0)
prec_dry_coord<-prec_coord
prec_dry_coord[!is.na(prec_dry_coord[,3]),3:9]<-0
prec_dry_coord[prec_cru_nonapixels,3]<-1
prec_dry_coord[prec_era5_nonapixels,4]<-1
prec_dry_coord[prec_gldas_nonapixels,5]<-1
prec_dry_coord[prec_berk_nonapixels,6]<-1
prec_dry_coord[prec_cams_nonapixels,7]<-1
prec_dry_coord[prec_ncep_nonapixels,8]<-1
prec_dry_coord[prec_av_nonapixels,9]<-1
colnames(prec_dry_coord)[1:9]<-c('x','y',c(paste0('dry_',1:7)))
prec_dry <- rasterFromXYZ(prec_dry_coord[, c('x','y',c(paste0('dry_',1:7)))])
crs(prec_dry)<-'+proj=longlat +datum=WGS84 +no_defs'
prec_dry<-rast(prec_dry)
count_ones <- app(prec_dry[[1:6]], function(x) sum(x == 1, na.rm = F))
plot(count_ones)
writeRaster(count_ones,'hot_events_2324_countones.tif',overwrite=T)
plot(prec_dry[[7]])
writeRaster(prec_dry[[7]],'hot_events_2324_occ_weightedav.tif')

poly <- as.polygons(count_ones)
poly_sf <- st_as_sf(poly)
names(poly_sf)[1] <- "count"
av_poly <- as.polygons(prec_dry[[7]])
av_poly_sf <- st_as_sf(av_poly)
stipple_sf <- av_poly_sf %>% filter(dry_7 == 1)

poly_sf$fill_color <- alpha(scales::col_numeric(c("mistyrose", "darkred"), 
                                                domain = c(0,6))(poly_sf$count),
                            poly_sf$count / 6)
ggplot() +
  geom_sf(data = poly_sf, aes(fill = fill_color), color = NA) +
  scale_fill_identity() +
  geom_sf(data = stipple_sf, linewidth = 0.3, color = "black",fill='NA') +
  geom_sf(data = world,fill = NA, color = "grey40", linewidth = 0.1) +
  theme_minimal()
plot(raster(count_ones),col=colorRampPalette(c("mistyrose", "darkred"))(6),zlim=c(0,6))

#Hot-dry events ####
prec_3imerg_dry1_df<-read.csv('df.3IMERG.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_chirps_dry1_df<-read.csv('df.chirps.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_cru_dry1_df<-read.csv('df.CRU.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_era5_dry1_df<-read.csv('df.ERA5.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_gldas_dry1_df<-read.csv('df.GLDAS.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_gpcc_dry1_df<-read.csv('df.GPCC.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_mswep_dry1_df<-read.csv('df.MSWEP.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_ncep_dry1_df<-read.csv('df.NCEP.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_berk_dry1_df<-read.csv('df.Berk.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_cams_dry1_df<-read.csv('df.CAMS.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_av_dry1_df<-read.csv('df.Tropics.climate.rspd.0024.mean.weightedR2.hotdry.1SD.removeEN.csv')[,-1]

#only look at 2023-2024
prec_3imerg_dry1_df<-prec_3imerg_dry1_df[,277:300]
prec_chirps_dry1_df<-prec_chirps_dry1_df[,277:300]
prec_cru_dry1_df<-prec_cru_dry1_df[,277:300]
prec_era5_dry1_df<-prec_era5_dry1_df[,277:300]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,277:300]
prec_gpcc_dry1_df<-prec_gpcc_dry1_df[,277:300]
prec_mswep_dry1_df<-prec_mswep_dry1_df[,277:300]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,277:300]
prec_berk_dry1_df<-prec_berk_dry1_df[,277:300]
prec_cams_dry1_df<-prec_cams_dry1_df[,277:300]
prec_av_dry1_df<-prec_av_dry1_df[,277:300]
#only look at 2015-2016
prec_3imerg_dry1_df<-prec_3imerg_dry1_df[,169:192]
prec_chirps_dry1_df<-prec_chirps_dry1_df[,169:192]
prec_cru_dry1_df<-prec_cru_dry1_df[,169:192]
prec_era5_dry1_df<-prec_era5_dry1_df[,169:192]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,169:192]
prec_gpcc_dry1_df<-prec_gpcc_dry1_df[,169:192]
prec_mswep_dry1_df<-prec_mswep_dry1_df[,169:192]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,169:192]
prec_berk_dry1_df<-prec_berk_dry1_df[,169:192]
prec_cams_dry1_df<-prec_cams_dry1_df[,169:192]
prec_av_dry1_df<-prec_av_dry1_df[,169:192]
#find area
prec_3imerg_nonapixels<-which(rowSums(!is.na(prec_3imerg_dry1_df))>0)
prec_chirps_nonapixels<-which(rowSums(!is.na(prec_chirps_dry1_df))>0)
prec_cru_nonapixels<-which(rowSums(!is.na(prec_cru_dry1_df))>0)
prec_era5_nonapixels<-which(rowSums(!is.na(prec_era5_dry1_df))>0)
prec_gldas_nonapixels<-which(rowSums(!is.na(prec_gldas_dry1_df))>0)
prec_gpcc_nonapixels<-which(rowSums(!is.na(prec_gpcc_dry1_df))>0)
prec_mswep_nonapixels<-which(rowSums(!is.na(prec_mswep_dry1_df))>0)
prec_ncep_nonapixels<-which(rowSums(!is.na(prec_ncep_dry1_df))>0)
prec_berk_nonapixels<-which(rowSums(!is.na(prec_berk_dry1_df))>0)
prec_cams_nonapixels<-which(rowSums(!is.na(prec_cams_dry1_df))>0)
prec_av_nonapixels<-which(rowSums(!is.na(prec_av_dry1_df))>0)

prec_dry_coord<-prec_coord
prec_dry_coord[!is.na(prec_dry_coord[,3]),3:13]<-0
prec_dry_coord[prec_3imerg_nonapixels,3]<-1
prec_dry_coord[prec_chirps_nonapixels,4]<-1
prec_dry_coord[prec_cru_nonapixels,5]<-1
prec_dry_coord[prec_era5_nonapixels,6]<-1
prec_dry_coord[prec_gldas_nonapixels,7]<-1
prec_dry_coord[prec_gpcc_nonapixels,8]<-1
prec_dry_coord[prec_mswep_nonapixels,9]<-1
prec_dry_coord[prec_ncep_nonapixels,10]<-1
prec_dry_coord[prec_berk_nonapixels,11]<-1
prec_dry_coord[prec_cams_nonapixels,12]<-1
prec_dry_coord[prec_av_nonapixels,13]<-1
colnames(prec_dry_coord)[1:13]<-c('x','y',c(paste0('dry_',1:11)))

#convert to 0-1 raster to show where dry periods occurred during EN
prec_dry <- rasterFromXYZ(prec_dry_coord[, c('x','y',c(paste0('dry_',1:11)))])
crs(prec_dry)<-'+proj=longlat +datum=WGS84 +no_defs'
prec_dry<-rast(prec_dry)
count_ones <- app(prec_dry[[1:10]], function(x) sum(x == 1, na.rm = F))
plot(count_ones)
writeRaster(count_ones,'./RDSfiles/hotdry_events_1516_countones.tif',overwrite=T)
plot(prec_dry[[11]])
writeRaster(prec_dry[[11]],'./RDSfiles/hotdry_events_1516_occ_weightedav.tif')

poly <- as.polygons(count_ones)
poly_sf <- st_as_sf(poly)
names(poly_sf)[1] <- "count"
av_poly <- as.polygons(prec_dry[[11]])
av_poly_sf <- st_as_sf(av_poly)
stipple_sf <- av_poly_sf %>% filter(dry_11 == 1)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
poly_sf$fill_color <- alpha(scales::col_numeric(c("papayawhip", "darkorange2"), 
                                                domain = c(0,10))(poly_sf$count),
                            poly_sf$count / 10)
ggplot() +
  geom_sf(data = poly_sf, aes(fill = fill_color), color = NA) +
  scale_fill_identity() +
  geom_sf(data = stipple_sf, linewidth = 0.3, color = "black",fill='NA') +
  geom_sf(data = world,fill = NA, color = "grey40", linewidth = 0.1) +
  theme_minimal()
plot(count_ones,col=colorRampPalette(c("papayawhip", "darkorange2"))(10))

# 5-2) Make density distribution for 23-24 and 15-16 (INSETS FIGURE 2) ####
dry_1516<-rast('./RDSfiles/dry_events_1516_occ_weightedav.tif')
hot_1516<-rast('./RDSfiles/hot_events_1516_occ_weightedav.tif')
hotdry_1516<-rast('./RDSfiles/hotdry_events_1516_occ_weightedav.tif')
dry_2324<-rast('./RDSfiles/dry_events_2324_occ_weightedav.tif')
hot_2324<-rast('./RDSfiles/hot_events_2324_occ_weightedav.tif')
hotdry_2324<-rast('./RDSfiles/hotdry_events_2324_occ_weightedav.tif')

#crop to continent
extent_am<-ext(-120,-25,-30.5,30.5)
extent_af<-ext(-25,50,-30.5,30.5)
extent_as<-ext(50,180,-30.5,30.5)
extent<-extent_as

dry_1516<-crop(dry_1516,extent)
hot_1516<-crop(hot_1516,extent)
hotdry_1516<-crop(hotdry_1516,extent)
dry_2324<-crop(dry_2324,extent)
hot_2324<-crop(hot_2324,extent)
hotdry_2324<-crop(hotdry_2324,extent)

dry_1516_val<-as.data.frame(dry_1516,xy=T)
dry_2324_val<-as.data.frame(dry_2324,xy=T)
dry_comb<-c(dry_1516_val[,3],dry_2324_val[,3])
dry_comb<-data.frame(val=dry_comb,year=c(rep('1516',nrow(dry_1516_val)),rep('2324',nrow(dry_2324_val))))
dry_plot<-ggplot(dry_comb, aes(x=val,fill=year) ) +
  geom_bar(aes(y = (..count..)/sum(..count..)*200),position = 'dodge') +
  scale_fill_manual(values = c("1516" = "grey60",
                               "2324" = 'darkgreen')) +
  ylim(0,100)
sum(dry_comb$year=='1516' & dry_comb$val==1)/sum(dry_comb$year=='1516')*100
sum(dry_comb$year=='2324' & dry_comb$val==1)/sum(dry_comb$year=='2324')*100

hot_1516_val<-as.data.frame(hot_1516,xy=T)
hot_2324_val<-as.data.frame(hot_2324,xy=T)
hot_comb<-c(hot_1516_val[,3],hot_2324_val[,3])
hot_comb<-data.frame(val=hot_comb,year=c(rep('1516',nrow(hot_1516_val)),rep('2324',nrow(hot_2324_val))))
hot_plot<-ggplot(hot_comb, aes(x=val,fill=year) ) +
  geom_bar(aes(y = (..count..)/sum(..count..)*200),position = 'dodge') +
  scale_fill_manual(values = c("1516" = "grey60",
                               "2324" = 'darkred'))+
  ylim(0,100)
sum(hot_comb$year=='1516' & hot_comb$val==1)/sum(hot_comb$year=='1516')*100
sum(hot_comb$year=='2324' & hot_comb$val==1)/sum(hot_comb$year=='2324')*100

hotdry_1516_val<-as.data.frame(hotdry_1516,xy=T)
hotdry_2324_val<-as.data.frame(hotdry_2324,xy=T)
hotdry_comb<-c(hotdry_1516_val[,3],hotdry_2324_val[,3])
hotdry_comb<-data.frame(val=hotdry_comb,year=c(rep('1516',nrow(hotdry_1516_val)),rep('2324',nrow(hotdry_2324_val))))
hotdry_plot<-ggplot(hotdry_comb, aes(x=val,fill=year) ) +
  geom_bar(aes(y = (..count..)/sum(..count..)*200),position = 'dodge') +
  scale_fill_manual(values = c("1516" = "grey60",
                               "2324" = 'darkorange2'))+
  ylim(0,100)
sum(hotdry_comb$year=='1516' & hotdry_comb$val==1)/sum(hotdry_comb$year=='1516')*100
sum(hotdry_comb$year=='2324' & hotdry_comb$val==1)/sum(hotdry_comb$year=='2324')*100

ggarrange(dry_plot, hot_plot, hotdry_plot)

# 6) Make map of event duration during 2023-24 for all products (FIGURE 3) ####
prec_3imerg_an<-stack('df.3IMERG.Tropics.climate.rspd.0024.prec.st.anomalies.removeEN.tif')
prec_an_df<-as.data.frame(prec_3imerg_an,xy=T)
prec_coord<-prec_an_df[,1:12]

#Dry events ####
prec_3imerg_dry1_df<-read.csv('df.3IMERG.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_chirps_dry1_df<-read.csv('df.chirps.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_cru_dry1_df<-read.csv('df.CRU.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_era5_dry1_df<-read.csv('df.ERA5.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_gldas_dry1_df<-read.csv('df.GLDAS.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_gpcc_dry1_df<-read.csv('df.GPCC.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_merra2_dry1_df<-read.csv('df.MERRA2.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_mswep_dry1_df<-read.csv('df.MSWEP.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_ncep_dry1_df<-read.csv('df.NCEP.Tropics.climate.rspd.0024.dry.1SD.removeEN.csv')[,-1]
prec_av_dry1_df<-read.csv('df.Tropics.climate.rspd.0024.prec.mean.weightedR2.dry.1SD.removeEN.csv')[,-1]

#only look at 2023-2024
prec_3imerg_dry1_df<-prec_3imerg_dry1_df[,277:300]
prec_chirps_dry1_df<-prec_chirps_dry1_df[,277:300]
prec_cru_dry1_df<-prec_cru_dry1_df[,277:300]
prec_era5_dry1_df<-prec_era5_dry1_df[,277:300]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,277:300]
prec_gpcc_dry1_df<-prec_gpcc_dry1_df[,277:300]
prec_mswep_dry1_df<-prec_mswep_dry1_df[,277:300]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,277:300]
prec_av_dry1_df<-prec_av_dry1_df[,277:300]
#OR only look at 2015-2016
prec_3imerg_dry1_df<-prec_3imerg_dry1_df[,169:192]
prec_chirps_dry1_df<-prec_chirps_dry1_df[,169:192]
prec_cru_dry1_df<-prec_cru_dry1_df[,169:192]
prec_era5_dry1_df<-prec_era5_dry1_df[,169:192]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,169:192]
prec_gpcc_dry1_df<-prec_gpcc_dry1_df[,169:192]
prec_mswep_dry1_df<-prec_mswep_dry1_df[,169:192]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,169:192]
prec_av_dry1_df<-prec_av_dry1_df[,169:192]
#get duration
prec_3imerg_duration<-apply(prec_3imerg_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_chirps_duration<-apply(prec_chirps_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_cru_duration<-apply(prec_cru_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_era5_duration<-apply(prec_era5_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_gldas_duration<-apply(prec_gldas_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_gpcc_duration<-apply(prec_gpcc_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_mswep_duration<-apply(prec_mswep_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_ncep_duration<-apply(prec_ncep_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_av_duration<-apply(prec_av_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))

prec_dry_dur<-prec_coord
prec_dry_dur[,3]<-prec_3imerg_duration
prec_dry_dur[,4]<-prec_chirps_duration
prec_dry_dur[,5]<-prec_cru_duration
prec_dry_dur[,6]<-prec_era5_duration
prec_dry_dur[,7]<-prec_gldas_duration
prec_dry_dur[,8]<-prec_gpcc_duration
prec_dry_dur[,9]<-prec_mswep_duration
prec_dry_dur[,10]<-prec_ncep_duration
prec_dry_dur[is.na(prec_coord[,3]),3:10]<-NA
colnames(prec_dry_dur)<-c('x','y',c(paste0('dur',1:8)))
prec_dry_dur_r <- rasterFromXYZ(prec_dry_dur[, c('x','y',c(paste0('dur',1:8)))])
crs(prec_dry_dur_r)<-'+proj=longlat +datum=WGS84 +no_defs'
plot(prec_dry_dur_r)

#total mean duration per product to put in table
colMeans(prec_dry_dur[,3:11],na.rm=T)
prec_av_duration[is.na(prec_coord[,3])]<-NA
mean(prec_av_duration,na.rm=T)
#continental mean duration per product to put in table
dur_am<-prec_dry_dur[prec_dry_dur$x<(-25),]
dur_af<-prec_dry_dur[prec_dry_dur$x>(-25)&prec_dry_dur$x<50,]
dur_as<-prec_dry_dur[prec_dry_dur$x>50,]
colMeans(dur_am[,3:10],na.rm=T)
colMeans(dur_af[,3:10],na.rm=T)
colMeans(dur_as[,3:10],na.rm=T)
mean(prec_av_duration[prec_dry_dur$x<(-25)],na.rm=T)
mean(prec_av_duration[prec_dry_dur$x>(-25)&prec_dry_dur$x<50],na.rm=T)
mean(prec_av_duration[prec_dry_dur$x>50],na.rm=T)
#using weighted average
prec_dry_dur<-prec_coord
prec_dry_dur[,3]<-prec_av_duration
prec_dry_dur_mean <- rasterFromXYZ(prec_dry_dur[,1:3])
crs(prec_dry_dur_mean)<-'+proj=longlat +datum=WGS84 +no_defs'
prec_dry_dur_mean<-mask(prec_dry_dur_mean,prec_3imerg_an[[1]])
prec_dry_dur_mean<-rast(prec_dry_dur_mean)
plot(prec_dry_dur_mean)
writeRaster(prec_dry_dur_mean,'dry_events_1516_weighted_meandur.tif',overwrite=T)

poly <- as.polygons(prec_dry_dur_mean)
poly_sf <- st_as_sf(poly)
names(poly_sf)[1] <- "dur"
summary(poly_sf$dur)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
poly_sf$fill_color <- alpha(scales::col_numeric(c("honeydew2", "darkgreen"), 
                                                domain = c(0,24))(poly_sf$dur),
                            poly_sf$dur / 24)
pal <- colorRampPalette(c("honeydew2","darkgreen"))
ggplot() +
  geom_sf(data = poly_sf, aes(fill = fill_color), color = NA) +
  geom_sf(data = world,fill = NA, color = "grey40",linewidth=0.1) +
  scale_fill_identity() +
  theme_minimal()
plot(raster(prec_dry_dur_mean),col=pal(25),zlim=c(0,24))

#Hot events ####
prec_cru_dry1_df<-read.csv('df.CRU.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_era5_dry1_df<-read.csv('df.ERA5.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_gldas_dry1_df<-read.csv('df.GLDAS.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_merra2_dry1_df<-read.csv('df.MERRA2.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_berk_dry1_df<-read.csv('df.Berk.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_cams_dry1_df<-read.csv('df.CAMS.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_ncep_dry1_df<-read.csv('df.NCEP.Tropics.climate.rspd.0024.hot.1SD.removeEN.csv')[,-1]
prec_av_dry1_df<-read.csv('df.Tropics.climate.rspd.0024.tav.mean.weightedR2.hot.1SD.removeEN.csv')[,-1]
#only look at 2023-2024
prec_cru_dry1_df<-prec_cru_dry1_df[,277:300]
prec_era5_dry1_df<-prec_era5_dry1_df[,277:300]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,277:300]
prec_berk_dry1_df<-prec_berk_dry1_df[,277:300]
prec_cams_dry1_df<-prec_cams_dry1_df[,277:300]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,277:300]
prec_av_dry1_df<-prec_av_dry1_df[,277:300]
#OR only look at 2015-2016
prec_cru_dry1_df<-prec_cru_dry1_df[,169:192]
prec_era5_dry1_df<-prec_era5_dry1_df[,169:192]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,169:192]
prec_berk_dry1_df<-prec_berk_dry1_df[,169:192]
prec_cams_dry1_df<-prec_cams_dry1_df[,169:192]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,169:192]
prec_av_dry1_df<-prec_av_dry1_df[,169:192]

prec_cru_duration<-apply(prec_cru_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_era5_duration<-apply(prec_era5_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_gldas_duration<-apply(prec_gldas_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_berk_duration<-apply(prec_berk_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_cams_duration<-apply(prec_cams_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_ncep_duration<-apply(prec_ncep_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_av_duration<-apply(prec_av_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))

prec_dry_dur<-prec_coord
prec_dry_dur[,3]<-prec_cru_duration
prec_dry_dur[,4]<-prec_era5_duration
prec_dry_dur[,5]<-prec_gldas_duration
prec_dry_dur[,6]<-prec_berk_duration
prec_dry_dur[,7]<-prec_cams_duration
prec_dry_dur[,8]<-prec_ncep_duration
prec_dry_dur[is.na(prec_coord[,3]),3:8]<-NA
colnames(prec_dry_dur)[1:8]<-c('x','y',c(paste0('dur',1:6)))
prec_dry_dur_r <- rasterFromXYZ(prec_dry_dur[, c('x','y',c(paste0('dur',1:6)))])
crs(prec_dry_dur_r)<-'+proj=longlat +datum=WGS84 +no_defs'
plot(prec_dry_dur_r)
#total mean duration per product to put in table
colMeans(prec_dry_dur[,3:8],na.rm=T)
prec_av_duration[is.na(prec_coord[,3])]<-NA
mean(prec_av_duration,na.rm=T)
#continental mean duration per product to put in table
dur_am<-prec_dry_dur[prec_dry_dur$x<(-25),]
dur_af<-prec_dry_dur[prec_dry_dur$x>(-25)&prec_dry_dur$x<50,]
dur_as<-prec_dry_dur[prec_dry_dur$x>50,]
colMeans(dur_am[,3:8],na.rm=T)
colMeans(dur_af[,3:8],na.rm=T)
colMeans(dur_as[,3:8],na.rm=T)
mean(prec_av_duration[prec_dry_dur$x<(-25)],na.rm=T)
mean(prec_av_duration[prec_dry_dur$x>(-25)&prec_dry_dur$x<50],na.rm=T)
mean(prec_av_duration[prec_dry_dur$x>50],na.rm=T)
#using weighted average
prec_dry_dur<-prec_coord
prec_dry_dur[,3]<-prec_av_duration
prec_dry_dur_mean <- rasterFromXYZ(prec_dry_dur[,1:3])
crs(prec_dry_dur_mean)<-'+proj=longlat +datum=WGS84 +no_defs'
prec_dry_dur_mean<-mask(prec_dry_dur_mean,prec_3imerg_an[[1]])
prec_dry_dur_mean<-rast(prec_dry_dur_mean)
plot(prec_dry_dur_mean)
writeRaster(prec_dry_dur_mean,'hot_events_1516_weighted_meandur.tif',overwrite=T)

poly <- as.polygons(prec_dry_dur_mean)
poly_sf <- st_as_sf(poly)
names(poly_sf)[1] <- "dur"
summary(poly_sf$dur)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
poly_sf$fill_color <- alpha(scales::col_numeric(c("mistyrose", "darkred"), 
                                                domain = c(0,24))(poly_sf$dur),
                            poly_sf$count / 24)
pal <- colorRampPalette(c("mistyrose","darkred"))
ggplot() +
  geom_sf(data = poly_sf, aes(fill = fill_color), color = NA) +
  geom_sf(data = world,fill = NA, color = "grey40",linewidth=0.1) +
  scale_fill_identity() +
  theme_minimal()
plot(raster(prec_dry_dur_mean),col=pal(25),zlim=c(0,24))

#Hot-dry events ####
prec_3imerg_dry1_df<-read.csv('df.3IMERG.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_chirps_dry1_df<-read.csv('df.chirps.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_cru_dry1_df<-read.csv('df.CRU.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_era5_dry1_df<-read.csv('df.ERA5.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_gldas_dry1_df<-read.csv('df.GLDAS.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_gpcc_dry1_df<-read.csv('df.GPCC.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_mswep_dry1_df<-read.csv('df.MSWEP.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_ncep_dry1_df<-read.csv('df.NCEP.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_berk_dry1_df<-read.csv('df.Berk.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_cams_dry1_df<-read.csv('df.CAMS.Tropics.climate.rspd.0024.hotdry.1SD.removeEN.csv')[,-1]
prec_av_dry1_df<-read.csv('df.Tropics.climate.rspd.0024.mean.weightedR2.hotdry.1SD.removeEN.csv')[,-1]
#only look at 2023-2024
prec_3imerg_dry1_df<-prec_3imerg_dry1_df[,277:300]
prec_chirps_dry1_df<-prec_chirps_dry1_df[,277:300]
prec_cru_dry1_df<-prec_cru_dry1_df[,277:300]
prec_era5_dry1_df<-prec_era5_dry1_df[,277:300]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,277:300]
prec_gpcc_dry1_df<-prec_gpcc_dry1_df[,277:300]
prec_mswep_dry1_df<-prec_mswep_dry1_df[,277:300]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,277:300]
prec_berk_dry1_df<-prec_berk_dry1_df[,277:300]
prec_cams_dry1_df<-prec_cams_dry1_df[,277:300]
prec_av_dry1_df<-prec_av_dry1_df[,277:300]
#OR only look at 2015-2016
prec_3imerg_dry1_df<-prec_3imerg_dry1_df[,169:192]
prec_chirps_dry1_df<-prec_chirps_dry1_df[,169:192]
prec_cru_dry1_df<-prec_cru_dry1_df[,169:192]
prec_era5_dry1_df<-prec_era5_dry1_df[,169:192]
prec_gldas_dry1_df<-prec_gldas_dry1_df[,169:192]
prec_gpcc_dry1_df<-prec_gpcc_dry1_df[,169:192]
prec_mswep_dry1_df<-prec_mswep_dry1_df[,169:192]
prec_ncep_dry1_df<-prec_ncep_dry1_df[,169:192]
prec_berk_dry1_df<-prec_berk_dry1_df[,169:192]
prec_cams_dry1_df<-prec_cams_dry1_df[,169:192]
prec_av_dry1_df<-prec_av_dry1_df[,169:192]
#get duration
prec_3imerg_duration<-apply(prec_3imerg_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_chirps_duration<-apply(prec_chirps_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_cru_duration<-apply(prec_cru_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_era5_duration<-apply(prec_era5_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_gldas_duration<-apply(prec_gldas_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_gpcc_duration<-apply(prec_gpcc_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_mswep_duration<-apply(prec_mswep_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_ncep_duration<-apply(prec_ncep_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_berk_duration<-apply(prec_berk_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_cams_duration<-apply(prec_cams_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))
prec_av_duration<-apply(prec_av_dry1_df,MARGIN = 1,FUN = function(x) length(which(!is.na(x))))

prec_dry_dur<-prec_coord
prec_dry_dur[,3]<-prec_3imerg_duration
prec_dry_dur[,4]<-prec_chirps_duration
prec_dry_dur[,5]<-prec_cru_duration
prec_dry_dur[,6]<-prec_era5_duration
prec_dry_dur[,7]<-prec_gldas_duration
prec_dry_dur[,8]<-prec_gpcc_duration
prec_dry_dur[,9]<-prec_mswep_duration
prec_dry_dur[,10]<-prec_ncep_duration
prec_dry_dur[,11]<-prec_berk_duration
prec_dry_dur[,12]<-prec_cams_duration
prec_dry_dur[is.na(prec_coord[,3]),3:12]<-NA
colnames(prec_dry_dur)<-c('x','y',c(paste0('dur',1:10)))
prec_dry_dur_r <- rasterFromXYZ(prec_dry_dur[, c('x','y',c(paste0('dur',1:10)))])
crs(prec_dry_dur_r)<-'+proj=longlat +datum=WGS84 +no_defs'
plot(prec_dry_dur_r)
#total mean duration per product to put in table
colMeans(prec_dry_dur[,3:12],na.rm=T)
prec_av_duration[is.na(prec_coord[,3])]<-NA
mean(prec_av_duration,na.rm=T)
#continental mean duration per product to put in table
dur_am<-prec_dry_dur[prec_dry_dur$x<(-25),]
dur_af<-prec_dry_dur[prec_dry_dur$x>(-25)&prec_dry_dur$x<50,]
dur_as<-prec_dry_dur[prec_dry_dur$x>50,]
colMeans(dur_am[,3:12],na.rm=T)
colMeans(dur_af[,3:12],na.rm=T)
colMeans(dur_as[,3:12],na.rm=T)
mean(prec_av_duration[prec_dry_dur$x<(-25)],na.rm=T)
mean(prec_av_duration[prec_dry_dur$x>(-25)&prec_dry_dur$x<50],na.rm=T)
mean(prec_av_duration[prec_dry_dur$x>50],na.rm=T)
#using weighted average
prec_dry_dur<-prec_coord
prec_dry_dur[,3]<-prec_av_duration
prec_dry_dur_mean <- rasterFromXYZ(prec_dry_dur[,1:3])
crs(prec_dry_dur_mean)<-'+proj=longlat +datum=WGS84 +no_defs'
prec_dry_dur_mean<-mask(prec_dry_dur_mean,prec_3imerg_an[[1]])
prec_dry_dur_mean<-rast(prec_dry_dur_mean)
plot(prec_dry_dur_mean)
writeRaster(prec_dry_dur_mean,'hotdry_events_1516_weighted_meandur.tif',overwrite=T)

poly <- as.polygons(prec_dry_dur_mean)
poly_sf <- st_as_sf(poly)
names(poly_sf)[1] <- "dur"
summary(poly_sf$dur)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
poly_sf$fill_color <- alpha(scales::col_numeric(c("papayawhip", "darkorange2"), 
                                                domain = c(0,24))(poly_sf$dur),
                            poly_sf$count / 24)
pal <- colorRampPalette(c("peachpuff","darkorange2"))
ggplot() +
  geom_sf(data = poly_sf, aes(fill = fill_color), color = NA) +
  geom_sf(data = world,fill = NA, color = "grey40",linewidth=0.1) +
  scale_fill_identity() +
  theme_minimal()
plot(raster(prec_dry_dur_mean),col=pal(25),zlim=c(0,24))

# 6-2) Make duration distribution per continent for 23-24 and 15-16 (INSETS FIGURE 3)####
dry_dur_1516<-rast('dry_events_1516_weighted_meandur.tif')
hot_dur_1516<-rast('hot_events_1516_weighted_meandur.tif')
hotdry_dur_1516<-rast('hotdry_events_1516_weighted_meandur.tif')
dry_dur_2324<-rast('dry_events_2324_weighted_meandur.tif')
hot_dur_2324<-rast('hot_events_2324_weighted_meandur.tif')
hotdry_dur_2324<-rast('hotdry_events_2324_weighted_meandur.tif')

#crop to continent
extent_am<-ext(-120,-25,-30.5,30.5)
extent_af<-ext(-25,50,-30.5,30.5)
extent_as<-ext(50,180,-30.5,30.5)
extent<-extent_as

dry_dur_1516<-crop(dry_dur_1516,extent)
hot_dur_1516<-crop(hot_dur_1516,extent)
hotdry_dur_1516<-crop(hotdry_dur_1516,extent)
dry_dur_2324<-crop(dry_dur_2324,extent)
hot_dur_2324<-crop(hot_dur_2324,extent)
hotdry_dur_2324<-crop(hotdry_dur_2324,extent)

#create violin plots
dry_1516_val<-as.data.frame(dry_dur_1516,xy=T)
dry_2324_val<-as.data.frame(dry_dur_2324,xy=T)
dry_comb<-data.frame(var1516=dry_1516_val[,3],var2324=dry_2324_val[,3])
dryplot<-ggplot(dry_comb, aes(x=x) ) +
  geom_density( aes(x = var2324, y = -..density..), fill= "darkgreen") +
  geom_density( aes(x = var1516, y = ..density..), fill="grey60" )
mean(dry_comb$var1516)
mean(dry_comb$var2324)

hot_1516_val<-as.data.frame(hot_dur_1516,xy=T)
hot_2324_val<-as.data.frame(hot_dur_2324,xy=T)
hot_comb<-data.frame(var1516=hot_1516_val[,3],var2324=hot_2324_val[,3])
hotplot<-ggplot(hot_comb, aes(x=x) ) +
  geom_density( aes(x = var2324, y = -..density..), fill= "darkred") +
  geom_density( aes(x = var1516, y = ..density..), fill="grey60" )
mean(hot_comb$var1516)
mean(hot_comb$var2324)

hotdry_1516_val<-as.data.frame(hotdry_dur_1516,xy=T)
hotdry_2324_val<-as.data.frame(hotdry_dur_2324,xy=T)
hotdry_comb<-data.frame(var1516=hotdry_1516_val[,3],var2324=hotdry_2324_val[,3])
hotdryplot<-ggplot(hotdry_comb, aes(x=x) ) +
  geom_density( aes(x = var2324, y = -..density..), fill= "darkorange2") +
  geom_density( aes(x = var1516, y = ..density..), fill="grey60" )
mean(hotdry_comb$var1516)
mean(hotdry_comb$var2324)

ggarrange(dryplot, hotplot, hotdryplot)
 ()


