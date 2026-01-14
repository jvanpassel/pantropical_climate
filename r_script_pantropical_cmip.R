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
library(pheatmap)
library(svglite)
library(RcppRoll)
library(RRphylo)
library(ragg)

#1) Choose which IPCC regions to include ####
ipcc<-shapefile('referenceRegions.shp')
ipcc_tropics<-crop(ipcc,extent(-180,230,-30.5,30.5))
#remove 4,5,8,9,15,18,21 because only tiny regions included
#remove 13,20,23,25,26 because no forest included
#remove 3,10,24 because less than 10 pixels included
ipcc_tropics$NAME
ipcc_tropics_updated<-ipcc_tropics[-c(3,4,5,8,9,10,13,15,18,20,21,23,24,25,26),]
plot(ipcc_tropics_updated)

#get number of non-NA pixels and total area per region:
prec_3imerg<-stack('df.3IMERG.Tropics.climate.rspd.0024.prec.mask50.tif')
ipcc_prec<-crop(prec_3imerg[[1]],ipcc_tropics)

pixel_area <- function(lat_deg) {
  R <- 6371 # Earth radius in km
  deg_lat_km <- (2 * pi * R) / 360
  dy <- 0.5 * deg_lat_km
  dx <- 0.5 * deg_lat_km * cos(lat_deg * pi / 180)
  area <- dx * dy
  return(area)
}
prec_area<-prec_3imerg[[1]]
latitudes <- yFromCell(prec_area, 1:ncell(prec_area))
areas <- pixel_area(latitudes)
values(prec_area) <- areas
prec_area<-mask(prec_area,prec_3imerg[[1]])

ipcc_prec_list<-list()
ipcc_n_area<-as.data.frame(matrix(nrow=26,ncol=3))
colnames(ipcc_n_area)<-c('name','n_pixels','area')
for(i in 1:length(ipcc_tropics)){
  ipcc_prec_i<-crop(ipcc_prec,ipcc_tropics[i,])
  ipcc_prec_i<-mask(ipcc_prec_i,ipcc_tropics[i,])
  prec_area_i<-crop(prec_area,ipcc_tropics[i,])
  prec_area_i<-mask(prec_area_i,ipcc_tropics[i,])
  prec_area_i_val<-as.data.frame(prec_area_i)
  ipcc_prec_list[[i]]<-ipcc_prec_i
  ipcc_n_area$name[i]<-ipcc_tropics[i,]$NAME
  ipcc_n_area$n_pixels[i]<-sum(!is.na(prec_area_i_val[,1]))
  ipcc_n_area$area[i]<-sum(prec_area_i_val[,1],na.rm=T) #in km²
}

#2) Convert bias corrected NC files into probability ratios ####
tas_brick<-rast('df.Tropics.climate.rspd.0024.tav.mean.weightedR2.tif')
tas_brick<-tas_brick-273.15
pr_brick<-rast('df.Tropics.climate.rspd.0024.prec.mean.weightedR2.tif')

models<-c('Timeseries.MCWD.ACCESS-CM2','Timeseries.MCWD.ACCESS-ESM1-5',
          'Timeseries.MCWD.AWI-CM-1-1-MR','Timeseries.MCWD.AWI-ESM-1-REcoM',
          'Timeseries.MCWD.BCC-CSM2-MR','Timeseries.MCWD.CAMS-CSM1-0',
          'Timeseries.MCWD.CanESM5','Timeseries.MCWD.CanESM5-1',
          'Timeseries.MCWD.CAS-ESM2-0','Timeseries.MCWD.CESM2-WACCM',
          'Timeseries.MCWD.CMCC-CM2-SR5','Timeseries.MCWD.CMCC-ESM2',
          'Timeseries.MCWD.EC-Earth3','Timeseries.MCWD.EC-Earth3-Veg',
          'Timeseries.MCWD.EC-Earth3-Veg-LR','Timeseries.MCWD.FGOALS-f3-L',
          'Timeseries.MCWD.FGOALS-g3','Timeseries.MCWD.FIO-ESM-2-0',
          'Timeseries.MCWD.INM-CM4-8','Timeseries.MCWD.INM-CM5-0',
          'Timeseries.MCWD.IPSL-CM6A-LR','Timeseries.MCWD.KACE-1-0-G',
          'Timeseries.MCWD.MIROC6','Timeseries.MCWD.MPI-ESM1-2-HR',
          'Timeseries.MCWD.MPI-ESM1-2-LR','Timeseries.MCWD.MRI-ESM2-0',
          'Timeseries.MCWD.NESM3','Timeseries.MCWD.NorESM2-LM',
          'Timeseries.MCWD.NorESM2-MM','Timeseries.MCWD.TaiESM1')
model_nrs_370<-c(1:3,5:8,10:13,15,17,19:26,28:30)
ssp<-c(126,245,370,585)
compute_PR <- function(events_control, years_control,
                       events_scenario, years_scenario) {
  P_control  <- events_control / years_control
  P_scenario <- events_scenario / years_scenario
  
  # Avoid division by zero
  if (P_control == 0) {
    return(NA)
  }
  
  PR <- P_scenario / P_control
  return(PR)
}
ipcc<-shapefile('./referenceRegions_ipcc/referenceRegions.shp')
ipcc_tropics<-crop(ipcc,extent(-180,230,-30.5,30.5))
ipcc_tropics_updated<-ipcc_tropics[-c(3,4,5,8,9,10,13,15,18,20,21,23,24,25,26),]
ipcc_tropics_updated<-vect(ipcc_tropics_updated)

# 2.1 Precipitation - rolling window max TS and rolling window PR for all IPCC regions ####
pr_regions_ssp_rollpr_list<-list()
pr_regions_ssp_ts_list<-list()
for(j in 1:length(ipcc_tropics_updated)){
  model_pr_ssp126_list<-list()
  model_pr_ssp245_list<-list()
  model_pr_ssp370_list<-list()
  model_pr_ssp585_list<-list()
  for(i in 1:length(model_nrs_370)){
    model_nrs_i<-model_nrs_370[i]
    ssp126_rasti<-rast(paste0('./RDS_CMIP6/ISIMIP_bias_adj/',models[model_nrs_i],'.ssp',ssp[1],'.prec.mask50_bias_adj.nc'))
    ssp126_rasti<-crop(ssp126_rasti,ipcc_tropics_updated[j,])
    model_pr_ssp126_list[i]<-mask(ssp126_rasti,ipcc_tropics_updated[j,])
    ssp245_rasti<-rast(paste0('./RDS_CMIP6/ISIMIP_bias_adj/',models[model_nrs_i],'.ssp',ssp[2],'.prec.mask50_bias_adj.nc'))
    ssp245_rasti<-crop(ssp245_rasti,ipcc_tropics_updated[j,])
    model_pr_ssp245_list[i]<-mask(ssp245_rasti,ipcc_tropics_updated[j,])
    ssp370_rasti<-rast(paste0('./RDS_CMIP6/ISIMIP_bias_adj/',models[model_nrs_i],'.ssp',ssp[3],'.prec.mask50_bias_adj.nc'))
    ssp370_rasti<-crop(ssp370_rasti,ipcc_tropics_updated[j,])
    model_pr_ssp370_list[i]<-mask(ssp370_rasti,ipcc_tropics_updated[j,])
    ssp585_rasti<-rast(paste0('./RDS_CMIP6/ISIMIP_bias_adj/',models[model_nrs_i],'.ssp',ssp[4],'.prec.mask50_bias_adj.nc'))
    ssp585_rasti<-crop(ssp585_rasti,ipcc_tropics_updated[j,])
    model_pr_ssp585_list[i]<-mask(ssp585_rasti,ipcc_tropics_updated[j,])
  }
  
  pr_ipcc<-mask(pr_brick,ipcc_tropics_updated[j,])
  pr_ipcc<-crop(pr_ipcc,ipcc_tropics_updated[j,])
  years<-rep(1:25,each=12)
  years_ssp<-rep(1:76,each=12)
  
  pr_ssp126_annualmean_list<-lapply(model_pr_ssp126_list,function(x) roll(x,n=12,fun = sum,type='to'))
  pr_ssp126_annualmean_1_list<-lapply(pr_ssp126_annualmean_list,function(x) global(x, fun=mean, na.rm=TRUE))
  pr_ssp126_annualmean_1_min<-lapply(pr_ssp126_annualmean_1_list,function(x) tapply(x$mean,years_ssp,min,na.rm=T))
  pr_ssp245_annualmean_list<-lapply(model_pr_ssp245_list,function(x) roll(x,n=12,fun = sum,type='to'))
  pr_ssp245_annualmean_1_list<-lapply(pr_ssp245_annualmean_list,function(x) global(x, fun=mean, na.rm=TRUE))
  pr_ssp245_annualmean_1_min<-lapply(pr_ssp245_annualmean_1_list,function(x) tapply(x$mean,years_ssp,min,na.rm=T))
  pr_ssp370_annualmean_list<-lapply(model_pr_ssp370_list,function(x) roll(x,n=12,fun = sum,type='to'))
  pr_ssp370_annualmean_1_list<-lapply(pr_ssp370_annualmean_list,function(x) global(x, fun=mean, na.rm=TRUE))
  pr_ssp370_annualmean_1_min<-lapply(pr_ssp370_annualmean_1_list,function(x) tapply(x$mean,years_ssp,min,na.rm=T))
  pr_ssp585_annualmean_list<-lapply(model_pr_ssp585_list,function(x) roll(x,n=12,fun = sum,type='to'))
  pr_ssp585_annualmean_1_list<-lapply(pr_ssp585_annualmean_list,function(x) global(x, fun=mean, na.rm=TRUE))
  pr_ssp585_annualmean_1_min<-lapply(pr_ssp585_annualmean_1_list,function(x) tapply(x$mean,years_ssp,min,na.rm=T))
  pr_anualmean<-roll(pr_ipcc,n=12,fun = sum,type='to')
  pr_annualmean_1<-global(pr_anualmean, fun=mean, na.rm=TRUE)
  pr_annualmean_min<-tapply(pr_annualmean_1$mean,years,min,na.rm=T)
  
  timing<-c(24:25)
  events_control  <- sum(pr_annualmean_min<=min(pr_annualmean_min[timing]),na.rm=T)
  years_control   <- 25
  
  #do rolling window
  ssp126_results_list<-list()
  for (i in seq_along(pr_ssp126_annualmean_1_min)) {
    mean_vec <- pr_ssp126_annualmean_1_min[[i]]
    n <- length(mean_vec)
    events_future <- numeric(n - 25 + 1)
    years_future  <- rep(25, n - 25 + 1)
    for (k in 1:(n - 25 + 1)) {
      window_k <- mean_vec[k:(k + 25 - 1)]
      events_future[k] <- sum(window_k<=min(pr_annualmean_min[timing]),na.rm=T)
    }
    ssp126_results_list[[i]] <- mapply(compute_PR,events_control,years_control,events_future,years_future)
  }
  ssp245_results_list<-list()
  for (i in seq_along(pr_ssp245_annualmean_1_min)) {
    mean_vec <- pr_ssp245_annualmean_1_min[[i]]
    n <- length(mean_vec)
    events_future <- numeric(n - 25 + 1)
    years_future  <- rep(25, n - 25 + 1)
    for (k in 1:(n - 25 + 1)) {
      window_k <- mean_vec[k:(k + 25 - 1)]
      events_future[k] <- sum(window_k<=min(pr_annualmean_min[timing]),na.rm=T)
    }
    ssp245_results_list[[i]] <- mapply(compute_PR,events_control,years_control,events_future,years_future)
  }
  ssp370_results_list<-list()
  for (i in seq_along(pr_ssp370_annualmean_1_min)) {
    mean_vec <- pr_ssp370_annualmean_1_min[[i]]
    n <- length(mean_vec)
    events_future <- numeric(n - 25 + 1)
    years_future  <- rep(25, n - 25 + 1)
    for (k in 1:(n - 25 + 1)) {
      window_k <- mean_vec[k:(k + 25 - 1)]
      events_future[k] <- sum(window_k<=min(pr_annualmean_min[timing]),na.rm=T)
    }
    ssp370_results_list[[i]] <- mapply(compute_PR,events_control,years_control,events_future,years_future)
  }
  ssp585_results_list<-list()
  for (i in seq_along(pr_ssp585_annualmean_1_min)) {
    mean_vec <- pr_ssp585_annualmean_1_min[[i]]
    n <- length(mean_vec)
    events_future <- numeric(n - 25 + 1)
    years_future  <- rep(25, n - 25 + 1)
    for (k in 1:(n - 25 + 1)) {
      window_k <- mean_vec[k:(k + 25 - 1)]
      events_future[k] <- sum(window_k<=min(pr_annualmean_min[timing]),na.rm=T)
    }
    ssp585_results_list[[i]] <- mapply(compute_PR,events_control,years_control,events_future,years_future)
  }
  ssp126_results<-do.call('rbind',ssp126_results_list)
  ssp245_results<-do.call('rbind',ssp245_results_list)
  ssp370_results<-do.call('rbind',ssp370_results_list)
  ssp585_results<-do.call('rbind',ssp585_results_list)
  
  #get mean and min/max PR time series per scenario in df
  ssp_results_pr<-as.data.frame(matrix(nrow=12,ncol=52))
  ssp_results_pr[1,]<-colMeans(ssp126_results)
  ssp_results_pr[2,]<-apply(ssp126_results, MARGIN=2, min)
  ssp_results_pr[3,]<-apply(ssp126_results, MARGIN=2, max)
  ssp_results_pr[4,]<-colMeans(ssp245_results)
  ssp_results_pr[5,]<-apply(ssp245_results, MARGIN=2, min)
  ssp_results_pr[6,]<-apply(ssp245_results, MARGIN=2, max)
  ssp_results_pr[7,]<-colMeans(ssp370_results)
  ssp_results_pr[8,]<-apply(ssp370_results, MARGIN=2, min)
  ssp_results_pr[9,]<-apply(ssp370_results, MARGIN=2, max)
  ssp_results_pr[10,]<-colMeans(ssp585_results)
  ssp_results_pr[11,]<-apply(ssp585_results, MARGIN=2, min)
  ssp_results_pr[12,]<-apply(ssp585_results, MARGIN=2, max)
  pr_regions_ssp_rollpr_list[[j]]<-ssp_results_pr
  
  #get mean and min/max time series per scenario in df
  pr_annualmean_df<-as.data.frame(matrix(nrow=13,ncol=101))
  pr_annualmean_df[1,1:25]<-pr_annualmean_min
  pr_ssp126_annualmean_1<-do.call('cbind',pr_ssp126_annualmean_1_min)
  pr_annualmean_df[2,26:101]<-rowMeans(pr_ssp126_annualmean_1)
  pr_annualmean_df[3,26:101]<-apply(pr_ssp126_annualmean_1, MARGIN=1, min)
  pr_annualmean_df[4,26:101]<-apply(pr_ssp126_annualmean_1, MARGIN=1, max)
  pr_ssp245_annualmean_1<-do.call('cbind',pr_ssp245_annualmean_1_min)
  pr_annualmean_df[5,26:101]<-rowMeans(pr_ssp245_annualmean_1)
  pr_annualmean_df[6,26:101]<-apply(pr_ssp245_annualmean_1, MARGIN=1, min)
  pr_annualmean_df[7,26:101]<-apply(pr_ssp245_annualmean_1, MARGIN=1, max)
  pr_ssp370_annualmean_1<-do.call('cbind',pr_ssp370_annualmean_1_min)
  pr_annualmean_df[8,26:101]<-rowMeans(pr_ssp370_annualmean_1)
  pr_annualmean_df[9,26:101]<-apply(pr_ssp370_annualmean_1, MARGIN=1, min)
  pr_annualmean_df[10,26:101]<-apply(pr_ssp370_annualmean_1, MARGIN=1, max)
  pr_ssp585_annualmean_1<-do.call('cbind',pr_ssp585_annualmean_1_min)
  pr_annualmean_df[11,26:101]<-rowMeans(pr_ssp585_annualmean_1)
  pr_annualmean_df[12,26:101]<-apply(pr_ssp585_annualmean_1, MARGIN=1, min)
  pr_annualmean_df[13,26:101]<-apply(pr_ssp585_annualmean_1, MARGIN=1, max)
  pr_regions_ssp_ts_list[[j]]<-pr_annualmean_df
  print(j)
}
saveRDS(pr_regions_ssp_rollpr_list,file='pr_24models_2324_ipcc_rollingts_min_rollingwindow_probabilityratios_list.RData')
saveRDS(pr_regions_ssp_ts_list,file='pr_24models_ipcc_rollingts_min_meantimeseries_list.RData')

# 2.2 Precipitation - plot probability ratios per region (FIGURE 4) #####
pr_regions_ssp_rollpr_list<-readRDS('pr_24models_2324_ipcc_rollingts_min_rollingwindow_probabilityratios_list.RData')
ipcc_names<-c('Amazon','Central_America_Mexico','East_Africa','East_Asia','North-East_Brazil',
              'Southern_Africa','South_Asia','Southeast_Asia','Southeastern_South_America',
              'West_Africa','West_Coast_South_America')
for(j in 1:length(ipcc_tropics_updated)){
  pr_regions_ssp_rollpr_j<-pr_regions_ssp_rollpr_list[[j]]
  colors <- c(rgb(0,52/255,102/255), "white", "darkorange2")
  break_points <- c(0, 1, 10)
  values_scaled <- (break_points - min(break_points)) / (max(break_points) - min(break_points))
  mat <- as.matrix(pr_regions_ssp_rollpr_j[c(1,4,7,10), ])
  df_long <- as.data.frame(as.table(mat))
  names(df_long) <- c("y", "x", "fill")
  df_long$x <- as.numeric(df_long$x)
  df_long$y <- as.numeric(df_long$y)
  heatmapgg<-ggplot(df_long, aes(x, rev(y), fill= fill)) + 
    geom_tile() +
    scale_fill_gradientn(
      colours = colors,
      values = values_scaled,
      limits = c(0,10),oob = scales::squish) +
    theme_void() + theme(legend.position="none")
  ggsave(heatmapgg,filename=paste0("prec_24models_rollingts_min_010_2324_pr_",ipcc_names[j],".png"))
}

#get values for text description
ipcc_names<-c('Amazon','Central_America_Mexico','East_Africa','East_Asia','North-East_Brazil',
              'Southern_Africa','South_Asia','Southeast_Asia','Southeastern_South_America',
              'West_Africa','West_Coast_South_America')
pr_regions_ssp_rollpr<-pr_regions_ssp_rollpr_list[[3]]
plot(as.numeric(pr_regions_ssp_rollpr[1,]),type='l')
plot(as.numeric(pr_regions_ssp_rollpr[4,]),type='l')
plot(as.numeric(pr_regions_ssp_rollpr[7,]),type='l')
plot(as.numeric(pr_regions_ssp_rollpr[10,]),type='l')
max(pr_regions_ssp_rollpr[1,])
max(pr_regions_ssp_rollpr[4,])
max(pr_regions_ssp_rollpr[7,])
max(pr_regions_ssp_rollpr[10,])

pr_regions_ssp_ts_list<-readRDS('pr_24models_ipcc_rollingts_min_meantimeseries_list.RData')
dates<-seq(as.Date('2000-01-01'),as.Date('2100-12-01'),by='year')
ipcc_names<-c('Amazon','Central_America_Mexico','East_Africa','East_Asia','North-East_Brazil',
              'Southern_Africa','South_Asia','Southeast_Asia','Southeastern_South_America',
              'West_Africa','West_Coast_South_America')
for(j in 1:length(ipcc_tropics_updated)){
  pr_regions_ssp_ts_j<-pr_regions_ssp_ts_list[[j]]
  ippc_region_plot <- ggplot()+
    geom_line(aes(x = dates,y = as.numeric(pr_regions_ssp_ts_j[1,]))) +
    geom_ribbon(aes(x = dates,ymin = as.numeric(pr_regions_ssp_ts_j[12,]), 
                    ymax = as.numeric(pr_regions_ssp_ts_j[13,])),fill=rgb(153/255,0,2/255),
                alpha=0.2)+
    geom_ribbon(aes(x = dates,ymin = as.numeric(pr_regions_ssp_ts_j[9,]), 
                    ymax = as.numeric(pr_regions_ssp_ts_j[10,])),fill=rgb(196/255,121/255,0),
                alpha=0.2)+
    geom_ribbon(aes(x = dates,ymin = as.numeric(pr_regions_ssp_ts_j[6,]), 
                    ymax = as.numeric(pr_regions_ssp_ts_j[7,])),fill=rgb(112/255,160/255,205/255),
                alpha=0.2)+
    geom_ribbon(aes(x = dates,ymin = as.numeric(pr_regions_ssp_ts_j[3,]), 
                    ymax = as.numeric(pr_regions_ssp_ts_j[4,])),fill=rgb(0,52/255,102/255),
                alpha=0.2)+
    geom_line(aes(x = dates,y = as.numeric(pr_regions_ssp_ts_j[1,]))) +
    geom_line(aes(x = dates,y = as.numeric(pr_regions_ssp_ts_j[11,])),col=rgb(153/255,0,2/255)) +
    geom_line(aes(x = dates,y = as.numeric(pr_regions_ssp_ts_j[8,])),col=rgb(196/255,121/255,0)) +
    geom_line(aes(x = dates,y = as.numeric(pr_regions_ssp_ts_j[5,])),col=rgb(112/255,160/255,205/255)) +
    geom_line(aes(x = dates,y = as.numeric(pr_regions_ssp_ts_j[2,])),col=rgb(0,52/255,102/255)) +
    geom_hline(yintercept=min(as.numeric(pr_regions_ssp_ts_j[1,24:25])),linetype='dashed') +
    ggtitle(ipcc_names[j]) +
    xlab('Time')+
    ylab('Annual precipitation') +
    ylim(min(pr_regions_ssp_ts_j,na.rm = T)-50,max(pr_regions_ssp_ts_j,na.rm = T)+50) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste0('prec_24models_rollingts_min_',ipcc_names[j],'.svg'),ippc_region_plot,
         width=5,height = 5,units='in')
}

# 2.3 Temperature - rolling window max TS and rolling window PR for all IPCC regions ####
tas_regions_ssp_rollpr_list<-list()
tas_regions_ssp_rollts_list<-list()
for(j in 1:length(ipcc_tropics_updated)){
  model_tas_ssp126_list<-list()
  model_tas_ssp245_list<-list()
  model_tas_ssp370_list<-list()
  model_tas_ssp585_list<-list()
  for(i in 1:length(model_nrs_370)){
    model_nrs_i<-model_nrs_370[i]
    ssp126_rasti<-rast(paste0('./RDS_CMIP6/ISIMIP_bias_adj/',models[model_nrs_i],'.ssp',ssp[1],'.tas.mask50_bias_adj.nc'))
    ssp126_rasti<-crop(ssp126_rasti,ipcc_tropics_updated[j,])
    model_tas_ssp126_list[i]<-mask(ssp126_rasti,ipcc_tropics_updated[j,])
    ssp245_rasti<-rast(paste0('./RDS_CMIP6/ISIMIP_bias_adj/',models[model_nrs_i],'.ssp',ssp[2],'.tas.mask50_bias_adj.nc'))
    ssp245_rasti<-crop(ssp245_rasti,ipcc_tropics_updated[j,])
    model_tas_ssp245_list[i]<-mask(ssp245_rasti,ipcc_tropics_updated[j,])
    ssp370_rasti<-rast(paste0('./RDS_CMIP6/ISIMIP_bias_adj/',models[model_nrs_i],'.ssp',ssp[3],'.tas.mask50_bias_adj.nc'))
    ssp370_rasti<-crop(ssp370_rasti,ipcc_tropics_updated[j,])
    model_tas_ssp370_list[i]<-mask(ssp370_rasti,ipcc_tropics_updated[j,])
    ssp585_rasti<-rast(paste0('./RDS_CMIP6/ISIMIP_bias_adj/',models[model_nrs_i],'.ssp',ssp[4],'.tas.mask50_bias_adj.nc'))
    ssp585_rasti<-crop(ssp585_rasti,ipcc_tropics_updated[j,])
    model_tas_ssp585_list[i]<-mask(ssp585_rasti,ipcc_tropics_updated[j,])
  }
  
  tas_ipcc<-mask(tas_brick,ipcc_tropics_updated[j,])
  tas_ipcc<-crop(tas_ipcc,ipcc_tropics_updated[j,])
  
  tas_anualmean<-roll(tas_ipcc,n=12,fun = mean,type='to')
  tas_annualmean_1<-global(tas_anualmean, fun=mean, na.rm=TRUE)
  years<-rep(1:25,each=12)
  tas_annualmean_max<-tapply(tas_annualmean_1$mean,years,max,na.rm=T)
  years_ssp<-rep(1:76,each=12)
  tas_ssp126_annualmean_list<-lapply(model_tas_ssp126_list,function(x) roll(x,n=12,fun = mean,type='to'))
  tas_ssp126_annualmean_1_list<-lapply(tas_ssp126_annualmean_list,function(x) global(x, fun=mean, na.rm=TRUE))
  tas_ssp126_annualmean_1_max<-lapply(tas_ssp126_annualmean_1_list,function(x) tapply(x$mean,years_ssp,max,na.rm=T))
  tas_ssp245_annualmean_list<-lapply(model_tas_ssp245_list,function(x) roll(x,n=12,fun = mean,type='to'))
  tas_ssp245_annualmean_1_list<-lapply(tas_ssp245_annualmean_list,function(x) global(x, fun=mean, na.rm=TRUE))
  tas_ssp245_annualmean_1_max<-lapply(tas_ssp245_annualmean_1_list,function(x) tapply(x$mean,years_ssp,max,na.rm=T))
  tas_ssp370_annualmean_list<-lapply(model_tas_ssp370_list,function(x) roll(x,n=12,fun = mean,type='to'))
  tas_ssp370_annualmean_1_list<-lapply(tas_ssp370_annualmean_list,function(x) global(x, fun=mean, na.rm=TRUE))
  tas_ssp370_annualmean_1_max<-lapply(tas_ssp370_annualmean_1_list,function(x) tapply(x$mean,years_ssp,max,na.rm=T))
  tas_ssp585_annualmean_list<-lapply(model_tas_ssp585_list,function(x) roll(x,n=12,fun = mean,type='to'))
  tas_ssp585_annualmean_1_list<-lapply(tas_ssp585_annualmean_list,function(x) global(x, fun=mean, na.rm=TRUE))
  tas_ssp585_annualmean_1_max<-lapply(tas_ssp585_annualmean_1_list,function(x) tapply(x$mean,years_ssp,max,na.rm=T))
  
  timing<-c(24:25)
  events_control <- sum(tas_annualmean_max>=max(tas_annualmean_max[timing]),na.rm=T)
  years_control <- 25
  
  #do rolling window
  ssp126_results_list<-list()
  for (i in seq_along(tas_ssp126_annualmean_1_max)) {
    mean_vec <- tas_ssp126_annualmean_1_max[[i]]
    n <- length(mean_vec)
    events_future <- numeric(n - 25 + 1)
    years_future  <- rep(25, n - 25 + 1)
    for (k in 1:(n - 25 + 1)) {
      window_k <- mean_vec[k:(k + 25 - 1)]
      events_future[k] <- sum(window_k >= max(tas_annualmean_max[timing]),na.rm=T)
    }
    ssp126_results_list[[i]] <- mapply(compute_PR,events_control,years_control,events_future,years_future)
  }
  ssp245_results_list<-list()
  for (i in seq_along(tas_ssp245_annualmean_1_max)) {
    mean_vec <- tas_ssp245_annualmean_1_max[[i]]
    n <- length(mean_vec)
    events_future <- numeric(n - 25 + 1)
    years_future  <- rep(25, n - 25 + 1)
    for (k in 1:(n - 25 + 1)) {
      window_k <- mean_vec[k:(k + 25 - 1)]
      events_future[k] <- sum(window_k >= max(tas_annualmean_max[timing]),na.rm=T)
    }
    ssp245_results_list[[i]] <- mapply(compute_PR,events_control,years_control,events_future,years_future)
  }
  ssp370_results_list<-list()
  for (i in seq_along(tas_ssp370_annualmean_1_max)) {
    mean_vec <- tas_ssp370_annualmean_1_max[[i]]
    n <- length(mean_vec)
    events_future <- numeric(n - 25 + 1)
    years_future  <- rep(25, n - 25 + 1)
    for (k in 1:(n - 25 + 1)) {
      window_k <- mean_vec[k:(k + 25 - 1)]
      events_future[k] <- sum(window_k >= max(tas_annualmean_max[timing]),na.rm=T)
    }
    ssp370_results_list[[i]] <- mapply(compute_PR,events_control,years_control,events_future,years_future)
  }
  ssp585_results_list<-list()
  for (i in seq_along(tas_ssp585_annualmean_1_max)) {
    mean_vec <- tas_ssp585_annualmean_1_max[[i]]
    n <- length(mean_vec)
    events_future <- numeric(n - 25 + 1)
    years_future  <- rep(25, n - 25 + 1)
    for (k in 1:(n - 25 + 1)) {
      window_k <- mean_vec[k:(k + 25 - 1)]
      events_future[k] <- sum(window_k >= max(tas_annualmean_max[timing]),na.rm=T)
    }
    ssp585_results_list[[i]] <- mapply(compute_PR,events_control,years_control,events_future,years_future)
  }
  ssp126_results<-do.call('rbind',ssp126_results_list)
  ssp245_results<-do.call('rbind',ssp245_results_list)
  ssp370_results<-do.call('rbind',ssp370_results_list)
  ssp585_results<-do.call('rbind',ssp585_results_list)
  
  #get mean and min/max PR time series per scenario in df
  ssp_results_pr<-as.data.frame(matrix(nrow=12,ncol=52))
  ssp_results_pr[1,]<-colMeans(ssp126_results)
  ssp_results_pr[2,]<-apply(ssp126_results, MARGIN=2, min)
  ssp_results_pr[3,]<-apply(ssp126_results, MARGIN=2, max)
  ssp_results_pr[4,]<-colMeans(ssp245_results)
  ssp_results_pr[5,]<-apply(ssp245_results, MARGIN=2, min)
  ssp_results_pr[6,]<-apply(ssp245_results, MARGIN=2, max)
  ssp_results_pr[7,]<-colMeans(ssp370_results)
  ssp_results_pr[8,]<-apply(ssp370_results, MARGIN=2, min)
  ssp_results_pr[9,]<-apply(ssp370_results, MARGIN=2, max)
  ssp_results_pr[10,]<-colMeans(ssp585_results)
  ssp_results_pr[11,]<-apply(ssp585_results, MARGIN=2, min)
  ssp_results_pr[12,]<-apply(ssp585_results, MARGIN=2, max)
  tas_regions_ssp_rollpr_list[[j]]<-ssp_results_pr
  
  #get mean and min/max time series per scenario in df
  tas_annualmean_df<-as.data.frame(matrix(nrow=13,ncol=101))
  tas_annualmean_df[1,1:25]<-tas_annualmean_max
  tas_ssp126_annualmean_1<-do.call('cbind',tas_ssp126_annualmean_1_max)
  tas_annualmean_df[2,26:101]<-rowMeans(tas_ssp126_annualmean_1)
  tas_annualmean_df[3,26:101]<-apply(tas_ssp126_annualmean_1, MARGIN=1, min)
  tas_annualmean_df[4,26:101]<-apply(tas_ssp126_annualmean_1, MARGIN=1, max)
  tas_ssp245_annualmean_1<-do.call('cbind',tas_ssp245_annualmean_1_max)
  tas_annualmean_df[5,26:101]<-rowMeans(tas_ssp245_annualmean_1)
  tas_annualmean_df[6,26:101]<-apply(tas_ssp245_annualmean_1, MARGIN=1, min)
  tas_annualmean_df[7,26:101]<-apply(tas_ssp245_annualmean_1, MARGIN=1, max)
  tas_ssp370_annualmean_1<-do.call('cbind',tas_ssp370_annualmean_1_max)
  tas_annualmean_df[8,26:101]<-rowMeans(tas_ssp370_annualmean_1)
  tas_annualmean_df[9,26:101]<-apply(tas_ssp370_annualmean_1, MARGIN=1, min)
  tas_annualmean_df[10,26:101]<-apply(tas_ssp370_annualmean_1, MARGIN=1, max)
  tas_ssp585_annualmean_1<-do.call('cbind',tas_ssp585_annualmean_1_max)
  tas_annualmean_df[11,26:101]<-rowMeans(tas_ssp585_annualmean_1)
  tas_annualmean_df[12,26:101]<-apply(tas_ssp585_annualmean_1, MARGIN=1, min)
  tas_annualmean_df[13,26:101]<-apply(tas_ssp585_annualmean_1, MARGIN=1, max)
  tas_regions_ssp_rollts_list[[j]]<-tas_annualmean_df
  print(j)
}
saveRDS(tas_regions_ssp_rollpr_list,file='tas_24models_2324_ipcc_rollingts_max_rollingwindow_probabilityratios_list.RData')
saveRDS(tas_regions_ssp_rollts_list,file='tas_24models_ipcc_rollingts_max_meantimeseries_list.RData')

# 2.4 Temperature - plot probability ratios per region (FIGURE 4) ####
ttas_regions_ssp_rollpr_list<-readRDS('tas_24models_2324_ipcc_rollingts_max_rollingwindow_probabilityratios_list.RData')
dates<-seq(as.Date('2049-01-01'),as.Date('2100-12-01'),by='year')
ipcc_names<-c('Amazon','Central_America_Mexico','East_Africa','East_Asia','North-East_Brazil',
              'Southern_Africa','South_Asia','Southeast_Asia','Southeastern_South_America',
              'West_Africa','West_Coast_South_America')
for(j in 1:length(ipcc_tropics_updated)){
  tas_regions_ssp_rollpr_j<-tas_regions_ssp_rollpr_list[[j]]
  colors <- c("white", "darkorange2")
  break_points <- c(1, 25)
  values_scaled <- (break_points - min(break_points)) / (max(break_points) - min(break_points))
  mat <- as.matrix(tas_regions_ssp_rollpr_j[c(1,4,7,10), ])
  df_long <- as.data.frame(as.table(mat))
  names(df_long) <- c("y", "x", "fill")
  df_long$x <- as.numeric(df_long$x)
  df_long$y <- as.numeric(df_long$y)
  heatmapgg<-ggplot(df_long, aes(x, rev(y), fill= fill)) + 
    geom_tile() +
    scale_fill_gradientn(
      colours = colors,
      #values = values_scaled,
      limits = c(1,25),oob = scales::squish) +
    theme_void() + theme(legend.position="none")
  ggsave(heatmapgg,filename=paste0("tas_24models_rollingts_max_125_2324_pr_",ipcc_names[j],".png"))
}

ipcc_names<-c('Amazon','Central_America_Mexico','East_Africa','East_Asia','North-East_Brazil',
              'Southern_Africa','South_Asia','Southeast_Asia','Southeastern_South_America',
              'West_Africa','West_Coast_South_America')
tas_regions_ssp_rollpr<-tas_regions_ssp_rollpr_list[[11]]
plot(as.numeric(tas_regions_ssp_rollpr[1,]),type='l')
plot(as.numeric(tas_regions_ssp_rollpr[4,]),type='l')
plot(as.numeric(tas_regions_ssp_rollpr[7,]),type='l')
plot(as.numeric(tas_regions_ssp_rollpr[10,]),type='l')
which(as.numeric(tas_regions_ssp_rollpr[4,])>24.5)
max(tas_regions_ssp_rollpr[1,])
max(tas_regions_ssp_rollpr[4,])
max(tas_regions_ssp_rollpr[7,])
max(tas_regions_ssp_rollpr[10,])

plot(as.numeric(tas_regions_ssp_rollpr_list[[1]][1,]),type='l',ylim=c(0,25)) #SSP126
for(j in 1:length(ipcc_tropics_updated)){
  tas_regions_ssp_rollpr_j<-tas_regions_ssp_rollpr_list[[j]]
  lines(as.numeric(tas_regions_ssp_rollpr_list[[j]][1,]))
  print(min(as.numeric(tas_regions_ssp_rollpr_list[[j]][1,])))
}
plot(as.numeric(tas_regions_ssp_rollpr_list[[1]][4,]),type='l',ylim=c(0,25))#SSP245
for(j in 2:length(ipcc_tropics_updated)){
  tas_regions_ssp_rollpr_j<-tas_regions_ssp_rollpr_list[[j]]
  lines(as.numeric(tas_regions_ssp_rollpr_list[[j]][4,]))
  print(max(as.numeric(tas_regions_ssp_rollpr_list[[j]][4,])))
}
plot(as.numeric(tas_regions_ssp_rollpr_list[[1]][7,]),type='l',ylim=c(0,25))#SSP370
for(j in 2:length(ipcc_tropics_updated)){
  tas_regions_ssp_rollpr_j<-tas_regions_ssp_rollpr_list[[j]]
  lines(as.numeric(tas_regions_ssp_rollpr_list[[j]][7,]))
  print(max(as.numeric(tas_regions_ssp_rollpr_list[[j]][7,])))
}
plot(as.numeric(tas_regions_ssp_rollpr_list[[1]][10,]),type='l',ylim=c(0,25))#SSP585
for(j in 2:length(ipcc_tropics_updated)){
  tas_regions_ssp_rollpr_j<-tas_regions_ssp_rollpr_list[[j]]
  lines(as.numeric(tas_regions_ssp_rollpr_list[[j]][10,]))
  print(max(as.numeric(tas_regions_ssp_rollpr_list[[j]][10,])))
}

tas_regions_ssp_ts_list<-readRDS('tas_24models_ipcc_rollingts_max_meantimeseries_list.RData')
dates<-seq(as.Date('2000-01-01'),as.Date('2100-12-01'),by='year')
ipcc_names<-c('Amazon','Central_America_Mexico','East_Africa','East_Asia','North-East_Brazil',
              'Southern_Africa','South_Asia','Southeast_Asia','Southeastern_South_America',
              'West_Africa','West_Coast_South_America')
for(j in 1:length(ipcc_tropics_updated)){
  tas_regions_ssp_ts_j<-tas_regions_ssp_ts_list[[j]]
  ippc_region_plot <- ggplot()+
    geom_ribbon(aes(x = dates,ymin = as.numeric(tas_regions_ssp_ts_j[12,]), 
                    ymax = as.numeric(tas_regions_ssp_ts_j[13,])),fill=rgb(153/255,0,2/255),
                alpha=0.2)+
    geom_ribbon(aes(x = dates,ymin = as.numeric(tas_regions_ssp_ts_j[9,]), 
                    ymax = as.numeric(tas_regions_ssp_ts_j[10,])),fill=rgb(196/255,121/255,0),
                alpha=0.2)+
    geom_ribbon(aes(x = dates,ymin = as.numeric(tas_regions_ssp_ts_j[6,]), 
                    ymax = as.numeric(tas_regions_ssp_ts_j[7,])),fill=rgb(112/255,160/255,205/255),
                alpha=0.2)+
    geom_ribbon(aes(x = dates,ymin = as.numeric(tas_regions_ssp_ts_j[3,]), 
                    ymax = as.numeric(tas_regions_ssp_ts_j[4,])),fill=rgb(0,52/255,102/255),
                alpha=0.2)+
    geom_line(aes(x = dates,y = as.numeric(tas_regions_ssp_ts_j[1,]))) +
    geom_line(aes(x = dates,y = as.numeric(tas_regions_ssp_ts_j[11,])),col=rgb(153/255,0,2/255)) +
    geom_line(aes(x = dates,y = as.numeric(tas_regions_ssp_ts_j[8,])),col=rgb(196/255,121/255,0)) +
    geom_line(aes(x = dates,y = as.numeric(tas_regions_ssp_ts_j[5,])),col=rgb(112/255,160/255,205/255)) +
    geom_line(aes(x = dates,y = as.numeric(tas_regions_ssp_ts_j[2,])),col=rgb(0,52/255,102/255)) +
    geom_hline(yintercept=max(as.numeric(tas_regions_ssp_ts_j[1,24:25])),linetype='dashed') +
    ggtitle(ipcc_names[j]) +
    xlab('Time')+
    ylab('Mean annual temperature') +
    ylim(min(tas_regions_ssp_ts_j,na.rm = T)-1,max(tas_regions_ssp_ts_j,na.rm = T)+1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(paste0('tas_24models_rollingts_max_',ipcc_names[j],'.svg'),ippc_region_plot,
         width=5,height = 5,units='in')
}


