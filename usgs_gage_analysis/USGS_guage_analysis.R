remove(list = ls()) # clear all workspace variables
cat("\014")         # clear command line

library(tidyverse)
library(waterData)
library(dataRetrieval)
library(lubridate)
library(ncdf4)
library(doParallel)
library(foreach)

root = "C:/Users/awlostowski/Documents/national-water-center/t-route"

# load RouteLink data
rl_dat <- nc_open(paste0(root,"/test/input/geo/Channels/","RouteLink_NHDPLUS.nwm.v2.0.4.nc"))


link = ncvar_get(rl_dat, 'link')
gages = ncvar_get(rl_dat, 'gages')

idx = which(gages != "               ")
gage_link <- data.frame(link = link[idx], gages = trimws(gages[idx],"l"))

# set core usaage
cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

print("starting parallel loop")
record_peaks <- foreach(i = 1:length(gage_link$gages),
                  .combine=rbind,
                  .packages = c('dataRetrieval','dplyr')) %dopar% {
                    
  # call masking function
  dat = readNWISpeak(siteNumbers = gage_link$gages[i])
  
  if (nrow(dat) > 0){
    
    dat <- dat %>%
      filter(peak_va == max(peak_va, na.rm = T)) %>%
      mutate(link = gage_link$link[i]) %>%
      select(link, site_no, peak_dt, peak_va)
   
    dat 
  }
  
}
stopCluster(cl)

# rename some columns more intuitively
record_peaks_rename <- record_peaks %>%
  mutate(peak_cms = round(peak_va * 0.0283168,2)) %>%
  select(-peak_va)
  
# write out data as csv
write.csv(record_peaks_rename,paste0(root,"/test/","usgs_record_peaks.csv"))



