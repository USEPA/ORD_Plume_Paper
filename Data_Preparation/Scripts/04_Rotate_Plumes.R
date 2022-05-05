library(tidyverse)
library(sf)
#library(spatialEco)
#library(spdep)
library(vroom)
library(geosphere)
library(here)

# Import Monitoring Well Data
timeDist <- vroom(here("Data_Preparation/Analysis_Data/Benzene_Time_Dist.tsv"))

# Import mean linear direction calculations
mld_files <- list.files(here("Data_Preparation/Analysis_Data"),
                        pattern = "PlumeDirections", full.names = TRUE)
mld <- vroom(mld_files)
start <- Sys.time()
for(n in 1:nrow(mld)){
  
  bearing <- mld[n,]
  pb <- bearing$PlumeBearing
  
  # move coordinates to 0,0 and create rotated convex hull of plume
  allWells <- timeDist%>%
    filter(GLOBAL_ID == bearing$GLOBAL_ID &
             LOGDATE > bearing$OriginTestDate - 28 & LOGDATE < bearing$OriginTestDate + 28)
  allWells$plumeBearing <- pb
  origin <- allWells%>%
    filter(POINT_TYPE == "ORIGIN")
  origin <- origin[1,]
  
  mws <- allWells%>%
    filter(POINT_TYPE == "MW")%>%
    mutate(newBearing = if_else(plumeBearing < 0,Bearing + abs(plumeBearing),99,
                                if_else(plumeBearing >= 0, Bearing - plumeBearing, 999)),
           newBearing = if_else(newBearing > 180, -180 + (newBearing - 180),
                               if_else(newBearing < -180, 180 - (newBearing - 180),newBearing)))
  
  rotate <- destPoint(c(0,0),mws$newBearing,mws$Distance_m)%>%
    as.data.frame()%>%
    st_as_sf(coords = c('lon', 'lat'), crs = 3310)%>%
    cbind(mws)
  
  originSf <- origin%>%
    mutate(shiftLat = 0,
           shiftLon = 0)%>%
    st_as_sf(coords = c("shiftLon","shiftLat"), crs = 3310)
    
  merge <- bind_rows(rotate,originSf)%>%
    mutate(dateID = paste0(GLOBAL_ID,"_",bearing$OriginTestDate))
  
  
  mwRotate <- rbind(mwRotate,merge)
  
  now <- Sys.time()
  
  print(paste0(round(100*(n/nrow(mld))),"% Complete ... ", round(as.numeric(difftime(now,start,units = "mins")))," minutes elapsed."))
}

