library(tidyverse)
library(sf)
library(vroom)
library(geosphere)
library(doParallel)
library(here)

# Import Monitoring Well Data
timeDist <- vroom(here("Data_Preparation/Analysis_Data/Benzene_Time_Dist.tsv"))

# Import mean linear direction calculations
mld <- vroom(here("Data_Preparation/Analysis_Data/PlumeDirections_bySite.csv"))
start <- Sys.time()
mwRotate <- st_sf(st_sfc(),crs=3310)# Create empty sf object to write points to.

# Create Cluster
cores <- detectCores()-2
cl <- makeCluster(cores)
registerDoParallel(cl)

sf <- foreach(n = 1:nrow(mld),
        .combine=rbind,
        .packages = c("tidyverse","sf","vroom","geosphere")) %dopar%{
  bearing <- mld[n,]
pb <- bearing$PlumeBearing

# move coordinates to 0,0 and create rotated convex hull of plume
allWells <- timeDist%>%
  filter(GLOBAL_ID == bearing$GLOBAL_ID)
allWells$plumeBearing <- pb
origin <- allWells%>%
  filter(POINT_TYPE == "ORIGIN")
origin <- origin[1,]

mws <- allWells%>%
  filter(POINT_TYPE == "MW")%>%
  mutate(newBearing = ifelse(plumeBearing < 0,Bearing + abs(plumeBearing),
                              ifelse(plumeBearing >= 0, Bearing - plumeBearing, NA)),
         newBearing = ifelse(newBearing > 180, -180 + (newBearing - 180),
                              ifelse(newBearing < -180, 180 - (abs(newBearing) - 180),newBearing)))

rotate <- destPoint(c(-120,38.01637),mws$newBearing,mws$Distance_m)%>%
  as.data.frame()%>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)%>%
  st_transform(3310)%>%
  cbind(mws)

originSf <- origin%>%
  mutate(shiftLat = 0,
         shiftLon = 0)%>%
  st_as_sf(coords = c("shiftLon","shiftLat"), crs = 3310)

bind_rows(rotate,originSf)%>%
  mutate(dateID = paste0(GLOBAL_ID,"_",bearing$OriginTestDate))
}

st_write(sf, here("Data_Preparation/Analysis_Data/Rotated_by_Site.gpkg"), layer = "All_Points")
