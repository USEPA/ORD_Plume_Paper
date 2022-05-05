# In this script, we will loop through every LUST site
# and every date that the point closest to the leak
# was sampled and determine the mean linear direction of the plume
# based on that sampling campaign. Then, a convex hull will be drawn
# to represent the plume extent. The plume will then be rotated so 
# that the mean direction is facing north. Then the origin point will be shifted
# to 0,0 so that we can overlay every plume.

## THIS SCRIPT WILL TAKE SEVERAL DAYS TO RUN

library(tidyverse)
library(sf)
library(spatialEco)
library(spdep)
library(vroom)
library(here)

source(here("Functions", "mean_linear_direction.R"))

# Import Benzene data
timeDist <- vroom(here("Data_Preparation/Analysis_Data/Benzene_Time_Dist.tsv"))

# List sites already run, then skip those
# run <- st_layers(here("Data_Preparation/Analysis_Data/rotatedPlumes.gpkg"))
# 
# run_GID <- data.frame(File = run$name)%>%
#   separate(File, into = "GLOBAL_ID", extra = "drop", sep = "_")
# 
# timeDist <- timeDist%>%
#   filter(!GLOBAL_ID %in% run_GID$GLOBAL_ID)

# Where is 0,0 in California Albers?
# test <- data.frame(Lat = 0, Lon = 0)%>%
#   st_as_sf(coords = c("Lon","Lat"),crs = 3310)%>%
#   st_transform(4326)


# Write a loop for every Global_ID

# Random sample for testing
# samp <- sample(timeDist$GLOBAL_ID,100, replace = FALSE)
# timeDist <- timeDist%>%
#   filter(GLOBAL_ID %in% samp)


i <- 0
mwRotate <- st_sf(st_sfc(),crs=3310)# Create empty sf object to write points to.
for (site in unique(timeDist$GLOBAL_ID)) {
  tryCatch({
    # Filter to the site
    siteDf <- timeDist%>%
      filter(GLOBAL_ID == site)
    
    # Filter to the origin points
    origin <- siteDf%>%
      filter(POINT_TYPE == "ORIGIN")
    
    # List the unique dates the origin was sampled
    originDates <- unique(origin$LOGDATE)
    
    # Loop through each date and filter monitoring well measurements based on
    # if measurements were taken within 28 days of the origin sample
    # Only proceed if there are at least 3 samples >= 5 ug/L in the sampling campaign
    for (date in unique(as.list(originDates))) {
      mws <- siteDf%>%
        filter(POINT_TYPE == "MW" & LOGDATE > date - 28 & LOGDATE < date + 28 & PARVAL >= 5)
      
      # If the conditions are satisfied, we'll map the convex hull
      # of both the area > 5 ug/L and the convex hull of the entire area
      # measured. Then we'll determine the mean linear direction and rotate
      # both polygons to north, move their origins to 0,0 and save them.
      if(nrow(mws)>2){
        
        # Run the weighted mean direction
        bearing <- mld(mws,coords = c("LONGITUDE","LATITUDE"),start_x = origin$LONGITUDE[1],
            start_y = origin$LATITUDE[1], date = str_replace_all(as.character(date),"-","_"),crs = 3310,
            weight = mws$PARVAL,filename = here("Data_Preparation/Analysis_Data/meanDirection.gpkg"))
        
        # move coordinates to 0,0 and create rotated convex hull of plume
        allWells <- siteDf%>%
          filter(LOGDATE > date - 28 & LOGDATE < date + 28)
        
        allWellsSf <- allWells%>%
          st_as_sf(coords = c("LONGITUDE","LATITUDE"), crs = 4326)%>%
          st_transform(3310)
        
        mwCoords <- as.data.frame(st_coordinates(allWellsSf))
        
        originSf <- st_as_sf(origin[1,], coords = c("LONGITUDE","LATITUDE"),
                             crs = 4326)%>%
          st_transform(3310)
        originCoords <- as.data.frame(st_coordinates(originSf))
        
        # Calculate shifted coordinates
        mwsShift <- allWells%>%
          mutate(shiftLat = mwCoords$Y - originCoords$Y,
                 shiftLon = mwCoords$X - originCoords$X)
        
        # convex hull for measurements > MCL
        # plumePoly <- mwsShift%>%
        #   st_as_sf(coords = c("shiftLon","shiftLat"), crs = 3310)%>%
        #   st_transform(4326)%>%
        #   st_union()%>%
        #   st_convex_hull()%>%
        #   st_cast("LINESTRING")%>%
        #   st_sf()
        # plumePoly$GLOBAL_ID = site
        # plumePoly$LOGDATE = date
        
        # st_write(plumePoly, here("projects/PetroleumHydrocarbons/data/meanDirection/shiftedPlumes.gpkg"),
        #          layer = paste0(site,"_",str_replace_all(date,"-","_")),append = FALSE)
        # 
        # Rotate all of the points so that mean direction is north from the origin
        mwsOnly <- mwsShift%>%
          filter(POINT_TYPE=="MW")
        
        #mwRotate <- st_sf(st_sfc(),crs=3310)
        
        for(n in 1:nrow(mwsOnly)){
          # Use trigonometry to rotate points around the axis of the origin point
          
          ## Find the length of the line between the monitoring well and it's rotated
          ## location by using the law of SINs, given the mean linear direction of the plume.
          ## Angles are A,B,C, side lengths are a,b,c. A = angle from north at the origin.
          
          pt <- mwsOnly[n,]
          
          A <- bearing
          
          B <- (180-A)/2
          
          b <- pt$Distance_m
          
          # sin(B)/b = sin(A)/a
          
          sinB <- sin(B*(pi/180))
          sinA <- sin(A*(pi/180))
          
          a <- b/sinB * sinA
          
          # Draw monitoring well in space (original bearing from 0,0)
          ptSf <- pt%>%
            st_as_sf(coords = c("shiftLon","shiftLat"), crs = 3310)
          
          # Draw origin point at (0,0)
          originSfShift <- mwsShift%>%
            filter(POINT_TYPE == "ORIGIN")%>%
            st_as_sf(coords = c("shiftLon","shiftLat"), crs = 3310)
          
          # Buffer the origin point based on given side length (b)
          originBuf <- originSfShift %>%
            st_buffer(pt$Distance_m)%>%
            st_cast("LINESTRING")
          
          # Buffer the monitoring well point based on the calculated length (a)
          ptBuf <- pt%>%
            st_as_sf(coords = c("shiftLon","shiftLat"), crs = 3310)%>%
            st_buffer(abs(a))%>%
            st_cast("LINESTRING")
          
          # Triangulate the two candidate points for the rotated point location
          intrsctPts <- st_intersection(originBuf,ptBuf)%>%
            st_cast("POINT")%>%
            select(COUNTY.1:Time.1)
          colnames(intrsctPts) <- colnames(originBuf)
          intrsctCoords <- as.data.frame(st_coordinates(intrsctPts))
          
          # Determine which point is the correct one to use
          # by calculating the bearing from the shifted origin
          # and comparing the result to the original bearing - dominant site bearing
          
          if(pt$Bearing - bearing >= -180 &pt$Bearing - bearing <= 180){
            newBearing <-  pt$Bearing - bearing
          }
          
          if(pt$Bearing - bearing < -180){
            newbearing <- 180 - (abs(pt$Bearing - bearing))
          }
          
          if(pt$Bearing - bearing > 180){
            newbearing <- -180 + (abs(pt$Bearing - bearing))
          }
          
          
          ptsA <- matrix(c(0,
                           intrsctCoords$X[1],
                           0,
                           intrsctCoords$Y[1]),,2)
          
          lineA = st_linestring(ptsA)%>%
            st_sfc(crs = 3310)%>%
            st_sf()%>%
            st_transform(4326)%>%
            mutate(Line = "A")
          
          ptsB <- matrix(c(0,
                           intrsctCoords$X[2],
                           0,
                           intrsctCoords$Y[2]),,2)
          
          lineB = st_linestring(ptsB)%>%
            st_sfc(crs = 3310)%>%
            st_sf()%>%
            st_transform(4326)%>%
            mutate(Line = "B")
          
          
          # Calculate new bearings
          candidateBearing <- bind_rows(lineA,lineB)
          candidateBearing$Bearing <- stplanr::line_bearing(candidateBearing)
          candidateBearing$NewBearing <- newBearing
          candidateBearing$absDif <- abs(candidateBearing$Bearing - candidateBearing$NewBearing)
          
          
          rotatedBearing <- candidateBearing%>%
            filter(absDif == min(absDif))%>%
            st_transform(3310)
          
          if(rotatedBearing$Line[1] == "A"){
            rotatedMw <- intrsctPts[1,]
          }
          if(rotatedBearing$Line[1] == "B"){
            rotatedMw <- intrsctPts[2,]
          }
          
          rotatedPts <- rotatedMw%>%
            rbind(originSfShift)
          
          mwRotate <- rbind(mwRotate,rotatedMw)
        }
        # Save rotated points
        # st_write(mwRotate,
        #          here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedPoints.gpkg"),
        #          layer = paste0(site,"_",str_replace_all(date,"-","_")),append = FALSE)
        
        # Calculate rotated convex hull > mcl (origin does not move)
        # gtMCLpts <- rbind(mwRotate,originSfShift)%>%
        #   filter(PARVAL >= 5)
        # 
        # plume <- gtMCLpts%>%
        #   st_union()%>%
        #   st_convex_hull()
        # 
        # st_write(plume,
        #          here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedPlumes.gpkg"),
        #          layer = paste0(site,"_",str_replace_all(date,"-","_")),append = FALSE)
        
        # Calculate rotated convex hull of study area
        # studyArea <- rbind(mwRotate,originSfShift)%>%
        #   st_union()%>%
        #   st_convex_hull()
        # 
        # st_write(studyArea,
        #          here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedStudyArea.gpkg"),
        #          layer = paste0(site,"_",str_replace_all(date,"-","_")),append = FALSE)
      }
    }
    print(paste0("Completed ",round(100*(i/length(unique(timeDist$GLOBAL_ID))),1),"% at ",Sys.time()))
    i <- i+1
  }, error=function(e){})
}

st_write(mwRotate,
         here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedPoints.gpkg"),
         layer = "Rotated_Points_All",append = FALSE)

#STARTED: 4:17 PM Friday

# Merge all plume layers into one
plumeLayers <- st_layers(here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedPlumes.gpkg"))
plDf <- data.frame(Layer = plumeLayers$name)%>%
  separate(Layer,into = c("GLOBAL_ID","Year","Month","Day"), sep = "_",remove = FALSE)%>%
  drop_na()

plumeSf <- st_read(here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedPlumes.gpkg"),
                   layer = plDf$Layer[1])%>%
  mutate(GLOBAL_ID = plDf$GLOBAL_ID[1],
         Date = lubridate::mdy(paste0(plDf$Month[1],"/",plDf$Day[1],"/",plDf$Year[1])))
for (n in 2:nrow(plDf)) {
  newPlume <- st_read(here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedPlumes.gpkg"),
                      layer = plDf$Layer[n])%>%
    mutate(GLOBAL_ID = plDf$GLOBAL_ID[n],
           Date = lubridate::mdy(paste0(plDf$Month[n],"/",plDf$Day[n],"/",plDf$Year[n])))
  
  plumeSf <- rbind(plumeSf,newPlume)
  
  print(paste0(round(100*(n/nrow(plDf)),2),"% complete at: ", Sys.time()))
}

st_write(plumeSf,here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotated_Polygons.gpkg"),
         layer = "Plumes")

# Merge study area layers into one
saLayers <- st_layers(here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedStudyArea.gpkg"))
saDf <- data.frame(Layer = saLayers$name)%>%
  separate(Layer,into = c("GLOBAL_ID","Year","Month","Day"), sep = "_",remove = FALSE)%>%
  drop_na()

saSf <- st_read(here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedStudyArea.gpkg"),
                layer = saDf$Layer[1])%>%
  mutate(GLOBAL_ID = saDf$GLOBAL_ID[1],
         Date = lubridate::mdy(paste0(saDf$Month[1],"/",saDf$Day[1],"/",saDf$Year[1])))
for (n in 2:nrow(saDf)) {
  newSA <- st_read(here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedStudyArea.gpkg"),
                   layer = saDf$Layer[n])%>%
    mutate(GLOBAL_ID = saDf$GLOBAL_ID[n],
           Date = lubridate::mdy(paste0(saDf$Month[n],"/",saDf$Day[n],"/",saDf$Year[n])))
  
  saSf <- rbind(saSf,newSA)
  
  print(paste0(round(100*(n/nrow(saDf)),2),"% complete at: ", Sys.time()))
}


st_write(saSf,here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotated_Polygons.gpkg"),
         layer = "Study_Area")



# Merge rotated points into a single layer and save
ptLayers <- st_layers(here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedPoints.gpkg"))
ptLayersdf <- data.frame("Layer" = ptLayers$name)

rotPts <- st_read(here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedPoints.gpkg"),
                  layer = ptLayersdf$Layer[1])
start <- Sys.time()
for(n in 2:nrow(ptLayersdf)){
  pts <- st_read(here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotatedPoints.gpkg"),
                 layer = ptLayersdf$Layer[n])
  
  rotPts <- rbind(rotPts,pts)
  
  print(paste0(round(100*(n/nrow(ptLayersdf)),2),"% Complete --- ",
               round(as.numeric(difftime(Sys.time(), start, units='mins')),2),
               " minutes elapsed"))
}

st_write(rotPts,here("CA_Benzene_Plumes/Data_Preparation/Analysis_Data/rotated_Polygons.gpkg"),
         layer = "Roints")




