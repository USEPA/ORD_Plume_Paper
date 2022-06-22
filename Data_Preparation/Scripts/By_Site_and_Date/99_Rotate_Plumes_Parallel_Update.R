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
library(stplanr)
library(doParallel)
library(here)

# Import Benzene data
timeDist <- vroom(here("Data_Preparation/Analysis_Data/Benzene_Time_Dist.tsv"))

source(here("Functions", "mean_linear_direction.R"))
print(paste0("Loading Data: ",Sys.time()))


# Create df to write new coordinates
rotatedTimeDist <- timeDist%>%
  mutate(LATITUDE = NA,
         LONGITUDE = NA)

# Write a loop for every Global_ID

# Random sample for testing
samp <- sample(timeDist$GLOBAL_ID,5, replace = FALSE)
timeDist <- timeDist%>%
  filter(GLOBAL_ID %in% samp)

# Test parallel
# Create df to write new coordinates
rotatedTimeDist <- timeDist%>%
  mutate(LATITUDE = NA,
         LONGITUDE = NA)

cores <- detectCores()-1
cl <- makeCluster(cores)
registerDoParallel(cl)


# testing
start <- Sys.time()
x <-foreach(c=unique(mtcars$cyl), .combine='rbind', .packages = c("tidyverse")) %:%
  foreach(g=unique(mtcars$gear), .combine='rbind') %dopar% {
    filt <- mtcars%>%
      filter(cyl == c, gear == g)
    data.frame("Cyl" = c, "Gear" = g, "Maxhp" = max(filt$hp),
               "Minhp" = min(filt$hp))
  }
end <- Sys.time()
x

insideTime <- end-start

# Write a loop for every Global_ID




parStart <- Sys.time()

print(paste0("Starting loop: ",parStart))

rotate <- foreach(site = unique(timeDist$GLOBAL_ID),
                  .combine=rbind,
                  .packages = c("tidyverse","sf","spatialEco",
                                "spdep","vroom","here","stplanr")) %:% 
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
    foreach (date = unique(as.list(originDates))) %:%
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
      
      # move coordinates to 0,0 and rotate based on bearing
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
      
      # Rotate all of the points so that mean direction is north from the origin
      mwsOnly <- mwsShift%>%
        filter(POINT_TYPE=="MW")
      
      #mwRotate <- st_sf(st_sfc(),crs=3310)
      
      foreach(n = 1:nrow(mwsOnly), .combine = "rbind")%dopar%{
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
        
        
        # Index location of monitoring well
        index <-which(rotatedTimeDist$GISID == rotatedMw$GISID &
                        rotatedTimeDist$LOGDATE == rotatedMw$LOGDATE)
        
        rotatedTimeDist$LONGITUDE[index] <- rotatedMw$LONGITUDE
        rotatedTimeDist$LATITUDE[index] <- rotatedMw$LATITUDE
      }
    }
    },
      error=function(e){})

parEnd <- Sys.time()
