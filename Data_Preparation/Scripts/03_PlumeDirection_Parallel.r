## mean linear direction

## This script is an optional implementation to find plume directions
## It uses parallel processing. By breaking the input data into 4 separate
## data sets and using 71 cores, I was able to complete this in about 15 hours.

library(tidyverse)
library(sf)
library(spatialEco)
library(spdep)
library(vroom)
library(stplanr)
library(doParallel)


print(paste0("Loading Data: ",Sys.time()))


# Load mean linear direction function
source("/work/GRDVULN/github/Plume_Paper/Functions/mean_linear_direction.R")

# Load testing date and counts data and refine for site selection
selection <- vroom("/work/GRDVULN/github/Plume_Paper/data/testCountsAndDates.csv")%>%
  filter(Tests_Above_MCL > 2)%>%
  mutate(dateID = paste0(GLOBAL_ID,"_",OriginTestDate))

# To save time, you can break the dataset into pieces, here we break into 4
# and we will run them in different instances.
nSites <- length(unique(selection$GLOBAL_ID))
breaks <- seq(1,nSites,round(nSites/4))
break1 <- breaks[2]+1
break2 <- breaks[3]+1
break3 <- breaks[4]+1

sub1 <- unique(selection$GLOBAL_ID)[1:breaks[2]]
sub2 <- unique(selection$GLOBAL_ID)[break1:breaks[3]]
sub3 <- unique(selection$GLOBAL_ID)[break2:breaks[4]]
sub4 <- unique(selection$GLOBAL_ID)[break3:nSites]

selection <- selection%>%
  filter(GLOBAL_ID %in% sub1)

# Import Benzene data and refine to selection
timeDist <- vroom("/work/GRDVULN/github/Plume_Paper/data/Benzene_Time_Dist.tsv")%>%
  filter(GLOBAL_ID %in% sub1)


print(paste0("Started at: ", Sys.time()))

# Create Cluster
cores <- detectCores()-1
cl <- makeCluster(cores)
registerDoParallel(cl)

print(paste0("Starting Processing of ", length(selection$dateID), " plumes using ", cores, " cores at: ", Sys.time()))

direction <- foreach(siteDate = unique(selection$dateID),
                       .combine=rbind,
                       .packages = c("tidyverse","sf","spatialEco",
                                     "spdep","vroom","here","stplanr")) %dopar%
  tryCatch({
    # Filter to the site and date range
    row <- selection%>%
      filter(dateID == siteDate)
    
    date <- row$OriginTestDate
    
    origin <- timeDist%>%
      filter(GLOBAL_ID == row$GLOBAL_ID &
               LOGDATE == date &
               POINT_TYPE == "ORIGIN"
      )
    
    mws <- timeDist%>%
      filter(GLOBAL_ID == row$GLOBAL_ID &
               LOGDATE >= date - 28 &
               LOGDATE <= date +28 &
               POINT_TYPE == "MW"
      )
    
    # Run the weighted mean direction
    bearing <- mld(mws,coords = c("LONGITUDE","LATITUDE"),start_x = origin$LONGITUDE[1],
                   start_y = origin$LATITUDE[1], date = str_replace_all(as.character(date),"-","_"),crs = 3310,
                   weight = mws$PARVAL,filename = here("Data_Preparation/Analysis_Data/meanDirection.gpkg"))
    
    data.frame("GLOBAL_ID" = row$GLOBAL_ID,
               "OriginTestDate" = row$OriginTestDate,
               "PlumeBearing" = bearing)
    
  },
  error=function(e){})
  
    
  
print(paste0("Completed at: ", Sys.time()))
  
vroom_write(direction, "/work/GRDVULN/github/Plume_Paper/PlumeDirections_1.csv", delim = ",")



