library(tidyverse)
library(sf)
library(vroom)
library(here)

plumes <- st_read(here("Data_Preparation/Analysis_Data/Rotated.gpkg"),
         layer = "Plume_Lines")

areas <- st_read(here("Data_Preparation/Analysis_Data/Rotated.gpkg"),
                  layer = "Monitor_Area_Lines")


boundaryIntersection <- data.frame()

# Iterate detection by plume
for(n in 1:nrow(plumes)){
  #subset plume
  plume <- plumes[n,]
  
  #subset monitoring area
  area <- areas%>%
    filter(dateID == plume$dateID)
  
  #Calculate intersection
  intrsct <- as.numeric(st_intersects(plume,area))
  
  newRow <- data.frame(dateID = plume$dateID, BoundIntrsct = intrsct)

  boundaryIntersection <- rbind(boundaryIntersection,newRow)
  
  print(paste0(round(100*(n/nrow(plumes)),2),"% Complete at ",Sys.time()))
}

# Change NA to FALSE
output <- boundaryIntersection%>%
  mutate(BoundIntrsct = if_else(is.na(BoundIntrsct),FALSE,TRUE))%>%
  separate(dateID, into = c("GLOBAL_ID","OriginTestDate"), sep = "_")

# Add some relevant statistics about each site

td <- vroom(here("Data_Preparation/Analysis_Data/Benzene_Time_Dist.tsv"))

stats <- td%>%
  group_by(GLOBAL_ID)%>%
  mutate(test_length = difftime(max(LOGDATE),min(LOGDATE),units = "days"),
         max_Parval = max(PARVAL),
         mean_Parval = mean(PARVAL),
         median_Parval = median(PARVAL),
         tot_Tests = n(),
         max_Dist = max(Distance_m),
         mean_Dist = mean(Distance_m),
         median_Dist = median(Distance_m))%>%
  ungroup()%>%
  group_by(GLOBAL_ID, FIELD_PT_NAME)%>%
  mutate(n_Wells = n())%>%
  ungroup()%>%
  select(GLOBAL_ID,test_length,max_Parval,mean_Parval,median_Parval,
         tot_Tests,max_Dist,mean_Dist,median_Dist,n_Wells)%>%
  distinct()
  
join <- stats%>%
  left_join(output)%>%
  filter(GLOBAL_ID %in% output$GLOBAL_ID)%>%
  distinct()

vroom_write(join,here("Data_Preparation/Analysis_Data/Boundary_Detection.csv"),
            delim = ",")
