library(tidyverse)
library(sf)
library(vroom)
library(here)

# Import boundary detection results
bd <- vroom(here("Data_Preparation/Analysis_Data/Boundary_Detection_by_site.csv"))%>%
  filter(BoundIntrsct == TRUE)


# Import Monitoring area lines
bound <- st_read(here("Data_Preparation/Analysis_Data/Rotated.gpkg"),
                 layer = "Monitor_Area_Lines")%>%
  separate(dateID, into = c("GLOBAL_ID","Date"),sep = "_")%>%
  filter(GLOBAL_ID %in% bd$GLOBAL_ID)


# Import Monitoring well data
pts <- st_read(here("Data_Preparation/Analysis_Data/Rotated.gpkg"),
               layer = "All_Points")%>%
  filter(PARVAL >= 5 & GLOBAL_ID %in% bd$GLOBAL_ID)

# Subset points that intersect the testing boundary
ptsIntersect <- st_intersection(pts, bound)
                                
# Summary statistics of benzene measurements at testing boundary
summary(ptsIntersect$PARVAL)

# Histogram of benzene measurements on testing boundary
ggplot(ptsIntersect)+
  geom_histogram(aes(x = PARVAL))

# Number of unique wells that tested above MCL at testing boundary
wells <- ptsIntersect%>%
  select(UID)%>%
  distinct()


