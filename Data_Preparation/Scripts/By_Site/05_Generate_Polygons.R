library(tidyverse)
library(sf)
library(here)

sf <- st_read(here("Data_Preparation/Analysis_Data/Rotated_by_Site.gpkg"), layer = "All_Points")%>%
  filter(Distance_m<1000)

# Create convex hulls of plumes
plumes <- sf%>%
  filter(PARVAL >= 5)%>%
  group_by(GLOBAL_ID)%>% 
  summarise()%>% 
  st_convex_hull()

plumes$Area_m <- st_area(plumes)

ggplot(plumes)+
  geom_histogram(aes(x = as.numeric(Area_m)))
  xlim(1,25000)+
  ylim(0,4000)

plumeRefine <- plumes%>%
  filter(as.numeric(Area_m) > 0)

st_write(plumeRefine,here("Data_Preparation/Analysis_Data/Rotated_by_Site.gpkg"),
         layer = "Plume_Polygons")

# Convert to linestring
plumeLines <- plumeRefine%>%
  st_cast("LINESTRING")

st_write(plumeLines,here("Data_Preparation/Analysis_Data/Rotated_by_Site.gpkg"),
         layer = "Plume_Lines")

# Create convex hulls of monitored areas
areas <- sf%>%
  group_by(GLOBAL_ID) %>% 
  summarise() %>% 
  st_convex_hull()

areas$Area_m <- st_area(areas)

areaRefine <- areas%>%
  filter(as.numeric(Area_m) > 0 & GLOBAL_ID %in% plumeRefine$GLOBAL_ID)

st_write(areaRefine,here("Data_Preparation/Analysis_Data/Rotated_by_Site.gpkg"),
         layer = "Monitor_Area_Polygons", append = FALSE)

# Convert to linestring
areaLines <- areaRefine%>%
  st_cast("LINESTRING")

st_write(areaLines,here("Data_Preparation/Analysis_Data/Rotated_by_Site.gpkg"),
         layer = "Monitor_Area_Lines", append = FALSE)
