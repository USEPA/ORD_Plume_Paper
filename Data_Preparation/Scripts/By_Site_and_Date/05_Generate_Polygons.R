library(tidyverse)
library(sf)
library(here)

sf <- st_read(here("Data_Preparation/Analysis_Data/Rotated.gpkg"), layer = "All_Points")

# Create convex hulls of plumes
plumes <- sf%>%
  filter(PARVAL >= 5)%>%
  group_by(dateID) %>% 
  summarise() %>% 
  st_convex_hull()

plumes$Area_m <- st_area(plumes)

ggplot(plumes)+
  geom_histogram(aes(x = as.numeric(Area_m)))+
  xlim(0,10000)+
  ylim(0,8000)

plumeRefine <- plumes%>%
  filter(as.numeric(Area_m) > 0)

st_write(plumeRefine,here("Data_Preparation/Analysis_Data/Rotated.gpkg"),
         layer = "Plume_Polygons")

# Convert to linestring
plumeLines <- plumeRefine%>%
  st_cast("LINESTRING")

st_write(plumeLines,here("Data_Preparation/Analysis_Data/Rotated.gpkg"),
         layer = "Plume_Lines")

# Create convex hulls of monitored areas
areas <- sf%>%
  group_by(dateID) %>% 
  summarise() %>% 
  st_convex_hull()

areas$Area_m <- st_area(areas)

areaRefine <- areas%>%
  filter(as.numeric(Area_m) > 0 & dateID %in% plumeRefine$dateID)

st_write(areaRefine,here("Data_Preparation/Analysis_Data/Rotated.gpkg"),
         layer = "Monitor_Area_Polygons", append = FALSE)

# Convert to linestring
areaLines <- areaRefine%>%
  st_cast("LINESTRING")

st_write(areaLines,here("Data_Preparation/Analysis_Data/Rotated.gpkg"),
         layer = "Monitor_Area_Lines", append = FALSE)
