library(here)
library(vroom)
library(sf)
library(tidyverse)


# Load Boundary Detection
bd <- vroom(here("Data_Preparation/Analysis_Data/Boundary_Detection_by_site.csv"))%>%
  filter(max_Dist < 1000)

# Load plume convex hulls
plumeAreas <- st_read(here("Data_Preparation/Analysis_Data/Rotated_by_Site.gpkg"),layer = "Plume_Lines")%>%
  filter(GLOBAL_ID %in% bd$GLOBAL_ID)

# Load monitoring convex hulls
monitorAreas <- st_read(here("Data_Preparation/Analysis_Data/Rotated_by_Site.gpkg"),layer = "Monitor_Area_Lines")%>%
  filter(GLOBAL_ID %in% bd$GLOBAL_ID)

# Join boundary detection to plume lines
areasBound <- plumeAreas%>%
  left_join(bd, by = "GLOBAL_ID")

# Filter on TRUE or FALSE
plumeAreasPC <- areasBound%>%
  filter(BoundIntrsct == FALSE)

plumeAreasBE <- areasBound%>%
  filter(BoundIntrsct == TRUE)

# Join boundary detection to monitoring area lines
maBound <- monitorAreas%>%
  left_join(bd, by = "GLOBAL_ID")

#Filter on TRUE or FALSE
maPC <- maBound%>%
  filter(BoundIntrsct == FALSE)

maBE <- maBound%>%
  filter(BoundIntrsct == TRUE)

# Add info for colors and facets
paPCsel <- plumeAreasPC%>%
  select(GLOBAL_ID,Area_m)%>%
  mutate(facet = "Capture : Plumes",
         cat = "Plumes",
         N = length(facet))
paBEsel <- plumeAreasBE%>%
  select(GLOBAL_ID,Area_m)%>%
  mutate(facet = "Exceedance : Plumes",
         cat = "Plumes",
         N = length(facet))
maPCsel <- maPC%>%
  select(GLOBAL_ID,Area_m)%>%
  mutate(facet = "Capture : Monitor Areas",
         cat = "Monitor Areas",
         N = length(facet))
maBEsel <- maBE%>%
  select(GLOBAL_ID,Area_m)%>%
  mutate(facet = "Exceedance : Monitor Areas",
         cat = "Monitor Areas",
         N = length(facet))

allAreas <- rbind(maPCsel,paPCsel,maBEsel,paBEsel)

areaPlots <- ggplot(allAreas)+
  geom_sf(alpha=.05, aes(color = cat))+
  scale_color_manual(values = c("#1b9e77","#d95f02"))+
  xlim(-400,400)+
  ylim(-200,200)+
  coord_sf(datum = st_crs(3310))+
  geom_segment(x=0,xend=0,y=-100,yend=100,linetype=3)+
  geom_segment(x=-100,xend=100,y=0,yend=0,linetype=3)+
  geom_label(x = 275,y = -175, aes(label =  paste0("N=",format(N,big.mark =","))), size = 2)+
  theme(plot.margin = margin(0,0,0,0, unit = "cm"),
        legend.position="bottom")+
  facet_wrap(~facet)+
  theme(plot.margin=unit(c(0,0,0,0), "cm"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1),title=""))

areaPlots

ggsave(plot = areaPlots, filename = here("figures/figure_06.pdf"), device = "pdf",
       height = 5, width = 6, units = "in")
