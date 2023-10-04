library(tidyverse)
library(sf)
library(here)
library(vroom)

exceedance <- vroom(here("Data_Preparation/Analysis_Data/Boundary_Detection_by_site.csv"))%>%
  filter(BoundIntrsct == TRUE)

# Import monitoring well points
sf <- st_read(here("Data_Preparation/Analysis_Data/Rotated_by_Site.gpkg"), layer = "All_Points")%>%
  filter(GLOBAL_ID %in% exceedance$GLOBAL_ID)

# Get unique wells
sf.unique <- sf%>%
  select(GLOBAL_ID,UID)%>%
  distinct()

# Determine which points were at the monitoring boundary
ma.lines <- st_read(here("Data_Preparation/Analysis_Data/Rotated_by_Site.gpkg"), layer = "Monitor_Area_Lines")%>%
  filter(GLOBAL_ID %in% sf.unique$GLOBAL_ID)%>%
  st_buffer(1)
## Iterate through sites and intersect the points
library(doParallel)
library(foreach)
library(tictoc)

cl <- makeCluster(5)  
registerDoParallel(cl)  

# Time it
tic()
bound.wells <- foreach(site = unique(sf.unique$GLOBAL_ID),.combine='rbind', .packages = c("dplyr","sf"))%dopar% {
  pts.site <- sf.unique%>%
    filter(GLOBAL_ID == site)
  
  line <- ma.lines%>%
    filter(GLOBAL_ID == site)%>%
    summarise()
  
  if(nrow(pts.site)>0 & nrow(line)==1){
    pts.site$boundary <- st_intersects(pts.site,line,sparse = FALSE)
    
    out <- pts.site%>%
      st_drop_geometry()
    
    return(out)
  }
  
}
registerDoSEQ()


toc()


bound.df <- sf%>%
  st_drop_geometry()%>%
  left_join(bound.wells)%>%
  drop_na(boundary)


sf$boundary <- st_intersects(sf,ma.lines, sparse = FALSE)

# Filter to monitoring wells at the monitoring well boundary
bound.wells <- sf%>%
  st_drop_geometry()%>%
  filter(boundary == TRUE)


# What is the bearing of the distance calculation for exceedance points?
# Max distance by bearing
max.val <- bound.df%>%
  filter(PARVAL >= 5)%>%
  # group_by(GLOBAL_ID)%>%
  # filter(Distance_m == max(Distance_m))%>%
  group_by(UID)%>%
  filter(PARVAL == max(PARVAL))%>%
  ungroup()%>%
  select(GLOBAL_ID,FIELD_PT_NAME,Distance_m,newBearing,PARVAL)%>%
  distinct()%>%
  mutate(bearingBin = cut(newBearing, breaks = 72, labels = seq(0,71)),
         bearingBin = as.numeric(bearingBin),
         bearingBin = bearingBin - 36,
         UID = paste0(GLOBAL_ID,"_",FIELD_PT_NAME))

bins <- max.val%>%
  group_by(bearingBin)%>%
  summarise(nLines = n(),
            med_PARVAL = median(PARVAL))

ggplot(bins)+
  geom_segment(x = 8.6, y = 0, xend = 8.6, yend = 200, color = "#7d420a", linetype = "solid")+
  geom_segment(x = -9.3, y = 0, xend = -9.3, yend = 200, color = "#7d420a", linetype = "solid")+
  geom_hline(yintercept = c(500,1000), color = "grey20")+
  geom_hline(yintercept = c(50,150), color = "grey60")+
  geom_col(aes(x = bearingBin, y = nLines, fill = med_PARVAL), color = NA,width = 1)+
  geom_col(aes(x = bearingBin, y = nLines),fill = NA, color = "grey20", linewidth = 0.1, width = 1, alpha = 0.5)+
  coord_polar(start = 3.15)+
  scale_fill_distiller(type = "seq",palette = "Spectral")+
  scale_y_continuous(breaks = seq(0, 200, 100))+
  annotate("label", x = 8, y = c(500,1000), label = c("500","1000"), size = 3)+
  annotate("text", x = -7.6, y = 130, label = "Up Gradient", color = "#7d420a", size = 3)+
  annotate("text", x = -6.3, y = 130, label = "Down Gradient", color = "#7d420a", size = 3)+
  annotate("text", x = -6, y = 1600, label = "Boundary Exceedances", color = "black", size = 5)+
  labs(title = "",
       fill = "Number of Distance Measures", x = "", y = "")+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12),
        legend.position = c(0.5, 0.3), legend.direction = "horizontal",
        plot.margin = margin(1, -4, -3, -4, "cm"),
        legend.background = element_rect(fill = "grey90", color = "black"),
        legend.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.background = element_rect(fill = "grey90"))

ggsave(filename = here("Water_Article/figures/Polar_Distance.png"), device = "png", dpi = 1000, plot = p,
       width = 8, units = "in")