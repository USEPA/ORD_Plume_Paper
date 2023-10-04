############
# FIGURE 2 #
############

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(cowplot)

source(here("Functions/mean_linear_direction.R"))

# Left Plot: Example of site with monitoring wells and measured concentrations
# Example data:
exSf <- data.frame(ID = c("A","B","C","D","E"),
                   LATITUDE = c(25,50,30,75,10),
                   LONGITUDE = c(-20,-5,16,10,25),
                   PARVAL = c(20,55,40,15,0))%>%
  st_as_sf(coords = c("LONGITUDE","LATITUDE"),crs = 3310, remove = FALSE)
exOrigin <- data.frame(LATITUDE = 0, LONGITUDE = 0)%>%
  st_as_sf(coords = c("LONGITUDE","LATITUDE"),crs = 3310, remove = FALSE)

# Draw lines
exSf$startLat <- 0
exSf$startLon <- 0

exLines <- st_sf(st_sfc(),crs=3310)
for(h in 1:nrow(exSf)){
  exPts <- matrix(c(exSf$startLon,
                    exSf$LONGITUDE[h],
                    exSf$startLat,
                    exSf$LATITUDE[h]),,2)
  
  exLine = st_linestring(exPts)%>%
    st_sfc(crs = 3310)%>%
    st_sf()%>%
    bind_cols(st_drop_geometry(exSf[h,]))
  exLines <- rbind(exLines,exLine)
}

lp1 <- ggplot()+
  geom_sf(data = exLines, aes(color = "Drawn Lines"))+
  geom_sf(data = exSf, aes(color = "Monitoring Wells"), size = 2)+
  geom_sf(data = exOrigin, aes(color = "Leak"), size = 2)+
  scale_color_manual(name = "", values = c("Drawn Lines" = "black",
                                           "Monitoring Wells" = "orange",
                                           "Leak" = "red"))+
  geom_sf_label(data = exSf,aes(label=paste0(ID,": ",PARVAL," ug/L")),nudge_y = 6, nudge_x = 14, size = 2)+
  coord_sf(datum=st_crs(3310))+
  labs(x = "meters", y = "meters")+
  xlim(-35,60)+
  theme(legend.direction = "horizontal", legend.position = c(0.4,-0.4),
        plot.margin = margin(0,0,1.5,0, "cm"))

# Right Plot: Example of vector math to calculate mean linear direction
exMWs <- st_drop_geometry(exSf)
start_x <- 0
start_y <- 0
crs <- 3310
startLon <- start_x
startLat <- start_y

weightedLines <- st_sf(st_sfc(),crs=3310)

for (j in 1:nrow(exMWs)) {
  #print(paste0("Iteration #",j))
  # Select one monitoring well
  mw <- exMWs[j,]
  rowWeight <- exMWs$PARVAL[j]
  # Create points matrix
  bearingPts <- matrix(c(start_x,
                         mw$LONGITUDE,
                         start_y,
                         mw$LATITUDE),,2)
  
  # Create a buffer from the origin point, using the weighted length
  # then select the point on the buffer closest to the monitoring well
  # Then draw a new line from the origin point to the selected point on
  # the buffer line. This is your new weighted line.
  
  bearingLine = st_linestring(bearingPts)%>%
    st_sfc(crs = 3310)%>%
    st_sf()%>%
    bind_cols(mw)
  
  #angle <- stplanr::line_bearing(bearingLine)
  
  length <- bearingLine%>%
    st_transform(crs = crs)%>%
    st_length()%>%
    as.numeric()
  
  weightedLength <- as.numeric(length)*as.numeric(rowWeight)
  
  if(weightedLength > 0){
    # Project mw point
    mwSf <- data.frame("X"=mw$LONGITUDE, "Y" = mw$LATITUDE)%>%
      st_as_sf(coords = c("X","Y"), crs = 3310)%>%
      st_transform(crs = crs)
    
    # Buffer then convert to points
    buf <- exOrigin%>%
      st_buffer(length)%>%
      st_cast("POINT")
    
    # Select nearest point to the monitoring well
    near <- st_nearest_feature(mwSf,buf)
    weightedEnd <- buf[near,]
    # weightedEndWGS <- weightedEnd%>%
    #   st_transform(4326)
    weCoords <- as.data.frame(st_coordinates(weightedEnd))
    
    # Calculate new endpoint from previous start point
    # Do this by calculating delta X/Y from origin, then
    # adding it to previous end point
    deltLat <- weCoords$Y - start_y
    endLat <- startLat + deltLat
    deltLon <- weCoords$X - start_x
    endLon <- startLon + deltLon
    
    # Create points matrix
    pts <- matrix(c(startLon,
                    endLon,
                    startLat,
                    endLat),,2)
    
    newLine = st_linestring(pts)%>%
      st_sfc(crs = 3310)%>%
      st_sf()%>%
      bind_cols(mw)
    
    # Test plot
    # ggplot()+
    #   geom_sf(data = originSf, shape = 4)+
    #   geom_sf(data = newLine)+
    #   geom_sf(data = mwSf, fill = "red", shape = 23)+
    #   geom_sf(data = buf)+
    #   geom_sf(data = weightedEnd, fill = "blue",shape = 21, size = 3)+
    #   geom_sf_label(data = mwSf, aes(label = PARVAL), nudge_y = 25, nudge_x = 25)+
    #   labs(title = "Weighting Directional Lines", x = "Longitude",y="Latitude")+
    #   theme(axis.text.x = element_text(angle = 70, vjust = 1.3, hjust=1.3))
    
    
    
    
    weightedLines <- rbind(weightedLines,newLine)
    
    # Update starting point for next iteration
    startLat <- endLat
    startLon <- endLon
    
  }
  
  
  
  # Once we reach the end, output the last endpoint to determine
  # mean angle
  
  if(j == nrow(exMWs)){
    #print("Finding Mean Direction ...")
    endPt <- matrix(c(exOrigin$LONGITUDE,
                      endLon,
                      exOrigin$LATITUDE,
                      endLat),,2)
    
    meanDirectionLine <- st_linestring(endPt)%>%
      st_sfc(crs = 3310)%>%
      st_sf()%>%
      st_transform(4326)
    
    bearing <- stplanr::line_bearing(meanDirectionLine)
    
    #meanDirectionLine$Bearing <- bearing
  }
}

# Get centroids for labels
exCntrds <- st_centroid(weightedLines)
mldCntrd <- st_centroid(meanDirectionLine)

# Endpoints for plotting
weightedPts <- st_cast(weightedLines, "POINT")

rp1 <- ggplot(weightedLines)+
  geom_sf()+
  geom_sf(data= meanDirectionLine, aes(color = "Mean Plume Direction"), linetype = "dashed")+
  geom_sf(data = weightedPts,color = "orange", size = 1.5,show.legend = FALSE)+
  geom_sf_label(data = exCntrds, aes(label = ID),nudge_x = 0, size = 2)+
  geom_sf_label(data = mldCntrd, label = paste0(round(bearing,2),"°")  ,nudge_x = 16, size = 2)+
  geom_sf(data = exOrigin, color = "red", size = 2)+
  scale_color_manual(name = "", values = c("Mean Plume Direction" = "blue"))+
  coord_sf(datum=st_crs(3310))+
  labs(x = "meters", y = "meters")+
  xlim(-80,80)+
  theme(legend.direction = "horizontal", legend.position = c(0.1,-0.4),#c(0.1,-0.2)
        plot.margin = margin(0,0,1.5,0, "cm"))

plot_row <- plot_grid(lp1,rp1,nrow=1, align = "h",labels = c('A', 'B'),
                      label_x = .2,
                      label_y = .9)  

# now add the title
title <- ggdraw() + 
  draw_label(
    "Example of Plume Mean Linear Direction Calculation",
    fontface = 'bold',
    x = 0.1,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

fig2 <- plot_grid(
  title, plot_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

fig2

#ggsave(plot = fig2,filename = here("figures/figure_02.pdf"), device = "pdf", height = 4, width = 6, units = "in")