mld <- function(x,coords,start_x,start_y,date="All",crs,weight,filename){
  st <- Sys.time()
  # These will be the initial starting coordinates and will
  # be updated every time the loop loops.
  startLon <- start_x
  startLat <- start_y
  
  originSf <- data.frame("X"=start_x, "Y" = start_y)%>%
    st_as_sf(coords = c("X","Y"), crs = 4326)%>%
    st_transform(crs = crs)
  
  weightedLines <- st_sf(st_sfc(),crs=4326)
  
  for (n in 1:nrow(x)) {
    #print(paste0("Iteration #",n))
    # Select one monitoring well
    mw <- x[n,]
    rowWeight <- weight[n]
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
      st_sfc(crs = 4326)%>%
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
        st_as_sf(coords = c("X","Y"), crs = 4326)%>%
        st_transform(crs = crs)
      
      # Buffer then convert to points
      buf <- originSf%>%
        st_buffer(weightedLength)%>%
        st_cast("POINT")
      
      # Select nearest point to the monitoring well
      near <- st_nearest_feature(mwSf,buf)
      weightedEnd <- buf[near,]
      weightedEndWGS <- weightedEnd%>%
        st_transform(4326)
      weCoords <- as.data.frame(st_coordinates(weightedEndWGS))
      
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
        st_sfc(crs = 4326)%>%
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
    
    if(n == nrow(x)){
      #print("Finding Mean Direction ...")
      endPt <- matrix(c(origin$LONGITUDE,
                        endLon,
                        origin$LATITUDE,
                        endLat),,2)
      
      meanDirectionLine <- st_linestring(endPt)%>%
        st_sfc(crs = 4326)%>%
        st_sf()
      
      bearing <- stplanr::line_bearing(meanDirectionLine)
      
      #meanDirectionLine$Bearing <- bearing
    }
  }
  et <- Sys.time()
  elap <- round(as.numeric(difftime(et,st,units = "secs")),1)
  print(paste0("Calculated mean linear direction of ",nrow(x)," points in ", elap," seconds ... bearing = ",bearing))
  # st_write(weightedLines,filename,
  #          layer = "Weighted_Lines",append = TRUE)
  # st_write(meanDirectionLine, filename,
  #          layer = "Mean_Linear_Direction",append = TRUE)
  
  return(bearing)
}
