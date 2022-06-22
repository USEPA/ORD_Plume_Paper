library(tidyverse)
library(sf)
library(vroom)
library(here)
library(stplanr)

# Determine which geotracker sites are LUST cleanup sites
lustSites <- vroom(here("Data_Preparation/Input_Data/Geotracker/sites.txt"),
                                   delim = "\t")%>%
  filter(CASE_TYPE == "LUST Cleanup Site")


# Find all field points at LUST cleanup sites which are monitoring wells with XY coordinates
xyFiles <- list.files(here("Data_Preparation/Input_Data/Geotracker/"),
                      pattern = "GeoXY.txt", full.names = TRUE)
xyAll <- vroom(xyFiles)%>%
  mutate(UID = paste0(GLOBAL_ID,"_",FIELD_PT_NAME))

# List all EDF files
edfFiles <- list.files(here("Data_Preparation/Input_Data/Geotracker/"),
                       full.names = TRUE, pattern = "EDF")

# Loop through EDF files for each county and extract benzene
# Options: c("MTBE","BZ","EBZ","BZME","XYLENES") - codes for btex
bz <- data.frame()
unitsCount <- data.frame()

for (n in 1:length(edfFiles)) {
  # Import a county
  cnty <- vroom(edfFiles[n])%>%
    filter(GLOBAL_ID %in% lustSites$GLOBAL_ID)%>%
    mutate(UID = paste0(GLOBAL_ID,"_",FIELD_PT_NAME))%>%
    filter(UID %in% xyAll$UID)%>% # Make sure we have XY data for sample
    filter(PARVQ %in% c("ND","<","="))%>%
    filter(PARLABEL == "BZ")%>%
    mutate(PARVAL = if_else(UNITS == "MG/L", PARVAL * 1000,PARVAL),
           UNITS = if_else(UNITS == "MG/L","UG/L",UNITS))%>%
    filter(UNITS == "UG/L")
  
  # Append county to statewide dataset
  bz <- bind_rows(bz,cnty)
  print(paste0(round(100*(n/length(edfFiles)),2),"% Complete at: ",Sys.time()))
}

# Import xy coordinates associated with sample points
xyFiles <- list.files(here("Data_Preparation/Input_Data/Geotracker/"),
                      full.names = TRUE,
                      pattern = "GeoXY.txt")

# Import sample point locations
# Some monitoring wells have the same exact location but different names,
# We need to combine these observations to avoid geometry issues in
# future analyses. We will create new IDs based on site ID and gps coordinates
geoXY <- vroom(xyFiles)%>%
  mutate(UID = paste0(GLOBAL_ID,"_",FIELD_PT_NAME),
         GISID = paste0(GLOBAL_ID,"_",LATITUDE,"_",LONGITUDE))

# Join xy data to benzene samples and filter out missing locations
bzXY <- bz%>%
  left_join(geoXY)%>%
  filter(!is.na(LATITUDE) & !is.na(LONGITUDE))

#Filter out results with "<" value greater than 5 and change <5 & "ND" to zero
filt <- bzXY%>%
  mutate(DROP = ifelse(PARVQ == "<" & PARVAL >5,"DROP","KEEP"))%>%
  filter(DROP == "KEEP")%>%
  select(!DROP)%>%
  mutate(PARVAL = ifelse(PARVQ == "<" | PARVQ == "ND",0,PARVAL))

# Map out maximum contaminant location of every site, using Benzene to start
benzeneOrigins <- filt%>%
  group_by(GLOBAL_ID)%>%
  filter(PARVAL == max(PARVAL))%>%
  filter(LOGDATE == min(LOGDATE))%>%
  ungroup()%>%
  select(GLOBAL_ID,GISID,PARVAL)%>%
  distinct()%>%
  filter(PARVAL > 5) # Only keep origin points where the MCL was broken

# Add point type to benzene xy data
bzXYtype <- bzXY%>%
  mutate(POINT_TYPE = if_else(GISID %in% benzeneOrigins$GISID,"ORIGIN","MW"))%>%
  drop_na(LATITUDE,LONGITUDE)%>%
  filter(GLOBAL_ID %in% benzeneOrigins$GLOBAL_ID)


# Iterate through sites and calculate distance and time from origin points
# to monitoring wells

mwDists <- data.frame()
n <- 0
for (site in unique(bzXYtype$GLOBAL_ID)) {
  origin <- bzXYtype%>%
    filter(GLOBAL_ID == site & POINT_TYPE == "ORIGIN")%>%
    filter(LOGDATE == min(LOGDATE))
  origin <- origin[1,]
  
  originSf <- origin%>%
    st_as_sf(coords = c("LONGITUDE","LATITUDE"),crs = 4326)%>%
    st_transform(3310)
  
  mws <- bzXYtype%>%
    filter(GLOBAL_ID == site & POINT_TYPE == "MW")%>%
    select(GISID,LATITUDE,LONGITUDE)%>%
    distinct()
  
  if(nrow(mws>0)){
    mwsSf <- mws%>%
      st_as_sf(coords = c("LONGITUDE","LATITUDE"),crs = 4326)%>%
      st_transform(3310)
    
    mws$Distance_m <- as.numeric(st_distance(originSf,mwsSf))
    
    # Calculate bearing
    
    ## Create lines
    bearings <- data.frame()
    for (l in 1:nrow(mwsSf)) {
      row <- mwsSf[l,]
      
      pts <- matrix(c(as.numeric(st_coordinates(originSf)[1]),
                      as.numeric(st_coordinates(row)[1]),
                      as.numeric(st_coordinates(originSf)[2]),
                      as.numeric(st_coordinates(row)[2])),,2)
      
      line = st_linestring(pts)%>%
        st_sfc(crs = 3310)%>%
        st_sf()%>%
        st_transform(4326)
      
      bearing <- line_bearing(line)
      
      newRow <- data.frame("GISID" = row$GISID, "Bearing_deg" = bearing)
      
      bearings <- bind_rows(bearings,newRow)
    }
    
    # Write bearing to table
    mws$Bearing <- bearings$Bearing_deg
    
    # combine data
    mwDists <- bind_rows(mwDists,mws)
    
    n <- n+1
    
    # Print status of loop
    if(n%%100==0){
      print(paste0(round(100*(n/length(unique(bzXYtype$GLOBAL_ID))),2),"% Complete at: ",Sys.time()))
    }
  }
}

# Add distance & Bearing back to benzene data & assign 0 distance to origin points.
bzD <- bzXYtype%>%
  left_join(mwDists)%>%
  mutate(Distance_m = if_else(POINT_TYPE == "ORIGIN",0,Distance_m))

# Calculate time differences
bzTD <- data.frame()
n <- 0
for (site in unique(bzD$GLOBAL_ID)) {
  originTime <- bzD%>%
    filter(GLOBAL_ID == site & PARVAL > 0)%>%
    filter(LOGDATE == min(lubridate::ymd(LOGDATE)))
  startDate <- lubridate::ymd(originTime$LOGDATE)[1]
  
  times <- bzD%>%
    filter(GLOBAL_ID == site)%>%
    mutate(Time = lubridate::ymd(LOGDATE) - startDate)
  
  bzTD <- bind_rows(bzTD,times)
  
  n <- n+1
  
  if(n%%100==0){
    print(paste0(round(100*(n/length(unique(bzD$GLOBAL_ID))),2),"% Complete at: ",Sys.time()))
  }
  
}

vroom_write(bzTD, here("Data_Preparation/Analysis_Data/Benzene_Time_Dist.tsv"),
            delim = "\t")
