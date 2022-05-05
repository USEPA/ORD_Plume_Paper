library(tidyverse)
library(vroom)
library(here)

# Import Benzene data
timeDist <- vroom(here("Data_Preparation/Analysis_Data/Benzene_Time_Dist.tsv"))

### STEP 1: Identify sites and dates with > 2 Monitoring Wells > 5 ug/L
testDates <- data.frame()
n=0
nSites <- length(unique(timeDist$GLOBAL_ID))
for (site in unique(timeDist$GLOBAL_ID)) {
  # Filter records for site
  sitePts <- timeDist%>%
    filter(GLOBAL_ID == site)
  
  # Look at origin point only
  originPt <- sitePts%>%
    filter(POINT_TYPE == "ORIGIN")
  
  # List dates that the origin was tested
  originDates <- unique(originPt$LOGDATE)
  
  # Loop through dates and figure out how many monitoring wells were tested within 28 days
  for(i in 1:length(originDates)){
    date <- lubridate::ymd(originDates[i])
    #print(date)
    wellsTested <- sitePts%>%
      filter(POINT_TYPE == "MW" & LOGDATE > date - 28 & LOGDATE < date + 28)
    overMCL <- wellsTested%>%
      filter(PARVAL >= 5)
    
    newRow <- data.frame("GLOBAL_ID" = site,
                         "OriginTestDate" = date,
                         "Total_Tests" = nrow(wellsTested),
                         "Tests_Above_MCL" = nrow(overMCL))
    testDates <- rbind(testDates,newRow)
  }
  
  n <- n+1
  
  if(n %% 100 ==0){
    print(paste0(round(100*(n/nSites),2),"% Complete --- ",Sys.time()))
  }
}

vroom_write(testDates, here("Data_Preparation/Analysis_Data/testCountsAndDates.csv"))

# Plot percent of tests above MCL
pct <- testDates%>%
  mutate(PctMCL = 100*(Tests_Above_MCL/Total_Tests))

ggplot(pct)+
  geom_histogram(aes(x = PctMCL),binwidth = 5)