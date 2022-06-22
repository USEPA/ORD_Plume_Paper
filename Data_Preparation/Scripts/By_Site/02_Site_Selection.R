library(tidyverse)
library(vroom)
library(here)

# Import Benzene data
timeDist <- vroom(here("Data_Preparation/Analysis_Data/Benzene_Time_Dist.tsv"))

### STEP 1: Identify sites with > 2 Monitoring Wells > 5 ug/L
candidates <- data.frame()
nSites <- length(unique(timeDist$GLOBAL_ID))
n <- 0
for (site in unique(timeDist$GLOBAL_ID)) {
  siteFilt <- timeDist%>%
    filter(GLOBAL_ID == site)%>%
    group_by(GISID)%>%
    mutate(maxPARVAL = max(PARVAL))%>%
    select(GLOBAL_ID,GISID,maxPARVAL)%>%
    filter(maxPARVAL >= 5)%>%
    distinct()
  
  newRow <- data.frame(GLOBAL_ID = site, nWellsGT5 = nrow(siteFilt))
  
  candidates <- rbind(candidates,newRow)
  n <- n+1
  print(paste0(round(100*(n/nSites),2),"% Complete at: ",Sys.time()))
  
}

# Plot results
ggplot(candidates)+
  geom_histogram(aes(x = nWellsGT5))


vroom_write(candidates, here("Data_Preparation/Analysis_Data/testCounts.csv"),
            delim = ",")

