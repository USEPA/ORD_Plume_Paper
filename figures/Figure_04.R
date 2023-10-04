library(here)
library(vroom)
library(ggplot2)
library(dplyr)

# Load Boundary Detection
bd <- vroom(here("Data_Preparation/Analysis_Data/Boundary_Detection_by_site.csv"))%>%
  filter(max_Dist < 1000)

# Load time and distance
td <- vroom(here("Data_Preparation/Analysis_Data/Benzene_Time_Dist.tsv"))

uniqueMW <- td%>%
  filter(POINT_TYPE == "MW" & GLOBAL_ID %in% bd$GLOBAL_ID)%>%
  select(GISID,Distance_m)%>%
  distinct()

distQuantiles <- as.data.frame(quantile(uniqueMW$Distance_m,c(.25,.5,.75,.9)))%>%
  mutate(Distance = round(.[[1]],2))%>%
  select(Distance)

uniqueMWPdf <- uniqueMW%>%
  mutate(Distance_m = ifelse(Distance_m > 200,200,Distance_m))

fig4 <- ggplot(uniqueMWPdf)+
  geom_histogram(aes(x = Distance_m),binwidth = 5,
                 fill = "white",color = "black")+
  geom_segment(data = distQuantiles, aes ( x = Distance,xend = Distance,y = 0, yend = 6300, color = rownames(distQuantiles)),size = 1.5)+
  guides(color=guide_legend(title="Percentile"))+
  scale_color_manual(values = c("#fbb4b9","#f768a1","#c51b8a","#7a0177"))+
  scale_x_continuous(breaks=c(0,50,100,150,200),
                     labels=c("0","50", "100","150",">200"))+
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  labs(title = "Monitoring Well Distances",
       x = "Distance [m]",y = "Count of Monitoring Wells")

fig4

ggsave(plot = fig4, filename = here("figures/figure_04.pdf"), device = "pdf",
       height = 4, width = 6, units = "in")