library(here)
library(vroom)
library(tidyverse)
library(cowplot)
# Load time and distance
td <- vroom(here("Data_Preparation/Analysis_Data/Benzene_Time_Dist.tsv"))

# 10 ug/L
maxDist10 <- td%>%
  filter(PARVAL >= 10)%>%
  group_by(GLOBAL_ID)%>%
  mutate(maxDist = max(Distance_m))%>%
  ungroup()%>%
  select(GLOBAL_ID, maxDist)%>%
  distinct()%>%
  filter(maxDist < 1000)

maxDist10Pdf <- maxDist10%>%
  mutate(maxDist = if_else(maxDist > 600,600, maxDist))%>%
  filter(maxDist>1)

# 5 ug/L
maxDist5 <- td%>%
  filter(PARVAL >= 5)%>%
  group_by(GLOBAL_ID)%>%
  mutate(maxDist = max(Distance_m))%>%
  ungroup()%>%
  select(GLOBAL_ID, maxDist)%>%
  distinct()%>%
  filter(maxDist < 1000)

maxDist5Pdf <- maxDist5%>%
  mutate(maxDist = if_else(maxDist > 600,600, maxDist))%>%
  filter(maxDist > 1)


# Get ninetieth percentile for 10 ug/L
ninety10 <- quantile(maxDist10$maxDist,.9)
# Convert from meters to feet
ft10 <- ninety10*3.28084

# Get ninetieth percentile for 5 ug/L
ninety5 <- quantile(maxDist5$maxDist,.9)

# --- DISTANCE PLOTS ---
max10Plot <- ggplot(maxDist10Pdf)+
  geom_histogram(aes(x = maxDist),binwidth = 10,color = "black",fill = "white")+
  xlim(0,640)+
  ylim(0,1000)+
  geom_segment(x=91.44 ,y =0, xend=91.44, yend =800)+
  annotate("text", x = 200, y = 780, label = "Reisinger et al. (2000)")+
  geom_segment(x=137 ,y =0, xend=137, yend =700 )+
  annotate("text", x = 240, y = 680, label = "Mace & Choi (1998)")+
  geom_segment(x=524 ,y =0, xend=524, yend =600 )+
  annotate("text", x = 435, y = 570, label = "Rice et al. (1995)")+
  labs(title = "Maximum Distance of Benzene > 10 ug/L (A) and 5 ug/L (B)",
       x = "Meters",y = "Count of Sites")+
  scale_x_continuous(breaks=c(0,200,400,600),
                     labels=c("0", "200", "400","> 600"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = margin(0, 0, -.1, 0, "cm"))

max5Plot <- ggplot(maxDist5Pdf)+
  geom_histogram(aes(x = maxDist),binwidth = 10,color = "black",fill = "white")+
  xlim(0,640)+
  ylim(0,1000)+
  geom_segment(x=168 ,y =0, xend=168, yend =700 )+
  annotate("text", x = 250, y = 680, label = "Shi et al. (2004)")+
  geom_segment(x=305 ,y =0, xend=305, yend =600)+
  annotate("text", x = 405, y = 580, label = "Riffai & Rixey (2004)")+
  geom_segment(x=502 ,y =0, xend=502, yend =450 )+
  annotate("text", x = 435, y = 380, label = str_wrap("Kamath et al. (2012)",13))+
  labs(x = "Meters",y = "Count of Sites")+
  scale_x_continuous(breaks=c(0,200,400,600),
                     labels=c("0", "200", "400","> 600"))+
  theme(plot.margin = margin(-.1, 0, 0, 0, "cm"))

fig5 <- plot_grid(max10Plot,max5Plot, align = "v",
                  nrow = 2,
                  labels = c("A","B"),
                  label_x = .1,
                  label_y = c(.85,.95))

ggsave(plot = fig5, filename = here("figures/figure_05.pdf"), device = "pdf",
       height = 5, width = 6, units = "in")