library(tidyverse)
library(vroom)
library(here)

# Import exceedance / capture data
bd <- vroom(here("Data_Preparation/Analysis_Data/Boundary_Detection_by_site.csv"))%>%
  select(GLOBAL_ID, BoundIntrsct)
td <- vroom(here("Data_Preparation/Analysis_Data/Benzene_Time_Dist.tsv"))

# Calculate plume lengths by site

plumeDist <- td%>%
  filter(PARVAL >= 5 & Distance_m < 1000)%>%
  group_by(GLOBAL_ID)%>%
  filter(Distance_m == max(Distance_m))%>%
  ungroup()%>%
  select(GLOBAL_ID,Distance_m)%>%
  left_join(bd)%>%
  drop_na()%>%
  distinct()

fig3 <- ggplot(plumeDist)+
  geom_boxplot(aes(x = BoundIntrsct, y = Distance_m))+
  labs(title = "Plume Lengths by Capture / Exceedance",
       x = "", y = "Plume Length [m]")+
  scale_x_discrete(labels=c("FALSE" = "Plume Capture", "TRUE" = "Boundary Exceedance"))+
  coord_cartesian(ylim = c(0,150))

ggsave(plot = fig3, filename = here("figures/figure_03.pdf"), device = "pdf")  