#Goal - plot all the hatchery and wild data for puyallup coho

library(ggplot2)
library(here)
library(zoo)
library(MASS)
library(GGally)
library(tidyverse)
library(broom)
library(ggpubr)

data <- read.csv(here("puyallup", "data","puyallup_final.csv"),
                 na.strings = c("NA",""))

data$coho1_hatchery_perhour_night[dim(data)[1]] <- 0
data$coho1_hatchery_perhour_day[dim(data)[1]] <- 0

puyallup_coho <- arrange(data, doy) %>% 
  filter(doy > 90, doy < 160) %>% 
  mutate(wild_day = coho1_wild_perhour_day, 
         wild_night = coho1_wild_perhour_night,
         hatchery_day = coho1_hatchery_perhour_day,
         hatchery_night = coho1_hatchery_perhour_night
         ) %>% 
  select(year, doy, wild_day, wild_night, hatchery_day, hatchery_night) %>% 
  pivot_wider(names_from = year, values_from = c(wild_day, wild_night,
                                                 hatchery_day,
                                                 hatchery_night))

#scale the data

for(i in 2:(dim(puyallup_coho)[2])){
  puyallup_coho[,i] = scale(puyallup_coho[,i])[,1]
  print(i)
  print(colnames(puyallup_coho)[i])
}

puyallup_coho_long <- puyallup_coho %>%
  pivot_longer(cols = -doy,
               names_to = c("origin", "daytime_category","year"),
               names_sep = "_",
               values_to = "coho_per_hour")

ggplot(puyallup_coho_long %>% 
         filter(daytime_category == "day")) +
  geom_line(aes(x = doy, y = coho_per_hour, color = origin), alpha = 0.5, linewidth = 1)+
  facet_wrap(~year, ncol = 5)+
  scale_color_manual(values = c("wild" = "salmon", "hatchery" = "cadetblue"))+
  theme_classic()+
  labs(x = "Day of Year", y = "Scaled Coho per Hour", color = "Origin")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white", color = "white")
  )   

ggsave(here("puyallup", "output", "puyallup_coho_day.png"), width = 12, height = 8, dpi = 300)

ggplot(puyallup_coho_long %>% 
         filter(daytime_category == "night")) +
  geom_line(aes(x = doy, y = coho_per_hour, color = origin), alpha = 0.5, linewidth = 1)+
  facet_wrap(~year, ncol = 5)+
  scale_color_manual(values = c("wild" = "salmon", "hatchery" = "cadetblue"))+
  theme_classic()+
  labs(x = "Day of Year", y = "Scaled Coho per Hour", color = "Origin")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white", color = "white")
  )

ggsave(here("puyallup", "output", "puyallup_coho_night.png"), width = 12, height = 8, dpi = 300)
