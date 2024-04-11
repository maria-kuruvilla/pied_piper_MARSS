#Goal -make a figure for each river to visualise the correlation
#between all the environmental variables

#Load libraries
library(ggplot2)
library(tidyverse)
library(here)
library(GGally)
library(modelr)

#read the data

data <- read.csv(here("data", "dungeness_aggregated_w_trap_efficiency.csv"))

GGally::ggpairs(data %>% 
                  filter(doy > 120 & doy < 200) %>%
                  select(temp, flow, photoperiod,
                         atu_april, lunar_phase, resid, temp_diff,
                         flow_diff, photo_diff
                         ),
                  
                  aes(alpha = 0.2)) +
  # geom_point(alpha = 0.5)
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(n.breaks = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)
        )+
  theme_classic()

ggsave(here("output", "correlation_matrix_dungeness_subset.png"),
       width = 10, height = 10, units = "in", dpi = 300)

#do same for puyallup river

puyallup_data <- read.csv(here("puyallup",  "data", 
                               "puyallup_final.csv"))


GGally::ggpairs(puyallup_data %>%
                  # filter(doy > 120 & doy < 200) %>%
                  select(temp_night, flow_night, photoperiod_night,
                         atu_night, lunar_phase_night, resid_night, 
                         temp_diff_night, flow_diff_night, photo_diff_night,
                         secchi_depth_night
                  ) %>% 
                  #remove"night" from the column names
                  rename_with(~sub("_night", "", .x)),
                  
                  aes(alpha = 0.2)) +
  # geom_point(alpha = 0.5)
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(n.breaks = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)
        )+
  theme_classic()

ggsave(here("output", "correlation_matrix_puyallup.png"),
       width = 10, height = 10, units = "in", dpi = 300)



#do same for skagit

skagit_data <- read.csv(here("skagit","data", "skagit_2010-2022_w_covariates.csv"))
colnames(skagit_data)

lm_temp_day <- lm(temp ~ 1+photoperiod, data = skagit_data)

skagit_data <- skagit_data %>% add_residuals(lm_temp_day)

GGally::ggpairs(skagit_data %>%
                  filter(trap == "screw", daytime_category == "night") %>%
                  select(temp, flow, photoperiod,
                         atu_solstice, lunar_phase, resid, 
                         temp_diff, flow_diff, photo_diff
                  ),aes(alpha = 0.2)) +
  # geom_point(alpha = 0.5)
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(n.breaks = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)
        )+
  theme_classic()

ggsave(here("output", "correlation_matrix_skagit.png"),
       width = 10, height = 10, units = "in", dpi = 300)
