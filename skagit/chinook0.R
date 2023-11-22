#goal - to try MARSS analysis for chinook0 Skagit River data

#load libraries

library(MARSS)

library(here)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)
library(GGally)
library(ggplot2)
library(tidyverse)

#load data

data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2019_w_covariates.csv"))

#using ggplot plot chinook0_wild_perhour as a function of doy, facet wrap for years
#add a x limit from 150 to 200

ggplot(data_day_night, aes(x = doy, y = chinook0_wild_num)) +
  geom_point() +
  xlim(150,200)+
  facet_wrap(~year, ncol = 2, scales = "free_y")
  

#do same for chinook_hatchery_num

ggplot(data_day_night, aes(x = doy, y = chinook0_hatchery_num)) +
  geom_point() +
  xlim(150,200)+
  facet_wrap(~year, ncol = 2, scales = "free_y")

#doy limits should be 150 to 190

#making new columns for chinook0_wild_perhour and chinook0_hatchery_perhour and
#coho1_wild_perhour and coho1_hatchery_perhour and chinook1_wild_perhour and
#chinook1_hatchery_perhour

data_day_night$chinook0_wild_perhour <- data_day_night$chinook0_wild_num/data_day_night$In
data_day_night$chinook0_hatchery_perhour <- data_day_night$chinook0_hatchery_num/data_day_night$In
data_day_night$coho1_wild_perhour <- data_day_night$coho1_wild_num/data_day_night$In
data_day_night$coho1_hatchery_perhour <- data_day_night$coho1_hatchery_num/data_day_night$In
data_day_night$chinook1_wild_perhour <- data_day_night$chinook1_wild_num/data_day_night$In
data_day_night$chinook1_hatchery_perhour <- data_day_night$chinook1_hatchery_num/data_day_night$In

#plotting the perhour data for coho1

ggplot(data_day_night, aes(x = doy, y = coho1_wild_perhour)) +
  geom_point() +
  # xlim(150,190)+
  facet_wrap(~year, ncol = 2, scales = "free_y")

ggplot(data_day_night, aes(x = doy, y = coho1_hatchery_perhour)) +
  geom_point() +
  # xlim(150,190)+
  facet_wrap(~year, ncol = 2, scales = "free_y")
  
#doy limits should be 100 and 150

#plotting chinook1

ggplot(data_day_night, aes(x = doy, y = chinook1_wild_perhour)) +
  geom_point() +
  # xlim(150,190)+
  facet_wrap(~year, ncol = 2, scales = "free_y")

ggplot(data_day_night, aes(x = doy, y = chinook1_hatchery_perhour)) +
  geom_point() +
  # xlim(150,190)+
  facet_wrap(~year, ncol = 2, scales = "free_y")

#okay only doing for chinook and coho


#make new columns of interpolated values for chinook0_hatchery_perhour 
# and coho1_hatchery_perhour

data_day_night$chinook0_hatchery_perhour_inp <- na.approx(data_day_night$chinook0_hatchery_perhour, 
                                                             na.rm = FALSE)
data_day_night$coho1_hatchery_perhour_inp <- na.approx(data_day_night$coho1_hatchery_perhour, 
                                                             na.rm = FALSE)

lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day_night)

data_day_night <- data_day_night %>% add_residuals(lm_temp_day)


covariates_chinook0_skagit <- arrange(data_day_night,doy) %>%
  filter(doy >150 & doy <= 190) %>%
  dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_solstice,
                lunar_phase, resid, temp_diff, flow_diff, photo_diff,
                chinook0_hatchery_perhour_inp, trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = c(
    temp, flow, photoperiod, atu_solstice, lunar_phase,
    resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_perhour_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

#plot the correlation between temp,flow,photoperiod,resid,temp_diff,flow_diff,photo_diff

data_day_night %>%
  dplyr::select(flow,photoperiod, lunar_phase, temp, temp_diff, atu_solstice, resid, 
                photo_diff,flow_diff) %>%
  GGally::ggpairs(upper=list(continuous='points'),
                  lower=list(continuous='cor'))
#group 1 - photoperiod, temp, atu, photo_diff are highly correlated

ggsave(here("skagit", "output", "covariates_correlation_skagit.png"), width = 10, height = 10)




