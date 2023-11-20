#goal - to try MARSS analysis for chinook0 green river data

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

data_day <- read.csv(here("green", "data", "green_covariates_day.csv"))
data_night <- read.csv(here("green", "data", "green_covariates_night.csv"))

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"


#make columns for difference in temp, flow, photoperiod and call it 
#temp_diff, flow_diff, photo_diff

data_day$temp_diff <- c(NA,diff(data_day$temp))
data_day$flow_diff <- c(NA,diff(data_day$flow))
data_day$photo_diff <- c(NA,diff(data_day$photoperiod))

data_night$temp_diff <- c(NA,diff(data_night$temp))
data_night$flow_diff <- c(NA,diff(data_night$flow))
data_night$photo_diff <- c(NA,diff(data_night$photoperiod))

lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day)

data_day <- data_day %>% add_residuals(lm_temp_day)

lm_temp_night <- lm(temp ~ 1+photoperiod, data = data_night)

data_night <- data_night %>% add_residuals(lm_temp_night)

#need to go back to the python code and add atu to the dataframe - done


data_day_night <- rbind(data_day,data_night)

#plot chinook0_wild_perhour as a function of doy

ggplot(data_day_night, aes(x = doy, y = chinook0_wild_perhour)) + geom_point() + facet_wrap(~daytime_category)

#do same for chinook hatchery

ggplot(data_day_night, aes(x = doy, y = chinook0_hatchery_perhour)) + geom_point() + facet_wrap(~daytime_category)

#all the dates are messed up in data_night


# covariates_chinook0 <- arrange(data_day_night,doy) %>%
#   filter(doy >130 & doy <= 200) %>%
#   dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_solstice, 
#                 lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
#                 chinook0_hatchery_perhour_interpolate) %>%
#   pivot_wider(names_from = c(year, daytime_category), values_from = c(
#     temp, flow, photoperiod, atu_april, lunar_phase, 
#     resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_perhour_interpolate)) %>%
#   column_to_rownames(var = "doy") %>%
#   as.matrix() %>%
#   t()

