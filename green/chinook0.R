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


lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day)

data_day <- data_day %>% add_residuals(lm_temp_day)

lm_temp_night <- lm(temp ~ 1+photoperiod, data = data_night)

data_night <- data_night %>% add_residuals(lm_temp_night)

#need to go back to the python code and add atu to the dataframe