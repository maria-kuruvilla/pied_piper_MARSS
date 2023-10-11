#Goal - to try MARSS package in a subset for chinook data with day and night separate
# data files

#steps

#interpolate the hatchery covariate values
#add residuals of photperiod and temperature to the dataset
#save the dataset
#create subset
#run MARSS

#Load libraries


library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)
library(MASS)
library(modelr)

#check system

if(.Platform$OS.type == "unix") {
  data_string = here("..","..","data","pied_piper","dungeness")
} else {
  data_string = here("..","..","..","OneDrive","Documents","data","pied_piper","dungeness")
}


#read data
data_day <- read.csv(here(data_string,"dungeness_unaggregated_day.csv"), header = TRUE)

data_night <- read.csv(here(data_string,"dungeness_unaggregated_night.csv"), header = TRUE)

#interpolate the hatchery covariate values
data_day$chinook0_hatchery_perhour_interpolate <- na.approx(data_day$chinook0_hatchery_perhour, na.rm = FALSE)
data_night$chinook0_hatchery_perhour_interpolate <- na.approx(data_night$chinook0_hatchery_perhour, na.rm = FALSE)

#add residuals of photoperiod and temperature to the dataset

lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day)

data_day <- data_day %>% add_residuals(lm_temp_day)

lm_temp_night <- lm(temp ~ 1+photoperiod, data = data_night)

data_night <- data_night %>% add_residuals(lm_temp_night)

#save the dataset

write.csv(data_day,here("data","dungeness_unaggregated_day_w_covariates.csv"))
write.csv(data_night,here("data","dungeness_unaggregated_night_w_covariates.csv"))


write.csv(data_day,here(data_string,"dungeness_unaggregated_day_w_covariates.csv"))
write.csv(data_night,here(data_string,"dungeness_unaggregated_night_w_covariates.csv"))

#subset the data



