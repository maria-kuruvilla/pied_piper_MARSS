#goal - run marss analysis for chinook0 in Puyallup
#with temp data for years >2010

#load libraries
library(MARSS)
library(ggplot2)
library(here)
library(GGally)
library(tidyverse)
library(zoo)
library(modelr)


#load data
data <- read.csv(here("puyallup", "data","puyallup_2004-2021_all_days_w_covariates_w_temp.csv"),
                 na.strings = c("NA",""))


data$chinook0_hatchery_perhour_night[dim(data)[1]] <- 0
data$chinook0_hatchery_perhour_day[dim(data)[1]] <- 0

#when secchi_depth_day is NA, replace with secchi_depth values
#when secchi_depth_night is NA, replace with secchi_depth values

data$secchi_depth_day[is.na(data$secchi_depth_day)] <- data$secchi_depth[is.na(data$secchi_depth_day)]

data$secchi_depth_night[is.na(data$secchi_depth_night)] <- data$secchi_depth[is.na(data$secchi_depth_night)]

#interpolating na values for flow_day, flow_night, sechhi_depth_day, secchi_depth_night, hatchery

data <- data %>% 
  mutate(flow_day_inp = na.approx(flow_day, na.rm = FALSE),
         flow_night_inp = na.approx(flow_night, na.rm = FALSE),
         temp_day_inp = na.approx(temp_day, na.rm = FALSE),
         temp_night_inp = na.approx(temp_night, na.rm = FALSE),
         secchi_depth_day_inp = na.approx(secchi_depth_day, na.rm = FALSE),
         secchi_depth_night_inp = na.approx(secchi_depth_night, na.rm = FALSE),
         chinook0_hatchery_perhour_night_inp = na.approx(chinook0_hatchery_perhour_night, na.rm = FALSE),
         chinook0_hatchery_perhour_day_inp = na.approx(chinook0_hatchery_perhour_day, na.rm = FALSE),
         photoperiod_day = photoperiod,
         photoperiod_night = photoperiod,
         lunar_phase_day = lunar_phase,
         lunar_phase_night = lunar_phase)


#making columns for flow diff, photoperiod diff
data <- data %>%
  mutate(flow_diff_day = c(NA,diff(flow_day_inp)),
         flow_diff_night = c(NA,diff(flow_night_inp)),
         photo_diff_day = c(NA,diff(photoperiod_day)),
         photo_diff_night = c(NA,diff(photoperiod_night)),
         temp_diff_day = c(NA,diff(temp_day_inp)),
         temp_diff_night = c(NA,diff(temp_night_inp))
  )


lm_temp_day <- lm(temp_day ~ 1+photoperiod_day, data = data)
summary(lm_temp_day)
lm_temp_night <- lm(temp_night ~ 1+photoperiod_night, data = data)

data <- data %>% add_residuals(lm_temp_day, var = "resid_day") %>% 
  add_residuals(lm_temp_night, var = "resid_night")


data %>%
  dplyr::select(flow_day,photoperiod_day, lunar_phase_day, secchi_depth_day_inp, 
                photo_diff_day,flow_diff_day, temp_day, temp_diff_day, resid_day) %>%
  GGally::ggpairs(upper=list(continuous='points'),
                  lower=list(continuous='cor'))
#temp is correlated with secchi depth, photoperiod and photo diff


data$Date = as.Date(data$Date, format = "%m/%d/%Y")
data$year = year(data$Date)





# I think the doy limits should be changed from 200 to 218 
covariates_chinook0_puyallup <- arrange(data,doy) %>%
  filter(doy >130 & doy <= 218) %>%
  dplyr::select(year,doy, flow_day_inp, flow_night_inp, 
                photoperiod_day, photoperiod_night, secchi_depth_day_inp, secchi_depth_night_inp,
                lunar_phase_day, lunar_phase_night, flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
                temp_day_inp, temp_night_inp, temp_diff_day, temp_diff_night, 
                resid_day, resid_night,
                chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_day_inp, flow_night_inp, secchi_depth_day_inp, secchi_depth_night_inp,
    photoperiod_day, photoperiod_night, lunar_phase_day, 
    lunar_phase_night, flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
    temp_day_inp, temp_night_inp, temp_diff_day, temp_diff_night,
    resid_day, resid_night,
    chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

