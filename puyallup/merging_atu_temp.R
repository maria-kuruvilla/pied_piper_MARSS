# Goal - to read the file with air temperature data and atu data
# read file with puyallup catch data
# merge the two data sets


library(ggplot2)
library(here)
library(GGally)
library(tidyverse)
library(zoo)
library(modelr)


#read in the air temperature data

atu_temp <- read_csv(here("puyallup","data", "puyallup_airtemp_atu.csv"))

#read in the puyallup catch data

puyallup_catch <- read_csv(here("puyallup","data", "puyallup_2004-2021_all_days_w_covariates_w_airtemp_coho_new.csv"))


#make new columns in atu_temp - doy and year, then only select those and temp, atu_solstice

atu_temp_subset <- atu_temp %>%
  mutate(doy = as.numeric(format(Date, format = "%j")),
         year = as.numeric(format(Date, format = "%Y"))) %>%
  select(doy, year, temp, atu_solstice)

#drop temp columns from puyallup_catch and them merge with atu_temp_subset on doy and year


puyallup_catch$chinook0_hatchery_perhour_night[dim(puyallup_catch)[1]] <- 0
puyallup_catch$chinook0_hatchery_perhour_day[dim(puyallup_catch)[1]] <- 0

#when secchi_depth_day is NA, replace with secchi_depth values
#when secchi_depth_night is NA, replace with secchi_depth values

puyallup_catch$secchi_depth_day[is.na(puyallup_catch$secchi_depth_day)] <- puyallup_catch$secchi_depth[is.na(puyallup_catch$secchi_depth_day)]

puyallup_catch$secchi_depth_night[is.na(puyallup_catch$secchi_depth_night)] <- puyallup_catch$secchi_depth[is.na(puyallup_catch$secchi_depth_night)]


puyallup_catch_edited <- puyallup_catch %>%
  select(-temp) %>%
  left_join(atu_temp_subset, by = c("doy", "year")) %>% 
  mutate(temp_day = temp, temp_night = temp,
         atu_day = atu_solstice, atu_night = atu_solstice,
         flow_day = na.approx(flow_day, na.rm = FALSE),
         flow_night = na.approx(flow_night, na.rm = FALSE),
         secchi_depth_day = na.approx(secchi_depth_day, na.rm = FALSE),
         secchi_depth_night = na.approx(secchi_depth_night, na.rm = FALSE),
         chinook0_hatchery_perhour_night = na.approx(chinook0_hatchery_perhour_night, na.rm = FALSE),
         chinook0_hatchery_perhour_day = na.approx(chinook0_hatchery_perhour_day, na.rm = FALSE),
         photoperiod_day = photoperiod,
         photoperiod_night = photoperiod,
         lunar_phase_day = lunar_phase,
         lunar_phase_night = lunar_phase,
         flow_diff_day = c(NA,diff(flow_day)),
         flow_diff_night = c(NA,diff(flow_night)),
         photo_diff_day = c(NA,diff(photoperiod_day)),
         photo_diff_night = c(NA,diff(photoperiod_night)),
         temp_diff_day = c(NA,diff(temp_day)),
         temp_diff_night = c(NA,diff(temp_night))
         )

lm_temp_day <- lm(temp_day ~ 1+photoperiod_day, data = puyallup_catch_edited)
summary(lm_temp_day)
lm_temp_night <- lm(temp_night ~ 1+photoperiod_night, data =puyallup_catch_edited)
summary(lm_temp_night)

puyallup_catch_edited <- puyallup_catch_edited %>%
  add_residuals(lm_temp_day, var = "resid_day") %>% 
  add_residuals(lm_temp_night, var = "resid_night")


#plotting all the data before saving

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = temp_day, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = temp_day, x = doy), alpha = 0.5) +
  facet_wrap(~year)

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = flow_day, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = flow_day, x = doy), alpha = 0.5) +
  facet_wrap(~year)

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = flow_diff_day, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = flow_diff_day, x = doy), alpha = 0.5) +
  facet_wrap(~year)

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = atu_day, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = atu_day, x = doy), alpha = 0.5) +
  facet_wrap(~year)

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = coho1_wild_num, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = coho1_wild_num, x = doy), alpha = 0.5) +
  facet_wrap(~year, scales = "free_y")

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = coho1_wild_num_day, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = coho1_wild_num_day, x = doy), alpha = 0.5) +
  facet_wrap(~year, scales = "free_y")

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = coho1_wild_num_night, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = coho1_wild_num_night, x = doy), alpha = 0.5) +
  facet_wrap(~year, scales = "free_y")

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = coho1_wild_perhour_night, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = coho1_wild_perhour_night, x = doy), alpha = 0.5) +
  facet_wrap(~year, scales = "free_y")

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = coho1_hatchery_num_night, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = coho1_hatchery_num_night, x = doy), alpha = 0.5) +
  facet_wrap(~year, scales = "free_y")

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = coho1_wild_perhour_estimate, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = coho1_wild_perhour_estimate, x = doy), alpha = 0.5) +
  facet_wrap(~year, scales = "free_y")

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = coho1_wild_num_night_estimate, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = coho1_wild_num_night_estimate, x = doy), alpha = 0.5) +
  facet_wrap(~year, scales = "free_y")

ggplot(puyallup_catch_edited) +
  geom_point(aes(y = trap_efficiency, x = doy), alpha = 0.5, size = 1) +
  geom_line(aes(y = trap_efficiency, x = doy), alpha = 0.5) +
  facet_wrap(~year, scales = "free_y")

#save puyallup_catch_edited

write.csv(puyallup_catch_edited, here("puyallup","data","puyallup_final.csv"))
