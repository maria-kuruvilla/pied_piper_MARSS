# Goal 

# try MARSS analysis with Puyallup data fo chinook0

#load packages

library(here)
library(MARSS)
library(ggplot2)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)
library(tidyverse)
library(GGally)


# load data

data <- read.csv(here("puyallup", "data","puyallup_2004-2018_all_days_w_covariates.csv"),
                 na.strings = c("NA",""))

#plotting chinook0_wild_perhour with doy to see timing of migration

ggplot(data, aes(x=doy, y=chinook0_wild_perhour_night)) +
  geom_point(alpha = 0.5) +
  labs(x = "Day of Year", y = "Chinook0 Wild per Hour") +
  theme_bw() + 
  ylim(0, 20)

#plotting chinook0_hatchery_perhour with doy to see timing of migration

ggplot(data, aes(x=doy, y=chinook0_hatchery_perhour_night)) +
  geom_point(alpha = 0.5) +
  labs(x = "Day of Year", y = "Chinook0 hatchery per Hour") +
  theme_bw()  + 
  ylim(0, 500)

#before interpolating, I need to add 0 to the end of 2018 chinook hatchery column
#if not, there are nans in the 2018 column even after interpolation

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
  mutate(flow_day_diff = c(NA,diff(flow_day_inp)),
         flow_night_diff = c(NA,diff(flow_night_inp)),
         photo_diff_day = c(NA,diff(photoperiod_day)),
         photo_diff_night = c(NA,diff(photoperiod_night))
         )

covariates_chinook0 <- arrange(data,doy) %>%
  filter(doy >130 & doy <= 200) %>%
  dplyr::select(year,doy, flow_day_inp, flow_night_inp, 
                photoperiod_day, photoperiod_night, secchi_depth_day_inp, secchi_depth_night_inp,
                lunar_phase_day, lunar_phase_night, flow_day_diff, flow_night_diff, photo_diff_day, photo_diff_night, 
                chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_day_inp, flow_night_inp, secchi_depth_day_inp, secchi_depth_night_inp,
    photoperiod_day, photoperiod_night, lunar_phase_day, 
    lunar_phase_night, flow_day_diff, flow_night_diff, photo_diff_day, photo_diff_night, 
    chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

#scaling the variables

num_years = 2018-2004+1
num_rows = num_years*2
num_covariates = 14
total_covariates = dim(covariates_chinook0)[1]

#scale and center the covariates

for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_chinook0[i,] = scale(covariates_chinook0[i,])[,1]
}

#just scale

for(i in (total_covariates/2 + 1):(total_covariates)){
  covariates_chinook0[i,] = scale(covariates_chinook0[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_chinook_summer_perhour <- arrange(data,doy) %>%
  filter(doy > 130 & doy <= 200) %>%
  mutate(log.value_day = log(chinook0_wild_perhour_day + 1), 
         log.value_night = log(chinook0_wild_perhour_night + 1)) %>%
  dplyr::select(log.value_day, log.value_night ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_day, log.value_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_chinook_summer_perhour)[1]){
  subset_chinook_summer_perhour[i,] = scale(subset_chinook_summer_perhour[i,])[,1]
}



#model list
#####
mod.list_0_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal"

)

mod.list_0_0_unequalq <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = matrix(c("q_d",rep(list(0),29),
               rep(list(0),1),"q_d",rep(list(0),28),
               rep(list(0),2),"q_d",rep(list(0),27),
               rep(list(0),3),"q_d",rep(list(0),26),
               rep(list(0),4),"q_d",rep(list(0),25),
               rep(list(0),5),"q_d",rep(list(0),24),
               rep(list(0),6),"q_d",rep(list(0),23),
               rep(list(0),7),"q_d",rep(list(0),22),
               rep(list(0),8),"q_d",rep(list(0),21),
               rep(list(0),9),"q_d",rep(list(0),20),
               rep(list(0),10),"q_d",rep(list(0),19),
               rep(list(0),11),"q_d",rep(list(0),18),
               rep(list(0),12),"q_d",rep(list(0),17),
               rep(list(0),13),"q_d",rep(list(0),16),
               rep(list(0),14),"q_d",rep(list(0),15),
               rep(list(0),15),"q_n",rep(list(0),14),
               rep(list(0),16),"q_n",rep(list(0),13),
               rep(list(0),17),"q_n",rep(list(0),12),
               rep(list(0),18),"q_n",rep(list(0),11),
               rep(list(0),19),"q_n",rep(list(0),10),
               rep(list(0),20),"q_n",rep(list(0),9),
               rep(list(0),21),"q_n",rep(list(0),8),
               rep(list(0),22),"q_n",rep(list(0),7),
               rep(list(0),23),"q_n",rep(list(0),6),
               rep(list(0),24),"q_n",rep(list(0),5),
               rep(list(0),25),"q_n",rep(list(0),4),
               rep(list(0),26),"q_n",rep(list(0),3),
               rep(list(0),27),"q_n",rep(list(0),2),
               rep(list(0),28),"q_n",rep(list(0),1),
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE)
  
)

mod.list_0_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal"
)

mod.list_1_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = "diagonal and equal"
)

mod.list_2_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(c("a",rep(list(0),29),"b",rep(list(0),29),
               0,"a",rep(list(0),29),"b",rep(list(0),28),
               0,0,"a",rep(list(0),29),"b",rep(list(0),27),
               0,0,0,"a",rep(list(0),29),"b",rep(list(0),26),
               0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),25),
               0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),24),
               0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),23),
               0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),22),
               0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),21),
               0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),20),
               0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),19),
               0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),18),
               0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),17),
               0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),16),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),15),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),14),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),13),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),12),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),11),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),10),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),9),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),8),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),7),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),6),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),5),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),4),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),3),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),2),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b",rep(list(0),1),
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",rep(list(0),29),"b"),
             
             30,60, byrow = TRUE)
  
)

mod.list_3_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(c("a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),
               rep(list(0),1),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),28),
               rep(list(0),2),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),27),
               rep(list(0),3),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),26),
               rep(list(0),4),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),25),
               rep(list(0),5),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),24),
               rep(list(0),6),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),23),
               rep(list(0),7),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),22),
               rep(list(0),8),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),21),
               rep(list(0),9),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),20),
               rep(list(0),10),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),19),
               rep(list(0),11),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),18),
               rep(list(0),12),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),17),
               rep(list(0),13),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),16),
               rep(list(0),14),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),15),
               rep(list(0),15),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),
               rep(list(0),16),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),13),
               rep(list(0),17),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),12),
               rep(list(0),18),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),11),
               rep(list(0),19),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),10),
               rep(list(0),20),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),9),
               rep(list(0),21),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),8),
               rep(list(0),22),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),7),
               rep(list(0),23),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),6),
               rep(list(0),24),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),5),
               rep(list(0),25),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),4),
               rep(list(0),26),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),3),
               rep(list(0),27),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),2),
               rep(list(0),28),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),1),
               rep(list(0),29),"a",rep(list(0),29),"b",rep(list(0),29),"c"),30,90, byrow = TRUE)
  
  
)


mod.list_4_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(c("a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),
               rep(list(0),1),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),28),
               rep(list(0),2),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),27),
               rep(list(0),3),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),26),
               rep(list(0),4),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),25),
               rep(list(0),5),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),24),
               rep(list(0),6),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),23),
               rep(list(0),7),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),22),
               rep(list(0),8),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),21),
               rep(list(0),9),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),20),
               rep(list(0),10),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),19),
               rep(list(0),11),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),18),
               rep(list(0),12),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),17),
               rep(list(0),13),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),16),
               rep(list(0),14),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),15),
               rep(list(0),15),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),14),
               rep(list(0),16),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),13),
               rep(list(0),17),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),12),
               rep(list(0),18),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),11),
               rep(list(0),19),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),10),
               rep(list(0),20),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),9),
               rep(list(0),21),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),8),
               rep(list(0),22),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),7),
               rep(list(0),23),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),6),
               rep(list(0),24),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),5),
               rep(list(0),25),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),4),
               rep(list(0),26),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),3),
               rep(list(0),27),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),2),
               rep(list(0),28),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),1),
               rep(list(0),29),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d"), 
             30, 120, byrow = TRUE)
)


mod.list_5_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(c("a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),
               rep(list(0),1),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),28),
               rep(list(0),2),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),27),
               rep(list(0),3),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),26),
               rep(list(0),4),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),25),
               rep(list(0),5),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),24),
               rep(list(0),6),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),23),
               rep(list(0),7),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),22),
               rep(list(0),8),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),21),
               rep(list(0),9),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),20),
               rep(list(0),10),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),19),
               rep(list(0),11),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),18),
               rep(list(0),12),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),17),
               rep(list(0),13),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),16),
               rep(list(0),14),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),15),
               rep(list(0),15),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),14),
               rep(list(0),16),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),13),
               rep(list(0),17),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),12),
               rep(list(0),18),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),11),
               rep(list(0),19),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),10),
               rep(list(0),20),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),9),
               rep(list(0),21),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),8),
               rep(list(0),22),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),7),
               rep(list(0),23),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),6),
               rep(list(0),24),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),5),
               rep(list(0),25),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),4),
               rep(list(0),26),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),3),
               rep(list(0),27),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),2),
               rep(list(0),28),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),1),
               rep(list(0),29),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e"),30,150, byrow = TRUE)
)

####

# let's make mod lists for hatchery
#########
mod.list_1_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(c("d_on_d",rep(list(0),14),"n_on_d",rep(list(0),14),
               rep(list(0),1),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),13),
               rep(list(0),2),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),12),
               rep(list(0),3),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),11),
               rep(list(0),4),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),10),
               rep(list(0),5),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),9),
               rep(list(0),6),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),8),
               rep(list(0),7),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),7),
               rep(list(0),8),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),6),
               rep(list(0),9),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),5),
               rep(list(0),10),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),4),
               rep(list(0),11),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),3),
               rep(list(0),12),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),2),
               rep(list(0),13),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),1),
               rep(list(0),14),"d_on_d",rep(list(0),14),"n_on_d",
               rep(list(0),15),"n_on_n",rep(list(0),14),
               rep(list(0),16),"n_on_n",rep(list(0),13),
               rep(list(0),17),"n_on_n",rep(list(0),12),
               rep(list(0),18),"n_on_n",rep(list(0),11),
               rep(list(0),19),"n_on_n",rep(list(0),10),
               rep(list(0),20),"n_on_n",rep(list(0),9),
               rep(list(0),21),"n_on_n",rep(list(0),8),
               rep(list(0),22),"n_on_n",rep(list(0),7),
               rep(list(0),23),"n_on_n",rep(list(0),6),
               rep(list(0),24),"n_on_n",rep(list(0),5),
               rep(list(0),25),"n_on_n",rep(list(0),4),
               rep(list(0),26),"n_on_n",rep(list(0),3),
               rep(list(0),27),"n_on_n",rep(list(0),2),
               rep(list(0),28),"n_on_n",rep(list(0),1),
               rep(list(0),29),"n_on_n"),30,30,byrow=TRUE)
  
  
)



mod.list_2_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(c("a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),14),
               rep(list(0),1),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),13),
               rep(list(0),2),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),12),
               rep(list(0),3),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),11),
               rep(list(0),4),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),10),
               rep(list(0),5),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),9),
               rep(list(0),6),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),8),
               rep(list(0),7),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),7),
               rep(list(0),8),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),6),
               rep(list(0),9),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),5),
               rep(list(0),10),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),4),
               rep(list(0),11),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),3),
               rep(list(0),12),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),2),
               rep(list(0),13),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),1),
               rep(list(0),14),"a",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",
               rep(list(0),15),"a",rep(list(0),29),"n_on_n",rep(list(0),14),
               rep(list(0),16),"a",rep(list(0),29),"n_on_n",rep(list(0),13),
               rep(list(0),17),"a",rep(list(0),29),"n_on_n",rep(list(0),12),
               rep(list(0),18),"a",rep(list(0),29),"n_on_n",rep(list(0),11),
               rep(list(0),19),"a",rep(list(0),29),"n_on_n",rep(list(0),10),
               rep(list(0),20),"a",rep(list(0),29),"n_on_n",rep(list(0),9),
               rep(list(0),21),"a",rep(list(0),29),"n_on_n",rep(list(0),8),
               rep(list(0),22),"a",rep(list(0),29),"n_on_n",rep(list(0),7),
               rep(list(0),23),"a",rep(list(0),29),"n_on_n",rep(list(0),6),
               rep(list(0),24),"a",rep(list(0),29),"n_on_n",rep(list(0),5),
               rep(list(0),25),"a",rep(list(0),29),"n_on_n",rep(list(0),4),
               rep(list(0),26),"a",rep(list(0),29),"n_on_n",rep(list(0),3),
               rep(list(0),27),"a",rep(list(0),29),"n_on_n",rep(list(0),2),
               rep(list(0),28),"a",rep(list(0),29),"n_on_n",rep(list(0),1),
               rep(list(0),29),"a",rep(list(0),29),"n_on_n"),30,60,byrow=TRUE)
  
)


mod.list_3_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(c("a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),14),
               rep(list(0),1),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),13),
               rep(list(0),2),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),12),
               rep(list(0),3),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),11),
               rep(list(0),4),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),10),
               rep(list(0),5),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),9),
               rep(list(0),6),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),8),
               rep(list(0),7),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),7),
               rep(list(0),8),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),6),
               rep(list(0),9),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),5),
               rep(list(0),10),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),4),
               rep(list(0),11),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),3),
               rep(list(0),12),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),2),
               rep(list(0),13),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),1),
               rep(list(0),14),"a",rep(list(0),29),"b",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",
               rep(list(0),15),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),14),
               rep(list(0),16),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),13),
               rep(list(0),17),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),12),
               rep(list(0),18),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),11),
               rep(list(0),19),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),10),
               rep(list(0),20),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),9),
               rep(list(0),21),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),8),
               rep(list(0),22),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),7),
               rep(list(0),23),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),6),
               rep(list(0),24),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),5),
               rep(list(0),25),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),4),
               rep(list(0),26),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),3),
               rep(list(0),27),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),2),
               rep(list(0),28),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n",rep(list(0),1),
               rep(list(0),29),"a",rep(list(0),29),"b",rep(list(0),29),"n_on_n"),30,90,byrow=TRUE)
  
)

mod.list_4_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(c("a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),14),
               rep(list(0),1),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),13),
               rep(list(0),2),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),12),
               rep(list(0),3),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),11),
               rep(list(0),4),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),10),
               rep(list(0),5),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),9),
               rep(list(0),6),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),8),
               rep(list(0),7),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),7),
               rep(list(0),8),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),6),
               rep(list(0),9),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),5),
               rep(list(0),10),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),4),
               rep(list(0),11),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),3),
               rep(list(0),12),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),2),
               rep(list(0),13),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),1),
               rep(list(0),14),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",
               rep(list(0),15),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),14),
               rep(list(0),16),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),13),
               rep(list(0),17),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),12),
               rep(list(0),18),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),11),
               rep(list(0),19),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),10),
               rep(list(0),20),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),9),
               rep(list(0),21),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),8),
               rep(list(0),22),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),7),
               rep(list(0),23),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),6),
               rep(list(0),24),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),5),
               rep(list(0),25),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),4),
               rep(list(0),26),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),3),
               rep(list(0),27),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),2),
               rep(list(0),28),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n",rep(list(0),1),
               rep(list(0),29),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"n_on_n"),30,120,byrow=TRUE)
  
)

mod.list_5_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(c("a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),14),
               rep(list(0),1),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),13),
               rep(list(0),2),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),12),
               rep(list(0),3),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),11),
               rep(list(0),4),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),10),
               rep(list(0),5),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),9),
               rep(list(0),6),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),8),
               rep(list(0),7),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),7),
               rep(list(0),8),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),6),
               rep(list(0),9),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),5),
               rep(list(0),10),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),4),
               rep(list(0),11),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),3),
               rep(list(0),12),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),2),
               rep(list(0),13),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),1),
               rep(list(0),14),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",
               rep(list(0),15),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),14),
               rep(list(0),16),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),13),
               rep(list(0),17),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),12),
               rep(list(0),18),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),11),
               rep(list(0),19),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),10),
               rep(list(0),20),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),9),
               rep(list(0),21),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),8),
               rep(list(0),22),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),7),
               rep(list(0),23),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),6),
               rep(list(0),24),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),5),
               rep(list(0),25),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),4),
               rep(list(0),26),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),3),
               rep(list(0),27),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),2),
               rep(list(0),28),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n",rep(list(0),1),
               rep(list(0),29),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"n_on_n"),30,150,byrow=TRUE)
  
)

mod.list_6_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(c("a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),14),
               rep(list(0),1),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),13),
               rep(list(0),2),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),12),
               rep(list(0),3),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),11),
               rep(list(0),4),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),10),
               rep(list(0),5),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),9),
               rep(list(0),6),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),8),
               rep(list(0),7),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),7),
               rep(list(0),8),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),6),
               rep(list(0),9),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),5),
               rep(list(0),10),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),4),
               rep(list(0),11),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),3),
               rep(list(0),12),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),2),
               rep(list(0),13),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",rep(list(0),1),
               rep(list(0),14),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"d_on_d",rep(list(0),14),"n_on_d",
               rep(list(0),15),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),14),
               rep(list(0),16),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),13),
               rep(list(0),17),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),12),
               rep(list(0),18),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),11),
               rep(list(0),19),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),10),
               rep(list(0),20),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),9),
               rep(list(0),21),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),8),
               rep(list(0),22),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),7),
               rep(list(0),23),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),6),
               rep(list(0),24),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),5),
               rep(list(0),25),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),4),
               rep(list(0),26),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),3),
               rep(list(0),27),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),2),
               rep(list(0),28),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n",rep(list(0),1),
               rep(list(0),29),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"n_on_n"),30,180,byrow=TRUE)
  
)

#####

#function for C matrix
#######


Cmat <- function(nyears,ncov,hatchery=0){
  vars = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")
  
  if(hatchery == 1){
    C <- matrix(list(0),nrow = nyears,ncol = ncov*nyears, byrow = TRUE)
    for(i in 1:nyears){
      for(j in 1:(ncov*nyears)){
        
        for(k in 1:ncov){
          # print(((k-1)*nyears+i))
          if(j == ((k-1)*nyears+i)){
            # print(i)
            # print(j)
            # print(vars[k])
            if(k==ncov){
              if(i<=nyears/2){
                C[i,j] <- "day"
              }
              else{
                C[i,j] <- "night"
              }
            }
            else{
              C[i,j] <- vars[k]
            }
            
            
          }
          
        }
      }
      
    }
    
  }
  else{
    C <- matrix(list(0),nrow = nyears,ncol = ncov*nyears, byrow = TRUE)
    for(i in 1:nyears){
      for(j in 1:(ncov*nyears)){
        
        for(k in 1:ncov){
          # print(((k-1)*nyears+i))
          if(j == ((k-1)*nyears+i)){
            # print(i)
            # print(j)
            # print(vars[k])
            C[i,j] <- vars[k]
            
          }
          
        }
      }
      
    }
  }
  
  return(C)
}


#######

#function for mod list

######

mod_list <- function(nyears,ncov,hatchery=0){
  
  if(ncov == 0){
    mod.list = list(
      B = "identity",
      U = "zero",
      Z = "identity",
      A = "zero",
      R = "diagonal and equal",
      Q = "diagonal and equal"
    )
  }
  else{
    if((ncov == 1) & (hatchery == 0)){
      C = "diagonal and equal"
    }
    else{
      C = Cmat(nyears,ncov,hatchery)
    }
    mod.list = list(
      B = "identity",
      U = "zero",
      Z = "identity",
      A = "zero",
      R = "diagonal and equal",
      Q = "diagonal and equal",
      C = C
    )
  }
  
  return(mod.list)
}

######

fit <- MARSS(subset_chinook_summer_perhour, model=mod.list_0_0, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
autoplot(fit)
fit$AICc #4319.137


fit_unequalq <- MARSS(subset_chinook_summer_perhour, model=mod.list_0_0_unequalq, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
autoplot(fit_unequalq)
fit_unequalq$AICc #4319.239
#not that different from equal q. I will just go with equal q for now

#try MARSS analysis with just one covariate - hatchery
c = covariates_chinook0[(12*15 + 1):(14*15),]
fit.model = c(list(c= c), mod.list_1_0_h)
fit_hatchery <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                      control=list(maxit=2000))
fit_hatchery$AICc #4318.47

ci_hatchery <- tidy(fit_hatchery)
ci_hatchery
#not a big effect of hatchery, but need to add other covariates

#look at correlation between all the covariates with one plot

#looking at correlated covariates
data %>%
  dplyr::select(flow_day,photoperiod_day, lunar_phase_day, secchi_depth_day_inp, 
                photo_diff_day,flow_day_diff) %>%
  GGally::ggpairs(upper=list(continuous='points'),
                  lower=list(continuous='cor'))
#photoperiod, secchi_depth and photoperiod_diff are correlated

get_covariate_combinations <- function(covariates) {
  n <- length(covariates)
  combinations <- lapply(1:n, function(x) combn(covariates, x, simplify = FALSE))
  unlist(combinations, recursive = FALSE)
}


list_combinations <- get_covariate_combinations(1:5)
out.tab.hatchery<- NULL
fits.hatchery <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  if(covariates[1] == 1){
    
    for(kk in c(2,3,6)){
      c = NULL
      name = NULL
      for(j in covariates){
        if(j == 1){
          k = kk
        }
        else if(j==2){
          k = 1
        }
        else if(j==3){
          k = 4
        }
        else if(j==4){
          k = 5
        }
        else if(j==5){
          k = 7
        }
        
        c = rbind(c,covariates_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
        name_long = rownames(covariates_chinook0)[1+(k-1)*num_rows]
        name = paste(name, substr(name_long,1,nchar(name_long)-5))
        
      }
      # print(c)
      print(name)
      c_num <- length(covariates)
      if(k==7){
        has_hatchery = TRUE
        c_num <- length(covariates)
        if(c_num == 1){

          fit.model = c(list(c= c), mod.list_1_0_h)
        }

        else if(c_num == 2){

          fit.model = c(list(c= c), mod.list_2_0_h)
        }
        else if(c_num == 3){

          fit.model = c(list(c= c), mod.list_3_0_h)
        }
        else if(c_num == 4){

          fit.model = c(list(c= c), mod.list_4_0_h)
        }
        else if(c_num == 5){

          fit.model = c(list(c= c), mod.list_5_0_h)
        }
      }
      else{
        has_hatchery = FALSE
        c_num <- length(covariates)
        if(c_num == 1){

          fit.model = c(list(c= c), mod.list_1_0)
        }

        else if(c_num == 2){

          fit.model = c(list(c= c), mod.list_2_0)
        }
        else if(c_num == 3){

          fit.model = c(list(c= c), mod.list_3_0)
        }
        else if(c_num == 4){

          fit.model = c(list(c= c), mod.list_4_0)
        }
        else if(c_num == 5){

          fit.model = c(list(c= c), mod.list_5_0)
        }
      }

      fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                   control=list(maxit=2000))


      out=data.frame(c=name, d = "None",
                     logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                     num.iter=fit$numIter, converged=!fit$convergence,
                     stringsAsFactors = FALSE)
      out.tab.hatchery=rbind(out.tab.hatchery,out)
      fits.hatchery=c(fits.hatchery,list(fit))
    }
    
  }
  else{
    for(j in covariates){
      if(j == 1){
        k = 0
      }
      else if(j==2){
        k = 1
      }
      else if(j==3){
        k = 4
      }
      else if(j==4){
        k = 5
      }
      else if(j==5){
        k = 7
      }
      
      c = rbind(c,covariates_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
      name_long = rownames(covariates_chinook0)[1+(k-1)*num_rows]
      num_individual = substr(name_long,1,nchar(name_long)-9)
      name = paste(name, substr(name_long,1,nchar(name_long)-9))
      
    }
    # print(c)
    print(name)
    c_num <- length(covariates)
    if(k==7){
      has_hatchery = TRUE
      c_num <- length(covariates)
      if(c_num == 1){

        fit.model = c(list(c= c), mod.list_1_0_h)
      }

      else if(c_num == 2){

        fit.model = c(list(c= c), mod.list_2_0_h)
      }
      else if(c_num == 3){

        fit.model = c(list(c= c), mod.list_3_0_h)
      }
      else if(c_num == 4){

        fit.model = c(list(c= c), mod.list_4_0_h)
      }
      else if(c_num == 5){

        fit.model = c(list(c= c), mod.list_5_0_h)
      }
    }
    else{
      has_hatchery = FALSE
      c_num <- length(covariates)
      if(c_num == 1){

        fit.model = c(list(c= c), mod.list_1_0)
      }

      else if(c_num == 2){

        fit.model = c(list(c= c), mod.list_2_0)
      }
      else if(c_num == 3){

        fit.model = c(list(c= c), mod.list_3_0)
      }
      else if(c_num == 4){

        fit.model = c(list(c= c), mod.list_4_0)
      }
      else if(c_num == 5){

        fit.model = c(list(c= c), mod.list_5_0)
      }
    }


    fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                 control=list(maxit=2000))


    out=data.frame(c=name, d = "None",
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence,
                   stringsAsFactors = FALSE)
    out.tab.hatchery=rbind(out.tab.hatchery,out)
    fits.hatchery=c(fits.hatchery,list(fit))

  }
  
  
}
out.tab.hatchery$deltaAICc <- out.tab.hatchery$AICc - min(out.tab.hatchery$AICc)
min.AICc <- order(out.tab.hatchery$AICc)
out.tab.hatchery.ordered <- out.tab.hatchery[min.AICc, ]
out.tab.hatchery.ordered

fits.hatchery[[55]]
ci_best_model <- tidy(fits.hatchery[[55]])

#plot the effects

ggplot(ci_best_model[c(33:36,38),], 
       aes(x = c("Photoperiod", "Flow", "Flow\n difference",
                 "Hatchery,\n day", "Hatchery,\n night"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  theme(axis.text.x=element_text(size=20),axis.title.y=element_text(size=24))

#plotting the fit

autoplot(fits.hatchery[[55]])


#first using the functions to see which of the correlated covariates best fit
#the data
out.tab <- NULL
fits <- NULL
nyears = 30
for(kk in c(2,3,6)){
  c = covariates_chinook0[((1+(kk-1)*num_rows):(kk*num_rows)),]
  name_long = rownames(covariates_chinook0)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0))
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab=rbind(out.tab,out)
  fits=c(fits,list(fit))
}

out.tab
#photoperiod is best

#now using the functions to see which of the uncorrelated covariates best fit



list_combinations <- get_covariate_combinations(1:5)
out.tab.hatchery<- NULL
fits.hatchery <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  for(j in covariates){
        if(j == 1){
          k = 3
        }
        else if(j==2){
          k = 1
        }
        else if(j==3){
          k = 4
        }
        else if(j==4){
          k = 5
        }
        else if(j==5){
          k = 7
        }
        
        c = rbind(c,covariates_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
        name_long = rownames(covariates_chinook0)[1+(k-1)*num_rows]
        name = paste(name, substr(name_long,1,nchar(name_long)-5))
        
      }
      # print(c)
      print(name)
      c_num <- length(covariates)
      if(k==7){
        has_hatchery = 1
        c_num <- length(covariates)
        fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery))
      }
      else{
        has_hatchery = 0
        c_num <- length(covariates)
        fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery))
      }
      fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                   control=list(maxit=2000))
      
      
      out=data.frame(c=name, d = "None",
                     logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                     num.iter=fit$numIter, converged=!fit$convergence,
                     stringsAsFactors = FALSE)
      out.tab.hatchery=rbind(out.tab.hatchery,out)
      fits.hatchery=c(fits.hatchery,list(fit))
    
  
}
out.tab.hatchery$deltaAICc <- out.tab.hatchery$AICc - min(out.tab.hatchery$AICc)
min.AICc <- order(out.tab.hatchery$AICc)
out.tab.hatchery.ordered <- out.tab.hatchery[min.AICc, ]
out.tab.hatchery.ordered

fits.hatchery[[28]]
ci_best_model <- tidy(fits.hatchery[[28]])


#plot the effects

ggplot(ci_best_model[c(33:37),], 
       aes(x = c("Photoperiod", "Flow", "Flow\n difference",
                 "Hatchery,\n day", "Hatchery,\n night"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  theme(axis.text.x=element_text(size=20),axis.title.y=element_text(size=24))


#need to change to day on night effect

