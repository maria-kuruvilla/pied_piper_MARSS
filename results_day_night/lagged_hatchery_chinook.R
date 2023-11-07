#goal - to try hatchery covariate with different lags

#load packages

library(here)
library(MARSS)
library(ggplot2)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)
library(GGally)
library(broom)
library(tidyverse)


if(.Platform$OS.type == "unix") {
  data_string = here("..","..","data","pied_piper","dungeness")
} else {
  data_string = here("..","..","..","OneDrive","Documents","data","pied_piper","dungeness")
}


#load data
data_day <- read.csv(here("data","dungeness_unaggregated_day_w_covariates.csv"))
data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

#lag for day
data_day$chinook0_hatchery_lag1 <- lag(data_day$chinook0_hatchery_perhour_interpolate, 1)
data_day$chinook0_hatchery_lag2 <- lag(data_day$chinook0_hatchery_perhour_interpolate, 2)
data_day$chinook0_hatchery_lag3 <- lag(data_day$chinook0_hatchery_perhour_interpolate, 3)

#lag for night
data_night$chinook0_hatchery_lag1 <- lag(data_night$chinook0_hatchery_perhour_interpolate, 1)
data_night$chinook0_hatchery_lag2 <- lag(data_night$chinook0_hatchery_perhour_interpolate, 2)
data_night$chinook0_hatchery_lag3 <- lag(data_night$chinook0_hatchery_perhour_interpolate, 3)


data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"

data_day_night <- rbind(data_day,data_night)

covariates_chinook0_lag1 <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy >130 & doy <= 200) %>%
  dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_april, 
                lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
                chinook0_hatchery_lag1) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = c(
    temp, flow, photoperiod, atu_april, lunar_phase, 
    resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_lag1)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


num_years = 2020-2005
num_rows = num_years*2

dim(covariates_chinook0_lag1)



for(i in 1:(dim(covariates_chinook0_lag1)[1]/2)){ # everything except resid, diffs and hatchery
  covariates_chinook0_lag1[i,] = scale(covariates_chinook0_lag1[i,])[,1]
}

#skip 2005 day
for(i in (dim(covariates_chinook0_lag1)[1]/2 +1):(dim(covariates_chinook0_lag1)[1] - num_rows)){
  covariates_chinook0_lag1[i,] = scale(covariates_chinook0_lag1[i,], center = FALSE, scale= TRUE)[,1]
}

#skip 2005 night
for(i in (dim(covariates_chinook0_lag1)[1] - num_rows + 2):(dim(covariates_chinook0_lag1)[1] - num_rows/2)){
  covariates_chinook0_lag1[i,] = scale(covariates_chinook0_lag1[i,], center = FALSE, scale= TRUE)[,1]
}


for(i in (dim(covariates_chinook0_lag1)[1] - num_rows/2 + 2):(dim(covariates_chinook0_lag1)[1])){
  covariates_chinook0_lag1[i,] = scale(covariates_chinook0_lag1[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_chinook_summer_perhour <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy > 130 & doy <= 200) %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
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




mod.list_1_0 <- list(
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
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
  C = "diagonal and equal"
)

mod.list_2_0 <- list(
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
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
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
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
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
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
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
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
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

mod.list_6_0 <- list(
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
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
  C = matrix(c("a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),29),
               rep(list(0),1),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),28),
               rep(list(0),2),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),27),
               rep(list(0),3),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),26),
               rep(list(0),4),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),25),
               rep(list(0),5),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),24),
               rep(list(0),6),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),23),
               rep(list(0),7),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),22),
               rep(list(0),8),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),21),
               rep(list(0),9),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),20),
               rep(list(0),10),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),19),
               rep(list(0),11),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),18),
               rep(list(0),12),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),17),
               rep(list(0),13),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),16),
               rep(list(0),14),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),15),
               rep(list(0),15),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),14),
               rep(list(0),16),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),13),
               rep(list(0),17),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),12),
               rep(list(0),18),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),11),
               rep(list(0),19),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),10),
               rep(list(0),20),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),9),
               rep(list(0),21),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),8),
               rep(list(0),22),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),7),
               rep(list(0),23),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),6),
               rep(list(0),24),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),5),
               rep(list(0),25),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),4),
               rep(list(0),26),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),3),
               rep(list(0),27),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),2),
               rep(list(0),28),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f",rep(list(0),1),
               rep(list(0),29),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),29),"d",rep(list(0),29),"e",rep(list(0),29),"f"),30,180, byrow = TRUE)
)



mod.list_1_0_h <- list(
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
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
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
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
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
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
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
  Q =matrix(c("q_d",rep(list(0),29),
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
              rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
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


#adding effect of day on night d_on_n variable

mod.list_4_0_h_lag <- list(
  U = "zero",
  R = "diagonal and equal",
  Q =matrix(c("q_d",rep(list(0),29),
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
              rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
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
               rep(list(0),15),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),14),
               rep(list(0),16),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),13),
               rep(list(0),17),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),12),
               rep(list(0),18),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),11),
               rep(list(0),19),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),10),
               rep(list(0),20),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),9),
               rep(list(0),21),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),8),
               rep(list(0),22),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),7),
               rep(list(0),23),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),6),
               rep(list(0),24),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),5),
               rep(list(0),25),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),4),
               rep(list(0),26),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),3),
               rep(list(0),27),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),2),
               rep(list(0),28),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n",rep(list(0),1),
               rep(list(0),29),"a",rep(list(0),29),"b",rep(list(0),29),"c",rep(list(0),14),"d_on_n",rep(list(0),14),"n_on_n"),30,120,byrow=TRUE)
  
)

mod.list_5_0_h <- list(
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
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
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
               rep(list(0),29),"q_n"),nrow = 30, ncol = 30, byrow = TRUE),
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

#######




num_years = 15 
num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:6)
out.tab_photo_diff<- NULL
fits_photo_diff <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  photoperiod_difference = 0
  temperature_difference = 0
  flow = 0
  lunar_phase = 0
  hatchery = 0
  flow_difference = 0
  for(j in covariates){
    if(j == 1){
      k = 9
      photoperiod_difference =1
    }
    else if(j==2){
      k = 7
      temperature_difference = 1
    }
    else if(j==3){
      k = 2
      flow = 1
    }
    else if(j==4){
      k = 5
      lunar_phase = 1
    }
    else if(j==5){
      k = 8
      flow_difference = 1
    }
    else if(j==6){
      k = 10
      hatchery = 1
    }
    
    c = rbind(c,covariates_chinook0_lag1[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_chinook0_lag1)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-9))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  
  if(hatchery == 1){
    
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
    else if(c_num == 6){
      
      fit.model = c(list(c= c), mod.list_6_0_h)
    }
    
  }
  else{
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
    else if(c_num == 6){
      
      fit.model = c(list(c= c), mod.list_6_0)
    }
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, photoperiod_difference = photoperiod_difference,
                 temperature_difference = temperature_difference, flow = flow,
                 flow_difference = flow_difference,
                 lunar_phase = lunar_phase, hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_photo_diff=rbind(out.tab_photo_diff,out)
  fits_photo_diff=c(fits_photo_diff,list(fit))
  
}

out.tab_photo_diff


fit <- MARSS(subset_chinook_summer_perhour, model=mod.list_0_0, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
photoperiod_difference = 0
temperature_difference = 0
flow = 0
lunar_phase = 0
hatchery = 0
flow_difference = 0
out=data.frame(c=name, photoperiod_difference = photoperiod_difference,
               temperature_difference = temperature_difference, flow = flow,
               flow_difference = flow_difference,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)


out.tab_photo_diff=rbind(out.tab_photo_diff,out)
fits_photo_diff=c(fits_photo_diff,list(fit))

out.tab_photo_diff$deltaAICc <- NULL
out.tab_photo_diff$rel.LL <- NULL
out.tab_photo_diff$weights <- NULL


weights <- akaike.weights(out.tab_photo_diff$AICc)

out.tab_photo_diff$deltaAICc <- weights$deltaAIC
out.tab_photo_diff$rel.LL <- weights$rel.LL
out.tab_photo_diff$weights <- weights$weights


min.AICc <- order(out.tab_photo_diff$AICc)
out.tab_photo_diff_lag1.ordered <- out.tab_photo_diff[min.AICc, ]
out.tab_photo_diff_lag1.ordered


#so hatchery lag 1 is not a variable in the best model

# but maybe what I need to do is to take the best model from the no lag 

#then add various lag to it and look at the estimates


#the best model had flow dff, photo diff, temp diff, and hatchery

#this will use mod.list_4_0_h




num_years = 15 
num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:6)[[47]]
out <- NULL
fits_photo_diff_lag1 <- list()
covariate_number <- length(list_combinations)
# print(paste("covariate number",covariate_number))
covariates <- list_combinations
# print(paste("covariates",covariates))
c = NULL
name = NULL
photoperiod_difference = 0
temperature_difference = 0
flow = 0
lunar_phase = 0
hatchery = 0
flow_difference = 0
for(j in covariates){
  if(j == 1){
    k = 9
    photoperiod_difference =1
  }
  else if(j==2){
    k = 7
    temperature_difference = 1
  }
  else if(j==3){
    k = 2
    flow = 1
  }
  else if(j==4){
    k = 5
    lunar_phase = 1
  }
  else if(j==5){
    k = 8
    flow_difference = 1
  }
  else if(j==6){
    k = 10
    hatchery = 1
  }
  
  c = rbind(c,covariates_chinook0_lag1[((1+(k-1)*num_rows):(k*num_rows)),])
  name_long = rownames(covariates_chinook0_lag1)[1+(k-1)*num_rows]
  name = paste(name, substr(name_long,1,nchar(name_long)-9))
  
}
# print(c)
print(paste("name",name))
c_num <- length(covariates)
print(paste("c_num",c_num))
print(rownames(c))
fit.model = c(list(c= c), mod.list_4_0_h)

fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

lag1_ci <- tidy(fit)

out_lag1=data.frame(lag=1, day_effect = lag1_ci$estimate[37],day_effect_std = lag1_ci$std.error[37],
               day_effect_low = lag1_ci$conf.low[37],day_effect_high = lag1_ci$conf.up[37],
               night_effect = lag1_ci$estimate[39],night_effect_std = lag1_ci$std.error[39],
               night_effect_low = lag1_ci$conf.low[39],night_effect_high = lag1_ci$conf.up[39],
               night_effect_on_day = lag1_ci$estimate[38],
               night_effect_on_day_std = lag1_ci$std.error[38],
               night_effect_on_day_low = lag1_ci$conf.low[38],
               night_effect_on_day_high = lag1_ci$conf.up[38])

out = rbind(out,out_lag1)



#writing for loop to make dataframe with lag n and run marss analysis
out=NULL
for(n in 1:3){
  c = NULL
  cov_name = paste0("chinook0_hatchery_lag",n)
  covariates_chinook0_lag <- arrange(data_day_night,doy) %>%
    filter(year != 2015 & doy >130 & doy <= 200) %>%
    dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_april, 
                  lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
                  cov_name) %>%
    pivot_wider(names_from = c(year, daytime_category), values_from = c(
      temp, flow, photoperiod, atu_april, lunar_phase, 
      resid, temp_diff, flow_diff, photo_diff, cov_name)) %>%
    column_to_rownames(var = "doy") %>%
    as.matrix() %>%
    t()
  
  
  
  for(i in 1:(dim(covariates_chinook0_lag)[1]/2)){ # everything except resid, diffs and hatchery
    covariates_chinook0_lag[i,] = scale(covariates_chinook0_lag[i,])[,1]
  }
  
  #skip 2005 day
  for(i in (dim(covariates_chinook0_lag)[1]/2 +1):(dim(covariates_chinook0_lag)[1] - num_rows)){
    covariates_chinook0_lag[i,] = scale(covariates_chinook0_lag[i,], center = FALSE, scale= TRUE)[,1]
  }
  
  #skip 2005 night
  for(i in (dim(covariates_chinook0_lag)[1] - num_rows + 2):(dim(covariates_chinook0_lag)[1] - num_rows/2)){
    covariates_chinook0_lag[i,] = scale(covariates_chinook0_lag[i,], center = FALSE, scale= TRUE)[,1]
  }
  
  
  for(i in (dim(covariates_chinook0_lag)[1] - num_rows/2 + 2):(dim(covariates_chinook0_lag)[1])){
    covariates_chinook0_lag[i,] = scale(covariates_chinook0_lag[i,], center = FALSE, scale= TRUE)[,1]
  }
  
  for(j in covariates){
    if(j == 1){
      k = 9
      photoperiod_difference =1
    }
    else if(j==2){
      k = 7
      temperature_difference = 1
    }
    else if(j==3){
      k = 2
      flow = 1
    }
    else if(j==4){
      k = 5
      lunar_phase = 1
    }
    else if(j==5){
      k = 8
      flow_difference = 1
    }
    else if(j==6){
      k = 10
      hatchery = 1
    }
    
    c = rbind(c,covariates_chinook0_lag[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_chinook0_lag)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-9))
    
  }
  # print(c)
  print(paste("name",name))
  c_num <- length(covariates)
  print(paste("c_num",c_num))
  # print(rownames(c))
  fit.model = c(list(c= c), mod.list_4_0_h_lag)
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  lag_ci <- tidy(fit)
  
  out_lag =data.frame(lag=n, day_effect = lag_ci$estimate[37],day_effect_std = lag_ci$std.error[37],
                      day_effect_low = lag_ci$conf.low[37],day_effect_high = lag_ci$conf.up[37],
                      night_effect = lag_ci$estimate[40],night_effect_std = lag_ci$std.error[40],
                      night_effect_low = lag_ci$conf.low[40],night_effect_high = lag_ci$conf.up[40],
                      night_effect_on_day = lag_ci$estimate[39],
                      night_effect_on_day_std = lag_ci$std.error[39],
                      night_effect_on_day_low = lag_ci$conf.low[39],
                      night_effect_on_day_high = lag_ci$conf.up[39],
                      day_effect_on_night = lag_ci$estimate[38],
                      day_effect_on_night_std = lag_ci$std.error[38],
                      day_effect_on_night_low = lag_ci$conf.low[38],
                      day_effect_on_night_high = lag_ci$conf.up[38])
                      
  
  out = rbind(out,out_lag)
  
}
out

#for lag = 0
name = NULL
c = NULL

covariates_chinook0 <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy >130 & doy <= 200) %>%
  dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_april, 
                lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
                chinook0_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = c(
    temp, flow, photoperiod, atu_april, lunar_phase, 
    resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

num_years = 2020-2005
num_rows = num_years*2

dim(covariates_chinook0)

#scale the covariates

for(i in 1:(dim(covariates_chinook0)[1]/2)){ # everything except resid, diffs and hatchery
  covariates_chinook0[i,] = scale(covariates_chinook0[i,])[,1]
}

#skip 2005 day
for(i in (dim(covariates_chinook0)[1]/2 +1):(dim(covariates_chinook0)[1] - num_rows)){
  covariates_chinook0[i,] = scale(covariates_chinook0[i,], center = FALSE, scale= TRUE)[,1]
}

#skip 2005 night
for(i in (dim(covariates_chinook0)[1] - num_rows + 2):(dim(covariates_chinook0)[1] - num_rows/2)){
  covariates_chinook0[i,] = scale(covariates_chinook0[i,], center = FALSE, scale= TRUE)[,1]
}


for(i in (dim(covariates_chinook0)[1] - num_rows/2 + 2):(dim(covariates_chinook0)[1])){
  covariates_chinook0[i,] = scale(covariates_chinook0[i,], center = FALSE, scale= TRUE)[,1]
}



for(j in covariates){
  if(j == 1){
    k = 9
    photoperiod_difference =1
  }
  else if(j==2){
    k = 7
    temperature_difference = 1
  }
  else if(j==3){
    k = 2
    flow = 1
  }
  else if(j==4){
    k = 5
    lunar_phase = 1
  }
  else if(j==5){
    k = 8
    flow_difference = 1
  }
  else if(j==6){
    k = 10
    hatchery = 1
  }
  
  c = rbind(c,covariates_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
  name_long = rownames(covariates_chinook0)[1+(k-1)*num_rows]
  name = paste(name, substr(name_long,1,nchar(name_long)-9))
  
}
# print(c)
print(paste("name",name))
c_num <- length(covariates)
print(paste("c_num",c_num))
# print(rownames(c))
fit.model = c(list(c= c), mod.list_4_0_h)

fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

lag_ci <- tidy(fit)

out_lag =data.frame(lag=0, day_effect = lag_ci$estimate[37],day_effect_std = lag_ci$std.error[37],
                    day_effect_low = lag_ci$conf.low[37],day_effect_high = lag_ci$conf.up[37],
                    night_effect = lag_ci$estimate[39],night_effect_std = lag_ci$std.error[39],
                    night_effect_low = lag_ci$conf.low[39],night_effect_high = lag_ci$conf.up[39],
                    night_effect_on_day = lag_ci$estimate[38],
                    night_effect_on_day_std = lag_ci$std.error[38],
                    night_effect_on_day_low = lag_ci$conf.low[38],
                    night_effect_on_day_high = lag_ci$conf.up[38],
                    day_effect_on_night = NA,
                    day_effect_on_night_std = NA,
                    day_effect_on_night_low = NA,
                    day_effect_on_night_high = NA)


out = rbind(out,out_lag)
out

#plot the lagged effects of hatchery

ggplot(out, 
       aes(x = lag,
           y = day_effect, ymin = day_effect_low, ymax = day_effect_high)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Lag", y = "Effect of hatchery during day") +
  theme(axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24))

ggsave(here("output","lagged_day_effect.jpeg"), width = 10, height = 8)


ggplot(out, 
       aes(x = lag,
           y = night_effect, ymin = night_effect_low, ymax = night_effect_high)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Lag", y = "Effect of hatchery at night") +
  theme(axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24))

ggsave(here("output","lagged_night_effect.jpeg"), width = 10, height = 8)


ggplot(out, 
       aes(x = lag,
           y = night_effect_on_day, ymin = night_effect_on_day_low, ymax = night_effect_on_day_high)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Lag", y = "Effect of night hatchery during day") +
  theme(axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24))

ggsave(here("output","lagged_night_on_day_effect.jpeg"), width = 10, height = 8)


ggplot(out, 
       aes(x = lag,
           y = day_effect_on_night, ymin = day_effect_on_night_low, ymax = day_effect_on_night_high)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Lag", y = "Effect of day hatchery at night") +
  theme(axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24))

ggsave(here("output","lagged_day_on_night_effect.jpeg"), width = 10, height = 8)

out

#make dataframe with lagged effects in order
#day
#0 is the effect of day on day with lag 0
#0.5 is effect of night on day with lag 0
#1 is effect of day with lag 1
#1.5 is the effect of night on day with lag 1
#2 is the effect of day with lag 2
#2.5 is the effect of night on day with lag 2
#3 is the effect of day with lag 3
#3.5 is the effect of night on day with lag 3

out_day <- data.frame(lag = c(0,0.5,1,1.5,2,2.5,3,3.5), estimate = c(out$day_effect[out$lag==0],
                                                                     out$night_effect_on_day[out$lag==0],
                                                                     out$day_effect[out$lag==1],
                                                                     out$night_effect_on_day[out$lag==1],
                                                                     out$day_effect[out$lag==2],
                                                                     out$night_effect_on_day[out$lag==2],
                                                                     out$day_effect[out$lag==3],
                                                                     out$night_effect_on_day[out$lag==3]
                                                                     ),
                      low=c(out$day_effect_low[out$lag==0],
                            out$night_effect_on_day_low[out$lag==0],
                            out$day_effect_low[out$lag==1],
                            out$night_effect_on_day_low[out$lag==1],
                            out$day_effect_low[out$lag==2],
                            out$night_effect_on_day_low[out$lag==2],
                            out$day_effect_low[out$lag==3],
                            out$night_effect_on_day_low[out$lag==3]),
                      high=c(out$day_effect_high[out$lag==0],
                            out$night_effect_on_day_high[out$lag==0],
                            out$day_effect_high[out$lag==1],
                            out$night_effect_on_day_high[out$lag==1],
                            out$day_effect_high[out$lag==2],
                            out$night_effect_on_day_high[out$lag==2],
                            out$day_effect_high[out$lag==3],
                            out$night_effect_on_day_high[out$lag==3]))


ggplot(out_day, 
       aes(x = lag,
           y = estimate, ymin = low, ymax = high)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Lag", y = "Effect of hatchery during day") +
  theme(axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24))

ggsave(here("output","lagged_day_effect_combined.jpeg"), width = 10, height = 8)


out_night <- data.frame(lag = c(0,0.5,1,1.5,2,2.5,3), estimate = c(out$night_effect[out$lag==0],
                                                                     out$day_effect_on_night[out$lag==1],
                                                                     out$night_effect[out$lag==1],
                                                                     out$day_effect_on_night[out$lag==2],
                                                                     out$night_effect[out$lag==2],
                                                                     out$day_effect_on_night[out$lag==3],
                                                                     out$night_effect[out$lag==3]
                                                                     ),
                      low=c(out$night_effect_low[out$lag==0],
                            out$day_effect_on_night_low[out$lag==1],
                            out$night_effect_low[out$lag==1],
                            out$day_effect_on_night_low[out$lag==2],
                            out$night_effect_low[out$lag==2],
                            out$day_effect_on_night_low[out$lag==3],
                            out$night_effect_low[out$lag==3]),
                      high=c(out$night_effect_high[out$lag==0],
                            out$day_effect_on_night_high[out$lag==1],
                            out$night_effect_high[out$lag==1],
                            out$day_effect_on_night_high[out$lag==2],
                            out$night_effect_high[out$lag==2],
                            out$day_effect_on_night_high[out$lag==3],
                            out$night_effect_high[out$lag==3]))

out_night

ggplot(out_night, 
       aes(x = lag,
           y = estimate, ymin = low, ymax = high)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Lag", y = "Effect of hatchery at night") +
  theme(axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24))

ggsave(here("output","lagged_night_effect_combined.jpeg"), width = 10, height = 8)
