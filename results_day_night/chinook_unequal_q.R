

library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)
library(GGally)

if(.Platform$OS.type == "unix") {
  data_string = here("..","..","data","pied_piper","dungeness")
} else {
  data_string = here("..","..","..","OneDrive","Documents","data","pied_piper","dungeness")
}

data_day <- read.csv(here("data","dungeness_unaggregated_day_w_covariates.csv"))
data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"

data_day_night <- rbind(data_day,data_night)

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

fit <- MARSS(subset_chinook_summer_perhour, model=mod.list_0_0, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
autoplot(fit)




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
    else if(j==6){
      k = 8
      flow_difference = 1
    }
    else if(j==5){
      k = 10
      hatchery = 1
    }
    
    c = rbind(c,covariates_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_chinook0)[1+(k-1)*num_rows]
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
out.tab_photo_diff.ordered <- out.tab_photo_diff[min.AICc, ]
out.tab_photo_diff.ordered

out.tab_photo_diff.ordered$cumulative_weights <- cumsum(out.tab_photo_diff.ordered$weights)

relative_importance_photo_diff <- sum(out.tab_photo_diff$weights[out.tab_photo_diff$photoperiod_difference==1])
relative_importance_temperature_difference <- sum(out.tab_photo_diff$weights[out.tab_photo_diff$temperature_difference==1])
relative_importance_flow <- sum(out.tab_photo_diff$weights[out.tab_photo_diff$flow==1])
relative_importance_lunar_phase <- sum(out.tab_photo_diff$weights[out.tab_photo_diff$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab_photo_diff$weights[out.tab_photo_diff$hatchery==1])
relative_importance_flow_diff <- sum(out.tab_photo_diff$weights[out.tab_photo_diff$flow_difference==1])

riv_photo_diff <- data.frame(variable = c("photoperiod difference",
                                          "temperature difference",
                                          "flow",
                                          "lunar phase",
                                          "flow difference",
                                          "hatchery"),
                             relative_importance = c(relative_importance_photo_diff,
                                                     relative_importance_temperature_difference,
                                                     relative_importance_flow,
                                                     relative_importance_lunar_phase,
                                                     relative_importance_flow_diff,
                                                     relative_importance_hatchery))


#sort according to relative importance
riv_photo_diff[order(riv_photo_diff$relative_importance, decreasing = TRUE),]

riv_photo_diff <- riv_photo_diff[order(riv_photo_diff$relative_importance, decreasing = TRUE),]
riv_photo_diff$relative_importance <- round(riv_photo_diff$relative_importance,1)

riv_photo_diff

autoplot(fits_photo_diff[[47]], plot.type = "fitted.ytT")
ggsave(here("output","fitted_y_chinook_best_model_day_night.jpeg"), width = 10, height = 8)


ci_best_photo_diff <- tidy(fits_photo_diff[[47]])

ggplot(ci_best_photo_diff[c(34:37,39),], 
       aes(x = c("Photoperiod\n difference", "Temperature\n difference", "Flow\n difference",
                 "Hatchery,\n day", "Hatchery,\n night"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  ggtitle("Chinook subyearlings") + 
  theme(plot.title = element_text(size = 24))+
  theme(axis.text.x=element_text(size=24),axis.title.y=element_text(size=24))

ggsave(here("output","chinook0_day_night_effect_final_photo_diff.jpeg"), width = 14, height = 12)

write.csv(riv_photo_diff, here("output","riv_photo_diff_chinook0_day_night.csv"))
