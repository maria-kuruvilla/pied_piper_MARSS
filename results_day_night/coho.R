# Goal

# to repeat the MARSS day ngiht analysis for coho

#load packages
library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)


#check system
if(.Platform$OS.type == "unix") {
  data_string = here("..","..","data","pied_piper","dungeness")
} else {
  data_string = here("..","..","..","OneDrive","Documents","data","pied_piper","dungeness")
}


#load data

data_day <- read.csv(here("data","dungeness_unaggregated_day_w_covariates.csv"))
data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"

#interpolate the hatchery covariate values
data_day$coho1_hatchery_perhour_interpolate <- na.approx(data_day$coho1_hatchery_perhour, na.rm = FALSE)
data_night$coho1_hatchery_perhour_interpolate <- na.approx(data_night$coho1_hatchery_perhour, na.rm = FALSE)


data_day_night <- rbind(data_day,data_night)

covariates_coho1 <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy >120 & doy <= 160) %>%
  dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_april, 
                lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
                coho1_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = c(
    temp, flow, photoperiod, atu_april, lunar_phase, 
    resid, temp_diff, flow_diff, photo_diff, coho1_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

num_years = 2020-2005
num_rows = num_years*2


dim(covariates_coho1)

#scale the covariates

for(i in 1:(dim(covariates_coho1)[1]/2)){ # everything except resid, diffs and hatchery
  covariates_coho1[i,] = scale(covariates_coho1[i,])[,1]
}

#only scale the rest, do not center
for(i in (dim(covariates_coho1)[1]/2 +1):(dim(covariates_coho1)[1])){
  covariates_coho1[i,] = scale(covariates_coho1[i,], center = FALSE, scale= TRUE)[,1]
}


#subset response variable
subset_coho_summer_perhour <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy > 120 & doy <= 160) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour)[1]){
  subset_coho_summer_perhour[i,] = scale(subset_coho_summer_perhour[i,])[,1]
}

#model list
#####
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


######

#trial 
#####
fit.model = mod.list_0_0
fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
autoplot(fit)

#trying if mod.list_2_0 with rep in the matrix works

c = covariates_coho1[1:60,]
fit.model = c(list(c= c), mod.list_2_0)
fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
autoplot(fit)

######

#running marss in loop

#making combination of covariates

get_covariate_combinations <- function(covariates) {
  n <- length(covariates)
  combinations <- lapply(1:n, function(x) combn(covariates, x, simplify = FALSE))
  unlist(combinations, recursive = FALSE)
}

num_years = 15
num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:5)
out.tab.coho <- NULL
fits.coho <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  if(covariates[1] == 1){
    print("looking at temp type variables")
    for(kk in c(1,3,4,6,9)){
      c = NULL
      name = NULL
      for(j in covariates){
        if(j == 1){
          k = kk
        }
        else if(j==2){
          k = 7
        }
        else if(j==3){
          k = 2
        }
        else if(j==4){
          k = 5
        }
        else if(j==5){
          k = 10
        }
        
        c = rbind(c,covariates_coho1[((1+(k-1)*num_rows):(k*num_rows)),])
        name_long = rownames(covariates_coho1)[1+(k-1)*num_rows]
        name = paste(name, substr(name_long,1,nchar(name_long)-9))
        
      }
      # print(c)
      print(name)
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
      
      print(paste("fitting marss with ", name))
      fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                   control=list(maxit=2000))
      
      
      out=data.frame(c=name, d = "None",
                     logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                     num.iter=fit$numIter, converged=!fit$convergence,
                     stringsAsFactors = FALSE)
      out.tab.coho=rbind(out.tab.coho,out)
      fits.coho=c(fits.coho,list(fit))
    }
  }
  
  else{
    for(j in covariates){
      if(j == 1){
        k = 1
      }
      else if(j==2){
        k = 7
      }
      else if(j==3){
        k = 2
      }
      else if(j==4){
        k = 5
      }
      else if(j==5){
        k = 10
      }
      
      c = rbind(c,covariates_coho1[((1+(k-1)*num_rows):(k*num_rows)),])
      name_long = rownames(covariates_coho1)[1+(k-1)*num_rows]
      num_individual = substr(name_long,1,nchar(name_long)-9)
      name = paste(name, substr(name_long,1,nchar(name_long)-9))
      
    }
    # print(c)
    print(name)
    c_num <- length(covariates)
    
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
    
    
    print(paste("fitting marss with ", name))
    fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                 control=list(maxit=2000))
    
    
    out=data.frame(c=name, d = "None",
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence,
                   stringsAsFactors = FALSE)
    out.tab.coho=rbind(out.tab.coho,out)
    fits.coho=c(fits.coho,list(fit))
  }
}



out.tab.coho$deltaAICc <- out.tab.coho$AICc - min(out.tab.coho$AICc)
min.AICc <- order(out.tab.coho$AICc)
out.tab.coho.ordered <- out.tab.coho[min.AICc, ]
out.tab.coho.ordered

fits.coho[[16]]
tidy(fits.coho[[16]])
autoplot(fits.coho[[16]])

###models with hatchery

######

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



########


# run marss in loop with hatchery variables in mod lists

########
num_years = 15
num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:5)
out.tab.coho.hatchery<- NULL
fits.coho.hatchery <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  if(covariates[1] == 1){
    
    for(kk in c(1,3,4,6,9)){
      c = NULL
      name = NULL
      for(j in covariates){
        if(j == 1){
          k = kk
        }
        else if(j==2){
          k = 7
        }
        else if(j==3){
          k = 2
        }
        else if(j==4){
          k = 5
        }
        else if(j==5){
          k = 10
        }
        
        c = rbind(c,covariates_coho1[((1+(k-1)*num_rows):(k*num_rows)),])
        name_long = rownames(covariates_coho1)[1+(k-1)*num_rows]
        name = paste(name, substr(name_long,1,nchar(name_long)-9))
        
      }
      # print(c)
      print(name)
      c_num <- length(covariates)
      if(k==10){
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
      
      fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                   control=list(maxit=2000))
      
      
      out=data.frame(c=name, d = "None",
                     logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                     num.iter=fit$numIter, converged=!fit$convergence,
                     stringsAsFactors = FALSE)
      out.tab.coho.hatchery=rbind(out.tab.coho.hatchery,out)
      fits.coho.hatchery=c(fits.coho.hatchery,list(fit))
    }
    
  }
  else{
    for(j in covariates){
      if(j == 1){
        k = 1
      }
      else if(j==2){
        k = 7
      }
      else if(j==3){
        k = 2
      }
      else if(j==4){
        k = 5
      }
      else if(j==5){
        k = 10
      }
      
      c = rbind(c,covariates_coho1[((1+(k-1)*num_rows):(k*num_rows)),])
      name_long = rownames(covariates_coho1)[1+(k-1)*num_rows]
      num_individual = substr(name_long,1,nchar(name_long)-9)
      name = paste(name, substr(name_long,1,nchar(name_long)-9))
      
    }
    # print(c)
    print(name)
    c_num <- length(covariates)
    if(k==10){
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
    
    
    fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                 control=list(maxit=2000))
    
    
    out=data.frame(c=name, d = "None",
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence,
                   stringsAsFactors = FALSE)
    out.tab.coho.hatchery=rbind(out.tab.coho.hatchery,out)
    fits.coho.hatchery=c(fits.coho.hatchery,list(fit))
    
  }
  
  
}
out.tab.coho.hatchery$deltaAICc <- out.tab.coho.hatchery$AICc - min(out.tab.coho.hatchery$AICc)
min.AICc <- order(out.tab.coho.hatchery$AICc)
out.tab.coho.hatchery.ordered <- out.tab.coho.hatchery[min.AICc, ]
out.tab.coho.hatchery.ordered

######

fits.coho.hatchery[[57]]
ci_coho <- tidy(fits.coho.hatchery[[57]])
autoplot(fits.coho.hatchery[[57]])

fitted_coho <- fitted(fits.coho.hatchery[[57]])
fitted_coho$year <- substr(fitted_coho$.rownames,1,4)
