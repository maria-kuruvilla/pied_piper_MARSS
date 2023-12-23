#Goal - to write code to run MARSS analysis for chinook0 data from Dungeness,
#Puyallup and Skagit Rivers combined

#Load libraries

library(MARSS)
library(ggplot2)
library(here)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)
library(GGally)
library(tidyverse)
library(broom)

#load data
dungeness <- read.csv(here("data", "dungeness_subset_covariates_diff.csv"))

puyallup <- read.csv(here("puyallup","data", "puyallup_2004-2021_all_days_w_covariates_edited.csv"))

skagit <- read.csv(here("skagit","data", "skagit_2010-2022_w_covariates.csv"))


#make column for river in all data and then rbind all data

dungeness$river <- "dungeness"
puyallup$river <- "puyallup"
skagit$river <- "skagit"


#check unique secchi_depth values in puyallup
unique(puyallup$secchi_depth_night)

#In puyallup when secchi_depth is NA, average secchi_depth_day and secchi_depth_night
#unless both are NA, then secchi_depth is NA
#first comvert all '-' to NA
puyallup <- puyallup %>%
  mutate(secchi_depth_day = ifelse(secchi_depth_day == " - ", NA, as.numeric(secchi_depth_day)),
         secchi_depth_night = ifelse(secchi_depth_night == " - ", NA, as.numeric(secchi_depth_night))) %>%
  mutate(flag = ifelse(is.na(secchi_depth_day) & is.na(secchi_depth_night), 1, 0)) %>%
  mutate(secchi_depth = ifelse(is.na(secchi_depth) & flag == 0,
                               (secchi_depth_day + secchi_depth_night)/2, secchi_depth))

#combine data for all rivers into one data frame and have NA for columns in rivers
#that do not have that columns
#rename atu_april to atu
dungeness_subset <- dungeness %>% 
  filter(year != 2015) %>%
  select(c("Date","doy","year","chinook0_wild_perhour","chinook0_hatchery_perhour",
           "chinook0_wild_num","chinook0_hatchery_num",
           "chinook0_hatchery_perhour_interpolate","In",
           "temp", "flow","photoperiod","atu_april","lunar_phase",
           "resid","temp_diff","flow_diff","photo_diff","river")) %>% 
  mutate("secchi_depth" = NA) %>%
  rename("atu" = "atu_april")

#the last row of chinook0_hatchery_perhour in puyallup should be set to 0
puyallup$chinook0_hatchery_perhour[dim(puyallup)[1]] <- 0

puyallup_subset <- puyallup %>% 
  select(c("Date","doy","year","chinook0_wild_perhour","chinook0_hatchery_perhour",
           "chinook0_wild_num","chinook0_hatchery_num","In",
           "flow","photoperiod","lunar_phase","secchi_depth","river")) %>% 
  mutate("temp" = NA, "atu" = NA, "resid" = NA, "temp_diff" = NA,
         "flow_diff" = NA, "photo_diff" = NA, 
         "chinook0_hatchery_perhour_interpolate" = na.approx(chinook0_hatchery_perhour, na.rm = FALSE),
         "secchi_depth" = na.approx(secchi_depth, na.rm = FALSE))


skagit_subset <- skagit %>%
  select(c("Date","doy","year","chinook0_wild_num","chinook0_hatchery_num",
           "temp", "flow","photoperiod","atu_solstice","lunar_phase",
           "temp_diff","flow_diff","photo_diff","river", "In", "trap", "daytime_category")) %>% 
  
  rename("atu" = "atu_solstice") %>% 
  group_by(Date,doy,year) %>%
  summarise("In" = sum(In, na.rm = TRUE),
            "chinook0_wild_num" = if(In==0) NA else sum(chinook0_wild_num, na.rm = TRUE),
            "chinook0_hatchery_num" = if(In==0) NA else sum(chinook0_hatchery_num, na.rm = TRUE),
            "temp" = mean(temp, na.rm = TRUE),
            "flow" = mean(flow, na.rm = TRUE),
            "photoperiod" = mean(photoperiod, na.rm = TRUE),
            "atu" = mean(atu, na.rm = TRUE),
            "lunar_phase" = mean(lunar_phase, na.rm = TRUE),
            "temp_diff" = mean(temp_diff, na.rm = TRUE),
            "flow_diff" = mean(flow_diff, na.rm = TRUE),
            "photo_diff" = mean(photo_diff, na.rm = TRUE),
            "river" = first(river)) %>%
  ungroup() %>%
  mutate("secchi_depth" = NA, "chinook0_wild_perhour" = chinook0_wild_num/In,
         "chinook0_hatchery_perhour" = chinook0_hatchery_num/In)
         # "chinook0_hatchery_perhour_interpolate" = na.approx(chinook0_hatchery_perhour, na.rm = FALSE))
  
#set last line of chinook0_hatchery_perhour in skagit to 0
skagit_subset$chinook0_hatchery_perhour[dim(skagit_subset)[1]] <- 0

#calculate chinook0_hatchery_perhour_interpolate in skagit_subset

skagit_subset$chinook0_hatchery_perhour_interpolate <- na.approx(skagit_subset$chinook0_hatchery_perhour, na.rm = FALSE)

#add residuals to skagit subset
lm_temp_day <- lm(temp ~ 1+photoperiod, data = skagit_subset)

skagit_subset <- skagit_subset %>% add_residuals(lm_temp_day)

all_rivers <- rbind(dungeness_subset, puyallup_subset, skagit_subset)

#subset response variable
chinook0_wild <- arrange(all_rivers,doy) %>%
  filter(doy > 130 & doy <= 200) %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,river) %>%
  pivot_wider(names_from = c(year,river), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(chinook0_wild)[1]){
  chinook0_wild[i,] = scale(chinook0_wild[i,])[,1]
}

covariates_chinook0 <- arrange(all_rivers,doy) %>%
  filter(doy >130 & doy <= 200) %>%
  dplyr::select(year,doy, river, secchi_depth, temp, flow, photoperiod, atu,
                lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
                chinook0_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = c(year, river), values_from = c(secchi_depth,
    temp, flow, photoperiod, atu, lunar_phase, 
    resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

#scale and center all rows till 270
#only scale all rows from 271

for(i in 1:(45*6)){
  covariates_chinook0[i,] = scale(covariates_chinook0[i,])[,1]
}

#skip row 451
for(i in ((45*6)+1):dim(covariates_chinook0)[1]){
  if(i == 451){
    next
  }
  covariates_chinook0[i,] = scale(covariates_chinook0[i,], center = FALSE)[,1]
}

#check by taking average of each row

for(i in 1:dim(covariates_chinook0)[1]){
  print(rownames(covariates_chinook0)[i])
  print(mean(covariates_chinook0[i,],na.rm = TRUE))
}





Cmat_all_rivers <- function(nyears,ncov){
  vars = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")
  
  
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
  
  
  return(C)
}



mod_list <- function(nyears,ncov){
  
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
    if((ncov == 1)){
      C = "diagonal and equal"
    }
    else{
      C = Cmat_all_rivers(nyears,ncov)
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


#let's try MARSS analysis with no covariates

nyears = 45

fit <- MARSS(chinook0_wild, model=mod_list(nyears,0,0), silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
fit$convergence
autoplot(fit)

#check 2018 puyallup
#yeah lots of missing data

#let's try MARSS analysis with just hatchery

nyears = 45
for(kk in c(11)){
  c = covariates_chinook0[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_chinook0)[1+(kk-1)*nyears]
  print(name_long)
  
}

fit.model = c(list(c= c), mod_list(nyears,1))
fit_hatchery <- MARSS(chinook0_wild, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

#there seems to be either a NA or inf in a row in c. check which row

for(i in 1:dim(c)[1]){
  print(i)
  print(any(is.na(c[i,])))
  print(any(is.infinite(c[i,])))
}

##### should I do this?

#replace NA with forward fill in each row of covariates_chinook0
for(i in 1:dim(covariates_chinook0)[1]){
  covariates_chinook0[i,] <- na.locf(covariates_chinook0[i,])
}
