# Goal 

# try the MARSS analysis of best model with data for years 2016-2022
#that has only the corrected data

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

data <- read.csv(here("puyallup", "data","puyallup_2004-2021_all_days_w_covariates.csv"),
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
         photo_diff_night = c(NA,diff(photoperiod_night))
  )

#convert date$Date to as.Date and then get year
#first check where date string is not in the right format
#specify the format and then convert to as.Date


data$Date = as.Date(data$Date)#, format = "%m/%d/%Y")
data$year = year(data$Date)

# I think the doy limits should be changed from 200 to 218 
covariates_chinook0 <- arrange(data,doy) %>%
  filter(doy >130 & doy <= 218 & year > 2015) %>%
  dplyr::select(year,doy, flow_day_inp, flow_night_inp, 
                photoperiod_day, photoperiod_night, secchi_depth_day_inp, secchi_depth_night_inp,
                lunar_phase_day, lunar_phase_night, flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
                chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_day_inp, flow_night_inp, secchi_depth_day_inp, secchi_depth_night_inp,
    photoperiod_day, photoperiod_night, lunar_phase_day, 
    lunar_phase_night, flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
    chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2021-2016+1
num_rows = num_years*2
num_covariates = 14
total_covariates = dim(covariates_chinook0)[1]


for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_chinook0[i,] = scale(covariates_chinook0[i,])[,1]
}

#just scale

for(i in (total_covariates/2 + 1):(total_covariates)){
  covariates_chinook0[i,] = scale(covariates_chinook0[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_chinook_summer_perhour <- arrange(data,doy) %>%
  filter(doy > 130 & doy <= 218 & year > 2015) %>%
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


#function for C matrix
#######


Cmat <- function(nyears,ncov,hatchery=0, day_on_night = FALSE){
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
                if(day_on_night){
                  C[i+nyears/2,j] <- "day_on_night"
                }
                
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

#function for Q matrix
#######


Qmat <- function(nyears){
  Q <- matrix(list(0),nrow = nyears,ncol = nyears, byrow = TRUE)
  for(i in 1:nyears){
    for(j in 1:nyears){
      if(i==j){
        if(i <= nyears/2){
          Q[i,j] <- "q_d"
        }
        else{
          Q[i,j] <- "q_n"
        }
      }
    }
  }
  return(Q)
}


#######

#function for mod list

######

mod_list <- function(nyears,ncov,hatchery=0, day_on_night = FALSE, unequal_q = FALSE){
  
  if(unequal_q){
    Q = Qmat(nyears)
    
  }
  else{
    Q = "diagonal and equal"
  }
  
  if(ncov == 0){
    mod.list = list(
      B = "identity",
      U = "zero",
      Z = "identity",
      A = "zero",
      R = "diagonal and equal",
      Q = Q
    )
  }
  else{
    if((ncov == 1) & (hatchery == 0)){
      C = "diagonal and equal"
    }
    else{
      C = Cmat(nyears,ncov,hatchery, day_on_night)
    }
    mod.list = list(
      B = "identity",
      U = "zero",
      Z = "identity",
      A = "zero",
      R = "diagonal and equal",
      Q = Q,
      C = C
    )
  }
  
  return(mod.list)
}

######

#function to get covariate combinations
##### 
get_covariate_combinations <- function(covariates) {
  n <- length(covariates)
  combinations <- lapply(1:n, function(x) combn(covariates, x, simplify = FALSE))
  unlist(combinations, recursive = FALSE)
}
#####

#first using the functions to see which of the correlated covariates best fit
#the data


#function for C matrix
#######


Cmat <- function(nyears,ncov,hatchery=0, day_on_night = FALSE){
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
                if(day_on_night){
                  C[i+nyears/2,j] <- "day_on_night"
                }
                
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

#function for Q matrix
#######


Qmat <- function(nyears){
  Q <- matrix(list(0),nrow = nyears,ncol = nyears, byrow = TRUE)
  for(i in 1:nyears){
    for(j in 1:nyears){
      if(i==j){
        if(i <= nyears/2){
          Q[i,j] <- "q_d"
        }
        else{
          Q[i,j] <- "q_n"
        }
      }
    }
  }
  return(Q)
}


#######

#function for mod list

######

mod_list <- function(nyears,ncov,hatchery=0, day_on_night = FALSE, unequal_q = FALSE){
  
  if(unequal_q){
    Q = Qmat(nyears)
    
  }
  else{
    Q = "diagonal and equal"
  }
  
  if(ncov == 0){
    mod.list = list(
      B = "identity",
      U = "zero",
      Z = "identity",
      A = "zero",
      R = "diagonal and equal",
      Q = Q
    )
  }
  else{
    if((ncov == 1) & (hatchery == 0)){
      C = "diagonal and equal"
    }
    else{
      C = Cmat(nyears,ncov,hatchery, day_on_night)
    }
    mod.list = list(
      B = "identity",
      U = "zero",
      Z = "identity",
      A = "zero",
      R = "diagonal and equal",
      Q = Q,
      C = C
    )
  }
  
  return(mod.list)
}

######

#function to get covariate combinations
##### 
get_covariate_combinations <- function(covariates) {
  n <- length(covariates)
  combinations <- lapply(1:n, function(x) combn(covariates, x, simplify = FALSE))
  unlist(combinations, recursive = FALSE)
}
#####

#first using the functions to see which of the correlated covariates best fit
#the data


#only looking at best model with flow, photoperiod difference, flow difference, and hatchery




c <- NULL
name_individual <- NULL
for(kk in c(1,5,6,7)){
  c = rbind(c,covariates_chinook0[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_chinook0)[1+(kk-1)*num_rows]
  name_individual = paste(name_individual,substr(name_long,1,nchar(name_long)-9))
  print(name_individual)
  
 
}

fit.model = c(list(c= c), mod_list(num_rows,4,1))
fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))


