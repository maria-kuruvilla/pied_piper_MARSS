# goal - make prediction figure for coho1 from puyallup
#analysis


#load libraries
library(MARSS)
library(ggplot2)
library(here)
library(GGally)
library(tidyverse)
library(zoo)
library(modelr)
library(qpcR)
library(broom)


#load data
data <- read.csv(here("puyallup", "data","puyallup_final.csv"),
                 na.strings = c("NA",""))

data$coho1_hatchery_perhour_night[dim(data)[1]] <- 0
data$coho1_hatchery_perhour_day[dim(data)[1]] <- 0


data <- data %>% 
  mutate(coho1_hatchery_perhour_night = na.approx(coho1_hatchery_perhour_night, na.rm = FALSE),
         coho1_hatchery_perhour_day = na.approx(coho1_hatchery_perhour_day, na.rm = FALSE))



covariates_coho1_puyallup_w_temp <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160) %>%
  dplyr::select(year,doy, flow_day, flow_night, 
                photoperiod_day, photoperiod_night, secchi_depth_day, secchi_depth_night,
                lunar_phase_day, lunar_phase_night, temp_day, temp_night, atu_day, atu_night,
                flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
                temp_diff_day, temp_diff_night, 
                resid_day, resid_night,
                coho1_hatchery_perhour_day, coho1_hatchery_perhour_night) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_day, flow_night, secchi_depth_day, secchi_depth_night,
    photoperiod_day, photoperiod_night, lunar_phase_day, 
    lunar_phase_night, temp_day, temp_night, atu_day, atu_night,
    flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
    temp_diff_day, temp_diff_night,
    resid_day, resid_night,
    coho1_hatchery_perhour_day, coho1_hatchery_perhour_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2021-2004+1
num_rows = num_years*2
total_covariates = dim(covariates_coho1_puyallup_w_temp)[1]


for(i in 1:(num_rows*6)){ # everything except diffs and hatchery
  covariates_coho1_puyallup_w_temp[i,] = scale(covariates_coho1_puyallup_w_temp[i,])[,1]
}

#just scale

for(i in (num_rows*6 + 1):(total_covariates)){
  covariates_coho1_puyallup_w_temp[i,] = scale(covariates_coho1_puyallup_w_temp[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_coho_summer_perhour <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160) %>%
  mutate(log.value_day = log(coho1_wild_perhour_day + 1), 
         log.value_night = log(coho1_wild_perhour_night + 1)) %>%
  dplyr::select(log.value_day, log.value_night ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_day, log.value_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour)[1]){
  subset_coho_summer_perhour[i,] = scale(subset_coho_summer_perhour[i,])[,1]
}


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

c<-NULL
for(kk in c(1,3,6,7,11)){
  c = rbind(c,covariates_coho1_puyallup_w_temp[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_coho1_puyallup_w_temp)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  
}
fit.model_check = c(list(c= c), mod_list(nyears,5,1, FALSE, FALSE))
fit_check <- MARSS(subset_coho_summer_perhour, model=fit.model_check, silent = TRUE, method = "BFGS",
                   control=list(maxit=2000))

autoplot(fit_check, plot.type = "fitted.ytT")

ggsave(here("puyallup","output","coho1_puyallup_model_fit.png"), width = 12, height = 10, units = "in")

autoplot(fit_check)
