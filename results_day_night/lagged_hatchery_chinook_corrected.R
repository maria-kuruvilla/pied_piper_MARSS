#Goal - to create a lagged hatchery variable for the chinook data with lag 0.5
#remove correlated hatchery covariates from same model

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

#load data
data_day <- read.csv(here("data","dungeness_unaggregated_day_w_covariates.csv"))
data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

#lag for day
data_day$chinook0_hatchery_lag1 <- lag(data_day$chinook0_hatchery_perhour_interpolate, 1)
data_day$chinook0_hatchery_lag0.5 <- data_night$chinook0_hatchery_perhour_interpolate
data_day$chinook0_hatchery_lag2 <- lag(data_day$chinook0_hatchery_perhour_interpolate, 2)

#lag for night
data_night$chinook0_hatchery_lag1 <- lag(data_night$chinook0_hatchery_perhour_interpolate, 1)
data_night$chinook0_hatchery_lag0.5 <- lag(data_day$chinook0_hatchery_perhour_interpolate, 1)
data_night$chinook0_hatchery_lag2 <- lag(data_night$chinook0_hatchery_perhour_interpolate, 2)

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"

data_day_night <- rbind(data_day,data_night)


covariates_chinook0_lag0.5 <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy >130 & doy <= 200) %>%
  dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_april, 
                lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
                chinook0_hatchery_lag0.5) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = c(
    temp, flow, photoperiod, atu_april, lunar_phase, 
    resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_lag0.5)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

num_years = 2020-2005
num_rows = num_years*2

dim(covariates_chinook0_lag0.5)



for(i in 1:(dim(covariates_chinook0_lag0.5)[1]/2)){ # everything except resid, diffs and hatchery
  covariates_chinook0_lag0.5[i,] = scale(covariates_chinook0_lag0.5[i,])[,1]
}

#skip 2005 day
for(i in (dim(covariates_chinook0_lag0.5)[1]/2 +1):(dim(covariates_chinook0_lag0.5)[1] - num_rows)){
  covariates_chinook0_lag0.5[i,] = scale(covariates_chinook0_lag0.5[i,], center = FALSE, scale= TRUE)[,1]
}

#skip 2005 night
for(i in (dim(covariates_chinook0_lag0.5)[1] - num_rows + 2):(dim(covariates_chinook0_lag0.5)[1] - num_rows/2)){
  covariates_chinook0_lag0.5[i,] = scale(covariates_chinook0_lag0.5[i,], center = FALSE, scale= TRUE)[,1]
}


for(i in (dim(covariates_chinook0_lag0.5)[1] - num_rows/2 + 2):(dim(covariates_chinook0_lag0.5)[1])){
  covariates_chinook0_lag0.5[i,] = scale(covariates_chinook0_lag0.5[i,], center = FALSE, scale= TRUE)[,1]
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


out.tab.lag1 <- NULL
fits.lag1 <- NULL
nyears = num_years*2
c=NULL
name = NULL
for(k in c(9,7,8,10)){
  c = rbind(c,covariates_chinook0_lag1[((1+(k-1)*nyears):(k*nyears)),])
  name_long = rownames(covariates_chinook0_lag1)[1+(k-1)*nyears]
  name = paste(name, substr(name_long,1,nchar(name_long)-9))
}
fit.model = c(list(c= c), mod_list(nyears,4,1,FALSE,TRUE))
fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
out=data.frame(c=name_individual, d = "None",
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)
fit
ci <- tidy(fit)
ci


#writing for loop to make dataframe with lag n and run marss analysis
out_lag=NULL
nyears = num_years*2
list_combinations <- get_covariate_combinations(1:6)[[47]]
out <- NULL
covariate_number <- length(list_combinations)
covariates <- list_combinations
for(n in c(0.5,1,2)){
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
  for(i in (dim(covariates_chinook0_lag)[1]/2 +1):(dim(covariates_chinook0_lag)[1] - nyears)){
    covariates_chinook0_lag[i,] = scale(covariates_chinook0_lag[i,], center = FALSE, scale= TRUE)[,1]
  }
  
  #skip 2005 night
  for(i in (dim(covariates_chinook0_lag)[1] - nyears + 2):(dim(covariates_chinook0_lag)[1] - nyears/2)){
    covariates_chinook0_lag[i,] = scale(covariates_chinook0_lag[i,], center = FALSE, scale= TRUE)[,1]
  }
  
  
  for(i in (dim(covariates_chinook0_lag)[1] - nyears/2 + 2):(dim(covariates_chinook0_lag)[1])){
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
    
    c = rbind(c,covariates_chinook0_lag[((1+(k-1)*nyears):(k*nyears)),])
    name_long = rownames(covariates_chinook0_lag)[1+(k-1)*nyears]
    name = paste(name, substr(name_long,1,nchar(name_long)-9))
    
  }
  # print(c)
  print(paste("name",name))
  c_num <- length(covariates)
  print(paste("c_num",c_num))
  # print(rownames(c))
  fit.model = c(list(c= c), mod_list(nyears,c_num,1,FALSE,TRUE))
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  lag_ci <- tidy(fit)
  
  out_lag =data.frame(lag=n, day_effect = lag_ci$estimate[37],day_effect_std = lag_ci$std.error[37],
                      day_effect_low = lag_ci$conf.low[37],day_effect_high = lag_ci$conf.up[37],
                      night_effect = lag_ci$estimate[38],night_effect_std = lag_ci$std.error[38],
                      night_effect_low = lag_ci$conf.low[38],night_effect_high = lag_ci$conf.up[38]
                      )
  
  
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
fit.model = c(list(c= c), mod_list(num_rows,c_num,1,FALSE,TRUE))

fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

lag_ci <- tidy(fit)

out_lag =data.frame(lag=0, day_effect = lag_ci$estimate[37],day_effect_std = lag_ci$std.error[37],
                    day_effect_low = lag_ci$conf.low[37],day_effect_high = lag_ci$conf.up[37],
                    night_effect = lag_ci$estimate[38],night_effect_std = lag_ci$std.error[38],
                    night_effect_low = lag_ci$conf.low[38],night_effect_high = lag_ci$conf.up[38]
                    )


out = rbind(out,out_lag)
out



ggplot(out, 
       aes(x = lag,
           y = day_effect, ymin = day_effect_low, ymax = day_effect_high)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Lag", y = "Effect of hatchery during day") +
  theme(axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24))

ggsave(here("output","lagged_day_effect_corrected.jpeg"), width = 10, height = 8)


ggplot(out, 
       aes(x = lag,
           y = night_effect, ymin = night_effect_low, ymax = night_effect_high)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Lag", y = "Effect of hatchery at night") +
  theme(axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24))

ggsave(here("output","lagged_night_effect_corrected.jpeg"), width = 10, height = 8)

