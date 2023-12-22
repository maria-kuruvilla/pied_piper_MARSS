#goal - do model selection for chinook0 data
#calculate relative variable importance

#load libraries

library(MARSS)

library(here)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)
library(GGally)
library(ggplot2)
library(tidyverse)

#load data

data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))

data_day_night$chinook0_wild_perhour <- data_day_night$chinook0_wild_num/data_day_night$In
data_day_night$chinook0_hatchery_perhour <- data_day_night$chinook0_hatchery_num/data_day_night$In
data_day_night$coho1_wild_perhour <- data_day_night$coho1_wild_num/data_day_night$In
data_day_night$coho1_hatchery_perhour <- data_day_night$coho1_hatchery_num/data_day_night$In
data_day_night$chinook1_wild_perhour <- data_day_night$chinook1_wild_num/data_day_night$In
data_day_night$chinook1_hatchery_perhour <- data_day_night$chinook1_hatchery_num/data_day_night$In



data_day_night$chinook0_hatchery_perhour_inp <- na.approx(data_day_night$chinook0_hatchery_perhour, 
                                                          na.rm = FALSE)
data_day_night$coho1_hatchery_perhour_inp <- na.approx(data_day_night$coho1_hatchery_perhour, 
                                                       na.rm = FALSE)

lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day_night)

data_day_night <- data_day_night %>% add_residuals(lm_temp_day)


covariates_chinook0_skagit <- arrange(data_day_night,doy) %>%
  filter(doy >150 & doy < 190) %>%
  dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_solstice,
                lunar_phase, resid, temp_diff, flow_diff, photo_diff,
                chinook0_hatchery_perhour_inp, trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = c(
    temp, flow, photoperiod, atu_solstice, lunar_phase,
    resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_perhour_inp), names_vary = "fastest") %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

data_day_night %>%
  dplyr::select(flow,photoperiod, lunar_phase, temp, temp_diff, atu_solstice, resid, 
                photo_diff,flow_diff) %>%
  GGally::ggpairs(upper=list(continuous='points'),
                  lower=list(continuous='cor'))

#scaling the variables

num_years = 2022-2010
num_rows = num_years*2
num_covariates = 10
total_covariates = dim(covariates_chinook0_skagit)[1]

#scaling and centering
for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_chinook0_skagit[i,] = scale(covariates_chinook0_skagit[i,])[,1]
}


#checking
for(i in 1:(total_covariates)){
  print(sum(covariates_chinook0_skagit[i,]))
}

#just scale
for(i in (total_covariates/2 + 1):(total_covariates)){
  if(sum(covariates_chinook0_skagit[i,]) != 0){
    covariates_chinook0_skagit[i,] = scale(covariates_chinook0_skagit[i,], center = FALSE, scale= TRUE)[,1]
  }
}




#subset response variable
subset_chinook_summer_perhour <- arrange(data_day_night,doy) %>%
  filter(doy > 150 & doy < 190) %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category,trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_chinook_summer_perhour)[1]){
  subset_chinook_summer_perhour[i,] = scale(subset_chinook_summer_perhour[i,])[,1]
}

######


#function for C matrix
#######


Cmat_skagit <- function(nyears,ncov,hatchery=0, day_on_night = FALSE){
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
              if(i%%2 == 1){
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


Qmat_skagit <- function(nyears){
  Q <- matrix(list(0),nrow = nyears,ncol = nyears, byrow = TRUE)
  for(i in 1:nyears){
    for(j in 1:nyears){
      if(i==j){
        if(i%%2 == 1){
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



mod_list <- function(nyears,ncov,hatchery=0, day_on_night = FALSE, unequal_q = FALSE){
  
  if(unequal_q){
    Q = Qmat_skagit(nyears)
    
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
      C = Cmat_skagit(nyears,ncov,hatchery, day_on_night)
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

#first see if equal errors or unequal errors have lower aicc
nyears = num_years*2*2
fit.model <- mod_list(nyears, 0, 0, FALSE, FALSE)

fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
fit$convergence
fit$AICc #2542.9
nyears = num_years*2*2
fit.model_unequalq <- mod_list(nyears, 0, 0, FALSE, TRUE)
fit_unequal <- MARSS(subset_chinook_summer_perhour, model=fit.model_unequalq, 
                     silent = TRUE, method = "BFGS",
                     control=list(maxit=2000))
fit_unequal$convergence
fit_unequal$AICc #2541

#let's use unequal q


out.tab.all.years <- NULL
fits.all.years <- NULL
nyears = num_years*2*2
for(kk in c(1,3,4,6,9)){
  c = covariates_chinook0_skagit[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_chinook0_skagit)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-15)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0, FALSE, TRUE))
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years=rbind(out.tab.all.years,out)
  fits.all.years=c(fits.all.years,list(fit))
}

out.tab.all.years
#residuals have the least aicc

out.tab.all.years$DeltaAICc <- out.tab.all.years$AICc - min(out.tab.all.years$AICc)
#reorder by delta aicc
out.tab.all.years <- out.tab.all.years[order(out.tab.all.years$DeltaAICc),]

#round to 2 decimal places
out.tab.all.years$DeltaAICc <- round(out.tab.all.years$DeltaAICc,2)
out.tab.all.years$logLik <- round(out.tab.all.years$logLik,2)
out.tab.all.years$AICc <- round(out.tab.all.years$AICc,2)


#save the dataframe
write.csv(out.tab.all.years, 
          here("skagit", "output","model_selection_correlated_variables_chinook0.csv"))

list_combinations <- get_covariate_combinations(1:7)
out.tab.hatchery.all.years <- NULL
fits.hatchery.all.years <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  residuals = 0
  flow = 0
  flow_difference = 0
  hatchery = 0
  temperature_difference = 0
  photoperiod = 0
  lunar_phase = 0
  for(j in covariates){
    if(j == 1){
      k = 6
      residuals = 1
    }
    else if(j==2){
      k = 2
      flow = 1
    }
    else if(j==3){
      k = 3
      photoperiod = 1
      
    }
    else if(j==4){
      k = 7
      temperature_difference = 1
    }
    else if(j==5){
      k = 8
      flow_difference = 1
    }
    else if(j==6){
      k = 5
      lunar_phase = 1
    }
    else if(j==7){
      k = 10
      hatchery = 1
    }
    
    c = rbind(c,covariates_chinook0_skagit[((1+(k-1)*nyears):(k*nyears)),])
    name_long = rownames(covariates_chinook0_skagit)[1+(k-1)*nyears]
    name = paste(name, substr(name_long,1,nchar(name_long)-15))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(k==10){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery, FALSE ,TRUE))
  }
  else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery, FALSE, TRUE))
  }
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, residuals=residuals, flow=flow, flow_difference=flow_difference,
                 temperature_difference=temperature_difference, photoperiod = photoperiod,
                 lunar_phase = lunar_phase ,hatchery=hatchery, 
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.hatchery.all.years=rbind(out.tab.hatchery.all.years,out)
  fits.hatchery.all.years=c(fits.hatchery.all.years,list(fit))
  
  
}
out.tab.hatchery.all.years$deltaAICc <- out.tab.hatchery.all.years$AICc - min(out.tab.hatchery.all.years$AICc)
min.AICc <- order(out.tab.hatchery.all.years$AICc)
out.tab.hatchery.ordered <- out.tab.hatchery.all.years[min.AICc, ]
out.tab.hatchery.ordered


#add model with no covariates
out.tab.hatchery.ordered$deltaAICc <- NULL
model_wo_covariates <- mod_list(nyears,0,0, FALSE, TRUE)
fit <- MARSS(subset_chinook_summer_perhour, model=model_wo_covariates, 
             silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
fits.hatchery.all.years=c(fits.hatchery.all.years,list(fit))
out=data.frame(c="no covariates", residuals=0, flow=0, flow_difference=0,
               temperature_difference=0, photoperiod = 0, lunar_phase = 0 ,
               hatchery=0, 
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)
out.tab.hatchery.ordered=rbind(out.tab.hatchery.ordered,out)
out.tab.hatchery.ordered$deltaAICc <- out.tab.hatchery.ordered$AICc - min(out.tab.hatchery.ordered$AICc)
out.tab.hatchery.ordered <- out.tab.hatchery.ordered[order(out.tab.hatchery.ordered$deltaAICc),]
out.tab.hatchery.ordered

#round to 2 decimal places
out.tab.hatchery.ordered$deltaAICc <- round(out.tab.hatchery.ordered$deltaAICc,2)
out.tab.hatchery.ordered$logLik <- round(out.tab.hatchery.ordered$logLik,2)
out.tab.hatchery.ordered$AICc <- round(out.tab.hatchery.ordered$AICc,2)


#save the dataframe

write.csv(out.tab.hatchery.ordered, 
          here("skagit", "output","model_selection_chinook0.csv"))

#save the fits
save(fits.hatchery.all.years, file = here("skagit", "output","fits_hatchery_all_years_chinook0.RData"))

#calculate the relative variable importance


out.tab.hatchery.ordered$deltaAICc <- NULL
out.tab.hatchery.ordered$rel.LL <- NULL
out.tab.hatchery.ordered$weights <- NULL


weights <- akaike.weights(out.tab.hatchery.ordered$AICc)

out.tab.hatchery.ordered$deltaAICc <- weights$deltaAIC
out.tab.hatchery.ordered$rel.LL <- weights$rel.LL
out.tab.hatchery.ordered$weights <- weights$weights


min.AICc <- order(out.tab.hatchery.ordered$AICc)
out.tab.hatchery.ordered <- out.tab.hatchery.ordered[min.AICc, ]
out.tab.hatchery.ordered

out.tab.hatchery.ordered$cumulative_weights <- cumsum(out.tab.hatchery.ordered$weights)

relative_importance_residuals <- sum(out.tab.hatchery.ordered$weights[out.tab.hatchery.ordered$residuals==1])
relative_importance_photo <- sum(out.tab.hatchery.ordered$weights[out.tab.hatchery.ordered$photoperiod==1])
relative_importance_temperature_difference <- sum(out.tab.hatchery.ordered$weights[out.tab.hatchery.ordered$temperature_difference==1])
relative_importance_flow <- sum(out.tab.hatchery.ordered$weights[out.tab.hatchery.ordered$flow==1])
relative_importance_lunar_phase <- sum(out.tab.hatchery.ordered$weights[out.tab.hatchery.ordered$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab.hatchery.ordered$weights[out.tab.hatchery.ordered$hatchery==1])
relative_importance_flow_diff <- sum(out.tab.hatchery.ordered$weights[out.tab.hatchery.ordered$flow_difference==1])

riv_chinook_skagit <- data.frame(variable = c("photoperiod",
                                          "residuals",
                                          "temperature difference",
                                          "flow",
                                          "lunar phase",
                                          "flow difference",
                                          "hatchery"),
                             relative_importance = c(relative_importance_photo,
                                                     relative_importance_residuals,
                                                     relative_importance_temperature_difference,
                                                     relative_importance_flow,
                                                     relative_importance_lunar_phase,
                                                     relative_importance_flow_diff,
                                                     relative_importance_hatchery))


#sort according to relative importance
riv_chinook_skagit[order(riv_chinook_skagit$relative_importance, decreasing = TRUE),]

riv_chinook_skagit <- riv_chinook_skagit[order(riv_chinook_skagit$relative_importance, decreasing = TRUE),]
riv_chinook_skagit$relative_importance <- round(riv_chinook_skagit$relative_importance,1)

riv_chinook_skagit

#save the relative importance table
write.csv(riv_chinook_skagit, here("output","riv_chinook0_skagit.csv"), row.names = FALSE)

out.tab.hatchery.ordered
#clean up table by reducing number of digits to 2 and by only keeping
#the variables of interest


out.tab.hatchery.ordered$rel.LL <- round(
  out.tab.hatchery.ordered$rel.LL,2)
out.tab.hatchery.ordered$weights <- round(
  out.tab.hatchery.ordered$weights,2)
out.tab.hatchery.ordered$cumulative_weights <- round(
  out.tab.hatchery.ordered$cumulative_weights,2)

#save the table
write.csv(out.tab.hatchery.ordered, here("output","chinook0_model_selection_w_weights.csv"), row.names = FALSE)

fits.hatchery.all.years[[34]]
autoplot(fits.hatchery.all.years[[34]])

autoplot(fits.hatchery.all.years[[64]])
#the model is selecting no error for all the day estimates
#there are very few day values
#perhaps it is best to add the day and night and have a daily value


#perhaps it is best to add the day and night