#goal - run MARSS analysis with no night on day effect
#and w photoperiod in model selection
#load libraries

library(MARSS)
library(broom)
library(here)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)
library(GGally)
library(ggplot2)
library(tidyverse)


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




num_years = 15 
num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:7)
out.tab_photo_diff2<- NULL
fits_photo_diff2 <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  photoperiod = 0
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
      k = 3
      photoperiod = 1
    }
    else if(j==3){
      k = 7
      temperature_difference = 1
    }
    else if(j==4){
      k = 2
      flow = 1
    }
    else if(j==5){
      k = 5
      lunar_phase = 1
    }
    else if(j==6){
      k = 8
      flow_difference = 1
    }
    else if(j==7){
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
  
  if(k==10){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery, FALSE, TRUE))
  }
  else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery, FALSE, TRUE))
  }
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, photoperiod = photoperiod,
                 photoperiod_difference = photoperiod_difference,
                 temperature_difference = temperature_difference, flow = flow,
                 flow_difference = flow_difference,
                 lunar_phase = lunar_phase, hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_photo_diff2=rbind(out.tab_photo_diff2,out)
  fits_photo_diff2=c(fits_photo_diff2,list(fit))
  
  
}

out.tab_photo_diff2$deltaAICc <- out.tab_photo_diff2$AICc - min(out.tab_photo_diff2$AICc)
min.AICc <- order(out.tab_photo_diff2$AICc)
out.tab_photo_diff.ordered2 <- out.tab_photo_diff2[min.AICc, ]
out.tab_photo_diff.ordered2

ci_best3 <- tidy(fits_photo_diff2[[79]])
ci_best3


#calculate the relative variable importance

mod.list_0_0 <- mod_list(num_rows,0,0, FALSE, TRUE)
fit <- MARSS(subset_chinook_summer_perhour, model=mod.list_0_0, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
photoperiod_difference = 0
photoperiod = 0
temperature_difference = 0
flow = 0
lunar_phase = 0
hatchery = 0
flow_difference = 0
name = "none"
out=data.frame(c=name, photoperiod = photoperiod,
               photoperiod_difference = photoperiod_difference,
               temperature_difference = temperature_difference, flow = flow,
               flow_difference = flow_difference,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)



#remove deltaAICc column from out.tab_photo_diff2
out.tab_photo_diff2$deltaAICc <- NULL
out.tab_photo_diff2=rbind(out.tab_photo_diff2,out)
fits_photo_diff2=c(fits_photo_diff2,list(fit))

out.tab_photo_diff2$deltaAICc <- NULL
out.tab_photo_diff2$rel.LL <- NULL
out.tab_photo_diff2$weights <- NULL


weights <- akaike.weights(out.tab_photo_diff2$AICc)

out.tab_photo_diff2$deltaAICc <- weights$deltaAIC
out.tab_photo_diff2$rel.LL <- weights$rel.LL
out.tab_photo_diff2$weights <- weights$weights


min.AICc <- order(out.tab_photo_diff2$AICc)
out.tab_photo_diff2.ordered <- out.tab_photo_diff2[min.AICc, ]
out.tab_photo_diff2.ordered

out.tab_photo_diff2.ordered$cumulative_weights <- cumsum(out.tab_photo_diff2.ordered$weights)

relative_importance_photo_diff <- sum(out.tab_photo_diff2$weights[out.tab_photo_diff2$photoperiod_difference==1])
relative_importance_photo <- sum(out.tab_photo_diff2$weights[out.tab_photo_diff2$photoperiod==1])
relative_importance_temperature_difference <- sum(out.tab_photo_diff2$weights[out.tab_photo_diff2$temperature_difference==1])
relative_importance_flow <- sum(out.tab_photo_diff2$weights[out.tab_photo_diff2$flow==1])
relative_importance_lunar_phase <- sum(out.tab_photo_diff2$weights[out.tab_photo_diff2$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab_photo_diff2$weights[out.tab_photo_diff2$hatchery==1])
relative_importance_flow_diff <- sum(out.tab_photo_diff2$weights[out.tab_photo_diff2$flow_difference==1])

riv_photo_diff <- data.frame(variable = c("photoperiod difference",
                                          "photoperiod",
                                          "temperature difference",
                                          "flow",
                                          "lunar phase",
                                          "flow difference",
                                          "hatchery"),
                             relative_importance = c(relative_importance_photo_diff,
                                                     relative_importance_photo,
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

#save the relative importance table
write.csv(riv_photo_diff, here("output",
                               "riv_photo_diff_chinook0_corrected_w_photoperiod.csv"), row.names = FALSE)

out.tab_photo_diff2.ordered
#clean up table by reducing number of digits to 2 and by only keeping
#the variables of interest

out.tab_photo_diff2.ordered$deltaAICc <- round(
  out.tab_photo_diff2.ordered$deltaAICc,2)
out.tab_photo_diff2.ordered$rel.LL <- round(
  out.tab_photo_diff2.ordered$rel.LL,2)
out.tab_photo_diff2.ordered$weights <- round(
  out.tab_photo_diff2.ordered$weights,2)
out.tab_photo_diff2.ordered$cumulative_weights <- round(
  out.tab_photo_diff2.ordered$cumulative_weights,2)

#save the table
write.csv(out.tab_photo_diff2.ordered, here("output",
                                            "chinook0_corrected_model_selection_w_photoperiod.csv"), 
          row.names = FALSE)

