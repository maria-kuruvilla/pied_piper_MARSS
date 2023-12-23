#Goal - repeat MARSS analysis for coho data in skagit
#adding 2020,2021,2022 data
# Load libraries
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

data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))

data_day_night$coho1_wild_perhour <- data_day_night$coho1_wild_num/data_day_night$In
data_day_night$coho1_hatchery_perhour <- data_day_night$coho1_hatchery_num/data_day_night$In


data_day_night$coho1_hatchery_perhour_inp <- na.approx(data_day_night$coho1_hatchery_perhour, 
                                                       na.rm = FALSE)

lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day_night)

data_day_night <- data_day_night %>% add_residuals(lm_temp_day)

#changing it from <=150 to <150 to see if it changes model selection results
covariates_coho1_skagit <- arrange(data_day_night,doy) %>%
  filter(doy >100 & doy < 150) %>%
  dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_solstice,
                lunar_phase, resid, temp_diff, flow_diff, photo_diff,
                coho1_hatchery_perhour_inp, trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = c(
    temp, flow, photoperiod, atu_solstice, lunar_phase,
    resid, temp_diff, flow_diff, photo_diff, coho1_hatchery_perhour_inp), names_vary = "fastest") %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()




#scaling the variables

num_years = 2022-2010
num_rows = num_years*2
num_covariates = 10
total_covariates = dim(covariates_coho1_skagit)[1]


for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_coho1_skagit[i,] = scale(covariates_coho1_skagit[i,])[,1]
}

#just scale


for(i in (total_covariates/2 + 1):(total_covariates)){
  if(sum(covariates_coho1_skagit[i,]) != 0){
    covariates_coho1_skagit[i,] = scale(covariates_coho1_skagit[i,], center = FALSE, scale= TRUE)[,1]
  }
}



#subset response variable
subset_coho_summer_perhour <- arrange(data_day_night,doy) %>%
  filter(doy > 100 & doy < 150) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category,trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour)[1]){
  subset_coho_summer_perhour[i,] = scale(subset_coho_summer_perhour[i,])[,1]
}


#function for C matrix
#######


Cmat_skagit <- function(nyears,ncov,hatchery=0){
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


Qmat_skagit_daynight <- function(nyears){
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


Qmat_skagit_trap_daynight <- function(nyears){
  count = 0
  Q <- matrix(list(0),nrow = nyears,ncol = nyears, byrow = TRUE)
  for(i in 1:nyears){
    for(j in 1:nyears){
      if(i==j){
        if(i%%2 == 1){
          count <- count + 1
          if(count%%2 == 1){
            Q[i,j] <- "q_scoop_day"
          }
          else{
            Q[i,j] <- "q_screw_day"
          }
        }
        else{
          # count <- count + 1
          if(count%%2 == 1){
            Q[i,j] <- "q_scoop_night"
          }
          else{
            Q[i,j] <- "q_screw_night"
          }
        }
        
      }
    }
  }
  return(Q)
}


#function to make a square matrix of 0s of length nyears
#if the row is equal to the column, then change the value
#if the row is 1,2,5,6,9,10,... then change the value to q_scoop
#if the row is 3,4,7,8,11,12,... then change the value to q_screw

Qmat_skagit_trap <- function(nyears){
  count = 0
  Q <- matrix(list(0),nrow = nyears,ncol = nyears, byrow = TRUE)
  for(i in 1:nyears){
    for(j in 1:nyears){
      if(i==j){
        if(i%%2 == 1){
          count <- count + 1
          }
        
        if(count%%2 == 1){
          Q[i,j] <- "q_scoop"
        }
        else{
          Q[i,j] <- "q_screw"
        }
      }
    }
  }
  return(Q)
}




#######



#######

#function for mod list

######

mod_list <- function(nyears,ncov,hatchery=0, unequal_q = FALSE, trap = FALSE, daynight = FALSE){
  
  if(unequal_q){
    
    if(trap & daynight){
      Q = Qmat_skagit_trap_daynight(nyears)
    }
    
    else if(trap){
      Q = Qmat_skagit_trap(nyears)
    }
    else if(daynight){
      Q = Qmat_skagit_daynight(nyears)
    }
    else{
      Q = Qmat_skagit(nyears)
    }
    
    
    
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
      C = Cmat_skagit(nyears,ncov,hatchery)
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


out.tab.all.years <- NULL
fits.all.years <- NULL
nyears = num_years*2*2
for(kk in c(1,3,4,6,9)){
  c = covariates_coho1_skagit[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_coho1_skagit)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-15)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0))
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years=rbind(out.tab.all.years,out)
  fits.all.years=c(fits.all.years,list(fit))
}

out.tab.all.years
#atu_solstice have the least aicc


#now try with all covariates

list_combinations <- get_covariate_combinations(1:6)
out.tab.hatchery.all.years <- NULL
fits.hatchery.all.years <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  for(j in covariates){
    if(j == 1){
      k = 4
    }
    else if(j==2){
      k = 2
    }
    else if(j==3){
      k = 5
    }
    else if(j==4){
      k = 7
    }
    else if(j==5){
      k = 8
    }
    else if(j==6){
      k = 10
    }
    
    c = rbind(c,covariates_coho1_skagit[((1+(k-1)*nyears):(k*nyears)),])
    name_long = rownames(covariates_coho1_skagit)[1+(k-1)*nyears]
    name = paste(name, substr(name_long,1,nchar(name_long)-15))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(k==10){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery))
  }
  else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery))
  }
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, d = "None",
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

ci_coho_skagit <- tidy(fits.hatchery.all.years[[25]])
autoplot(fits.hatchery.all.years[[25]])

#plot the estimates of the best model

ggplot(ci_coho_skagit[c(51:54),], aes(x = c("ATU","Flow","Hatchery,\nday",
                                            "Hatchery,\nnight"),
                                      y = estimate, 
                                      ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  ggtitle("Coho yearlings") + 
  theme(plot.title = element_text(size = 20))+
  theme(axis.text.x=element_text(size=14),axis.title.y=element_text(size=14))

ggsave(here("skagit","output","coho_skagit_best_model_all_years.jpeg"), width = 10, height = 10, units = "in")


#checking if there is support for different errors for day and night

fit.model = mod_list(nyears,0,0,FALSE, FALSE, FALSE)
fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
fit$AICc
#2978

fit.model_unequal_q_day_night = mod_list(nyears,0,0,TRUE, FALSE, TRUE)
fit_unequal_q_day_night <- MARSS(subset_coho_summer_perhour, 
                                 model=fit.model_unequal_q_day_night, 
                                 silent = TRUE, method = "BFGS",
                                 control=list(maxit=2000))
fit_unequal_q_day_night$AICc
#2945


fit.model_unequal_q_trap = mod_list(nyears,0,0,TRUE, TRUE, TRUE)
fit_unequal_q_trap <- MARSS(subset_coho_summer_perhour, 
                                 model=fit.model_unequal_q_trap, 
                                 silent = TRUE, method = "BFGS",
                                 control=list(maxit=2000))
fit_unequal_q_trap$AICc
#2948

#just keep the unequal q for day and night





out.tab.all.years.unequalq <- NULL
fits.all.years.unequalq <- NULL
nyears = num_years*2*2
for(kk in c(1,3,4,6,9)){
  c = covariates_coho1_skagit[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_coho1_skagit)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-15)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0,TRUE,FALSE,TRUE))
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years.unequalq=rbind(out.tab.all.years.unequalq,out)
  fits.all.years.unequalq=c(fits.all.years.unequalq,list(fit))
}

out.tab.all.years.unequalq
#atu_solstice have the least aicc


#now try with all covariates

list_combinations <- get_covariate_combinations(1:6)
out.tab.hatchery.all.years.unequalq <- NULL
fits.hatchery.all.years.unequalq <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  for(j in covariates){
    if(j == 1){
      k = 4
      #atu
    }
    else if(j==2){
      k = 2
      #flow
    }
    else if(j==3){
      k = 5
      #lunar phase
    }
    else if(j==4){
      k = 7
      #temperature difference
    }
    else if(j==5){
      k = 8
      #flow difference
    }
    else if(j==6){
      k = 10
      #hatchery
    }
    
    c = rbind(c,covariates_coho1_skagit[((1+(k-1)*nyears):(k*nyears)),])
    name_long = rownames(covariates_coho1_skagit)[1+(k-1)*nyears]
    name = paste(name, substr(name_long,1,nchar(name_long)-15))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(k==10){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery, TRUE, FALSE, TRUE))
  }
  else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery))
  }
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.hatchery.all.years.unequalq=rbind(out.tab.hatchery.all.years.unequalq,out)
  fits.hatchery.all.years.unequalq=c(fits.hatchery.all.years.unequalq,list(fit))
  
  
}
out.tab.hatchery.all.years.unequalq$deltaAICc <- out.tab.hatchery.all.years.unequalq$AICc - min(out.tab.hatchery.all.years.unequalq$AICc)
min.AICc <- order(out.tab.hatchery.all.years.unequalq$AICc)
out.tab.hatchery.ordered.unequalq <- out.tab.hatchery.all.years.unequalq[min.AICc, ]
out.tab.hatchery.ordered.unequalq

ci <- tidy(fits.hatchery.all.years.unequalq[[25]])

ggplot(ci[c(52:55),], aes(x = c("ATU","Flow","Hatchery,\nday",
                                "Hatchery,\nnight"),
                          y = estimate, 
                          ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  ggtitle("Coho yearlings") + 
  theme(plot.title = element_text(size = 20))+
  theme(axis.text.x=element_text(size=14),axis.title.y=element_text(size=14))

ggsave(here("skagit","output","coho_skagit_best_model_all_years.jpeg"), width = 10, height = 10, units = "in")




#trying out just the best model with unequal q

nyears = num_years*2*2
name = NULL
c = NULL
for(k in c(4,2,10)){
  c = rbind(c,covariates_coho1_skagit[((1+(k-1)*nyears):(k*nyears)),])
  name_long = rownames(covariates_coho1_skagit)[1+(k-1)*nyears]
  name = paste(name, substr(name_long,1,nchar(name_long)-15))
  
}

fit.model = c(list(c= c), mod_list(nyears,3,1,TRUE,FALSE,TRUE))
fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
out=data.frame(c=name, d = "None",
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)

fit$AICc
#2870 was with equal q
#2830 with unequal q

ci <- tidy(fit)

ggplot(ci[c(52:55),], aes(x = c("ATU","Flow","Hatchery,\nday",
                                            "Hatchery,\nnight"),
                                      y = estimate, 
                                      ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  ggtitle("Coho yearlings") + 
  theme(plot.title = element_text(size = 20))+
  theme(axis.text.x=element_text(size=14),axis.title.y=element_text(size=14))

ggsave(here("skagit","output","coho_skagit_best_model_all_years.jpeg"), width = 10, height = 10, units = "in")


autoplot(fit, plot.type = "fitted.ytT")

ggsave(here("skagit","output","coho_skagit_best_model_all_years_fitted.jpeg"), width = 10, height = 10, units = "in")
