#goal - to try MARSS analysis for chinook0 Skagit River data

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

#using ggplot plot chinook0_wild_perhour as a function of doy, facet wrap for years
#add a x limit from 150 to 200

ggplot(data_day_night, aes(x = doy, y = chinook0_wild_perhour)) +
  geom_point() +
  # xlim(150,200)+
  facet_wrap(~year, ncol = 2, scales = "free_y")
  

#do same for chinook_hatchery_num

ggplot(data_day_night, aes(x = doy, y = chinook0_hatchery_perhour)) +
  geom_point() +
  # xlim(150,200)+
  facet_wrap(~year, ncol = 2, scales = "free_y")

#doy limits should be 150 to 189

#making new columns for chinook0_wild_perhour and chinook0_hatchery_perhour and
#coho1_wild_perhour and coho1_hatchery_perhour and chinook1_wild_perhour and
#chinook1_hatchery_perhour

data_day_night$chinook0_wild_perhour <- data_day_night$chinook0_wild_num/data_day_night$In
data_day_night$chinook0_hatchery_perhour <- data_day_night$chinook0_hatchery_num/data_day_night$In
data_day_night$coho1_wild_perhour <- data_day_night$coho1_wild_num/data_day_night$In
data_day_night$coho1_hatchery_perhour <- data_day_night$coho1_hatchery_num/data_day_night$In
data_day_night$chinook1_wild_perhour <- data_day_night$chinook1_wild_num/data_day_night$In
data_day_night$chinook1_hatchery_perhour <- data_day_night$chinook1_hatchery_num/data_day_night$In

#plotting the perhour data for coho1

ggplot(data_day_night, aes(x = doy, y = coho1_wild_perhour)) +
  geom_point() +
  # xlim(150,190)+
  facet_wrap(~year, ncol = 2, scales = "free_y")

ggplot(data_day_night, aes(x = doy, y = coho1_hatchery_perhour)) +
  geom_point() +
  # xlim(150,190)+
  facet_wrap(~year, ncol = 2, scales = "free_y")
  
#doy limits should be 100 and 150

#plotting chinook1

ggplot(data_day_night, aes(x = doy, y = chinook1_wild_perhour)) +
  geom_point() +
  # xlim(150,190)+
  facet_wrap(~year, ncol = 2, scales = "free_y")

ggplot(data_day_night, aes(x = doy, y = chinook1_hatchery_perhour)) +
  geom_point() +
  # xlim(150,190)+
  facet_wrap(~year, ncol = 2, scales = "free_y")

#okay only doing for chinook and coho


#make new columns of interpolated values for chinook0_hatchery_perhour 
# and coho1_hatchery_perhour

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

#plot the correlation between temp,flow,photoperiod,resid,temp_diff,flow_diff,photo_diff

data_day_night %>%
  dplyr::select(flow,photoperiod, lunar_phase, temp, temp_diff, atu_solstice, resid, 
                photo_diff,flow_diff) %>%
  GGally::ggpairs(upper=list(continuous='points'),
                  lower=list(continuous='cor'))
#group 1 - photoperiod, temp, atu, photo_diff, resid are highly correlated - 1,3,4,6,9

ggsave(here("skagit", "output", "covariates_correlation_skagit.png"), width = 10, height = 10)

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

# Zmatrix - making shortcut for Z matrix for MARSS analysis 

Z1 <- factor(c("day_2010","night_2010","day_2010","night_2010",
               "day_2012","night_2012","day_2012","night_2012",
               "day_2013","night_2013","day_2013","night_2013",
               "day_2014","night_2014","day_2014","night_2014",
               "day_2015","night_2015","day_2015","night_2015",
               "day_2016","night_2016","day_2016","night_2016",
               "day_2017","night_2017","day_2017","night_2017",
               "day_2018","night_2018","day_2018","night_2018",
               "day_2019","night_2019","day_2019","night_2019"))

######

mod_list_skagit <- function(nyears,ncov,hatchery=0, day_on_night = FALSE, unequal_q = FALSE){
  
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
      Z = Z1,
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
      Z = Z1,
      A = "zero",
      R = "diagonal and equal",
      Q = Q,
      C = C
    )
  }
  
  return(mod.list)
}

######

fit.model <- mod_list(num_years, 0, 0, FALSE, FALSE)
fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
fit$convergence
fit
autoplot(fit)


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


#######

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



out.tab.all.years <- NULL
fits.all.years <- NULL
nyears = num_years*2*2
for(kk in c(1,3,4,6,9)){
  c = covariates_chinook0_skagit[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_chinook0_skagit)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-15)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0))
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years=rbind(out.tab.all.years,out)
  fits.all.years=c(fits.all.years,list(fit))
}

out.tab.all.years
#residuals have the least aicc





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
      k = 6
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
    fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery))
  }
  else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery))
  }
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
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


ci_best <- fits.hatchery.all.years[[9]]
autoplot(fits.hatchery.all.years[[9]])

#need to try different q for day and night

fit.model_equal_q = mod_list(nyears,0,0, day_on_night = FALSE, unequal_q = FALSE)

fit_equal_q <- MARSS(subset_chinook_summer_perhour, model=fit.model_equal_q, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

fit_equal_q$AICc

fit.model_unequal_q = mod_list(nyears,0,0, day_on_night = FALSE, unequal_q = TRUE)

fit_unequal_q <- MARSS(subset_chinook_summer_perhour, 
                       model=fit.model_unequal_q, silent = TRUE, method = "BFGS",
                        control=list(maxit=2000))

fit_unequal_q$AICc

#maybe try model selection with unequal q


