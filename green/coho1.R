#goal - to try MARSS analysis for coho1 green river data

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

data_day <- read.csv(here("green", "data", "green_covariates_day.csv"))
data_night <- read.csv(here("green", "data", "green_covariates_night.csv"))

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"


#make columns for difference in temp, flow, photoperiod and call it 
#temp_diff, flow_diff, photo_diff

data_day$temp_diff <- c(NA,diff(data_day$temp_inp))
data_day$flow_diff <- c(NA,diff(data_day$flow))
data_day$photo_diff <- c(NA,diff(data_day$photoperiod))

data_night$temp_diff <- c(NA,diff(data_night$temp_inp))
data_night$flow_diff <- c(NA,diff(data_night$flow))
data_night$photo_diff <- c(NA,diff(data_night$photoperiod))

lm_temp_day <- lm(temp_inp ~ 1+photoperiod, data = data_day)

data_day <- data_day %>% add_residuals(lm_temp_day)

lm_temp_night <- lm(temp_inp ~ 1+photoperiod, data = data_night)

data_night <- data_night %>% add_residuals(lm_temp_night)

#need to go back to the python code and add atu to the dataframe - done

#adding coho1_hatchery_perhour=0 at the end of the data from
data_day$coho1_hatchery_perhour[dim(data_day)[1]] <- 0
data_night$coho1_hatchery_perhour[dim(data_night)[1]] <- 0

#adding column for interpolated values of coho1_hatchery_perhour to data_day

data_day$coho1_hatchery_perhour_inp <- na.approx(data_day$coho1_hatchery_perhour,na.rm = FALSE)

#doing same for night

data_night$coho1_hatchery_perhour_inp <- na.approx(data_night$coho1_hatchery_perhour,na.rm = FALSE)


data_day_night <- rbind(data_day,data_night)

#plot coho1_wild_perhour as a function of doy

ggplot(data_day_night, aes(x = doy, y = coho1_wild_perhour)) + 
  geom_point() + 
  facet_wrap(~year, scales = "free_y")

#do same for coho hatchery

ggplot(data_day_night, aes(x = doy, y = coho1_hatchery_perhour)) + 
  geom_point() + 
  facet_wrap(~year, scales = "free_y")

ggplot(data_day_night, aes(x = doy, y = coho1_mixed_perhour)) + 
  geom_point() + 
  facet_wrap(~year, scales = "free_y")


covariates_coho1_green <- arrange(data_day_night,doy) %>%
  filter((year == 2014 | year == 2017) & doy >70 & doy <= 150) %>%
  dplyr::select(year,doy, daytime_category, temp_inp, flow_inp, photoperiod, atu_solstice,
                lunar_phase, resid, temp_diff, flow_diff, photo_diff,
                coho1_hatchery_perhour_inp) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = c(
    temp_inp, flow_inp, photoperiod, atu_solstice, lunar_phase,
    resid, temp_diff, flow_diff, photo_diff, coho1_hatchery_perhour_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()



#scaling the variables

num_years = 2
num_rows = num_years*2
num_covariates = 10
total_covariates = dim(covariates_coho1_green)[1]


for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_coho1_green[i,] = scale(covariates_coho1_green[i,])[,1]
}

#just scale
#need to skip rows which are just zeros - 2017,2018,2021 day and 2018, 2021 night

for(i in (total_covariates/2 + 1):(total_covariates)){
  if(sum(covariates_coho1_green[i,]) != 0){
    covariates_coho1_green[i,] = scale(covariates_coho1_green[i,], center = FALSE, scale= TRUE)[,1]
  }
}




#subset response variable
subset_coho_summer_perhour <- arrange(data_day_night,doy) %>%
  filter((year == 2014 | year == 2017) & doy >70 & doy <= 150) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour)[1]){
  subset_coho_summer_perhour[i,] = scale(subset_coho_summer_perhour[i,])[,1]
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



out.tab.all.years <- NULL
fits.all.years <- NULL
nyears = num_years*2
for(kk in c(1,2,3,4,6)){
  c = covariates_coho1_green[((1+(kk-1)*num_rows):(kk*num_rows)),]
  name_long = rownames(covariates_coho1_green)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
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
#atu solstice has least aicc - 4


out.tab.all.years2 <- NULL
fits.all.years2 <- NULL
nyears = num_years*2
for(kk in c(7,9)){
  c = covariates_coho1_green[((1+(kk-1)*num_rows):(kk*num_rows)),]
  name_long = rownames(covariates_coho1_green)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0))
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years2=rbind(out.tab.all.years2,out)
  fits.all.years2=c(fits.all.years2,list(fit))
}

out.tab.all.years2

#temp diff has much lower aicc


list_combinations <- get_covariate_combinations(1:5)
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
      k = 7
    }
    else if(j==3){
      k = 5
    }
    else if(j==4){
      k = 8
    }
    else if(j==5){
      k = 10
    }
    
    c = rbind(c,covariates_coho1_green[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_coho1_green)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-5))
    
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

ci_best <- tidy(fits.hatchery.all.years[[20]])

ggplot(ci_best[c(7:10),], 
       aes(x = c("ATU", "Lunar phase","Hatchery\n day", "Hatchery\n night"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  theme(axis.text.x=element_text(size=20),axis.title.y=element_text(size=24))



ggsave(here("green","output","coho_all_years_effects_estimate_best_model.png"), width = 10, height = 8, units = "in", dpi = 300)



#maybe include turbidity??

