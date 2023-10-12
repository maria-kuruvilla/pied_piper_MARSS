#Goal - to try MARSS package in a subset for chinook data with day and night separate
# data files

#steps

#interpolate the hatchery covariate values
#add residuals of photperiod and temperature to the dataset
#save the dataset
#create subset
#run MARSS

#Load libraries


library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)
library(MASS)
library(modelr)

#check system

if(.Platform$OS.type == "unix") {
  data_string = here("..","..","data","pied_piper","dungeness")
} else {
  data_string = here("..","..","..","OneDrive","Documents","data","pied_piper","dungeness")
}


#read data
data_day <- read.csv(here(data_string,"dungeness_unaggregated_day.csv"), header = TRUE)

data_night <- read.csv(here(data_string,"dungeness_unaggregated_night.csv"), header = TRUE)

#interpolate the hatchery covariate values
data_day$chinook0_hatchery_perhour_interpolate <- na.approx(data_day$chinook0_hatchery_perhour, na.rm = FALSE)
data_night$chinook0_hatchery_perhour_interpolate <- na.approx(data_night$chinook0_hatchery_perhour, na.rm = FALSE)

#add residuals of photoperiod and temperature to the dataset

lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day)

data_day <- data_day %>% add_residuals(lm_temp_day)

lm_temp_night <- lm(temp ~ 1+photoperiod, data = data_night)

data_night <- data_night %>% add_residuals(lm_temp_night)

#save the dataset

write.csv(data_day,here("data","dungeness_unaggregated_day_w_covariates.csv"))
write.csv(data_night,here("data","dungeness_unaggregated_night_w_covariates.csv"))


write.csv(data_day,here(data_string,"dungeness_unaggregated_day_w_covariates.csv"))
write.csv(data_night,here(data_string,"dungeness_unaggregated_night_w_covariates.csv"))

#subset the data

#trying to row bind data_day and data_night before subsetting
data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"

data_day_night <- rbind(data_day,data_night)

#subset

covariates_2008_chinook0 <- arrange(data_day_night,doy) %>%
  filter(year == 2008 & doy >130 & doy <= 200) %>%
  dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_april, 
                lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
                chinook0_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = c(
    temp, flow, photoperiod, atu_april, lunar_phase, 
    resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

#scale variables

for(i in 1:(dim(covariates_2008_chinook0)[1]-10)){ # everything except resid, diffs and hatchery
  covariates_2008_chinook0[i,] = scale(covariates_2008_chinook0[i,])[,1]
}


for(i in 11:20){
  covariates_2008_chinook0[i,] = scale(covariates_2008_chinook0[i,], center = FALSE, scale= TRUE)[,1]
}

#subset response variable
subset_chinook_summer_perhour <- arrange(data_day_night,doy) %>%
  filter(year == 2008 & doy > 130 & doy <= 200) %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_chinook_summer_perhour)[1]){
  subset_chinook_summer_perhour[i,] = scale(subset_chinook_summer_perhour[i,])[,1]
}


#run MARSS
mod.list_0_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal"
)

fit.model = mod.list_0_0
fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
#it works

#try with covariates

#making model list

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
  C = matrix(list("a",0,"b",0,
                  
                  0,"a",0,"b"),
             2,4, byrow = TRUE)
)

mod.list_3_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,"b",0,"c",0,
                  
                  0,"a",0,"b",0,"c"),
             2,6, byrow = TRUE)
)

mod.list_4_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,"b",0,"c",0,"d",0,
                  
                  0,"a",0,"b",0,"c",0,"d"),
             2,8, byrow = TRUE)
)

mod.list_5_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,"b",0,"c",0,"d",0,"e",0,
                  
                  0,"a",0,"b",0,"c",0,"d",0,"e"),
             2,10, byrow = TRUE)
)




#making combination of covariates

get_covariate_combinations <- function(covariates) {
  n <- length(covariates)
  combinations <- lapply(1:n, function(x) combn(covariates, x, simplify = FALSE))
  unlist(combinations, recursive = FALSE)
}

list_combinations <- get_covariate_combinations(1:5)
out.tab<- NULL
fits <- list()
######
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
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
    
    c = rbind(c,covariates_2008_chinook0[((1+(k-1)*2):(k*2)),])
    name_long = rownames(covariates_2008_chinook0)[1+(k-1)*2]
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

  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))


  out=data.frame(c=name, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab=rbind(out.tab,out)
  fits=c(fits,list(fit))

}
out.tab
#####
#let's try with different hatchery effects for day and night

#####
mod.list_1_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = "diagonal and unequal"
)

mod.list_2_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,"b",0,
                  
                  0,"a",0,"b_night"),
             2,4, byrow = TRUE)
)

mod.list_3_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,"b",0,"c",0,
                  
                  0,"a",0,"b",0,"c_night"),
             2,6, byrow = TRUE)
)

mod.list_4_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,"b",0,"c",0,"d",0,
                  
                  0,"a",0,"b",0,"c",0,"d_night"),
             2,8, byrow = TRUE)
)

mod.list_5_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,"b",0,"c",0,"d",0,"e",0,
                  
                  0,"a",0,"b",0,"c",0,"d",0,"e_night"),
             2,10, byrow = TRUE)
)

out.tab.hatchery<- NULL
fits.hatchery <- list()

for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
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
    
    c = rbind(c,covariates_2008_chinook0[((1+(k-1)*2):(k*2)),])
    name_long = rownames(covariates_2008_chinook0)[1+(k-1)*2]
    name_individual = substr(name_long,1,nchar(name_long)-9)
    # print(name_individual)
    name = paste(name, substr(name_long,1,nchar(name_long)-9))
    
  }
  # print(c)
  print(name)
  if(name_individual == "chinook0_hatchery_perhour_interpolate"){
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
  
  # print(has_hatchery)
  

  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))


  out=data.frame(c=name, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.hatchery=rbind(out.tab.hatchery,out)
  fits.hatchery=c(fits.hatchery,list(fit))
  
}
out.tab.hatchery
out.tab.hatchery[order(out.tab.hatchery$AICc),]
out.tab[order(out.tab$AICc),]


#####


#trying for all temp type variables


#######

num_years = 1
num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:5)
out.tab<- NULL
fits <- list()
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
        
        c = rbind(c,covariates_2008_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
        name_long = rownames(covariates_2008_chinook0)[1+(k-1)*num_rows]
        name_individual = substr(name_long,1,nchar(name_long)-9)
        name = paste(name, substr(name_long,1,nchar(name_long)-9))
        
      }
      # print(c)
      print(name)
      c_num <- length(covariates)
      if(name_individual == "chinook0_hatchery_perhour_interpolate"){
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
      
      fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                   control=list(maxit=2000))
      
      
      out=data.frame(c=name, d = "None",
                     logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                     num.iter=fit$numIter, converged=!fit$convergence,
                     stringsAsFactors = FALSE)
      out.tab=rbind(out.tab,out)
      fits=c(fits,list(fit))
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
      
      c = rbind(c,covariates_2008_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
      name_long = rownames(covariates_2008_chinook0)[1+(k-1)*num_rows]
      num_individual = substr(name_long,1,nchar(name_long)-9)
      name = paste(name, substr(name_long,1,nchar(name_long)-9))
      
    }
    # print(c)
    print(name)
    c_num <- length(covariates)
    if(name_individual == "chinook0_hatchery_perhour_interpolate"){
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
    
    
    fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                 control=list(maxit=2000))
    
    
    out=data.frame(c=name, d = "None",
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence,
                   stringsAsFactors = FALSE)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
    
  }
  
  
}
out.tab$deltaAICc <- out.tab$AICc - min(out.tab$AICc)
min.AICc <- order(out.tab$AICc)
out.tab.ordered <- out.tab[min.AICc, ]
out.tab.ordered


