#goal - to combine day and night and scoop and trap variables into one column

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

#add coho1_hatchery_num for all values of trap and daytime_category
#add coho1_wild_num for all values of trap and daytime_category
#take mean of temp, flow, lunar_phase, photoperiod, photo_diff, temp_diff, flow_diff, atu_solstice
#for all values of trap and daytime_category
#drop trap and daytime_category
data <- data_day_night %>% 
  select(coho1_hatchery_num, coho1_wild_num, trap, Date, daytime_category, temp,
         flow, lunar_phase, photoperiod, photo_diff, temp_diff, flow_diff, atu_solstice,
         doy, year, In) %>% 
  group_by(Date) %>%
  summarize(coho1_hatchery_num = if(all(is.na(coho1_hatchery_num))) NA else sum(coho1_hatchery_num, na.rm = TRUE),
            coho1_wild_num = if(all(is.na(coho1_wild_num))) NA else sum(coho1_wild_num, na.rm = TRUE),
            In = if(all(is.na(In))) NA else sum(In, na.rm = TRUE),
            temp = mean(temp),
            flow = mean(flow),
            lunar_phase = mean(lunar_phase),
            photoperiod = mean(photoperiod),
            photo_diff = mean(photo_diff),
            temp_diff = mean(temp_diff),
            flow_diff = mean(flow_diff),
            atu_solstice = mean(atu_solstice),
            doy = first(doy),
            year = first(year)) %>%
  ungroup()

data %>% 
  ggplot(aes(x = doy, y = coho1_hatchery_num)) +
  geom_point() +
  geom_line() +
  facet_wrap(~year, scales = "free_y")  

data <- data %>% 
  mutate(coho1_hatchery_perhour = coho1_hatchery_num/In,
         coho1_wild_perhour = coho1_wild_num/In,
         coho1_hatchery_perhour_inp = na.approx(coho1_hatchery_perhour, na.rm = FALSE),
         temp = na.approx(temp),
         flow = na.approx(flow),
         lunar_phase = na.approx(lunar_phase),
         photoperiod = na.approx(photoperiod),
         photo_diff = na.approx(photo_diff),
         temp_diff = na.approx(temp_diff),
         flow_diff = na.approx(flow_diff),
         atu_solstice = na.approx(atu_solstice))

ggplot(data) +
  geom_line(aes(x = doy, y = coho1_hatchery_perhour), color = "cadetblue") +
  geom_line(aes(x = doy, y = coho1_wild_perhour), color = "salmon") +
  facet_wrap(~year, scales = "free_y") +
  #add dotted line for coho1_hatchery_perhour_inp
  geom_line(aes(x = doy, y = coho1_hatchery_perhour_inp), color = "cadetblue", linetype = "dotted")


lm_temp_day <- lm(temp ~ 1+photoperiod, data = data)

data <- data %>% add_residuals(lm_temp_day)



covariates_coho1_skagit_combined <- arrange(data,doy) %>%
  filter(doy >100 & doy < 150) %>%
  dplyr::select(year,doy,  temp, flow, photoperiod, atu_solstice,
                lunar_phase, resid, temp_diff, flow_diff, photo_diff,
                coho1_hatchery_perhour_inp) %>%
  pivot_wider(names_from = c(year), values_from = c(
    temp, flow, photoperiod, atu_solstice, lunar_phase,
    resid, temp_diff, flow_diff, photo_diff, coho1_hatchery_perhour_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2022-2010
num_rows = num_years
num_covariates = 10
total_covariates = dim(covariates_coho1_skagit_combined)[1]

#scaling and centering
for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_coho1_skagit_combined[i,] = scale(covariates_coho1_skagit_combined[i,])[,1]
}

#just scale
for(i in (total_covariates/2 + 1):(total_covariates)){
  if(sum(covariates_coho1_skagit_combined[i,]) != 0){
    covariates_coho1_skagit_combined[i,] = scale(covariates_coho1_skagit_combined[i,], center = FALSE, scale= TRUE)[,1]
  }
  
  
}


#subset response variable
subset_coho_summer_perhour_combined <- arrange(data,doy) %>%
  filter(doy > 100 & doy < 150) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour_combined)[1]){
  subset_coho_summer_perhour_combined[i,] = scale(subset_coho_summer_perhour_combined[i,])[,1]
}

Cmat_skagit_year <- function(nyears,ncov,hatchery=0, year_effect = TRUE){
  vars = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")
  #make a list of years from 2010 to 2022 but each year has to be repeated 4 times
  # years = c(rep("2010",4),rep("2012",4),rep("2013",4),rep("2014",4),rep("2015",4),
  # rep("2016",4),rep("2017",4),rep("2018",4),rep("2019",4),rep("2020",4),rep("2021",4),rep("2022",4
  # ))  
  years = c("2010","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022")
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
              if(year_effect == TRUE){
                C[i,j] <- years[i]
              }
              
              else{
                C[i,j] <- vars[k]
              }
              
              
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


mod_list_skagit <- function(nyears,ncov,hatchery=0, year_effect = TRUE){
  
  Q = "diagonal and equal"
  
  
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
      
      C = Cmat_skagit_year(nyears,ncov,hatchery, year_effect)
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


#fit MARSS model with just hatchery and year effect
kk=10
nyears = (2022-2010)
c = covariates_coho1_skagit_combined[((1+(kk-1)*nyears):(kk*nyears)),]
name_long = rownames(covariates_coho1_skagit_combined)[1+(kk-1)*nyears]
name_individual = substr(name_long,1,nchar(name_long)-9)
print(name_individual)
fit.model = c(list(c= c), mod_list_skagit(nyears,1,1,TRUE))
fit_year_effect_coho <- MARSS(subset_coho_summer_perhour_combined, model=fit.model, silent = TRUE, method = "BFGS",
                         control=list(maxit=2000))
ci_year_effect_coho <- tidy(fit_year_effect_coho)

fit_year_effect_coho$AICc #2550



out.tab.all.years.combined_coho <- NULL
fits.all.years.combined_coho <- NULL
c <- NULL
nyears = num_years*2*2
for(kk in c(1,2,3,4,5,6,7,8,9)){
  c = covariates_coho1_skagit_combined[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_coho1_skagit_combined)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-5)
  print(name_individual)
  fit.model = c(list(c= c), mod_list_skagit(nyears,1,0, FALSE))
  fit <- MARSS(subset_coho_summer_perhour_combined, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years.combined_coho=rbind(out.tab.all.years.combined_coho,out)
  fits.all.years.combined_coho=c(fits.all.years.combined_coho,list(fit))
}

out.tab.all.years.combined_coho
#calculate delta AICc and order by that
out.tab.all.years.combined_coho$delta.AICc = out.tab.all.years.combined_coho$AICc - min(out.tab.all.years.combined_coho$AICc)
out.tab.all.years.combined_coho = out.tab.all.years.combined_coho[order(out.tab.all.years.combined_coho$delta.AICc),]
#round to 2 decimal places
out.tab.all.years.combined_coho$delta.AICc = round(out.tab.all.years.combined_coho$delta.AICc,2)
out.tab.all.years.combined_coho

#save the dataframe
write.csv(out.tab.all.years.combined, file = here("skagit", "output","model_selection_correlated_covariates_coho1.csv"))


#try model selection


list_combinations <- get_covariate_combinations(1:6)
out.tab.hatchery.all.years.combined_coho <- NULL
fits.hatchery.all.years.combined_coho <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  atu = 0
  flow = 0
  flow_difference = 0
  hatchery = 0
  temperature_difference = 0
  lunar_phase = 0
  for(j in covariates){
    if(j == 1){
      k = 4
      atu = 1
    }
    else if(j==2){
      k = 2
      flow = 1
    }
    else if(j==3){
      k = 5
      lunar_phase = 1
      
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
      k = 10
      hatchery = 1
    }
    
    
    c = rbind(c,covariates_coho1_skagit_combined[((1+(k-1)*nyears):(k*nyears)),])
    name_long = rownames(covariates_coho1_skagit_combined)[1+(k-1)*nyears]
    name = paste(name, substr(name_long,1,nchar(name_long)-5))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(k==10){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list_skagit(nyears,c_num,has_hatchery, FALSE ))
  }
  else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list_skagit(nyears,c_num,has_hatchery, FALSE))
  }
  fit <- MARSS(subset_coho_summer_perhour_combined, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, atu=atu, 
                 flow=flow, flow_difference=flow_difference,
                 temperature_difference=temperature_difference,
                 lunar_phase = lunar_phase ,hatchery=hatchery, 
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.hatchery.all.years.combined_coho=rbind(out.tab.hatchery.all.years.combined_coho,out)
  fits.hatchery.all.years.combined_coho=c(fits.hatchery.all.years.combined_coho,list(fit))
  
  
}
out.tab.hatchery.all.years.combined_coho$deltaAICc <- out.tab.hatchery.all.years.combined_coho$AICc - min(out.tab.hatchery.all.years.combined_coho$AICc)
min.AICc <- order(out.tab.hatchery.all.years.combined_coho$AICc)
out.tab.hatchery.ordered.combined_coho <- out.tab.hatchery.all.years.combined_coho[min.AICc, ]
out.tab.hatchery.ordered.combined_coho

tidy(fits.hatchery.all.years.combined_coho[[6]])

tidy(fits.hatchery.all.years.combined_coho[[23]])

autoplot(fits.hatchery.all.years.combined_coho[[6]])

autoplot(fits.hatchery.all.years.combined_coho[[23]])
#i'm gonna change the doy to 100 to 150 to see if results change

#keep scoop and screw separate

#drop trap and daytime_category
data_trap <- data_day_night %>% 
  select(coho1_hatchery_num, coho1_wild_num, trap, Date, daytime_category, temp,
         flow, lunar_phase, photoperiod, photo_diff, temp_diff, flow_diff, atu_solstice,
         doy, year, In) %>% 
  group_by(Date, doy, year, trap) %>%
  summarize(coho1_hatchery_num = if(all(is.na(coho1_hatchery_num))) NA else sum(coho1_hatchery_num, na.rm = TRUE),
            coho1_wild_num = if(all(is.na(coho1_wild_num))) NA else sum(coho1_wild_num, na.rm = TRUE),
            In = if(all(is.na(In))) NA else sum(In, na.rm = TRUE),
            temp = mean(temp),
            flow = mean(flow),
            lunar_phase = mean(lunar_phase),
            photoperiod = mean(photoperiod),
            photo_diff = mean(photo_diff),
            temp_diff = mean(temp_diff),
            flow_diff = mean(flow_diff),
            atu_solstice = mean(atu_solstice)) %>%
  ungroup()


data_trap <- data_trap %>% 
  mutate(coho1_hatchery_perhour = coho1_hatchery_num/In,
         coho1_wild_perhour = coho1_wild_num/In,
         coho1_hatchery_perhour_inp = na.approx(coho1_hatchery_perhour, na.rm = FALSE),
         temp = na.approx(temp),
         flow = na.approx(flow),
         lunar_phase = na.approx(lunar_phase),
         photoperiod = na.approx(photoperiod),
         photo_diff = na.approx(photo_diff),
         temp_diff = na.approx(temp_diff),
         flow_diff = na.approx(flow_diff),
         atu_solstice = na.approx(atu_solstice))

ggplot(data_trap) +
  geom_line(aes(x = doy, y = coho1_hatchery_perhour), color = "cadetblue") +
  geom_line(aes(x = doy, y = coho1_wild_perhour), color = "salmon") +
  facet_wrap(~year+trap, scales = "free_y") +
  #add dotted line for coho1_hatchery_perhour_inp
  geom_line(aes(x = doy, y = coho1_hatchery_perhour_inp), color = "cadetblue", linetype = "dotted")


lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_trap)

data_trap <- data_trap %>% add_residuals(lm_temp_day)

#changing doy limit from 150 to 175
covariates_coho1_skagit_combined_trap <- arrange(data_trap,doy) %>%
  filter(doy >100 & doy < 175) %>%
  dplyr::select(year,doy, trap, temp, flow, photoperiod, atu_solstice,
                lunar_phase, resid, temp_diff, flow_diff, photo_diff,
                coho1_hatchery_perhour_inp) %>%
  pivot_wider(names_from = c(year,trap), values_from = c(
    temp, flow, photoperiod, atu_solstice, lunar_phase,
    resid, temp_diff, flow_diff, photo_diff, coho1_hatchery_perhour_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2022-2010
num_rows = num_years*2
num_covariates = 10
total_covariates = dim(covariates_coho1_skagit_combined_trap)[1]

#scaling and centering
for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_coho1_skagit_combined_trap[i,] = scale(covariates_coho1_skagit_combined_trap[i,])[,1]
}

#just scale
for(i in (total_covariates/2 + 1):(total_covariates)){
  if(sum(covariates_coho1_skagit_combined_trap[i,]) != 0){
    covariates_coho1_skagit_combined_trap[i,] = scale(covariates_coho1_skagit_combined_trap[i,], center = FALSE, scale= TRUE)[,1]
  }
  
  
}


#subset response variable
subset_coho_summer_perhour_combined_trap <- arrange(data_trap,doy) %>%
  filter(doy > 100 & doy < 175) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,trap) %>%
  pivot_wider(names_from = c(year,trap), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour_combined_trap)[1]){
  subset_coho_summer_perhour_combined_trap[i,] = scale(subset_coho_summer_perhour_combined_trap[i,])[,1]
}

Cmat_skagit_year <- function(nyears,ncov,hatchery=0, year_effect = TRUE){
  vars = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")
  #make a list of years from 2010 to 2022 but each year has to be repeated 4 times
  # years = c(rep("2010",4),rep("2012",4),rep("2013",4),rep("2014",4),rep("2015",4),
  # rep("2016",4),rep("2017",4),rep("2018",4),rep("2019",4),rep("2020",4),rep("2021",4),rep("2022",4
  # ))  
  years = c("2010","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022")
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
              if(year_effect == TRUE){
                C[i,j] <- years[i]
              }
              
              else{
                C[i,j] <- vars[k]
              }
              
              
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


mod_list_skagit <- function(nyears,ncov,hatchery=0, year_effect = TRUE){
  
  Q = "diagonal and equal"
  
  
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
      
      C = Cmat_skagit_year(nyears,ncov,hatchery, year_effect)
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



out.tab.all.years.combined_trap_coho <- NULL
fits.all.years.combined_trap_coho <- NULL
nyears = num_years*2
for(kk in c(1,2,3,4,5,6,7,8,9,10)){
  c = covariates_coho1_skagit_combined_trap[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_coho1_skagit_combined_trap)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-11)
  print(name_individual)
  fit.model = c(list(c= c), mod_list_skagit(nyears,1,0, FALSE))
  fit <- MARSS(subset_coho_summer_perhour_combined_trap, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years.combined_trap_coho=rbind(out.tab.all.years.combined_trap_coho,out)
  fits.all.years.combined_trap_coho=c(fits.all.years.combined_trap_coho,list(fit))
}

out.tab.all.years.combined_trap_coho
#calculate delta AICc and order by that
out.tab.all.years.combined_trap_coho$delta.AICc = out.tab.all.years.combined_trap_coho$AICc - min(out.tab.all.years.combined_trap_coho$AICc)
out.tab.all.years.combined_trap_coho = out.tab.all.years.combined_trap_coho[order(out.tab.all.years.combined_trap_coho$delta.AICc),]
#round to 2 decimal places
out.tab.all.years.combined_trap_coho$delta.AICc = round(out.tab.all.years.combined_trap$delta.AICc,2)
out.tab.all.years.combined_trap_coho


#try model selection


list_combinations <- get_covariate_combinations(1:6)
out.tab.hatchery.all.years.combined_trap_coho <- NULL
fits.hatchery.all.years.combined_trap_coho <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  temperature = 0
  flow = 0
  flow_difference = 0
  hatchery = 0
  temperature_difference = 0
  lunar_phase = 0
  for(j in covariates){
    if(j == 1){
      k = 1
      temperature = 1
    }
    else if(j==2){
      k = 2
      flow = 1
    }
    else if(j==3){
      k = 5
      lunar_phase = 1
      
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
      k = 10
      hatchery = 1
    }
    
    
    c = rbind(c,covariates_coho1_skagit_combined_trap[((1+(k-1)*nyears):(k*nyears)),])
    name_long = rownames(covariates_coho1_skagit_combined_trap)[1+(k-1)*nyears]
    name = paste(name, substr(name_long,1,nchar(name_long)-11))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(k==10){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list_skagit(nyears,c_num,has_hatchery, FALSE ))
  }
  else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list_skagit(nyears,c_num,has_hatchery, FALSE))
  }
  fit <- MARSS(subset_coho_summer_perhour_combined_trap, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, temperature=temperature,
                 flow=flow, flow_difference=flow_difference,
                 temperature_difference=temperature_difference,
                 lunar_phase = lunar_phase ,hatchery=hatchery, 
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.hatchery.all.years.combined_trap_coho=rbind(out.tab.hatchery.all.years.combined_trap_coho,out)
  fits.hatchery.all.years.combined_trap_coho=c(fits.hatchery.all.years.combined_trap_coho,list(fit))
  
  
}
out.tab.hatchery.all.years.combined_trap_coho$deltaAICc <- out.tab.hatchery.all.years.combined_trap_coho$AICc - min(out.tab.hatchery.all.years.combined_trap_coho$AICc)
min.AICc <- order(out.tab.hatchery.all.years.combined_trap_coho$AICc)
out.tab.hatchery.ordered.combined_trap_coho <- out.tab.hatchery.all.years.combined_trap_coho[min.AICc, ]
out.tab.hatchery.ordered.combined_trap_coho

#combine trap data but keep day night separate and plot the data
data_day <- data_day_night %>% 
  select(coho1_hatchery_num, coho1_wild_num, trap, Date, daytime_category, temp,
         flow, lunar_phase, photoperiod, photo_diff, temp_diff, flow_diff, atu_solstice,
         doy, year, In) %>% 
  group_by(Date, doy, year, daytime_category) %>%
  summarize(coho1_hatchery_num = if(all(is.na(coho1_hatchery_num))) NA else sum(coho1_hatchery_num, na.rm = TRUE),
            coho1_wild_num = if(all(is.na(coho1_wild_num))) NA else sum(coho1_wild_num, na.rm = TRUE),
            In = if(all(is.na(In))) NA else sum(In, na.rm = TRUE),
            temp = mean(temp),
            flow = mean(flow),
            lunar_phase = mean(lunar_phase),
            photoperiod = mean(photoperiod),
            photo_diff = mean(photo_diff),
            temp_diff = mean(temp_diff),
            flow_diff = mean(flow_diff),
            atu_solstice = mean(atu_solstice)) %>%
  ungroup()

data_day <- data_day %>% 
  mutate(coho1_hatchery_perhour = coho1_hatchery_num/In,
         coho1_wild_perhour = coho1_wild_num/In,
         coho1_hatchery_perhour_inp = na.approx(coho1_hatchery_perhour, na.rm = FALSE),
         temp = na.approx(temp),
         flow = na.approx(flow),
         lunar_phase = na.approx(lunar_phase),
         photoperiod = na.approx(photoperiod),
         photo_diff = na.approx(photo_diff),
         temp_diff = na.approx(temp_diff),
         flow_diff = na.approx(flow_diff),
         atu_solstice = na.approx(atu_solstice))

ggplot(data_day) +
  geom_line(aes(x = doy, y = coho1_hatchery_perhour), color = "cadetblue") +
  geom_line(aes(x = doy, y = coho1_wild_perhour), color = "salmon") +
  facet_wrap(~year+daytime_category, scales = "free_y") +
  #add dotted line for coho1_hatchery_perhour_inp
  geom_line(aes(x = doy, y = coho1_hatchery_perhour_inp), color = "cadetblue", linetype = "dotted")

#it does not make sense to have values by day and night separate, so combine them

#for now I want to keep daytime category and trap separate and for doy 100 to 150 
#just to see if I can recreate initial results or whether that was a mistake

data <- data_day_night %>% 
  select(coho1_hatchery_num, coho1_wild_num, trap, Date, daytime_category, temp,
         flow, lunar_phase, photoperiod, photo_diff, temp_diff, flow_diff, atu_solstice,
         doy, year, In) %>% 
  group_by(Date,  trap, daytime_category, doy, year) %>%
  summarize(coho1_hatchery_num = if(all(is.na(coho1_hatchery_num))) NA else sum(coho1_hatchery_num, na.rm = TRUE),
            coho1_wild_num = if(all(is.na(coho1_wild_num))) NA else sum(coho1_wild_num, na.rm = TRUE),
            In = if(all(is.na(In))) NA else sum(In, na.rm = TRUE),
            temp = mean(temp),
            flow = mean(flow),
            lunar_phase = mean(lunar_phase),
            photoperiod = mean(photoperiod),
            photo_diff = mean(photo_diff),
            temp_diff = mean(temp_diff),
            flow_diff = mean(flow_diff),
            atu_solstice = mean(atu_solstice)) %>%
  ungroup()

data %>% 
  ggplot(aes(x = doy, y = coho1_hatchery_num)) +
  geom_point() +
  geom_line() +
  facet_wrap(~year, scales = "free_y")  

data <- data %>% 
  mutate(coho1_hatchery_perhour = coho1_hatchery_num/In,
         coho1_wild_perhour = coho1_wild_num/In,
         coho1_hatchery_perhour_inp = na.approx(coho1_hatchery_perhour, na.rm = FALSE),
         temp = na.approx(temp),
         flow = na.approx(flow),
         lunar_phase = na.approx(lunar_phase),
         photoperiod = na.approx(photoperiod),
         photo_diff = na.approx(photo_diff),
         temp_diff = na.approx(temp_diff),
         flow_diff = na.approx(flow_diff),
         atu_solstice = na.approx(atu_solstice))


lm_temp_day <- lm(temp ~ 1+photoperiod, data = data)

data <- data %>% add_residuals(lm_temp_day)



covariates_coho1_skagit_combined <- arrange(data,doy) %>%
  filter(doy >100 & doy < 150) %>%
  dplyr::select(year,doy, trap, daytime_category, temp, flow, photoperiod, atu_solstice,
                lunar_phase, resid, temp_diff, flow_diff, photo_diff,
                coho1_hatchery_perhour_inp) %>%
  pivot_wider(names_from = c(year, trap, daytime_category), values_from = c(
    temp, flow, photoperiod, atu_solstice, lunar_phase,
    resid, temp_diff, flow_diff, photo_diff, coho1_hatchery_perhour_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2022-2010
num_rows = num_years*2*2
num_covariates = 10
total_covariates = dim(covariates_coho1_skagit_combined)[1]

#scaling and centering
for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_coho1_skagit_combined[i,] = scale(covariates_coho1_skagit_combined[i,])[,1]
}

#just scale
for(i in (total_covariates/2 + 1):(total_covariates)){
  if(sum(covariates_coho1_skagit_combined[i,]) != 0){
    covariates_coho1_skagit_combined[i,] = scale(covariates_coho1_skagit_combined[i,], center = FALSE, scale= TRUE)[,1]
  }
  
  
}


#subset response variable
subset_coho_summer_perhour_combined <- arrange(data,doy) %>%
  filter(doy > 100 & doy < 150) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy, trap, daytime_category) %>%
  pivot_wider(names_from = c(year, trap, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour_combined)[1]){
  subset_coho_summer_perhour_combined[i,] = scale(subset_coho_summer_perhour_combined[i,])[,1]
}

Cmat_skagit_year <- function(nyears,ncov,hatchery=0, year_effect = TRUE){
  vars = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")
  #make a list of years from 2010 to 2022 but each year has to be repeated 4 times
  # years = c(rep("2010",4),rep("2012",4),rep("2013",4),rep("2014",4),rep("2015",4),
  # rep("2016",4),rep("2017",4),rep("2018",4),rep("2019",4),rep("2020",4),rep("2021",4),rep("2022",4
  # ))  
  years = c("2010","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022")
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
              if(year_effect == TRUE){
                C[i,j] <- years[i]
              }
              
              else{
                C[i,j] <- vars[k]
              }
              
              
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


mod_list_skagit <- function(nyears,ncov,hatchery=0, year_effect = TRUE){
  
  Q = "diagonal and equal"
  
  
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
      
      C = Cmat_skagit_year(nyears,ncov,hatchery, year_effect)
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


out.tab.all.years.combined_coho <- NULL
fits.all.years.combined_coho <- NULL
nyears = num_years*2*2
for(kk in c(1,2,3,4,5,6,7,8,9)){
  c = covariates_coho1_skagit_combined[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_coho1_skagit_combined)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-5)
  print(name_individual)
  fit.model = c(list(c= c), mod_list_skagit(nyears,1,0, FALSE))
  fit <- MARSS(subset_coho_summer_perhour_combined, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years.combined_coho=rbind(out.tab.all.years.combined_coho,out)
  fits.all.years.combined_coho=c(fits.all.years.combined_coho,list(fit))
}


#calculate delta AICc and order by that
out.tab.all.years.combined_coho$delta.AICc = out.tab.all.years.combined_coho$AICc - min(out.tab.all.years.combined_coho$AICc)
out.tab.all.years.combined_coho = out.tab.all.years.combined_coho[order(out.tab.all.years.combined_coho$delta.AICc),]
#round to 2 decimal places
out.tab.all.years.combined_coho


list_combinations <- get_covariate_combinations(1:6)
out.tab.hatchery.all.years.combined_coho <- NULL
fits.hatchery.all.years.combined_coho <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  atu = 0
  flow = 0
  flow_difference = 0
  hatchery = 0
  temperature_difference = 0
  lunar_phase = 0
  for(j in covariates){
    if(j == 1){
      k = 4
      atu = 1
    }
    else if(j==2){
      k = 2
      flow = 1
    }
    else if(j==3){
      k = 5
      lunar_phase = 1
      
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
      k = 10
      hatchery = 1
    }
    
    
    c = rbind(c,covariates_coho1_skagit_combined[((1+(k-1)*nyears):(k*nyears)),])
    name_long = rownames(covariates_coho1_skagit_combined)[1+(k-1)*nyears]
    name = paste(name, substr(name_long,1,nchar(name_long)-5))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(k==10){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list_skagit(nyears,c_num,has_hatchery, FALSE ))
  }
  else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list_skagit(nyears,c_num,has_hatchery, FALSE))
  }
  fit <- MARSS(subset_coho_summer_perhour_combined, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, atu=atu,
                 flow=flow, flow_difference=flow_difference,
                 temperature_difference=temperature_difference,
                 lunar_phase = lunar_phase ,hatchery=hatchery, 
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.hatchery.all.years.combined_coho=rbind(out.tab.hatchery.all.years.combined_coho,out)
  fits.hatchery.all.years.combined_coho=c(fits.hatchery.all.years.combined_coho,list(fit))
  
  
}
out.tab.hatchery.all.years.combined_coho$deltaAICc <- out.tab.hatchery.all.years.combined_coho$AICc - min(out.tab.hatchery.all.years.combined_coho$AICc)
min.AICc <- order(out.tab.hatchery.all.years.combined_coho$AICc)
out.tab.hatchery.ordered.combined_coho <- out.tab.hatchery.all.years.combined_coho[min.AICc, ]
out.tab.hatchery.ordered.combined_coho
