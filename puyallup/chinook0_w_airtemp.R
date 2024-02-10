#goal - run marss analysis for chinook0 in Puyallup
#with air temperature data

#load libraries
library(MARSS)
library(ggplot2)
library(here)
library(GGally)
library(tidyverse)
library(zoo)
library(modelr)


#load data
data <- read.csv(here("puyallup", "data","puyallup_2004-2021_all_days_w_covariates_w_airtemp.csv"),
                 na.strings = c("NA",""))


data$chinook0_hatchery_perhour_night[dim(data)[1]] <- 0
data$chinook0_hatchery_perhour_day[dim(data)[1]] <- 0

#when secchi_depth_day is NA, replace with secchi_depth values
#when secchi_depth_night is NA, replace with secchi_depth values

data$secchi_depth_day[is.na(data$secchi_depth_day)] <- data$secchi_depth[is.na(data$secchi_depth_day)]

data$secchi_depth_night[is.na(data$secchi_depth_night)] <- data$secchi_depth[is.na(data$secchi_depth_night)]



data <- data %>% 
  mutate(flow_day_inp = na.approx(flow_day, na.rm = FALSE),
         flow_night_inp = na.approx(flow_night, na.rm = FALSE),
         temp_day_inp = na.approx(temp, na.rm = FALSE),
         temp_night_inp = na.approx(temp, na.rm = FALSE),
         secchi_depth_day_inp = na.approx(secchi_depth_day, na.rm = FALSE),
         secchi_depth_night_inp = na.approx(secchi_depth_night, na.rm = FALSE),
         chinook0_hatchery_perhour_night_inp = na.approx(chinook0_hatchery_perhour_night, na.rm = FALSE),
         chinook0_hatchery_perhour_day_inp = na.approx(chinook0_hatchery_perhour_day, na.rm = FALSE),
         photoperiod_day = photoperiod,
         photoperiod_night = photoperiod,
         lunar_phase_day = lunar_phase,
         lunar_phase_night = lunar_phase)

#making columns for flow diff, photoperiod diff, temp_diff
data <- data %>%
  mutate(flow_diff_day = c(NA,diff(flow_day_inp)),
         flow_diff_night = c(NA,diff(flow_night_inp)),
         photo_diff_day = c(NA,diff(photoperiod_day)),
         photo_diff_night = c(NA,diff(photoperiod_night)),
         temp_diff_day = c(NA,diff(temp_day_inp)),
         temp_diff_night = c(NA,diff(temp_night_inp))
  )


lm_temp_day <- lm(temp_day_inp ~ 1+photoperiod_day, data = data)
summary(lm_temp_day)
lm_temp_night <- lm(temp_night_inp ~ 1+photoperiod_night, data = data)

data <- data %>% add_residuals(lm_temp_day, var = "resid_day") %>% 
  add_residuals(lm_temp_night, var = "resid_night")


data %>%
  dplyr::select(flow_day,photoperiod_day, lunar_phase_day, secchi_depth_day_inp, 
                photo_diff_day,flow_diff_day, temp_day_inp, temp_diff_day, resid_day) %>%
  GGally::ggpairs(upper=list(continuous='points'),
                  lower=list(continuous='cor'))
#photoperiod and temp
#secchi_depth and photo diff
#secchi depth and temp
#photo diff and resid
#temp and resid

#do same with night values

data %>%
  dplyr::select(flow_night,photoperiod_night, lunar_phase_night, secchi_depth_night_inp, 
                photo_diff_night,flow_diff_night, temp_night_inp, temp_diff_night, resid_night) %>%
  GGally::ggpairs(upper=list(continuous='points'),
                  lower=list(continuous='cor'))
#it is the same

covariates_chinook0_puyallup_w_temp <- arrange(data,doy) %>%
  filter(doy >130 & doy <= 218) %>%
  dplyr::select(year,doy, flow_day_inp, flow_night_inp, 
                photoperiod_day, photoperiod_night, secchi_depth_day_inp, secchi_depth_night_inp,
                lunar_phase_day, lunar_phase_night, temp_day_inp, temp_night_inp, 
                flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
                temp_diff_day, temp_diff_night, 
                resid_day, resid_night,
                chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_day_inp, flow_night_inp, secchi_depth_day_inp, secchi_depth_night_inp,
    photoperiod_day, photoperiod_night, lunar_phase_day, 
    lunar_phase_night, temp_day_inp, temp_night_inp, 
    flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
     temp_diff_day, temp_diff_night,
    resid_day, resid_night,
    chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2021-2004+1
num_rows = num_years*2
num_covariates = 14
total_covariates = dim(covariates_chinook0_puyallup_w_temp)[1]


for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_chinook0_puyallup_w_temp[i,] = scale(covariates_chinook0_puyallup_w_temp[i,])[,1]
}

#just scale

for(i in (total_covariates/2 + 1):(total_covariates)){
  covariates_chinook0_puyallup_w_temp[i,] = scale(covariates_chinook0_puyallup_w_temp[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_chinook_summer_perhour <- arrange(data,doy) %>%
  filter(doy > 130 & doy <= 218) %>%
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


out.tab.all.years <- NULL
fits.all.years <- NULL
nyears = num_years*2
for(kk in c(2,3,5,7,9)){
  c = covariates_chinook0_puyallup_w_temp[((1+(kk-1)*num_rows):(kk*num_rows)),]
  name_long = rownames(covariates_chinook0_puyallup_w_temp)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0, FALSE, FALSE))
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years=rbind(out.tab.all.years,out)
  fits.all.years=c(fits.all.years,list(fit))
}

out.tab.all.years$delta_AICc = out.tab.all.years$AICc - min(out.tab.all.years$AICc)
out.tab.all.years_ordered = out.tab.all.years[order(out.tab.all.years$delta_AICc),]
out.tab.all.years_ordered
#make column for delta aicc
#order by delta aicc
out.tab.all.years_ordered$delta_AICc = out.tab.all.years_ordered$AICc - min(out.tab.all.years_ordered$AICc)
out.tab.all.years_ordered
#save table
write.csv(out.tab.all.years_ordered, file = here("puyallup","output",
                                                 "chinook_correlated_covariates_selection.csv"))




#photoperiod difference and temperature is best

photoperiod_difference = 0
temperature = 0
flow = 0
lunar_phase = 0
hatchery = 0
flow_difference = 0
temperature_difference = 0

list_combinations <- get_covariate_combinations(1:7)
out.tab.hatchery.all.years <- NULL
fits.hatchery.all.years <- list()
for(i in 1:length(list_combinations)){
  photoperiod_difference = 0
  temperature = 0
  flow = 0
  lunar_phase = 0
  hatchery = 0
  flow_difference = 0
  temperature_difference = 0
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  for(j in covariates){
    if(j == 1){
      k = 7
      photoperiod_difference = 1
    }
    else if(j==2){
      k = 1
      flow = 1
    }
    else if(j==3){
      k = 5
      temperature = 1
    }
    else if(j==4){
      k = 4
      lunar_phase = 1
    }
    else if(j==5){
      k = 6
      flow_difference = 1
    }
    else if(j==6){
      k = 8
      temperature_difference = 1
    }
    else if(j==7){
      k = 10
      hatchery = 1
    }
    
    c = rbind(c,covariates_chinook0_puyallup_w_temp[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_chinook0_puyallup_w_temp)[1+(k-1)*num_rows]
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
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, photoperiod_difference = photoperiod_difference,
                 temperature = temperature, temperature_difference = temperature_difference,
                 flow = flow, lunar_phase = lunar_phase, flow_difference = flow_difference,
                 hatchery = hatchery,
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

fits.hatchery.all.years[[121]]

ci_puyallup_chinook_airtemp <- tidy(fits.hatchery.all.years[[121]])
ci_puyallup_chinook_airtemp


ggplot(ci_puyallup_chinook_airtemp[c(39:45),], 
       aes(x = c("Photoperiod difference", "Flow", "Temperature", "Lunar phase",
                 "Flow\n difference",
                 "Hatchery,\n day", "Hatchery,\n night"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect", title = "Puyallup, Chinook sub yearlings") +
  theme_bw()+
  theme(axis.text.x=element_text(size=18),axis.title.y=element_text(size=18),
        title = element_text(size = 18))


                        