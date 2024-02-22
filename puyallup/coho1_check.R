# goal - check the results of the model with all the 
# covariates from the best model but with fewer days in the
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



covariates_coho1_puyallup_w_temp_short <- arrange(data,doy) %>%
  filter(doy > 110 & doy <= 160) %>%
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
total_covariates = dim(covariates_coho1_puyallup_w_temp_short)[1]


for(i in 1:(num_rows*6)){ # everything except diffs and hatchery
  covariates_coho1_puyallup_w_temp_short[i,] = scale(covariates_coho1_puyallup_w_temp_short[i,])[,1]
}

#just scale

for(i in (num_rows*6 + 1):(total_covariates)){
  covariates_coho1_puyallup_w_temp_short[i,] = scale(covariates_coho1_puyallup_w_temp_short[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_coho_summer_perhour_short <- arrange(data,doy) %>%
  filter(doy > 110 & doy <= 160) %>%
  mutate(log.value_day = log(coho1_wild_perhour_day + 1), 
         log.value_night = log(coho1_wild_perhour_night + 1)) %>%
  dplyr::select(log.value_day, log.value_night ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_day, log.value_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour_short)[1]){
  subset_coho_summer_perhour_short[i,] = scale(subset_coho_summer_perhour_short[i,])[,1]
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

fit.model_equal_q = mod_list(nyears,0,0, FALSE, FALSE)
fit_equal_q <- MARSS(subset_coho_summer_perhour_short, model=fit.model_equal_q, silent = TRUE, method = "BFGS",
                     control=list(maxit=2000))
fit_equal_q$AICc #4089
fit.model_unequal_q = mod_list(nyears,0,0, FALSE, TRUE)
fit_unequal_q <- MARSS(subset_coho_summer_perhour_short, model=fit.model_unequal_q, silent = TRUE, method = "BFGS",
                       control=list(maxit=2000))
fit_unequal_q$AICc #4090
#unequal is better or equal


autoplot(fit_equal_q)


out.tab.all.years.coho_short <- NULL
fits.all.years.coho_short <- NULL
nyears = num_years*2
for(kk in c(2,3,5,6,8,10)){
  c = covariates_coho1_puyallup_w_temp_short[((1+(kk-1)*num_rows):(kk*num_rows)),]
  name_long = rownames(covariates_coho1_puyallup_w_temp_short)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0, FALSE, FALSE))
  fit <- MARSS(subset_coho_summer_perhour_short, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years.coho_short=rbind(out.tab.all.years.coho_short,out)
  fits.all.years.coho_short=c(fits.all.years.coho_short,list(fit))
}

out.tab.all.years.coho_short$delta_AICc = out.tab.all.years.coho_short$AICc - min(out.tab.all.years.coho_short$AICc)
out.tab.all.years.coho_short_ordered = out.tab.all.years.coho_short[order(out.tab.all.years.coho_short$delta_AICc),]
out.tab.all.years.coho_short_ordered

#photoperiod and temp - photo
#secchi_depth and photo diff - secchi
#secchi depth and temp - secchi
#photo diff and resid - resid
#temp and resid - resid
#photo_diff and #atu - atu

c<-NULL
for(kk in c(1,3,6,7,11)){
  c = rbind(c,covariates_coho1_puyallup_w_temp_short[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_coho1_puyallup_w_temp_short)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  
}
fit.model_check = c(list(c= c), mod_list(nyears,5,1, FALSE, FALSE))
fit_check <- MARSS(subset_coho_summer_perhour_short, model=fit.model_check, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

autoplot(fit_check)

ci_check <- tidy(fit_check)
ci_check
#plot
ggplot(ci_check[39:44,], aes(x = c("ATU", "Flow", 
                                   "Photoperiod",
                                   "Flow\n difference",
                                   "Hatchery,\n day", "Hatchery,\n night"), 
                             y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Coho Puyallup",
       x = "Covariate",
       y = "Estimate") +
  theme_minimal()


## I think I should keep just the night values but for the whole duration



covariates_coho1_puyallup_w_temp_night <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160) %>%
  dplyr::select(year,doy, flow_night, 
                photoperiod_night, secchi_depth_night,
                lunar_phase_night, temp_night, atu_night,
                flow_diff_night, photo_diff_night, 
                temp_diff_night, 
                resid_night,
                coho1_hatchery_perhour_night) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_night, secchi_depth_night,
    photoperiod_night,
    lunar_phase_night, temp_night, atu_night,
    flow_diff_night, photo_diff_night, 
    temp_diff_night,
    resid_night,
    coho1_hatchery_perhour_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2021-2004+1
num_rows = num_years

total_covariates = dim(covariates_coho1_puyallup_w_temp_night)[1]


for(i in 1:(num_rows*6)){ # everything except diffs and hatchery
  covariates_coho1_puyallup_w_temp_night[i,] = scale(covariates_coho1_puyallup_w_temp_night[i,])[,1]
}

#just scale

for(i in (num_rows*6 + 1):(total_covariates)){
  covariates_coho1_puyallup_w_temp_night[i,] = scale(covariates_coho1_puyallup_w_temp_night[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_coho_summer_perhour_night <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160) %>%
  mutate(log.value_night = log(coho1_wild_perhour_night + 1)) %>%
  dplyr::select(log.value_night ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour_night)[1]){
  subset_coho_summer_perhour_night[i,] = scale(subset_coho_summer_perhour_night[i,])[,1]
}



c<-NULL
for(kk in c(1,3,6,7,11)){
  c = rbind(c,covariates_coho1_puyallup_w_temp_night[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_coho1_puyallup_w_temp_night)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  
}
fit.model_check_night = c(list(c= c), mod_list(num_rows,5,0, FALSE, FALSE))
fit_check_night <- MARSS(subset_coho_summer_perhour_night, model=fit.model_check_night, silent = TRUE, method = "BFGS",
                   control=list(maxit=2000))

autoplot(fit_check_night)

ci_check <- tidy(fit_check_night)
ci_check

#keepingjust values for a shorter duration



covariates_coho1_puyallup_w_temp_night <- arrange(data,doy) %>%
  filter(doy > 110 & doy <= 160) %>%
  dplyr::select(year,doy, flow_night, 
                photoperiod_night, secchi_depth_night,
                lunar_phase_night, temp_night, atu_night,
                flow_diff_night, photo_diff_night, 
                temp_diff_night, 
                resid_night,
                coho1_hatchery_perhour_night) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_night, secchi_depth_night,
    photoperiod_night,
    lunar_phase_night, temp_night, atu_night,
    flow_diff_night, photo_diff_night, 
    temp_diff_night,
    resid_night,
    coho1_hatchery_perhour_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2021-2004+1
num_rows = num_years

total_covariates = dim(covariates_coho1_puyallup_w_temp_night)[1]


for(i in 1:(num_rows*6)){ # everything except diffs and hatchery
  covariates_coho1_puyallup_w_temp_night[i,] = scale(covariates_coho1_puyallup_w_temp_night[i,])[,1]
}

#just scale

for(i in (num_rows*6 + 1):(total_covariates)){
  covariates_coho1_puyallup_w_temp_night[i,] = scale(covariates_coho1_puyallup_w_temp_night[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_coho_summer_perhour_night <- arrange(data,doy) %>%
  filter(doy > 110 & doy <= 160) %>%
  mutate(log.value_night = log(coho1_wild_perhour_night + 1)) %>%
  dplyr::select(log.value_night ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour_night)[1]){
  subset_coho_summer_perhour_night[i,] = scale(subset_coho_summer_perhour_night[i,])[,1]
}



c<-NULL
for(kk in c(1,3,6,7,11)){
  c = rbind(c,covariates_coho1_puyallup_w_temp_night[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_coho1_puyallup_w_temp_night)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  
}
fit.model_check_night = c(list(c= c), mod_list(num_rows,5,0, FALSE, FALSE))
fit_check_night <- MARSS(subset_coho_summer_perhour_night, model=fit.model_check_night, silent = TRUE, method = "BFGS",
                         control=list(maxit=2000))

autoplot(fit_check_night)

ci_check <- tidy(fit_check_night)
ci_check

#susbet fewer years and start at doy 90





covariates_coho1_puyallup_w_temp_years <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160, year>=2011) %>%
  dplyr::select(year,doy, flow_night, 
                photoperiod_night, secchi_depth_night,
                lunar_phase_night, temp_night, atu_night,
                flow_diff_night, photo_diff_night, 
                temp_diff_night, 
                resid_night,
                coho1_hatchery_perhour_night) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_night, secchi_depth_night,
    photoperiod_night,
    lunar_phase_night, temp_night, atu_night,
    flow_diff_night, photo_diff_night, 
    temp_diff_night,
    resid_night,
    coho1_hatchery_perhour_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()



num_years = 2021-2011+1
num_rows = num_years

total_covariates = dim(covariates_coho1_puyallup_w_temp_years)[1]


for(i in 1:(num_rows*6)){ # everything except diffs and hatchery
  covariates_coho1_puyallup_w_temp_years[i,] = scale(covariates_coho1_puyallup_w_temp_years[i,])[,1]
}

#just scale

for(i in (num_rows*6 + 1):(total_covariates)){
  covariates_coho1_puyallup_w_temp_years[i,] = scale(covariates_coho1_puyallup_w_temp_years[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_coho_summer_perhour_years <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160, year >=2011) %>%
  mutate(log.value_night = log(coho1_wild_perhour_night + 1)) %>%
  dplyr::select(log.value_night ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour_years)[1]){
  subset_coho_summer_perhour_years[i,] = scale(subset_coho_summer_perhour_years[i,])[,1]
}



c<-NULL
for(kk in c(1,3,6,7,11)){
  c = rbind(c,covariates_coho1_puyallup_w_temp_years[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_coho1_puyallup_w_temp_years)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  
}
fit.model_check_years = c(list(c= c), mod_list(num_rows,5,0, FALSE, FALSE))
fit_check_years <- MARSS(subset_coho_summer_perhour_years, model=fit.model_check_years, silent = TRUE, method = "BFGS",
                         control=list(maxit=2000))

autoplot(fit_check_years)

ci_check <- tidy(fit_check_years)
ci_check

#model residuals are better. I will try with the same years but both day
#and night now


covariates_coho1_puyallup_w_temp_years2 <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160, year > 2010) %>%
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


num_years = 2021-2011+1
num_rows = num_years*2 #day and night

total_covariates = dim(covariates_coho1_puyallup_w_temp_years2)[1]


for(i in 1:(num_rows*6)){ # everything except diffs and hatchery
  covariates_coho1_puyallup_w_temp_years2[i,] = scale(covariates_coho1_puyallup_w_temp_years2[i,])[,1]
}

#just scale

for(i in (num_rows*6 + 1):(total_covariates)){
  covariates_coho1_puyallup_w_temp_years2[i,] = scale(covariates_coho1_puyallup_w_temp_years2[i,], center = FALSE, scale= TRUE)[,1]
}




#subset response variable
subset_coho_summer_perhour_years2 <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160, year >=2011) %>%
  mutate(log.value_day = log(coho1_wild_perhour_day + 1),
         log.value_night = log(coho1_wild_perhour_night + 1)
         
         ) %>%
  dplyr::select(log.value_day,log.value_night , 
                year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_day,
                                                    log.value_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour_years2)[1]){
  subset_coho_summer_perhour_years2[i,] = scale(subset_coho_summer_perhour_years2[i,])[,1]
}


c<-NULL
for(kk in c(1,3,6,7,11)){
  c = rbind(c,covariates_coho1_puyallup_w_temp_years2[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_coho1_puyallup_w_temp_years2)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  
}
fit.model_check_years2 = c(list(c= c), mod_list(num_rows,5,1, FALSE, FALSE))
fit_check_years2 <- MARSS(subset_coho_summer_perhour_years2, model=fit.model_check_years2, silent = TRUE, method = "BFGS",
                         control=list(maxit=2000))


fit.model_check_years2 = c(list(c= c), mod_list(num_rows,5,1, FALSE, TRUE))
fit_check_years2 <- MARSS(subset_coho_summer_perhour_years2, model=fit.model_check_years2, silent = TRUE, method = "BFGS",
                          control=list(maxit=2000))



autoplot(fit_check_years2)

ci_check <- tidy(fit_check_years2)
ci_check

#just look at day



covariates_coho1_puyallup_w_temp_years <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160, year>=2011) %>%
  dplyr::select(year,doy, flow_day, 
                photoperiod_day, secchi_depth_day,
                lunar_phase_day, temp_day, atu_day,
                flow_diff_day, photo_diff_day, 
                temp_diff_day, 
                resid_day,
                coho1_hatchery_perhour_day) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_day, secchi_depth_day,
    photoperiod_day,
    lunar_phase_day, temp_day, atu_day,
    flow_diff_day, photo_diff_day, 
    temp_diff_day,
    resid_day,
    coho1_hatchery_perhour_day)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()



num_years = 2021-2011+1
num_rows = num_years

total_covariates = dim(covariates_coho1_puyallup_w_temp_years)[1]


for(i in 1:(num_rows*6)){ # everything except diffs and hatchery
  covariates_coho1_puyallup_w_temp_years[i,] = scale(covariates_coho1_puyallup_w_temp_years[i,])[,1]
}

#just scale

for(i in (num_rows*6 + 1):(total_covariates)){
  covariates_coho1_puyallup_w_temp_years[i,] = scale(covariates_coho1_puyallup_w_temp_years[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_coho_summer_perhour_years <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160, year >=2011) %>%
  mutate(log.value_day = log(coho1_wild_perhour_day + 1)) %>%
  dplyr::select(log.value_day ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_day)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour_years)[1]){
  subset_coho_summer_perhour_years[i,] = scale(subset_coho_summer_perhour_years[i,])[,1]
}



c<-NULL
for(kk in c(1,3,6,7,11)){
  c = rbind(c,covariates_coho1_puyallup_w_temp_years[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_coho1_puyallup_w_temp_years)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  
}
fit.model_check_years = c(list(c= c), mod_list(num_rows,5,0, FALSE, FALSE))
fit_check_years <- MARSS(subset_coho_summer_perhour_years, model=fit.model_check_years, silent = TRUE, method = "BFGS",
                         control=list(maxit=2000))

autoplot(fit_check_years)

ci_check <- tidy(fit_check_years)
ci_check
