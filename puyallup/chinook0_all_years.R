# Goal 

# try MARSS analysis with Puyallup data for chinook0
#all years from 2004-2021
#with and without separate error for day and night
#with and without day effect on night

#load packages

library(here)
library(MARSS)
library(ggplot2)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)
library(tidyverse)
library(GGally)



# load data

data <- read.csv(here("puyallup", "data","puyallup_2004-2021_all_days_w_covariates.csv"),
                 na.strings = c("NA",""))

#plotting chinook0_wild_perhour with doy to see timing of migration

ggplot(data, aes(x=doy, y=chinook0_wild_perhour_night)) +
  geom_point(alpha = 0.5) +
  labs(x = "Day of Year", y = "Chinook0 Wild per Hour") +
  theme_bw() #+ 
  #ylim(0, 20)

#plotting chinook0_hatchery_perhour with doy to see timing of migration
filter(data, chinook0_hatchery_perhour_night>1 & doy>200)
ggplot(data, aes(x=doy, y=chinook0_hatchery_perhour_night)) +
  geom_point(alpha = 0.5) +
  labs(x = "Day of Year", y = "Chinook0 hatchery per Hour") +
  theme_bw()  + 
  ylim(0, 500)

# I think the doy limits should be changed from 200 to 218 


#before interpolating, I need to add 0 to the end of 2018 chinook hatchery column
#if not, there are nans in the 2018 column even after interpolation

data$chinook0_hatchery_perhour_night[dim(data)[1]] <- 0
data$chinook0_hatchery_perhour_day[dim(data)[1]] <- 0

#when secchi_depth_day is NA, replace with secchi_depth values
#when secchi_depth_night is NA, replace with secchi_depth values

data$secchi_depth_day[is.na(data$secchi_depth_day)] <- data$secchi_depth[is.na(data$secchi_depth_day)]

data$secchi_depth_night[is.na(data$secchi_depth_night)] <- data$secchi_depth[is.na(data$secchi_depth_night)]

#interpolating na values for flow_day, flow_night, sechhi_depth_day, secchi_depth_night, hatchery

data <- data %>% 
  mutate(flow_day_inp = na.approx(flow_day, na.rm = FALSE),
         flow_night_inp = na.approx(flow_night, na.rm = FALSE),
         secchi_depth_day_inp = na.approx(secchi_depth_day, na.rm = FALSE),
         secchi_depth_night_inp = na.approx(secchi_depth_night, na.rm = FALSE),
         chinook0_hatchery_perhour_night_inp = na.approx(chinook0_hatchery_perhour_night, na.rm = FALSE),
         chinook0_hatchery_perhour_day_inp = na.approx(chinook0_hatchery_perhour_day, na.rm = FALSE),
         photoperiod_day = photoperiod,
         photoperiod_night = photoperiod,
         lunar_phase_day = lunar_phase,
         lunar_phase_night = lunar_phase)


#making columns for flow diff, photoperiod diff
data <- data %>%
  mutate(flow_diff_day = c(NA,diff(flow_day_inp)),
         flow_diff_night = c(NA,diff(flow_night_inp)),
         photo_diff_day = c(NA,diff(photoperiod_day)),
         photo_diff_night = c(NA,diff(photoperiod_night))
  )

#convert date$Date to as.Date and then get year
#first check where date string is not in the right format
#specify the format and then convert to as.Date


data$Date = as.Date(data$Date, format = "%m/%d/%Y")
data$year = year(data$Date)

# I think the doy limits should be changed from 200 to 218 
covariates_chinook0 <- arrange(data,doy) %>%
  filter(doy >130 & doy <= 218) %>%
  dplyr::select(year,doy, flow_day_inp, flow_night_inp, 
                photoperiod_day, photoperiod_night, secchi_depth_day_inp, secchi_depth_night_inp,
                lunar_phase_day, lunar_phase_night, flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
                chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_day_inp, flow_night_inp, secchi_depth_day_inp, secchi_depth_night_inp,
    photoperiod_day, photoperiod_night, lunar_phase_day, 
    lunar_phase_night, flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
    chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2021-2004+1
num_rows = num_years*2
num_covariates = 14
total_covariates = dim(covariates_chinook0)[1]


for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_chinook0[i,] = scale(covariates_chinook0[i,])[,1]
}

#just scale

for(i in (total_covariates/2 + 1):(total_covariates)){
  covariates_chinook0[i,] = scale(covariates_chinook0[i,], center = FALSE, scale= TRUE)[,1]
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
#####

#first using the functions to see which of the correlated covariates best fit
#the data

out.tab.all.years <- NULL
fits.all.years <- NULL
nyears = num_years*2
for(kk in c(2,3,6)){
  c = covariates_chinook0[((1+(kk-1)*num_rows):(kk*num_rows)),]
  name_long = rownames(covariates_chinook0)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
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
#photoperiod difference is best


#now using the functions to see which of the uncorrelated covariates best fit



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
      k = 6
    }
    else if(j==2){
      k = 1
    }
    else if(j==3){
      k = 4
    }
    else if(j==4){
      k = 5
    }
    else if(j==5){
      k = 7
    }
    
    c = rbind(c,covariates_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_chinook0)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-5))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(k==7){
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



fits.hatchery.all.years[[28]] #5252
ci_best_model <- tidy(fits.hatchery.all.years[[28]])


#plot the effects

ggplot(ci_best_model[c(39:43),], 
       aes(x = c("Photoperiod\ndifference", "Flow", "Flow\n difference",
                 "Hatchery,\n day", "Hatchery,\n night"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  theme(axis.text.x=element_text(size=20),axis.title.y=element_text(size=24))

#save this file
ggsave(here("puyallup","output","chinook_all_years_effects_estimate.png"), width = 10, height = 8, units = "in", dpi = 300)

autoplot(fits.hatchery.all.years[[28]], plot.type = "fitted.ytT")
ggsave(here("puyallup","output","chinook_all_years_fitted.png"), width = 10, height = 8, units = "in", dpi = 300)


#checking whether unequal q's are needed

#fitting model with no covariates
#Equal q
fit.model_equal_q = mod_list(nyears,0,hatchery = 0, unequal_q = FALSE)
fit_equal_q <- MARSS(subset_chinook_summer_perhour, model=fit.model_equal_q, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
fit_equal_q$convergence
fit_equal_q$AICc #5584

#unequal q for day and night
fit.model_unequal_q = mod_list(nyears,0,hatchery = 0, unequal_q = TRUE)
fit_unequal_q <- MARSS(subset_chinook_summer_perhour, model=fit.model_unequal_q, 
                       silent = TRUE, method = "BFGS", control=list(maxit=2000))
fit_unequal_q$convergence
fit_unequal_q$AICc #5585

#should have equal q

#now let's test day on night effect

c = NULL
name = NULL
fits.hatchery.all.years2 <- list()
out.tab.hatchery.all.years2 <- NULL
for(k in c(6,1,5,7)){
  c = rbind(c,covariates_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
  name_long = rownames(covariates_chinook0)[1+(k-1)*num_rows]
  name = paste(name, substr(name_long,1,nchar(name_long)-5))
}
c_num = 4
has_hatchery = 1
fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery,day_on_night = TRUE, unequal_q = FALSE))
fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

fit$convergence
fit$AICc #5246
ci <- tidy(fit)

ggplot(ci[c(39:44),], 
       aes(x = c("Photoperiod\ndifference", "Flow", "Flow\n difference",
                 "Hatchery,\n day", "Hatchery,\n day on night", "Hatchery,\n night"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  theme(axis.text.x=element_text(size=20),axis.title.y=element_text(size=24))
ggsave(here("puyallup","output","chinook_all_years_effects_estimate_day_on_night.png"), width = 10, height = 8, units = "in", dpi = 300)


#from data, plot the fraction of wild hatchery during day

data_subset <- data %>% 
  select(chinook0_wild_day_fraction,chinook0_hatchery_day_fraction)

data_subset_long <- data_subset %>% 
  pivot_longer(cols = c(chinook0_wild_day_fraction,chinook0_hatchery_day_fraction),
               names_to = c("species","origin"),
               names_pattern = "(.*)_(.*)_day_fraction",
               values_to = "proportion_day") %>%
  mutate(proportion_night = 1 - proportion_day) %>%
  pivot_longer(cols = c(proportion_day,proportion_night),
               names_to = "daytime_category",
               names_pattern = "proportion_(.*)",
               values_to = "proportion")

data_long_subset_avg <-   data_subset_long %>% 
  filter(!is.na(proportion)) %>%
  group_by(species,daytime_category, origin) %>%
  summarize(average_prop = mean(proportion, na.rm = T), n = n(), se = sd(proportion, na.rm = T)/sqrt(n))

ggplot(data = data_long_subset_avg, aes(y = average_prop, x = daytime_category, fill = origin)) +
  geom_col(position = "dodge")+
  facet_wrap(~species)+
  geom_errorbar(aes(ymin = average_prop-se, ymax = average_prop+se), 
                position = "dodge", alpha = 0.5)+
  theme(axis.text.x=element_text(size=24),axis.title.y=element_text(size=24))+
  labs(x = "",
       y = "Average proportion of \nfish caught per hour", 
       fill = "Origin")+
  scale_fill_manual(values = c("cadetblue","salmon"))

ggsave(here("puyallup","output","chinook_all_years_day_night_proportion.png"), width = 10, height = 8, units = "in", dpi = 300)
