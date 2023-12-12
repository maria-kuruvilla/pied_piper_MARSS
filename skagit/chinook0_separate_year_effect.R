#Goals - 
#try MARSS model with just hatchery
#try MARSS model with just hatchery but with separate year effects
#try MARSS model selection with separate year effects for hatchery


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

data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))

data_day_night$chinook0_wild_perhour <- data_day_night$chinook0_wild_num/data_day_night$In
data_day_night$chinook0_hatchery_perhour <- data_day_night$chinook0_hatchery_num/data_day_night$In
data_day_night$coho1_wild_perhour <- data_day_night$coho1_wild_num/data_day_night$In
data_day_night$coho1_hatchery_perhour <- data_day_night$coho1_hatchery_num/data_day_night$In
data_day_night$chinook1_wild_perhour <- data_day_night$chinook1_wild_num/data_day_night$In
data_day_night$chinook1_hatchery_perhour <- data_day_night$chinook1_hatchery_num/data_day_night$In




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



Cmat_skagit <- function(nyears,ncov,hatchery=0, year_effect = TRUE){
  vars = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")
  #make a list of years from 2010 to 2022 but each year has to be repeated 4 times
  years = c(rep("2010",4),rep("2012",4),rep("2013",4),rep("2014",4),rep("2015",4),
            rep("2016",4),rep("2017",4),rep("2018",4),rep("2019",4),rep("2020",4),rep("2021",4),rep("2022",4
            ))  
  # years = c("2010","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022")
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
                if(year_effect == TRUE){
                  C[i,j] <- paste("day",years[i], sep = "_")
                }
                else{
                  C[i,j] <- "day"
                }
                
                
              }
              else{
                
                if(year_effect == TRUE){
                  C[i,j] <- paste("night",years[i], sep = "_")
                }
                else{
                  C[i,j] <- "night"
                }
                
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

mod_list <- function(nyears,ncov,hatchery=0, day_on_night = FALSE, unequal_q = FALSE, year_effect = TRUE){
  
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
      
      C = Cmat_skagit(nyears,ncov,hatchery, year_effect)
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

#fit MARSS model with just hatchery and no year effect
kk=10
nyears = (2022-2010)*2*2
c = covariates_chinook0_skagit[((1+(kk-1)*nyears):(kk*nyears)),]
name_long = rownames(covariates_chinook0_skagit)[1+(kk-1)*nyears]
name_individual = substr(name_long,1,nchar(name_long)-15)
print(name_individual)
fit.model = c(list(c= c), mod_list(nyears,1,1,FALSE,FALSE,FALSE))
fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
ci_no_year_effect <- tidy(fit)
fit$AICc #2544.859

#fit MARSS model with just hatchery and year effect
kk=10
nyears = (2022-2010)*2*2
c = covariates_chinook0_skagit[((1+(kk-1)*nyears):(kk*nyears)),]
name_long = rownames(covariates_chinook0_skagit)[1+(kk-1)*nyears]
name_individual = substr(name_long,1,nchar(name_long)-15)
print(name_individual)
fit.model = c(list(c= c), mod_list(nyears,1,1,FALSE,FALSE,TRUE))
fit_year_effect <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
ci_year_effect <- tidy(fit_year_effect)

fit_year_effect$AICc #2550

#plot the covariates_chinook0_skagit[433:480,] by year
#make facet labels on one line
#plot only daytime_category = "night"

ggplot(data_day_night %>% 
         filter(daytime_category == "night")) +
  geom_line(aes(x = doy, y= chinook0_hatchery_perhour_inp), color = "cadetblue",
            linewidth = 1, alpha = 0.5) +
  geom_line(aes(x = doy, y = chinook0_wild_perhour), color = "salmon",
            linewidth = 1, alpha = 0.5) +
  facet_wrap(~year+trap, scales = "free", 
             labeller = label_wrap_gen(multi_line=FALSE)) +
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "day", y = "day")+
  xlim(100,200)


#in data_day_night, group by doy and year 