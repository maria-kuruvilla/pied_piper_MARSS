#Goal - to remake the figures with the estimates of the effects
#of the covariates from the MARSS analysis

# one figure for Chinook with 3 panels for the 3 rivers

# one figure for Coho with 2 panels for the 2 rivers

#load libraries

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
library(ggpubr)

#load data

#dungeness



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


covariates_chinook0_dungeness <- arrange(data_day_night,doy) %>%
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


dim(covariates_chinook0_dungeness)

#scale the covariates

for(i in 1:(dim(covariates_chinook0_dungeness)[1]/2)){ # everything except resid, diffs and hatchery
  covariates_chinook0_dungeness[i,] = scale(covariates_chinook0_dungeness[i,])[,1]
}

#skip 2005 day
for(i in (dim(covariates_chinook0_dungeness)[1]/2 +1):(dim(covariates_chinook0_dungeness)[1] - num_rows)){
  covariates_chinook0_dungeness[i,] = scale(covariates_chinook0_dungeness[i,], center = FALSE, scale= TRUE)[,1]
}

#skip 2005 night
for(i in (dim(covariates_chinook0_dungeness)[1] - num_rows + 2):(dim(covariates_chinook0_dungeness)[1] - num_rows/2)){
  covariates_chinook0_dungeness[i,] = scale(covariates_chinook0_dungeness[i,], center = FALSE, scale= TRUE)[,1]
}


for(i in (dim(covariates_chinook0_dungeness)[1] - num_rows/2 + 2):(dim(covariates_chinook0_dungeness)[1])){
  covariates_chinook0_dungeness[i,] = scale(covariates_chinook0_dungeness[i,], center = FALSE, scale= TRUE)[,1]
}




#subset response variable
chinook0_dungeness <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy > 130 & doy <= 200) %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(chinook0_dungeness)[1]){
  chinook0_dungeness[i,] = scale(chinook0_dungeness[i,])[,1]
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


#######

nyears = num_years*2
c <- NULL
name_individual <- NULL
for(kk in c(7,8,9,10)){
  c = rbind(c,covariates_chinook0_dungeness[((1+(kk-1)*nyears):(kk*nyears)),])
  name_long = rownames(covariates_chinook0_dungeness)[1+(kk-1)*nyears]
  name_individual = paste(name_individual,substr(name_long,1,nchar(name_long)-9))
  print(name_individual)
  
}

fit.model = c(list(c= c), mod_list(nyears,4,1,FALSE,TRUE))
fit_chinook0_dungeness <- MARSS(chinook0_dungeness, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
ci_chinook0_dungeness <- tidy(fit_chinook0_dungeness)


#plotting
#remove grid lines
#dodge the axis labels on the x axis
chinook0_dungeness_plot <- ggplot(ci_chinook0_dungeness[c(34:38),], 
                                  aes(x = c("Temperature\n difference", 
                                            "Flow\n difference",
                                            "Photoperiod\n difference", 
                                            "Hatchery,\n day", 
                                            "Hatchery,\n night"),
                                      y = estimate, 
                                      ymin = conf.low, 
                                      ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect"
       # title = "Dungeness River, Chinook sub yearlings"
  )+
  theme_classic() +
  theme(axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y = element_text(size = 14)) + 
    scale_x_discrete(guide = guide_axis(n.dodge=2))



#puyallup 
data <- read.csv(here("puyallup", "data","puyallup_2004-2021_all_days_w_covariates_edited.csv"),
                 na.strings = c("NA",""))

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
covariates_chinook0_puyallup <- arrange(data,doy) %>%
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
total_covariates = dim(covariates_chinook0_puyallup)[1]


for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_chinook0_puyallup[i,] = scale(covariates_chinook0_puyallup[i,])[,1]
}

#just scale

for(i in (total_covariates/2 + 1):(total_covariates)){
  covariates_chinook0_puyallup[i,] = scale(covariates_chinook0_puyallup[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
chinook0_puyallup <- arrange(data,doy) %>%
  filter(doy > 130 & doy <= 218) %>%
  mutate(log.value_day = log(chinook0_wild_perhour_day + 1), 
         log.value_night = log(chinook0_wild_perhour_night + 1)) %>%
  dplyr::select(log.value_day, log.value_night ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_day, log.value_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(chinook0_puyallup)[1]){
  chinook0_puyallup[i,] = scale(chinook0_puyallup[i,])[,1]
}


num_years = 2021-2004+1
num_rows = num_years*2
c <- NULL
name_individual = NULL
for(kk in c(1,5,6,7)){
  c = rbind(c,covariates_chinook0_puyallup[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_chinook0_puyallup)[1+(kk-1)*num_rows]
  name_individual = paste(name_individual,substr(name_long,1,nchar(name_long)-9))
  print(name_individual)
  
}

fit.model = c(list(c= c), mod_list(num_rows,4,1, FALSE, FALSE))
fit_chinook0_puyallup <- MARSS(chinook0_puyallup, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

ci_chinook0_puyallup <- tidy(fit_chinook0_puyallup)


#plotting
#remove grid lines
#dodge the axis labels on the x axis
chinook0_puyallup_plot <- ggplot(ci_chinook0_puyallup[c(39:43),], 
                                  aes(x = c("Flow", 
                                            "Flow\n difference",
                                            "Photoperiod\n difference", 
                                            "Hatchery,\n day", 
                                            "Hatchery,\n night"),
                                      y = estimate, 
                                      ymin = conf.low, 
                                      ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect"
       # title = "Dungeness River, Chinook sub yearlings"
  )+
  theme_classic() +
  theme(axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y = element_text(size = 14)) + 
  scale_x_discrete(guide = guide_axis(n.dodge=2),
                   limits = c("Flow\n difference",
                              "Photoperiod\n difference", 
                              "Hatchery,\n day", 
                              "Hatchery,\n night",
                              "Flow"))



#skagit

#load data
data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))

data_day_night$chinook0_wild_perhour <- data_day_night$chinook0_wild_num/data_day_night$In
data_day_night$chinook0_hatchery_perhour <- data_day_night$chinook0_hatchery_num/data_day_night$In


data_day_night$chinook0_hatchery_perhour_inp <- na.approx(data_day_night$chinook0_hatchery_perhour, 
                                                          na.rm = FALSE)


lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day_night)

data_day_night <- data_day_night %>% add_residuals(lm_temp_day)



covariates_chinook0_skagit_night <- arrange(data_day_night,doy) %>%
  filter(doy >150 & doy < 189 & daytime_category == 'night') %>%
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
total_covariates = dim(covariates_chinook0_skagit_night)[1]


for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_chinook0_skagit_night[i,] = scale(covariates_chinook0_skagit_night[i,])[,1]
}

#just scale


for(i in (total_covariates/2 + 1):(total_covariates)){
  if(sum(covariates_chinook0_skagit_night[i,]) != 0){
    covariates_chinook0_skagit_night[i,] = scale(covariates_chinook0_skagit_night[i,], 
                                                 center = FALSE, scale= TRUE)[,1]
  }
}

#subset response variable
chinook0_skagit_night <- arrange(data_day_night,doy) %>%
  filter(doy > 150 & doy < 189 & daytime_category == 'night') %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category,trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(chinook0_skagit_night)[1]){
  chinook0_skagit_night[i,] = scale(chinook0_skagit_night[i,])[,1]
}


#function for C matrix
#######


Cmat_skagit_night <- function(nyears,ncov,hatchery=0){
  vars = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")
  
  
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
  
  return(C)
}


#######





#function for mod list

######

mod_list <- function(nyears,ncov,hatchery=0){
  
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
      C = Cmat_skagit_night(nyears,ncov,hatchery)
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


c <- NULL
name_individual <- NULL
nyears = num_years*2
for(kk in c(3,6,7,10)){
  c = rbind(c,covariates_chinook0_skagit_night[((1+(kk-1)*nyears):(kk*nyears)),])
  name_long = rownames(covariates_chinook0_skagit_night)[1+(kk-1)*nyears]
  name_individual = paste(name_individual, substr(name_long,1,nchar(name_long)-15))
  print(name_individual)
  
}

fit.model = c(list(c= c), mod_list(nyears,4,1))
fit_chinook0_skagit <- MARSS(chinook0_skagit_night, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

ci_chinook0_skagit <- tidy(fit_chinook0_skagit)


chinook0_skagit_plot <- ggplot(ci_chinook0_skagit[c(27:30),], 
                                 aes(x = c("Photoperiod", 
                                           "Residuals",
                                           "Temperature\n difference", 
                                           "Hatchery"),
                                     y = estimate, 
                                     ymin = conf.low, 
                                     ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect"
       # title = "Dungeness River, Chinook sub yearlings"
  )+
  theme_classic() +
  theme(axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y = element_text(size = 14)) + 
  scale_x_discrete(guide = guide_axis(n.dodge=2),
                   limits = c("Hatchery",
                              "Photoperiod", 
                              "Residuals",
                              "Temperature\n difference"
                              ))

#put all the figures together and label the panels a,b,c

ggpubr::ggarrange(chinook0_dungeness_plot,
                  chinook0_puyallup_plot, 
                  chinook0_skagit_plot,
                  labels = c("a", "b", "c"), ncol = 3, nrow = 1, 
                  common.legend = TRUE, legend = "right")

#save the figure with width and height in inches

ggsave(here("output","chinook0_covariates_estimates_2.jpeg"), width = 18.5, height = 8)


ggpubr::ggarrange(chinook0_dungeness_plot,
                  chinook0_puyallup_plot, 
                  chinook0_skagit_plot,
                  labels = c("a", "b", "c"), ncol = 1, nrow = 3, 
                  common.legend = TRUE, legend = "right")


ggsave(here("output","chinook0_covariates_estimates_long.jpeg"), width = 8, height = 12)

