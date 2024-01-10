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
                                            "Hatchery, day", 
                                            "Hatchery, night"),
                                      y = estimate, 
                                      ymin = conf.low, 
                                      ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
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
                                            "Hatchery, day", 
                                            "Hatchery, night"),
                                      y = estimate, 
                                      ymin = conf.low, 
                                      ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
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
                              "Hatchery, day", 
                              "Hatchery, night",
                              "Photoperiod\n difference", 
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
                                 aes(x = c("Photoperiod\n", 
                                           "Residuals",
                                           "Temperature\n difference", 
                                           "Hatchery"),
                                     y = estimate, 
                                     ymin = conf.low, 
                                     ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
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
                              "Photoperiod\n", 
                              "Residuals",
                              "Temperature\n difference"
                              ))

chinook0_skagit_plot
#increase the y axis title font size in chinook0_skagit_plot,
#chinook0_puyallup_plot, and chinook0_dungeness_plot

#change y axis ticks in chinook0_skagit_plot to -1, -0.5, -0.25, 0, 0.25
chinook0_skagit_plot2 <- chinook0_skagit_plot + theme(axis.text.x=element_text(size=24),
                                                      axis.title.y=element_blank(),
                                                      axis.text.y = element_text(size = 18))+
  scale_y_continuous(breaks = c(-0.1,0,0.1,0.2), limits = c(-0.15,0.15))+
  scale_x_discrete(expand = expansion(add = 1.15), guide = guide_axis(n.dodge=3),
                   limits = c("Hatchery",
                              "Photoperiod\n", 
                              "Residuals",
                              "Temperature\n difference"
                              ))

chinook0_skagit_plot2
#remove y axis title in chinook0_puyallup_plot and chinook0_dungeness_plot
chinook0_puyallup_plot2 <- chinook0_puyallup_plot + theme(axis.title.y=element_blank(),
                                                          axis.text.y = element_text(size = 18),
                                                          axis.text.x=element_text(size=24)) +
  scale_x_discrete(guide = guide_axis(n.dodge=3),limits = c("Flow\n difference","Photoperiod\n difference","Hatchery, day", 
                                                            "Hatchery, night", "Flow")) + 
  scale_y_continuous(breaks = c(-0.25,0,0.25), limits = c(-0.32,0.32))

chinook0_puyallup_plot2

#reduce space between x axis ticks in chinook0_dungeness_plot
chinook0_dungeness_plot2 <- chinook0_dungeness_plot + theme(axis.title.y=element_text(size=28),
                                                            axis.title.x=element_text(size=28),
                                                            axis.text.y = element_text(size = 18),
                                                            axis.text.x=element_text(size=24)) +
  scale_x_discrete(expand = expansion(add = 1.1), #guide = guide_axis(n.dodge=3),
                   limits = c("Flow\n difference",
                              "Photoperiod\n difference", 
                              "Hatchery, day", 
                              "Hatchery, night",
                              "Temperature\n difference")
                   )+
  scale_y_continuous(breaks = c(-0.15,0,0.15), limits = c(-0.2,0.2))+
  coord_flip()
  

chinook0_dungeness_plot2


#put all the figures together and label the panels a,b,c

ggpubr::ggarrange(chinook0_dungeness_plot,
                  chinook0_puyallup_plot, 
                  chinook0_skagit_plot,
                  labels = c("a", "b", "c"), ncol = 3, nrow = 1, 
                  common.legend = TRUE, legend = "right")

ggsave(here("output","chinook0_covariates_estimates_2.jpeg"), width = 18.5, height = 8)


ggpubr::ggarrange(chinook0_dungeness_plot2,
                  chinook0_puyallup_plot2, 
                  chinook0_skagit_plot2,
                  labels = c("a", "b", "c"), ncol = 3, nrow = 1, font.label = list(size = 28),
                  common.legend = TRUE, legend = "right",
                  widths = c(1.2, 1,1.1))

#save the figure with width and height in inches

ggsave(here("output","chinook0_covariates_estimates_3.jpeg"), width = 18.5, height = 8)


ggpubr::ggarrange(chinook0_dungeness_plot,
                  chinook0_puyallup_plot, 
                  chinook0_skagit_plot,
                  labels = c("a", "b", "c"), ncol = 1, nrow = 3, 
                  common.legend = TRUE, legend = "right")


ggsave(here("output","chinook0_covariates_estimates_long_@.jpeg"), width = 6, height = 8)


#flip the order of the scale_x_discrete limits
chinook0_skagit_plot3 <- chinook0_skagit_plot + theme(axis.text.x=element_text(size=24, family = "Sans"),
                                                      axis.title.y=element_text(size=24, family = "Sans", 
                                                                                margin = margin(t = 10, r = 0, b = 0, l = 10)),
                                                      axis.title.x=element_text(size=24, family = "Sans", 
                                                                                margin = margin(t = 15, r = 20, b = 0, l = 0)),
                                                      axis.text.y = element_text(size = 24, family = "Sans"))+
  scale_y_continuous(breaks = c(-0.1,0,0.1,0.2), limits = c(-0.15,0.15))+
  scale_x_discrete(#expand = expansion(add = 1.15), #guide = guide_axis(n.dodge=3),
                   limits = c("Temperature\n difference",
                              "Residuals",
                              "Photoperiod\n", 
                              "Hatchery"
                   )) + 
  coord_flip()

chinook0_skagit_plot3
#remove y axis title in chinook0_puyallup_plot and chinook0_dungeness_plot

#flip the order of the scale_x_discrete limits

chinook0_puyallup_plot3 <- chinook0_puyallup_plot + theme(axis.title.y=element_text(size=24, family = "Sans", 
                                                                                    margin = margin(t = 10, r = 0, b = 0, l = 10)),
                                                          axis.title.x=element_text(size=24, family = "Sans", 
                                                                                    margin = margin(t = 15, r = 0, b = 0, l = 10)),
                                                          axis.text.y = element_text(size = 24, family = "Sans"),
                                                          axis.text.x=element_text(size=24, family = "Sans")) +
  scale_x_discrete(#guide = guide_axis(n.dodge=3),
                   limits = c("Flow","Hatchery, day", "Hatchery, night", "Photoperiod\n difference",
                              "Flow\n difference"
                              
                              )) + 
  scale_y_continuous(breaks = c(0,0.1, 0.2), limits = c(-0.05,0.32))+
  coord_flip()

chinook0_puyallup_plot3

#reduce space between x axis ticks in chinook0_dungeness_plot
chinook0_dungeness_plot3 <- chinook0_dungeness_plot + theme(axis.title.y=element_text(size=24, family = "Sans"),
                                                            axis.title.x=element_text(size=24, family = "Sans", 
                                                                                      margin = margin(t = 15, r = 0, b = 0, l = 10)),
                                                            axis.text.y = element_text(size = 24, family = "Sans"),
                                                            axis.text.x=element_text(size=24, family = "Sans")) +
  scale_x_discrete(#expand = expansion(add = 1.1), #guide = guide_axis(n.dodge=3),
                   limits = c("Temperature\n difference",
                              "Hatchery, day", 
                              "Hatchery, night",
                              "Photoperiod\n difference",
                              "Flow\n difference"
                              )
  )+
  scale_y_continuous(breaks = c(-0.15,0,0.15), limits = c(-0.2,0.2))+
  coord_flip()


chinook0_dungeness_plot3


ggpubr::ggarrange(chinook0_dungeness_plot3,
                  chinook0_puyallup_plot3, 
                  chinook0_skagit_plot3,
                  labels = c("a", "b", "c"), ncol = 3, nrow = 1, font.label = list(size = 28, family = "Sans"),
                  common.legend = TRUE, legend = "right",
                  widths = c(1.15, 1.25, 1), heights = c(0.5,0.5,0.5))

#save the figure with width and height in inches

ggsave(here("output","chinook0_covariates_estimates_4.jpeg"), width = 18.5, height = 8)


####### coho

# Dungeness


#check system
if(.Platform$OS.type == "unix") {
  data_string = here("..","..","data","pied_piper","dungeness")
} else {
  data_string = here("..","..","..","OneDrive","Documents","data","pied_piper","dungeness")
}


#load data

data_day <- read.csv(here("data","dungeness_unaggregated_day_w_covariates.csv"))
data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"

#interpolate the hatchery covariate values
data_day$coho1_hatchery_perhour_interpolate <- na.approx(data_day$coho1_hatchery_perhour, na.rm = FALSE)
data_night$coho1_hatchery_perhour_interpolate <- na.approx(data_night$coho1_hatchery_perhour, na.rm = FALSE)


data_day_night <- rbind(data_day,data_night)

covariates_coho1 <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy >120 & doy <= 160) %>%
  dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_april, 
                lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
                coho1_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = c(
    temp, flow, photoperiod, atu_april, lunar_phase, 
    resid, temp_diff, flow_diff, photo_diff, coho1_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

num_years = 2020-2005
num_rows = num_years*2


dim(covariates_coho1)

#scale the covariates

for(i in 1:(dim(covariates_coho1)[1]/2)){ # everything except resid, diffs and hatchery
  covariates_coho1[i,] = scale(covariates_coho1[i,])[,1]
}

#only scale the rest, do not center
for(i in (dim(covariates_coho1)[1]/2 +1):(dim(covariates_coho1)[1])){
  covariates_coho1[i,] = scale(covariates_coho1[i,], center = FALSE, scale= TRUE)[,1]
}


#subset response variable
subset_coho_summer_perhour <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy > 120 & doy <= 160) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour)[1]){
  subset_coho_summer_perhour[i,] = scale(subset_coho_summer_perhour[i,])[,1]
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
#function to get covariate combinations
##### 
get_covariate_combinations <- function(covariates) {
  n <- length(covariates)
  combinations <- lapply(1:n, function(x) combn(covariates, x, simplify = FALSE))
  unlist(combinations, recursive = FALSE)
}
#####


c <- NULL
name_individual = NULL
nyears = num_years*2
for(kk in c(2,3,7,8,10)){
  c = rbind(c,covariates_coho1[((1+(kk-1)*nyears):(kk*nyears)),])
  name_long = rownames(covariates_coho1)[1+(kk-1)*nyears]
  name_individual = paste(name_individual,substr(name_long,1,nchar(name_long)-9))
  print(name_individual)
  

}
fit.model = c(list(c= c), mod_list(nyears,5,1,FALSE,FALSE))
fit_coho_dungeness <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

ci_coho1_dungeness <- tidy(fit_coho_dungeness)
ci_coho1_dungeness

coho1_dungeness_plot <- ggplot(ci_coho1_dungeness[c(33:38),], 
       aes(x = c("Flow", "Photoperiod", "Temperature\n difference",  "Flow\n difference", "Hatchery,\nday", "Hatchery,\nnight"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")  +
  labs(x = "", y = "Estimate of effect",
  )+
  theme_classic()+
  theme(axis.title.y=element_text(size=28),
        axis.text.y = element_text(size = 18),
        axis.text.x=element_text(size=24)) +
  scale_y_continuous(breaks = seq(-0.1,0,0.1), limits = c(-0.15,0.15))+
  scale_x_discrete(limits = c("Flow\n difference", "Hatchery,\nday", "Hatchery,\nnight",
                              "Flow", "Photoperiod", "Temperature\n difference"),
                   guide = guide_axis(n.dodge = 2),
                   expand = c(0.15,0.15)
                   )

coho1_dungeness_plot
#skagit



data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))

data_day_night$coho1_wild_perhour <- data_day_night$coho1_wild_num/data_day_night$In
data_day_night$coho1_hatchery_perhour <- data_day_night$coho1_hatchery_num/data_day_night$In


data_day_night$coho1_hatchery_perhour_inp <- na.approx(data_day_night$coho1_hatchery_perhour, 
                                                       na.rm = FALSE)

lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day_night)

data_day_night <- data_day_night %>% add_residuals(lm_temp_day)


covariates_coho1_skagit_night <- arrange(data_day_night,doy) %>%
  filter(doy >100 & doy < 150 & daytime_category == 'night') %>%
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
total_covariates = dim(covariates_coho1_skagit_night)[1]


for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_coho1_skagit_night[i,] = scale(covariates_coho1_skagit_night[i,])[,1]
}

#just scale


for(i in (total_covariates/2 + 1):(total_covariates)){
  if(sum(covariates_coho1_skagit_night[i,]) != 0){
    covariates_coho1_skagit_night[i,] = scale(covariates_coho1_skagit_night[i,], 
                                              center = FALSE, scale= TRUE)[,1]
  }
}

#subset response variable
subset_coho_summer_perhour_night <- arrange(data_day_night,doy) %>%
  filter(doy > 100 & doy < 150 & daytime_category == 'night') %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category,trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour_night)[1]){
  subset_coho_summer_perhour_night[i,] = scale(subset_coho_summer_perhour_night[i,])[,1]
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
name_individual = NULL
nyears = num_years*2
for(kk in c(2,4,5,10)){
  c = rbind(c,covariates_coho1_skagit_night[((1+(kk-1)*nyears):(kk*nyears)),])
  name_long = rownames(covariates_coho1_skagit_night)[1+(kk-1)*nyears]
  name_individual = paste(name_individual, substr(name_long,1,nchar(name_long)-15))
  print(name_individual)
  
  
}

fit.model = c(list(c= c), mod_list(nyears,4,1))
fit_coho1_skagit <- MARSS(subset_coho_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))


ci_coho1_skagit <- tidy(fit_coho1_skagit)
ci_coho1_skagit

coho1_skagit_plot <- ggplot(ci_coho1_skagit[c(27:30),], 
                               aes(x = c("Flow", "ATU", "Lunar phase\n","Hatchery\n"),
                                   y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")  +
  labs(x = "", y = "",
  )+
  theme_classic()+
  theme(axis.title.y=element_text(size=28),
        axis.text.y = element_text(size = 18),
        axis.text.x=element_text(size=24)) +
  scale_y_continuous(breaks = seq(-0.1,0,0.1), limits = c(-0.1,0.1))+
  scale_x_discrete(limits = c("Hatchery\n","Lunar phase\n",
                              "Flow", "ATU"),
                   guide = guide_axis(n.dodge = 2))

coho1_skagit_plot


ggpubr::ggarrange(coho1_dungeness_plot,
                  coho1_skagit_plot,
                  labels = c("a", "b"), ncol = 2, nrow = 1, font.label = list(size = 28),
                  common.legend = TRUE, legend = "right",
                  widths = c(1.2, 1))

#save the figure with width and height in inches

ggsave(here("output","coho1_covariates_estimates.jpeg"), width = 18.5, height = 8)


#make coho1_dungness_plot2 in the same way as chinook0_dungeness_plot3


coho1_dungeness_plot2 <- ggplot(ci_coho1_dungeness[c(33:38),], 
                               aes(x = c("Flow", "Photoperiod", "Temperature\n difference", 
                                         "Flow\n difference", "Hatchery, day", "Hatchery, night"),
                                   y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")  +
  labs(x = "", y = "Estimate of effect",
  )+
  theme_classic()+
  theme(axis.title.y=element_text(size=24),
        axis.title.x=element_text(size=24),
        axis.text.y = element_text(size = 24),
        axis.text.x=element_text(size=24)) +
  scale_y_continuous(breaks = seq(-0.1,0.1,0.1), limits = c(-0.15,0.15))+
  scale_x_discrete(limits = c("Temperature\n difference", "Photoperiod", "Flow",
                              "Hatchery, night","Hatchery, day", 
                              "Flow\n difference"
                               ),
                   #guide = guide_axis(n.dodge = 2),
                   #expand = c(0.15,0.15)
  )+
  coord_flip()

coho1_dungeness_plot2




coho1_skagit_plot2 <- ggplot(ci_coho1_skagit[c(27:30),], 
                            aes(x = c("Flow", "ATU", "Lunar phase\n","Hatchery\n"),
                                y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")  +
  labs(x = "", y = "Estimate of effect",
  )+
  theme_classic()+
  theme(axis.title.y=element_text(size=24),
        axis.title.x=element_text(size=24),
        axis.text.y = element_text(size = 24),
        axis.text.x=element_text(size=24)) +
  scale_y_continuous(breaks = seq(-0.1,0.1,0.1), limits = c(-0.12,0.12))+
  scale_x_discrete(limits = c("ATU", "Flow","Lunar phase\n",
                              "Hatchery\n")) + 
  coord_flip()
                   #guide = guide_axis(n.dodge = 2))

coho1_skagit_plot2


ggpubr::ggarrange(coho1_dungeness_plot2,
                  coho1_skagit_plot2,
                  labels = c("a", "b"), ncol = 2, nrow = 1, font.label = list(size = 28),
                  common.legend = TRUE, legend = "right",
                  widths = c(1, 1))


ggsave(here("output","coho1_covariates_estimates_2.jpeg"), width = 18.5, height = 8)
