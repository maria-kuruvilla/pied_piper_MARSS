# Goal 

# try MARSS analysis with edited Puyallup data for chinook0
#all years from 2004-2021
#with and without separate error for day and night
#without day effect on night

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

data <- read.csv(here("puyallup", "data","puyallup_2004-2021_all_days_w_covariates_edited.csv"),
                 na.strings = c("NA",""))

#plotting chinook0_wild_perhour with doy to see timing of migration

ggplot(data, aes(x = doy, y = chinook0_wild_perhour_day)) +
  geom_point(alpha = 0.5) +
  labs(x = "Day of Year", y = "Chinook0 Wild per Hour") +
  theme_bw()+
  facet_wrap(~year, ncol = 4, scales = "free_y") +
  xlim(50, 250)

ggplot(data, aes(x = doy, y = chinook0_wild_perhour_night)) +
  geom_point(alpha = 0.5) +
  labs(x = "Day of Year", y = "Chinook0 Wild per Hour") +
  theme_bw()+
  facet_wrap(~year, ncol = 4, scales = "free_y") +
  xlim(50, 250)

#do same for hatchery

ggplot(data, aes(x = doy, y = chinook0_hatchery_perhour_day)) +
  geom_point(alpha = 0.5) +
  labs(x = "Day of Year", y = "Chinook0 Hatchery per Hour") +
  theme_bw()+
  facet_wrap(~year, ncol = 4, scales = "free_y") +
  xlim(50, 250)

ggplot(data, aes(x = doy, y = chinook0_hatchery_perhour_night)) +
  geom_point(alpha = 0.5) +
  labs(x = "Day of Year", y = "Chinook0 Hatchery per Hour") +
  theme_bw()+
  facet_wrap(~year, ncol = 4, scales = "free_y") +
  xlim(50, 250)



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

#first run marss analysis without any covaraites to see if equal errors or unequal errors 
#fit the data better

fit.model_equal_q = mod_list(nyears,0,0, FALSE, FALSE)
fit_equal_q <- MARSS(subset_chinook_summer_perhour, model=fit.model_equal_q, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
fit_equal_q$AICc #5584.9
fit.model_unequal_q = mod_list(nyears,0,0, FALSE, TRUE)
fit_unequal_q <- MARSS(subset_chinook_summer_perhour, model=fit.model_unequal_q, silent = TRUE, method = "BFGS",
                       control=list(maxit=2000))
fit_unequal_q$AICc #5585.8

#equal errors are better

data %>%
  dplyr::select(flow_day,photoperiod_day, lunar_phase_day, secchi_depth_day_inp, 
                photo_diff_day,flow_diff_day) %>%
  GGally::ggpairs(upper=list(continuous='points'),
                  lower=list(continuous='cor'))
#actually only photoperiod difference and secchi depth are correlated

#save this figure
ggsave(here("puyallup","output","covariate_correlations.png"), width = 10, height = 10, units = "in")

#first using the functions to see which of the correlated covariates best fit
#the data


out.tab.all.years <- NULL
fits.all.years <- NULL
nyears = num_years*2
for(kk in c(2,6)){
  c = covariates_chinook0[((1+(kk-1)*num_rows):(kk*num_rows)),]
  name_long = rownames(covariates_chinook0)[1+(kk-1)*num_rows]
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

#photoperiod difference is best

photoperiod_difference = 0
photoperiod = 0
flow = 0
lunar_phase = 0
hatchery = 0
flow_difference = 0

list_combinations <- get_covariate_combinations(1:6)
out.tab.hatchery.all.years <- NULL
fits.hatchery.all.years <- list()
for(i in 1:length(list_combinations)){
  photoperiod_difference = 0
  photoperiod = 0
  flow = 0
  lunar_phase = 0
  hatchery = 0
  flow_difference = 0
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  for(j in covariates){
    if(j == 1){
      k = 6
      photoperiod_difference = 1
    }
    else if(j==2){
      k = 1
      flow = 1
    }
    else if(j==3){
      k = 3
      photoperiod = 1
    }
    else if(j==4){
      k = 4
      lunar_phase = 1
    }
    else if(j==5){
      k = 5
      flow_difference = 1
    }
    else if(j==6){
      k = 7
      hatchery = 1
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
  
  
  out=data.frame(c=name, photoperiod_difference = photoperiod_difference,
                 photoperiod = photoperiod,
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

#calculate relative variable importance


fit <- MARSS(subset_chinook_summer_perhour, model=mod_list(nyears,0,0,FALSE,FALSE), 
             silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
photoperiod_difference = 0
photoperiod = 0
flow = 0
lunar_phase = 0
hatchery = 0
flow_difference = 0
out=data.frame(c=name, photoperiod_difference = photoperiod_difference, flow = flow,
               photoperiod = photoperiod,
               flow_difference = flow_difference,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)

#first drop the delta Aicc column from out.tab.hatchery.all.years
out.tab.hatchery.all.years$deltaAICc <- NULL

out.tab.hatchery.all.years=rbind(out.tab.hatchery.all.years,out)
fits.hatchery.all.years=c(fits.hatchery.all.years,list(fit))

# out.tab.hatchery.all.years$deltaAICc <- NULL
out.tab.hatchery.all.years$rel.LL <- NULL
out.tab.hatchery.all.years$weights <- NULL


weights <- akaike.weights(out.tab.hatchery.all.years$AICc)

out.tab.hatchery.all.years$deltaAICc <- weights$deltaAIC
out.tab.hatchery.all.years$rel.LL <- weights$rel.LL
out.tab.hatchery.all.years$weights <- weights$weights


min.AICc <- order(out.tab.hatchery.all.years$AICc)
out.tab.hatchery.all.years.ordered <- out.tab.hatchery.all.years[min.AICc, ]
out.tab.hatchery.all.years.ordered

out.tab.hatchery.all.years.ordered$cumulative_weights <- cumsum(out.tab.hatchery.all.years.ordered$weights)

relative_importance_photo_diff <- sum(out.tab.hatchery.all.years$weights[out.tab.hatchery.all.years$photoperiod_difference==1])
relative_importance_flow <- sum(out.tab.hatchery.all.years$weights[out.tab.hatchery.all.years$flow==1])
relative_importance_lunar_phase <- sum(out.tab.hatchery.all.years$weights[out.tab.hatchery.all.years$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab.hatchery.all.years$weights[out.tab.hatchery.all.years$hatchery==1])
relative_importance_flow_diff <- sum(out.tab.hatchery.all.years$weights[out.tab.hatchery.all.years$flow_difference==1])
relative_importance_photoperiod <- sum(out.tab.hatchery.all.years$weights[out.tab.hatchery.all.years$photoperiod==1])

riv_photo_diff_puyallup <- data.frame(variable = c("photoperiod difference",
                                          "photoperiod",
                                          "flow",
                                          "lunar phase",
                                          "flow difference",
                                          "hatchery"),
                             relative_importance = c(relative_importance_photo_diff,
                                                     relative_importance_photoperiod,
                                                     relative_importance_flow,
                                                     relative_importance_lunar_phase,
                                                     relative_importance_flow_diff,
                                                     relative_importance_hatchery))


#sort according to relative importance
riv_photo_diff_puyallup[order(riv_photo_diff_puyallup$relative_importance, decreasing = TRUE),]

riv_photo_diff_puyallup <- riv_photo_diff_puyallup[order(riv_photo_diff_puyallup$relative_importance, decreasing = TRUE),]
riv_photo_diff_puyallup$relative_importance <- round(riv_photo_diff_puyallup$relative_importance,1)

riv_photo_diff_puyallup

#save the dataframe riv_photo_diff_puyallup

write.csv(riv_photo_diff_puyallup, file = here("puyallup",
                                               "output","riv_photo_diff_puyallup.csv"))

out.tab.hatchery.all.years.ordered

#round the weights, deltaAICc, AICc, rel.LL, and logLik columns to 2 decimal places
out.tab.hatchery.all.years.ordered$weights <- round(out.tab.hatchery.all.years.ordered$weights,2)
out.tab.hatchery.all.years.ordered$deltaAICc <- round(out.tab.hatchery.all.years.ordered$deltaAICc,2)
out.tab.hatchery.all.years.ordered$AICc <- round(out.tab.hatchery.all.years.ordered$AICc,2)
out.tab.hatchery.all.years.ordered$rel.LL <- round(out.tab.hatchery.all.years.ordered$rel.LL,2)
out.tab.hatchery.all.years.ordered$logLik <- round(out.tab.hatchery.all.years.ordered$logLik,2)

out.tab.hatchery.all.years.ordered

#save the dataframe out.tab.hatchery.all.years.ordered

write.csv(out.tab.hatchery.all.years.ordered, file = here("puyallup",
                                                          "output",
                                                          "model_selection.csv"))

#plot the results from the best model
autoplot(fits.hatchery.all.years[[28]])

#make my make own plot by using predict() function in MARSS

#predict the model
pred <- predict(fits.hatchery.all.years[[28]], type = "ytT", interval = "confidence")
head(pred)

#make labels for all unique pred$pred$.rownames

#make a dataframe of the unique pred$pred$.rownames

pred.rownames <- substring(unique(pred$pred$.rownames), 11, 
                           length(unique(pred$pred$.rownames)))
names(pred.rownames) <- unique(pred$pred$.rownames)

#add a column in pred$pred for label and use pred.rownames to label the rows

pred$pred$label <- pred.rownames[pred$pred$.rownames]

#make a new column for daynight_category and year by parsing the label column

pred$pred$daynight_category <- ifelse(grepl("day", pred$pred$label), "day", "night")
#use last 4 charachters of label column to get year
pred$pred$year <- as.numeric(substr(pred$pred$label, nchar(pred$pred$label)-3, nchar(pred$pred$label)))
#make new columns for if year is less than 2012 and call it 1 else 0
pred$pred$year.2012 <- ifelse(pred$pred$year < 2012, 1, 0)

#reorder pred$pred according to year and then daynight_category
pred$pred <- pred$pred[order(pred$pred$year, pred$pred$daynight_category),]

#plot the predictions
#make points transparent by using alpha = 0.5
#add confidence interval by using geom_ribbon()
#change the number of columns in facet_grid() so that the figures are not too narrow
#remove legend
#make figure wider and shorter
ggplot(data = pred$pred)+
  geom_line(aes(x = t, y = estimate, color = "predicted"))+
  geom_point(aes(x = t, y = y, color = "observed"), alpha = 0.5)+
  geom_ribbon(aes(x = t, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  scale_color_manual(name = "Legend", values = c("predicted" = "#257085", "observed" = "black"))+
  labs(x = "Time", y = "log(CPUE + 1)", title = "Predicted vs Observed Chinook subyearlings, Puyallup River")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~year + daynight_category,  scales = "free_y")+
  theme_bw()+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid"
    ),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "none"
  )

#save figure with ggsave()

ggsave(here("puyallup", "output", "predicted_vs_observed.png"), width = 15, height = 5)



