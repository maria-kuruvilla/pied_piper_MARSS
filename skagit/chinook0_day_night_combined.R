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

#add chinook0_hatchery_num for all values of trap and daytime_category
#add chinook0_wild_num for all values of trap and daytime_category
#take mean of temp, flow, lunar_phase, photoperiod, photo_diff, temp_diff, flow_diff, atu_solstice
#for all values of trap and daytime_category
#drop trap and daytime_category
data <- data_day_night %>% 
  select(chinook0_hatchery_num, chinook0_wild_num, trap, Date, daytime_category, temp,
         flow, lunar_phase, photoperiod, photo_diff, temp_diff, flow_diff, atu_solstice,
         doy, year, In) %>% 
  group_by(Date) %>%
  summarize(chinook0_hatchery_num = if(all(is.na(chinook0_hatchery_num))) NA else sum(chinook0_hatchery_num, na.rm = TRUE),
         chinook0_wild_num = if(all(is.na(chinook0_wild_num))) NA else sum(chinook0_wild_num, na.rm = TRUE),
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
  
#plot to check that it worked

data_day_night %>% 
  ggplot(aes(x = doy, y = chinook0_hatchery_num)) +
  geom_point() +
  geom_line() +
  facet_wrap(~trap + year, scales = "free_y")

data %>% 
  ggplot(aes(x = doy, y = chinook0_hatchery_num)) +
  geom_point() +
  geom_line() +
  facet_wrap(~year, scales = "free_y")  

#in data, interpolate values of chinook0_hatchery_num, temp, flow, lunar_phase, 
#photoperiod, photo_diff, temp_diff, flow_diff, atu_solstice

#make new column for chinook0_hatchery_perhour and chinook0_wild_perhour

data <- data %>% 
  mutate(chinook0_hatchery_perhour = chinook0_hatchery_num/In,
         chinook0_wild_perhour = chinook0_wild_num/In,
         chinook0_hatchery_perhour_inp = na.approx(chinook0_hatchery_perhour, na.rm = FALSE),
         temp = na.approx(temp),
         flow = na.approx(flow),
         lunar_phase = na.approx(lunar_phase),
         photoperiod = na.approx(photoperiod),
         photo_diff = na.approx(photo_diff),
         temp_diff = na.approx(temp_diff),
         flow_diff = na.approx(flow_diff),
         atu_solstice = na.approx(atu_solstice))

ggplot(data) +
  geom_line(aes(x = doy, y = chinook0_hatchery_perhour), color = "cadetblue") +
  geom_line(aes(x = doy, y = chinook0_wild_perhour), color = "salmon") +
  facet_wrap(~year, scales = "free_y") +
  #add dotted line for chinook0_hatchery_perhour_inp
  geom_line(aes(x = doy, y = chinook0_hatchery_perhour_inp), color = "cadetblue", linetype = "dotted")

#get data in marss format
#then run marss analysis with separate year effect for hatchery

lm_temp_day <- lm(temp ~ 1+photoperiod, data = data)

data <- data %>% add_residuals(lm_temp_day)

covariates_chinook0_skagit_combined <- arrange(data,doy) %>%
  filter(doy >150 & doy < 190) %>%
  dplyr::select(year,doy,  temp, flow, photoperiod, atu_solstice,
                lunar_phase, resid, temp_diff, flow_diff, photo_diff,
                chinook0_hatchery_perhour_inp) %>%
  pivot_wider(names_from = c(year), values_from = c(
    temp, flow, photoperiod, atu_solstice, lunar_phase,
    resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_perhour_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2022-2010
num_rows = num_years
num_covariates = 10
total_covariates = dim(covariates_chinook0_skagit_combined)[1]

#scaling and centering
for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_chinook0_skagit_combined[i,] = scale(covariates_chinook0_skagit_combined[i,])[,1]
}


#checking
for(i in 1:(total_covariates)){
  print(sum(covariates_chinook0_skagit_combined[i,]))
}

#just scale
for(i in (total_covariates/2 + 1):(total_covariates)){
  if(sum(covariates_chinook0_skagit_combined[i,]) != 0){
    covariates_chinook0_skagit_combined[i,] = scale(covariates_chinook0_skagit_combined[i,], center = FALSE, scale= TRUE)[,1]
  }
}

#plot just the chinook0_hatchery_perhour_inp rows from covariates_chinook0_skagit_combined
#facet wrap by year
#see if it looks right

#scale the chinook0_hatchery_perhour_inp for every year
chinook0_scale_check <- arrange(data,doy) %>%
  filter(doy >150 & doy < 190) %>%
  dplyr::select(year,doy,chinook0_hatchery_perhour_inp) %>% 
  pivot_wider(names_from = year, values_from = chinook0_hatchery_perhour_inp) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix()

#scale each column
for(i in 1:dim(chinook0_scale_check)[2]){
  chinook0_scale_check[,i] = scale(chinook0_scale_check[,i], center = FALSE)[,1]
}

#plot
#make rownames the doy
data_long_scaled <- chinook0_scale_check %>% 
  as.data.frame() %>%
  rownames_to_column(var = "doy") %>%
  pivot_longer(cols = -doy,names_to = "year", values_to = "chinook0_hatchery_perhour_inp")

ggplot(data_long_scaled) +
  geom_line(aes(x = doy, y = chinook0_hatchery_perhour_inp, group = year, color = year))+
  facet_wrap(~year, scales = "free_y")

#subset response variable
subset_chinook_summer_perhour_combined <- arrange(data,doy) %>%
  filter(doy > 150 & doy < 190) %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_chinook_summer_perhour_combined)[1]){
  subset_chinook_summer_perhour_combined[i,] = scale(subset_chinook_summer_perhour_combined[i,])[,1]
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
c = covariates_chinook0_skagit_combined[((1+(kk-1)*nyears):(kk*nyears)),]
name_long = rownames(covariates_chinook0_skagit_combined)[1+(kk-1)*nyears]
name_individual = substr(name_long,1,nchar(name_long)-9)
print(name_individual)
fit.model = c(list(c= c), mod_list_skagit(nyears,1,1,TRUE))
fit_year_effect <- MARSS(subset_chinook_summer_perhour_combined, model=fit.model, silent = TRUE, method = "BFGS",
                         control=list(maxit=2000))
ci_year_effect <- tidy(fit_year_effect)

fit_year_effect$AICc #2550

#goal - plot the exact data that is going into the model
#try the same with coho

#plot all the rows of c

ggplot(as.data.frame(t(c))) +
  geom_line(aes(x = 151:189, y = chinook0_hatchery_perhour_inp_2022), color = "cadetblue")+
  geom_line(data = data_frame_scaled,aes(x = doy, y = year_2022), color = "salmon")


data_frame_scaled <- as.data.frame(t(subset_chinook_summer_perhour_combined))
#rename the columns to year_2010, year_2012, etc
colnames(data_frame_scaled) <- paste0("year_",colnames(data_frame_scaled))
data_frame_scaled$doy <- as.numeric(rownames(data_frame_scaled))

ggplot(data_frame_scaled)+
  geom_line(aes(x = doy, y = year_2010))

autoplot(fit_year_effect)




out.tab.all.years.combined <- NULL
fits.all.years.combined <- NULL
nyears = num_years
for(kk in c(1,2,3,4,5,6,7,8,9)){
  c = covariates_chinook0_skagit_combined[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_chinook0_skagit_combined)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-5)
  print(name_individual)
  fit.model = c(list(c= c), mod_list_skagit(nyears,1,0, FALSE))
  fit <- MARSS(subset_chinook_summer_perhour_combined, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years.combined=rbind(out.tab.all.years.combined,out)
  fits.all.years.combined=c(fits.all.years.combined,list(fit))
}

out.tab.all.years.combined
#calculate delta AICc and order by that
out.tab.all.years.combined$delta.AICc = out.tab.all.years.combined$AICc - min(out.tab.all.years.combined$AICc)
out.tab.all.years.combined = out.tab.all.years.combined[order(out.tab.all.years.combined$delta.AICc),]
#round to 2 decimal places
out.tab.all.years.combined$delta.AICc = round(out.tab.all.years.combined$delta.AICc,2)
out.tab.all.years.combined

#save the dataframe
write.csv(out.tab.all.years.combined, file = here("skagit", "output","model_selection_correlated_covariates_chinook0.csv"))


#try model with resid,photoperiod,flow,lunar_phase, temp_diff,flow_diff,hatchery



nyears = num_years
c <- NULL
name <- NULL
for(kk in c(2,3,5,6,7,8,10)){
  c = rbind(c,covariates_chinook0_skagit_combined[((1+(kk-1)*nyears):(kk*nyears)),])
  name_long = rownames(covariates_chinook0_skagit_combined)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-5)
  name = paste(name, name_individual)
  print(name_individual)
}
fit.model = c(list(c= c), mod_list_skagit(nyears,7,0, FALSE))
fit <- MARSS(subset_chinook_summer_perhour_combined, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
out=data.frame(c=name_individual,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)
# out.tab.all.years.combined=rbind(out.tab.all.years.combined,out)
# fits.all.years.combined=c(fits.all.years.combined,list(fit))
# 
# 
# out.tab.all.years.combined
fit
tidy(fit)

#try model selection


list_combinations <- get_covariate_combinations(1:7)
out.tab.hatchery.all.years.combined <- NULL
fits.hatchery.all.years.combined <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  residuals = 0
  flow = 0
  flow_difference = 0
  hatchery = 0
  temperature_difference = 0
  photoperiod = 0
  lunar_phase = 0
  for(j in covariates){
    if(j == 1){
      k = 6
      residuals = 1
    }
    else if(j==2){
      k = 2
      flow = 1
    }
    else if(j==3){
      k = 3
      photoperiod = 1
      
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
      k = 5
      lunar_phase = 1
    }
    else if(j==7){
      k = 10
      hatchery = 1
    }
    
    c = rbind(c,covariates_chinook0_skagit_combined[((1+(k-1)*nyears):(k*nyears)),])
    name_long = rownames(covariates_chinook0_skagit_combined)[1+(k-1)*nyears]
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
  fit <- MARSS(subset_chinook_summer_perhour_combined, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, residuals=residuals, flow=flow, flow_difference=flow_difference,
                 temperature_difference=temperature_difference, photoperiod = photoperiod,
                 lunar_phase = lunar_phase ,hatchery=hatchery, 
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.hatchery.all.years.combined=rbind(out.tab.hatchery.all.years.combined,out)
  fits.hatchery.all.years.combined=c(fits.hatchery.all.years.combined,list(fit))
  
  
}
out.tab.hatchery.all.years.combined$deltaAICc <- out.tab.hatchery.all.years.combined$AICc - min(out.tab.hatchery.all.years.combined$AICc)
min.AICc <- order(out.tab.hatchery.all.years.combined$AICc)
out.tab.hatchery.ordered.combined <- out.tab.hatchery.all.years.combined[min.AICc, ]
out.tab.hatchery.ordered.combined

autoplot(fits.hatchery.all.years.combined[[9]])
ci_best_chinook_skagit <- tidy(fits.hatchery.all.years.combined[[9]])

#plot the estimates of the best model

ggplot(ci_best_chinook_skagit[c(15:16),], 
       aes(x = c("Residuals", "Photoperiod"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect", title = "Skagit, Chinook sub yearlings") +
  theme_bw()+
  theme(axis.text.x=element_text(size=18),axis.title.y=element_text(size=18),
        title = element_text(size = 18))
ggsave(here("skagit","output","chinook_all_years_effects_estimate_combined.png"), 
       width = 10, height = 8, units = "in", dpi = 300)
  


#calculate relative variable importance



out.tab.hatchery.ordered.combined$deltaAICc <- NULL
out.tab.hatchery.ordered.combined$rel.LL <- NULL
out.tab.hatchery.ordered.combined$weights <- NULL


weights <- akaike.weights(out.tab.hatchery.ordered.combined$AICc)

out.tab.hatchery.ordered.combined$deltaAICc <- weights$deltaAIC
out.tab.hatchery.ordered.combined$rel.LL <- weights$rel.LL
out.tab.hatchery.ordered.combined$weights <- weights$weights


min.AICc <- order(out.tab.hatchery.ordered.combined$AICc)
out.tab.hatchery.ordered.combined <- out.tab.hatchery.ordered.combined[min.AICc, ]
out.tab.hatchery.ordered.combined

out.tab.hatchery.ordered.combined$cumulative_weights <- cumsum(out.tab.hatchery.ordered.combined$weights)

relative_importance_residuals <- sum(out.tab.hatchery.ordered.combined$weights[out.tab.hatchery.ordered.combined$residuals==1])
relative_importance_photo <- sum(out.tab.hatchery.ordered.combined$weights[out.tab.hatchery.ordered.combined$photoperiod==1])
relative_importance_temperature_difference <- sum(out.tab.hatchery.ordered.combined$weights[out.tab.hatchery.ordered.combined$temperature_difference==1])
relative_importance_flow <- sum(out.tab.hatchery.ordered.combined$weights[out.tab.hatchery.ordered.combined$flow==1])
relative_importance_lunar_phase <- sum(out.tab.hatchery.ordered.combined$weights[out.tab.hatchery.ordered.combined$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab.hatchery.ordered.combined$weights[out.tab.hatchery.ordered.combined$hatchery==1])
relative_importance_flow_diff <- sum(out.tab.hatchery.ordered.combined$weights[out.tab.hatchery.ordered.combined$flow_difference==1])

riv_chinook_skagit <- data.frame(variable = c("photoperiod",
                                              "residuals",
                                              "temperature difference",
                                              "flow",
                                              "lunar phase",
                                              "flow difference",
                                              "hatchery"),
                                 relative_importance = c(relative_importance_photo,
                                                         relative_importance_residuals,
                                                         relative_importance_temperature_difference,
                                                         relative_importance_flow,
                                                         relative_importance_lunar_phase,
                                                         relative_importance_flow_diff,
                                                         relative_importance_hatchery))


#sort according to relative importance
riv_chinook_skagit[order(riv_chinook_skagit$relative_importance, decreasing = TRUE),]

riv_chinook_skagit <- riv_chinook_skagit[order(riv_chinook_skagit$relative_importance, decreasing = TRUE),]
riv_chinook_skagit$relative_importance <- round(riv_chinook_skagit$relative_importance,1)

riv_chinook_skagit

#save the relative importance table
write.csv(riv_chinook_skagit, here("output","riv_chinook0_skagit.csv"), row.names = FALSE)

out.tab.hatchery.ordered.combined
#clean up table by reducing number of digits to 2 and by only keeping
#the variables of interest


out.tab.hatchery.ordered.combined$rel.LL <- round(
  out.tab.hatchery.ordered.combined$rel.LL,2)
out.tab.hatchery.ordered.combined$weights <- round(
  out.tab.hatchery.ordered.combined$weights,2)
out.tab.hatchery.ordered.combined$cumulative_weights <- round(
  out.tab.hatchery.ordered.combined$cumulative_weights,2)

#save the table
write.csv(out.tab.hatchery.ordered.combined, here("skagit","output","chinook0_model_selection_w_weights.csv"), row.names = FALSE)



#Maybe try all rivers combined model