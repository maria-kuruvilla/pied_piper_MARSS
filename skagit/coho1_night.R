#Goal - repeat MARSS analysis for coho data in skagit
#only the night values


# Load libraries


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




out.tab.all.years_night <- NULL
fits.all.years_night <- NULL
nyears = num_years*2
for(kk in c(1,3,4,6,9)){
  c = covariates_coho1_skagit_night[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_coho1_skagit_night)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-15)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0))
  fit <- MARSS(subset_coho_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years_night=rbind(out.tab.all.years_night,out)
  fits.all.years_night=c(fits.all.years_night,list(fit))
}

out.tab.all.years_night
#atu_solstice have the least aicc


out.tab.all.years_night$deltaAICc <- (out.tab.all.years_night$AICc - 
                                        min(out.tab.all.years_night$AICc))

out.tab.all.years_night <- out.tab.all.years_night[order(out.tab.all.years_night$AICc),]


#now try with all covariates

list_combinations <- get_covariate_combinations(1:6)
out.tab.hatchery.all.years_night <- NULL
fits.hatchery.all.years_night <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  atu = 0
  flow= 0
  lunar_phase = 0
  temperature_difference = 0
  flow_difference = 0
  hatchery = 0
  
  
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
    
    c = rbind(c,covariates_coho1_skagit_night[((1+(k-1)*nyears):(k*nyears)),])
    name_long = rownames(covariates_coho1_skagit_night)[1+(k-1)*nyears]
    name = paste(name, substr(name_long,1,nchar(name_long)-15))
    
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
  fit <- MARSS(subset_coho_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, atu = atu, flow = flow, lunar_phase = lunar_phase,
                 temperature_difference = temperature_difference, 
                 flow_difference = flow_difference,
                 hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.hatchery.all.years_night=rbind(out.tab.hatchery.all.years_night,out)
  fits.hatchery.all.years_night=c(fits.hatchery.all.years_night,list(fit))
  
  
}
out.tab.hatchery.all.years_night$deltaAICc <- out.tab.hatchery.all.years_night$AICc - min(out.tab.hatchery.all.years_night$AICc)
min.AICc <- order(out.tab.hatchery.all.years_night$AICc)
out.tab.hatchery.ordered_night <- out.tab.hatchery.all.years_night[min.AICc, ]
out.tab.hatchery.ordered_night

#save the dataframe
write.csv(out.tab.hatchery.ordered_night, 
          file = here("skagit","output","model_selection_skagit_coho_night.csv"))


ci_night <- tidy(fits.hatchery.all.years_night[[44]])


ggplot(ci_night[c(27:30),], aes(x = c("ATU","Flow","Lunar phase",
                                "Hatchery"),
                          y = estimate, 
                          ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  ggtitle("Coho yearlings, Skagit") + 
  theme_bw() +
  theme(plot.title = element_text(size = 20))+
  theme(axis.text.x=element_text(size=14),axis.title.y=element_text(size=14))
#save figure

ggsave(here("skagit","output","coho_night_estimates.png"), width = 10, height = 8, units = "in")


#relative importance of variables

out.tab.hatchery.ordered_night$deltaAICc <- NULL
out.tab.hatchery.ordered_night$rel.LL <- NULL
out.tab.hatchery.ordered_night$weights <- NULL


weights <- akaike.weights(out.tab.hatchery.ordered_night$AICc)

out.tab.hatchery.ordered_night$deltaAICc <- weights$deltaAIC
out.tab.hatchery.ordered_night$rel.LL <- weights$rel.LL
out.tab.hatchery.ordered_night$weights <- weights$weights

#round the deltaAICc, rel.LL, and weights to 2 decimal places

out.tab.hatchery.ordered_night$deltaAICc <- round(out.tab.hatchery.ordered_night$deltaAICc,2)
out.tab.hatchery.ordered_night$rel.LL <- round(out.tab.hatchery.ordered_night$rel.LL,2)
out.tab.hatchery.ordered_night$weights <- round(out.tab.hatchery.ordered_night$weights,2)



out.tab.hatchery.ordered_night$cumulative_weights <- cumsum(out.tab.hatchery.ordered_night$weights)

relative_importance_atu <- sum(out.tab.hatchery.ordered_night$weights[out.tab.hatchery.ordered_night$atu==1])
relative_importance_temperature_difference <- sum(out.tab.hatchery.ordered_night$weights[out.tab.hatchery.ordered_night$temperature_difference==1])
relative_importance_flow <- sum(out.tab.hatchery.ordered_night$weights[out.tab.hatchery.ordered_night$flow==1])
relative_importance_lunar_phase <- sum(out.tab.hatchery.ordered_night$weights[out.tab.hatchery.ordered_night$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab.hatchery.ordered_night$weights[out.tab.hatchery.ordered_night$hatchery==1])
relative_importance_flow_diff <- sum(out.tab.hatchery.ordered_night$weights[out.tab.hatchery.ordered_night$flow_difference==1])

riv_coho_skagit_night <- data.frame(variable = c("atu",
                                          "temperature difference",
                                          "flow",
                                          "lunar phase",
                                          "flow difference",
                                          "hatchery"),
                             relative_importance = c(relative_importance_atu,
                                                     relative_importance_temperature_difference,
                                                     relative_importance_flow,
                                                     relative_importance_lunar_phase,
                                                     relative_importance_flow_diff,
                                                     relative_importance_hatchery))


#sort according to relative importance

riv_coho_skagit_night <- riv_coho_skagit_night[order(riv_coho_skagit_night$relative_importance, decreasing = TRUE),]
riv_coho_skagit_night$relative_importance <- round(riv_photo_diff$relative_importance,1)
riv_coho_skagit_night

#save the dataframe

write.csv(riv_coho_skagit_night, 
          file = here("skagit","output","relative_importance_skagit_coho_night.csv"))
