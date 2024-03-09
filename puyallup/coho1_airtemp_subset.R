#goal - run marss analysis for coho1 in Puyallup
#with air temperature data but without 2009 and 2021

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



ggplot(data %>% 
         filter(year != 2009, year != 2021))+
  geom_point(aes(x = doy, y = coho1_wild_perhour_night), alpha = 0.5)+
  geom_point(aes(x = doy, y = coho1_hatchery_perhour_night), alpha = 0.5, color = "cadetblue")+
  facet_wrap(~year, scales = "free_y")+
  scale_x_continuous(limit = c(90,160))+
  scale_y_log10()



covariates_coho1_puyallup_w_temp <- arrange(data,doy) %>%
  filter(doy >90 & doy <= 160, year != 2009, year != 2021) %>%
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

num_years = 2021-2004+1 - 2
num_rows = num_years*2
total_covariates = dim(covariates_coho1_puyallup_w_temp)[1]
nyears = num_rows



for(i in 1:(num_rows*6)){ # everything except diffs and hatchery
  covariates_coho1_puyallup_w_temp[i,] = scale(covariates_coho1_puyallup_w_temp[i,])[,1]
}

#just scale

for(i in (num_rows*6 + 1):(total_covariates)){
  covariates_coho1_puyallup_w_temp[i,] = scale(covariates_coho1_puyallup_w_temp[i,], center = FALSE, scale= TRUE)[,1]
}


subset_coho_summer_perhour <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160, year != 2009, year != 2021) %>%
  mutate(log.value_day = log(coho1_wild_perhour_day + 1), 
         log.value_night = log(coho1_wild_perhour_night + 1)) %>%
  dplyr::select(log.value_day, log.value_night ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_day, log.value_night)) %>%
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

nyears = num_years*2
fit.model_equal_q = mod_list(nyears,0,0, FALSE, FALSE)
fit_equal_q <- MARSS(subset_coho_summer_perhour, model=fit.model_equal_q, silent = TRUE, method = "BFGS",
                     control=list(maxit=2000))
fit_equal_q$AICc #4829

fit.model_unequal_q = mod_list(nyears,0,0, FALSE, TRUE)
fit_unequal_q <- MARSS(subset_coho_summer_perhour, model=fit.model_unequal_q, silent = TRUE, method = "BFGS",
                       control=list(maxit=2000))
fit_unequal_q$AICc #4822


out.tab.all.years.coho <- NULL
fits.all.years.coho <- NULL
nyears = num_years*2
for(kk in c(2,3,5,6,8,10)){
  c = covariates_coho1_puyallup_w_temp[((1+(kk-1)*num_rows):(kk*num_rows)),]
  name_long = rownames(covariates_coho1_puyallup_w_temp)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0, FALSE, TRUE))
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  out=data.frame(c=name_individual, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.all.years.coho=rbind(out.tab.all.years.coho,out)
  fits.all.years.coho=c(fits.all.years.coho,list(fit))
}

out.tab.all.years.coho$delta_AICc = out.tab.all.years.coho$AICc - min(out.tab.all.years.coho$AICc)
out.tab.all.years.coho_ordered = out.tab.all.years.coho[order(out.tab.all.years.coho$delta_AICc),]
out.tab.all.years.coho_ordered
#make column for delta aicc
#order by delta aicc
out.tab.all.years.coho_ordered$delta_AICc = out.tab.all.years.coho_ordered$AICc - min(out.tab.all.years.coho_ordered$AICc)
out.tab.all.years.coho_ordered
#save table
write.csv(out.tab.all.years.coho_ordered, file = here("puyallup","output",
                                                      "coho_correlated_covariates_selection_atu_subset.csv"))

#checking just the best model
c<-NULL
for(kk in c(1,3,6,7,11)){
  c = rbind(c,covariates_coho1_puyallup_w_temp[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_coho1_puyallup_w_temp)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  
}
fit.model = c(list(c= c), mod_list(nyears,5,1, FALSE, TRUE))
fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
fit




checkpoint_file <- here("puyallup","output","checkpoint.txt")

# Function to read the last executed iteration from the file
get_last_iteration <- function() {
  if (file.exists(checkpoint_file)) {
    last_iteration <- as.numeric(readLines(checkpoint_file))
  } else {
    last_iteration <- 1
  }
  return(last_iteration)
}

# Function to update the last executed iteration in the file
update_checkpoint <- function(iteration) {
  writeLines(as.character(iteration), checkpoint_file)
}


atu = 0
photoperiod = 0
flow = 0
lunar_phase = 0
hatchery = 0
flow_difference = 0
temp_difference = 0
nyears = num_rows
list_combinations <- get_covariate_combinations(1:7)
out.tab.hatchery.all.years.coho <- NULL
fits.hatchery.all.years.coho <- list()
num_rows
start_iteration <- get_last_iteration()
for(i in start_iteration:length(list_combinations)){
  atu = 0
  photoperiod = 0
  flow = 0
  lunar_phase = 0
  hatchery = 0
  flow_difference = 0
  temp_difference = 0
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  for(j in covariates){
    if(j == 1){
      k = 6
      atu = 1
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
      k = 7
      flow_difference = 1
    }
    else if(j==6){
      k = 9
      temp_difference = 1
    }
    else if(j==7){
      k = 11
      hatchery = 1
    }
    
    c = rbind(c,covariates_coho1_puyallup_w_temp[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_coho1_puyallup_w_temp)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-5))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(k==11){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery,FALSE,TRUE))
  }
  else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(nyears,c_num,has_hatchery,FALSE,TRUE))
  }
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, atu = atu, photoperiod = photoperiod,
                 flow = flow, lunar_phase = lunar_phase, flow_difference = flow_difference,
                 temp_difference = temp_difference,
                 hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab.hatchery.all.years.coho=rbind(out.tab.hatchery.all.years.coho,out)
  fits.hatchery.all.years.coho=c(fits.hatchery.all.years.coho,list(fit))
  update_checkpoint(i + 1)
  
}
out.tab.hatchery.all.years.coho$deltaAICc <- out.tab.hatchery.all.years.coho$AICc - min(out.tab.hatchery.all.years.coho$AICc)
min.AICc <- order(out.tab.hatchery.all.years.coho$AICc)
out.tab.hatchery.ordered_coho <- out.tab.hatchery.all.years.coho[min.AICc, ]
out.tab.hatchery.ordered_coho

fits.hatchery.ordered_coho_best <- fits.hatchery.all.years.coho[[103]]

ci_fits.hatchery.ordered_coho_best <- tidy(fits.hatchery.ordered_coho_best)
ci_fits.hatchery.ordered_coho_best

ggplot(ci_fits.hatchery.ordered_coho_best[36:41,], 
       aes(x = c("ATU", "Flow", "Photoperiod",
                 "Flow\n difference",
                 "Hatchery,\n day", "Hatchery,\n night"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect", title = "Puyallup, Coho yearlings") +
  theme_bw()+
  theme(axis.text.x=element_text(size=18),axis.title.y=element_text(size=18),
        title = element_text(size = 18))


coho1_puyallup_plot <- ggplot(ci_fits.hatchery.ordered_coho_best[36:41,], 
                              aes(x = c("ATU", "Flow", "Photoperiod",
                                        "Flow\ndifference","Hatchery,\nday","Hatchery,\nnight"),
                                  y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")  +
  geom_rect(aes(xmin = 2.5, xmax = 4.5, ymin = -0.28, ymax = 0.28), col = "cadetblue", alpha = 0.0) +
  geom_text(aes(x = 6.35, y = -0.225, label = "Puyallup"), size = 10)+
  labs(x = "", y = "Estimate of effect",
  )+
  theme_classic()+
  theme(axis.title.y=element_text(size=24),
        axis.title.x=element_text(size=24),
        axis.text.y = element_text(size = 24),
        axis.text.x=element_text(size=24)) +
  scale_y_continuous(breaks = seq(-0.2,0.2,0.2), limits = c(-0.3,0.3), n.breaks = 3)+
  scale_x_discrete(limits = c("ATU", "Flow","Hatchery,\nnight","Hatchery,\nday", "Photoperiod",
                              "Flow\ndifference"
  )) +
  coord_flip()
#guide = guide_axis(n.dodge = 2))

coho1_puyallup_plot

ggsave(here("puyallup","output",
            "puyallup_coho_best_model_estimates_new.png"), width = 10, height = 5, units = "in")



fit.model = mod_list(nyears,0,0,FALSE,TRUE)
fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
atu = 0
photoperiod = 0
flow = 0
lunar_phase = 0
hatchery = 0
flow_difference = 0
temp_difference = 0

out=data.frame(c="No covariates", atu = atu, photoperiod = photoperiod,
               flow = flow, lunar_phase = lunar_phase, flow_difference = flow_difference,
               temp_difference = temp_difference,
               hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)
out.tab.hatchery.all.years.coho$deltaAICc <- NULL
out.tab.hatchery.all.years.coho=rbind(out.tab.hatchery.all.years.coho,out)
fits.hatchery.all.years.coho=c(fits.hatchery.all.years.coho,list(fit))

out.tab.hatchery.all.years.coho$deltaAICc <- out.tab.hatchery.all.years.coho$AICc - min(out.tab.hatchery.all.years.coho$AICc)

#calculate weights

weights <- akaike.weights(out.tab.hatchery.all.years.coho$AICc)

out.tab.hatchery.all.years.coho$deltaAICc <- weights$deltaAIC
out.tab.hatchery.all.years.coho$rel.LL <- weights$rel.LL
out.tab.hatchery.all.years.coho$weights <- weights$weights

min.AICc <- order(out.tab.hatchery.all.years.coho$AICc)
out.tab.hatchery.all.years.coho.ordered <- out.tab.hatchery.all.years.coho[min.AICc, ]
out.tab.hatchery.all.years.coho.ordered

out.tab.hatchery.all.years.coho.ordered$cumulative_weights <- cumsum(out.tab.hatchery.all.years.coho.ordered$weights)

relative_importance_atu <- sum(out.tab.hatchery.all.years.coho$weights[out.tab.hatchery.all.years.coho$atu==1])
relative_importance_flow <- sum(out.tab.hatchery.all.years.coho$weights[out.tab.hatchery.all.years.coho$flow==1])
relative_importance_lunar_phase <- sum(out.tab.hatchery.all.years.coho$weights[out.tab.hatchery.all.years.coho$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab.hatchery.all.years.coho$weights[out.tab.hatchery.all.years.coho$hatchery==1])
relative_importance_flow_diff <- sum(out.tab.hatchery.all.years.coho$weights[out.tab.hatchery.all.years.coho$flow_difference==1])
relative_importance_photoperiod <- sum(out.tab.hatchery.all.years.coho$weights[out.tab.hatchery.all.years.coho$photoperiod==1])
relative_importance_temp_diff <- sum(out.tab.hatchery.all.years.coho$weights[out.tab.hatchery.all.years.coho$temp_difference==1])

riv_puyallup_coho <- data.frame(variable = c("atu",
                                             "photoperiod",
                                             "lunar phase",
                                             "temp difference",
                                             "flow",
                                             "flow difference",
                                             "hatchery"),
                                relative_importance = c(relative_importance_atu,
                                                        relative_importance_photoperiod,
                                                        relative_importance_lunar_phase,
                                                        relative_importance_temp_diff,
                                                        relative_importance_flow,
                                                        relative_importance_flow_diff,
                                                        relative_importance_hatchery))

riv_puyallup_coho[order(riv_puyallup_coho$relative_importance, decreasing = TRUE),]

riv_puyallup_coho <- riv_puyallup_coho[order(riv_puyallup_coho$relative_importance, decreasing = TRUE),]
riv_puyallup_coho$relative_importance <- round(riv_puyallup_coho$relative_importance,1)

riv_puyallup_coho

#save the dataframe riv_puyallup_coho

write.csv(riv_puyallup_coho, file = here("puyallup",
                                         "output","riv_puyallup_coho_new.csv"))


out.tab.hatchery.all.years.coho.ordered$weights <- round(out.tab.hatchery.all.years.coho.ordered$weights,2)
out.tab.hatchery.all.years.coho.ordered$deltaAICc <- round(out.tab.hatchery.all.years.coho.ordered$deltaAICc,2)
out.tab.hatchery.all.years.coho.ordered$AICc <- round(out.tab.hatchery.all.years.coho.ordered$AICc,2)
out.tab.hatchery.all.years.coho.ordered$rel.LL <- round(out.tab.hatchery.all.years.coho.ordered$rel.LL,2)
out.tab.hatchery.all.years.coho.ordered$logLik <- round(out.tab.hatchery.all.years.coho.ordered$logLik,2)
out.tab.hatchery.all.years.coho.ordered$cumulative_weights <- round(out.tab.hatchery.all.years.coho.ordered$cumulative_weights,2)

out.tab.hatchery.all.years.coho.ordered

write.csv(out.tab.hatchery.all.years.coho.ordered, file = here("puyallup",
                                                          "output",
                                                          "model_selection_coho_new.csv"))
