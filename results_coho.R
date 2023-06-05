#going to try coho again with fish per hour data, hatchery in the D covariate and without converting
# 0 to na


#packages
library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)
library(broom)
library(ggfortify)
library(modelr)
library(MASS)


data <- read.csv(here("data","dungeness_subset_covariates.csv"), header = TRUE)

data$coho1_hatchery_perhour_interpolate <- na.approx(data$coho1_hatchery_perhour)

subset_perhour <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 120 & doy <= 160) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

#scale the values

for(i in 1:dim(subset_perhour)[1]){
  subset_perhour[i,] = scale(subset_perhour[i,])[,1]
}






covariates_all_years <- arrange(data,doy) %>%
  filter(year != 2015 & doy >120 & doy <= 160) %>%
  dplyr::select(year,doy, temp, atu_april, photoperiod, lunar_phase, flow, coho1_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = year, values_from = c(
    temp, atu_april, photoperiod, lunar_phase, flow, coho1_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 76:90){
  covariates_all_years[i,] = scale(covariates_all_years[i,], center = FALSE, scale= TRUE)[,1]
}







#base model with no covariates

mod.list_base <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal"
)


fit <- MARSS(subset_perhour, model=mod.list_base, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))


fit$AICc

autoplot(fit)


#hatchery as d covariate

mod.list_h_d <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  D = "diagonal and unequal",
  d = covariates_all_years[76:90,] 
)


fit_h_d <- MARSS(subset_perhour, model=mod.list_h_d, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))


fit_h_d$AICc

autoplot(fit_h_d)

ci_h_d <- tidy(fit_h_d)


ggplot(ci_h_d[18:32,], aes(x = term, y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Estimate of Hatchery effect") +
  ggtitle("Effect of hatchery Coho") + 
  theme(plot.title = element_text(size = 20))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))


#residuals of temp
lm_temp <- lm(temp ~ 1+photoperiod, data)
summary(lm_temp)
residuals(lm_temp)
length(residuals(lm_temp))
length(data$temp)


data <- data %>% add_residuals(lm_temp)


covariates_all_years <- arrange(data,doy) %>%
  filter(year != 2015 & doy >120 & doy <= 160) %>%
  dplyr::select(year,doy, temp, atu_april, photoperiod, lunar_phase, resid, flow, coho1_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = year, values_from = c(
    temp, atu_april, photoperiod, lunar_phase, resid, flow, coho1_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


for(i in 1:(dim(covariates_all_years)[1]-45)){
  covariates_all_years[i,] = scale(covariates_all_years[i,])[,1]
}


for(i in 76:90){
  covariates_all_years[i,] = scale(covariates_all_years[i,])[,1]
}

for(i in 91:105){
  covariates_all_years[i,] = scale(covariates_all_years[i,], center = FALSE, scale= TRUE)[,1]
}




#loop through and fit models with temp, atu, photoperiod


mod.list_base <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  D = "diagonal and unequal",
  d = covariates_all_years[91:105,] 
)



mod.list_h_d <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = "diagonal and equal",
  D = "diagonal and unequal",
  d = covariates_all_years[91:105,] 
)


#putting it into a loop

c.models <- list( 
  c1 = "zero", # No covariates
  c2 = covariates_all_years[1:15,], # Temperature
  c3 = covariates_all_years[16:30,], # ATU April
  c4 = covariates_all_years[31:45,]  #Photoperiod
  
)
names(c.models) <- c("No covariates", 
                     "Temperature", 
                     "ATU", 
                     "Photoperiod")


out.tab <- NULL
fits <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(names(c.models)[i] == "No covariates"){
    
    fit.model = mod.list_base
  }
  
  else{
    fit.model = c(list(c= c.models[[i]]), mod.list_h_d)
    
  }
  
  
    
  
  
  
  fit <- MARSS(subset_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d="Hatchery",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab=rbind(out.tab,out)
  fits=c(fits,list(fit))
  
  
}



min.AICc <- order(out.tab$AICc)
out.tab.1 <- out.tab[min.AICc, ]
out.tab.1$DeltaAICc <- out.tab.1$AICc - min(out.tab.1$AICc)
out.tab.1



# add other covariates to hatchery and photoperiod



c.models <- list( 
  
  
  c1 = covariates_all_years[c(46:60,31:45),],  #Lunar phase and photoperiod
  c2 = covariates_all_years[c(61:75,31:45),],  #Resid and photoperiod
  c3 = covariates_all_years[c(76:90,31:45),]  #Flow and photoperiod
)
names(c.models) <- c("Lunar Phase, Photoperiod", "Temperature residuals, Photoperiod", "Flow, Photoperiod")

mod.list_h_d_two <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p"),
             15,30, byrow = TRUE),
  D = "diagonal and unequal",
  d = covariates_all_years[91:105,] 
)



out.tab2 <- NULL
fits2 <- list()
for(i in 1:length(c.models)){
  
  fit.model = c(list(c= c.models[[i]]), mod.list_h_d_two)
  
  fit <- MARSS(subset_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d="Hatchery",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab2=rbind(out.tab2,out)
  fits2=c(fits2,list(fit))
  
  
}



min.AICc <- order(out.tab2$AICc)
out.tab.2 <- out.tab2[min.AICc, ]
out.tab.2$DeltaAICc <- out.tab.2$AICc - min(out.tab.2$AICc)
out.tab.2



#flow, photoperiod, hatchery are best. Now will add lunar phase and temp resid to it

c.models <- list( 
  
  
  c1 = covariates_all_years[c(46:60,76:90,31:45),],  #Lunar phase, flow, and photoperiod
  c2 = covariates_all_years[c(61:75,76:90,31:45),]  #Resid, flow, and photoperiod
  
)
names(c.models) <- c("Lunar Phase, Flow, Photoperiod", 
                     "Temperature residuals, Flow, Photoperiod")

mod.list_h_d_three <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,"p"),
             15,45, byrow = TRUE),
  D = "diagonal and unequal",
  d = covariates_all_years[91:105,] 
)




out.tab3 <- NULL
fits3 <- list()
for(i in 1:length(c.models)){
  
  fit.model = c(list(c= c.models[[i]]), mod.list_h_d_three)
  
  fit <- MARSS(subset_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d="Hatchery",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab3=rbind(out.tab3,out)
  fits3=c(fits3,list(fit))
  
  
}



min.AICc <- order(out.tab3$AICc)
out.tab.3 <- out.tab3[min.AICc, ]
out.tab.3$DeltaAICc <- out.tab.3$AICc - min(out.tab.3$AICc)
out.tab.3


#now that I know that only photoperiod and flow are needed, I will use 2015


subset_perhour_all <- arrange(data,doy) %>%
  filter(doy > 120 & doy <= 160) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

#scale the values

for(i in 1:dim(subset_perhour_all)[1]){
  subset_perhour_all[i,] = scale(subset_perhour_all[i,])[,1]
}

covariates_all_years_final <- arrange(data,doy) %>%
  filter(doy >120 & doy <= 160) %>%
  dplyr::select(year,doy, photoperiod, flow, coho1_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = year, values_from = c(
    photoperiod, flow, coho1_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


for(i in 1:(dim(covariates_all_years_final)[1]-16)){
covariates_all_years_final[i,] = scale(covariates_all_years_final[i,])[,1]
}



for(i in 33:48){
  covariates_all_years_final[i,] = scale(covariates_all_years_final[i,], center = FALSE, scale= TRUE)[,1]
}


mod.list_h_d_final <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"flow"),
             16,32, byrow = TRUE),
  c = covariates_all_years_final[1:32,],
  D = "diagonal and unequal",
  d = covariates_all_years_final[33:48,] 
)

final <- MARSS(subset_perhour_all, model=mod.list_h_d_final, silent = TRUE, method = "BFGS", 
               control=list(maxit=2000))

ci_final <- tidy(final)
ci_final

ggplot(ci_final[19:34,], aes(x = 2005:2020, y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Estimate of Hatchery effect") +
  ggtitle("Effect of hatchery Coho") + 
  theme(plot.title = element_text(size = 20))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

ggsave(here("output","coho_hatchery_effect.jpeg"), width = 10, height = 8)

autoplot(final)


autoplot(final, plot.type = "fitted.ytT")
ggsave(here("output","fitted_y_coho.jpeg"), width = 16, height = 10)

ggplot(ci_final[35:36,], aes(x = c("Photoperiod", "Flow"), y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of environmental effect") +
  ggtitle("Environmental effects") + 
  theme(plot.title = element_text(size = 20))+
  theme(axis.text.x=element_text(size=14),axis.title.y=element_text(size=14))

ggsave(here("output","coho_environmental_effect.jpeg"), width = 7500, height = 7000, units = "px")


#adding flow to observation covariate

