

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

data$steelheadsmolt_hatchery_perhour_interpolate <- na.approx(data$steelheadsmolt_hatchery_perhour)

lm_temp <- lm(temp ~ 1+photoperiod, data)


data <- data %>% add_residuals(lm_temp)

plot(data$doy[data$doy >120 & data$doy < 160], 
     data$steelheadsmolt_wild_perhour[data$doy >120 & data$doy < 160])

plot(data$doy[data$doy >120 & data$doy < 160], 
     data$steelheadsmolt_hatchery_perhour[data$doy >120 & data$doy < 160])


subset_perhour_steelhead <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 120 & doy <= 160) %>%
  mutate(log.value = log(steelheadsmolt_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

#scale the values

for(i in 1:dim(subset_perhour_steelhead)[1]){
  subset_perhour_steelhead[i,] = scale(subset_perhour_steelhead[i,])[,1]
}



covariates_all_years_steelhead <- arrange(data,doy) %>%
  filter(year != 2015 & doy >120 & doy <= 160) %>%
  dplyr::select(year,doy, temp, atu_april, photoperiod, lunar_phase, resid, flow, steelheadsmolt_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = year, values_from = c(
    temp, atu_april, photoperiod, lunar_phase, resid, flow, steelheadsmolt_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


for(i in 1:(dim(covariates_all_years_steelhead)[1]-45)){
  covariates_all_years[i,] = scale(covariates_all_years[i,])[,1]
}


for(i in 76:90){
  covariates_all_years_steelhead[i,] = scale(covariates_all_years_steelhead[i,])[,1]
}

for(i in 91:102){
  covariates_all_years_steelhead[i,] = scale(covariates_all_years_steelhead[i,], center = FALSE, scale= TRUE)[,1]
}


for(i in 104:105){
  covariates_all_years_steelhead[i,] = scale(covariates_all_years_steelhead[i,], center = FALSE, scale= TRUE)[,1]
}

make_poly_columns <- function(covariates_all_years){
  
  covariates_all_years_new <- matrix(NA,135,40)
  
  for(i in 1:15){
    covariates_all_years_new[i,] <- t(poly(covariates_all_years[i,],2)[,1])
    
    covariates_all_years_new[15+i,] <- t(poly(covariates_all_years[i,],2)[,2])
    print(rownames(covariates_all_years)[i])
    
    if(i==1){
      list1 <- rownames(covariates_all_years)[i]
      list2 <- paste0(rownames(covariates_all_years)[i], "_square")
    }
    else{
      list1 <- c(list1, rownames(covariates_all_years)[i])
      list2 <- c(list2, paste0(rownames(covariates_all_years)[i], "_square"))
    }
  }
  
  
  for(i in 31:45){
    covariates_all_years_new[i,] <- t(poly(covariates_all_years[i,],2)[,1])
    
    covariates_all_years_new[15+i,] <- t(poly(covariates_all_years[i,],2)[,2])
    print(rownames(covariates_all_years)[i])
    
    if(i==31){
      list3 <- rownames(covariates_all_years)[i]
      list4 <- paste0(rownames(covariates_all_years)[i], "_square")
    }
    else{
      list3 <- c(list3, rownames(covariates_all_years)[i])
      list4 <- c(list4, paste0(rownames(covariates_all_years)[i], "_square"))
    }
  }
  
  for(i in 76:90){
    covariates_all_years_new[i-15,] <- t(poly(covariates_all_years[i,],2)[,1])
    
    covariates_all_years_new[i,] <- t(poly(covariates_all_years[i,],2)[,2])
    print(rownames(covariates_all_years)[i])
    
    if(i==76){
      list5 <- rownames(covariates_all_years)[i]
      list6 <- paste0(rownames(covariates_all_years)[i], "_square")
    }
    else{
      list5 <- c(list5, rownames(covariates_all_years)[i])
      list6 <- c(list6, paste0(rownames(covariates_all_years)[i], "_square"))
    }
  }
  
  
  
  colnames(covariates_all_years_new) <- colnames(covariates_all_years)
  
  
  
  covariates_all_years_new[91:105,]<- covariates_all_years[16:30,]
  covariates_all_years_new[106:120,]<- covariates_all_years[61:75,]
  covariates_all_years_new[121:135,]<- covariates_all_years[91:105,]
  
  
  rownames(covariates_all_years_new) <- c(list1,list2, list3, list4, list5, list6,
                                          rownames(covariates_all_years)[16:30],
                                          rownames(covariates_all_years)[61:75],
                                          rownames(covariates_all_years)[91:105])
  
  return(covariates_all_years_new)
  
}


covariates_all_years_steelhead_new <- make_poly_columns(covariates_all_years_steelhead)

mod.list <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal"
)

mod.list_base <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  D = "diagonal and unequal",
  d = covariates_all_years_steelhead_new[121:135,] 
)

mod.list_h_d <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = "diagonal and equal",
  D = "diagonal and unequal",
  d = covariates_all_years_steelhead_new[121:135,] 
)

fit_steelhead <- MARSS(subset_perhour_steelhead, model = mod.list, method = "BFGS")
fit_steelhead$AICc #1509
autoplot(fit_steelhead)


fit_steelhead_h_d <- MARSS(subset_perhour_steelhead, model = mod.list_base, method = "BFGS")
fit_steelhead_h_d$AICc #944
autoplot(fit_steelhead_h_d)




c.models <- list( 
  c1 = "zero", # No covariates
  c2 = covariates_all_years_steelhead_new[1:15,], # Temperature
  c3 = covariates_all_years_steelhead_new[31:45,], # photoperiod
  c4 = covariates_all_years_steelhead_new[61:75,]  #flow
  
)
names(c.models) <- c("No covariates", 
                     "Temperature", 
                     "Photoperiod", 
                     "Flow")


out.tab_linear_steelhead <- NULL
fits_linear_steelhead <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(names(c.models)[i] == "No covariates"){
    
    fit.model = mod.list_base
  }
  
  else{
    fit.model = c(list(c= c.models[[i]]), mod.list_h_d)
    
  }
  
  
  
  
  
  
  fit <- MARSS(subset_perhour_steelhead, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d="Hatchery",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_linear_steelhead=rbind(out.tab_linear_steelhead,out)
  fits_linear_steelhead=c(fits_linear_steelhead,list(fit))
  
  
}



min.AICc <- order(out.tab_linear_steelhead$AICc)
out.tab.linear_steelhead <- out.tab_linear_steelhead[min.AICc, ]
out.tab.linear_steelhead$DeltaAICc <- out.tab.linear_steelhead$AICc - min(out.tab.linear_steelhead$AICc)
out.tab.linear_steelhead


ci_steelhead <- tidy(fits_linear_steelhead[[3]])

ci_steelhead


ggplot(ci_steelhead[18:32,], aes(x = term, y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Estimate of Hatchery effect") +
  ggtitle("Effect of hatchery Steelhead") + 
  theme(plot.title = element_text(size = 20))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

autoplot(fits_linear_steelhead[[3]])


mod.list_h_d_poly <- list(
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
  d = covariates_all_years_steelhead_new[121:135,] 
)

#putting it into a loop

c.models <- list(
  c2 = covariates_all_years_steelhead_new[1:30,], # Temperature
  c3 = covariates_all_years_steelhead_new[31:60,], # photoperiod
  c4 = covariates_all_years_steelhead_new[61:90,]  #flow
  
)
names(c.models) <- c("Temperature square", 
                     "Photoperiod square", 
                     "Flow square")


out.tab_poly_steelhead <- NULL
fits_poly_steelhead <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(names(c.models)[i] == "No covariates"){
    
    fit.model = mod.list_base
  }
  
  else{
    fit.model = c(list(c= c.models[[i]]), mod.list_h_d_poly)
    
  }
  
  
  
  
  
  
  fit <- MARSS(subset_perhour_steelhead, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d="Hatchery",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_poly_steelhead=rbind(out.tab_poly_steelhead,out)
  fits_poly_steelhead=c(fits_poly_steelhead,list(fit))
  
  
}



min.AICc <- order(out.tab_poly_steelhead$AICc)
out.tab.poly_steelhead <- out.tab_poly_steelhead[min.AICc, ]
out.tab.poly_steelhead$DeltaAICc <- out.tab.poly_steelhead$AICc - min(out.tab.poly_steelhead$AICc)
out.tab.poly_steelhead


#need to try photoperiod, temp and atu

#photoperiod square is best 
#add lunar phase, flow, resid