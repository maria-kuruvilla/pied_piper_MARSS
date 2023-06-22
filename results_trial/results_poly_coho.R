#trying poly of flow, temperature, photoperiod

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

lm_temp <- lm(temp ~ 1+photoperiod, data)


data <- data %>% add_residuals(lm_temp)


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

#making poly
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

mod.list_base <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  D = "diagonal and unequal",
  d = covariates_all_years_new[121:135,] 
)

mod.list_h_d <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = "diagonal and equal",
  D = "diagonal and unequal",
  d = covariates_all_years_new[121:135,] 
)


#putting it into a loop

c.models <- list( 
  c1 = "zero", # No covariates
  c2 = covariates_all_years_new[1:15,], # Temperature
  c3 = covariates_all_years_new[31:45,], # photoperiod
  c4 = covariates_all_years_new[61:75,]  #flow
  
)
names(c.models) <- c("No covariates", 
                     "Temperature", 
                     "Photoperiod", 
                     "Flow")


out.tab_linear <- NULL
fits_linear <- list()
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
  out.tab_linear=rbind(out.tab_linear,out)
  fits_linear=c(fits_linear,list(fit))
  
  
}



min.AICc <- order(out.tab_linear$AICc)
out.tab.linear <- out.tab_linear[min.AICc, ]
out.tab.linear$DeltaAICc <- out.tab.linear$AICc - min(out.tab.linear$AICc)
out.tab.linear



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
  d = covariates_all_years_new[121:135,] 
)

#putting it into a loop

c.models <- list(
  c2 = covariates_all_years_new[1:30,], # Temperature
  c3 = covariates_all_years_new[31:60,], # photoperiod
  c4 = covariates_all_years_new[61:90,]  #flow
  
)
names(c.models) <- c("Temperature square", 
                     "Photoperiod square", 
                     "Flow square")


out.tab_poly <- NULL
fits_poly <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(names(c.models)[i] == "No covariates"){
    
    fit.model = mod.list_base
  }
  
  else{
    fit.model = c(list(c= c.models[[i]]), mod.list_h_d_poly)
    
  }
  
  
  
  
  
  
  fit <- MARSS(subset_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d="Hatchery",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_poly=rbind(out.tab_poly,out)
  fits_poly=c(fits_poly,list(fit))
  
  
}



min.AICc <- order(out.tab_poly$AICc)
out.tab.poly <- out.tab_poly[min.AICc, ]
out.tab.poly$DeltaAICc <- out.tab.poly$AICc - min(out.tab.poly$AICc)
out.tab.poly

ci_poly_photo <- tidy(fits_poly[[2]])
ci_poly_photo

# I think poly is not much better than linear