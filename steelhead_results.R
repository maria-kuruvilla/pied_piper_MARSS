packages
library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)
library(broom)
library(ggfortify)
library(modelr)
library(MASS)

library(parallel)

# data wrangling
######### 
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
  covariates_all_years_steelhead[i,] = scale(covariates_all_years_steelhead[i,])[,1]
}


for(i in 76:90){
  covariates_all_years_steelhead[i,] = scale(covariates_all_years_steelhead[i,])[,1]
}

for(i in 91:102){
  covariates_all_years_steelhead[i,] = scale(covariates_all_years_steelhead[i,], center = FALSE, scale= TRUE)[,1]
}

#103 - 2018 does not have any hatchery

for(i in 104:105){
  covariates_all_years_steelhead[i,] = scale(covariates_all_years_steelhead[i,], center = FALSE, scale= TRUE)[,1]
}

make_poly_columns <- function(covariates_all_years){
  
  covariates_all_years_new <- matrix(NA,150,40)
  
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
  covariates_all_years_new[136:150,]<- covariates_all_years[46:60,]
  covariates_all_years_new[106:120,]<- covariates_all_years[61:75,]
  covariates_all_years_new[121:135,]<- covariates_all_years[91:105,]
  
  
  rownames(covariates_all_years_new) <- c(list1,list2, list3, list4, list5, list6,
                                          rownames(covariates_all_years)[16:30],
                                          rownames(covariates_all_years)[61:75],
                                          rownames(covariates_all_years)[91:105],
                                          rownames(covariates_all_years)[46:60])
  
  return(covariates_all_years_new)
  
}


covariates_all_years_steelhead_new <- make_poly_columns(covariates_all_years_steelhead)

#######

# make one for loop to try all models

# different model structures



# no covariates in c and d

mod.list_0_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal"
)

# no covariates in c but hatchery in d

mod.list_0_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  D = "diagonal and unequal",
  d = covariates_all_years_steelhead_new[121:135,] 
)


# photoperiod, atu, temperature in c and hatchery in d


mod.list_1_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = "diagonal and equal",
  D = "diagonal and unequal",
  d = covariates_all_years_steelhead_new[121:135,] 
)


# photo square, temp square, flow square in d


mod.list_2_1 <- list(
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


c.models <- list(
  c1 = "zero",
  c2 = "zero",
  c3 = covariates_all_years_steelhead_new[1:15,], # Temperature 
  c4 = covariates_all_years_steelhead_new[31:45,], # photoperiod 
  c5 = covariates_all_years_steelhead_new[91:105,],  #atu
  c6 = covariates_all_years_steelhead_new[1:30,], # Temperature sqaure
  c7 = covariates_all_years_steelhead_new[31:60,], # photoperiod sqaure
  c8 = covariates_all_years_steelhead_new[61:90,]  #flow sqaure 
  
)
names(c.models) <- c("No covariates", 
                     "No covariates in c, hatchery in d", 
                     "Temperature, hatchery in d", 
                     "Photoperiod, hatchery in d", 
                     "ATU, hatchery in d", 
                     "Temperature sqaure, hatchery in d", 
                     "Photoperiod square, hatchery in d", 
                     "Flow square, hatchery in d")

c_num = c(0,0,1,1,1,2,2,2)
d_num = c(0,1,1,1,1,1,1,1)

out.tab_steelhead1 <- NULL
fits_steelhead1 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 0 & d_num[i] == 0){
    
    fit.model = mod.list_0_0
  }
  
  else if(c_num[i] == 0 & d_num[i] == 1){
    
    fit.model = mod.list_0_1
  }
  
  else if(c_num[i] == 1 & d_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_1_1)
  }
  
  else if(c_num[i] == 2 & d_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_2_1)
  }
  
  
  
  
  
  
  fit <- MARSS(subset_perhour_steelhead, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i],
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_steelhead1=rbind(out.tab_steelhead1,out)
  fits_steelhead1=c(fits_steelhead1,list(fit))
  
  
}




min.AICc <- order(out.tab_steelhead1$AICc)
out.tab.steelhead1 <- out.tab_steelhead1[min.AICc, ]
out.tab.steelhead1


# then add lunar phase, resid, flow in c with hatchery in d

mod.list_3_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,"c"),
             15,45, byrow = TRUE),
  D = "diagonal and unequal",
  d = covariates_all_years_steelhead_new[121:135,] 
)


c.models <- list(
  
  c9 = covariates_all_years_steelhead_new[c(31:60,136:150),], # photo sqaure and lunar phase
  c10 = covariates_all_years_steelhead_new[c(31:60,106:120),],  #photo sqaure and resid
  c11 = covariates_all_years_steelhead_new[31:75,]  # photo square and flow
  
)

names(c.models) <- c(
                     "Photoperiod square, lunar phase, hatchery in d",
                     "Photoperiod square, residuals, hatchery in d",
                     "Photoperiod square, flow, hatchery in d"
                     )

out.tab_steelhead2 <- NULL
fits_steelhead2 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  
  fit.model = c(list(c= c.models[[i]]), mod.list_3_1)
  fit <- MARSS(subset_perhour_steelhead, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i],
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_steelhead2=rbind(out.tab_steelhead2,out)
  fits_steelhead2=c(fits_steelhead2,list(fit))
  
  
}




min.AICc <- order(out.tab_steelhead2$AICc)
out.tab.steelhead2 <- out.tab_steelhead2[min.AICc, ]
out.tab.steelhead2



# then add flow, resid and hatchery in d


mod.list_2_2 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  c = covariates_all_years_steelhead_new[31:60,],
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b"),
             15,30, byrow = TRUE),
  D = matrix(list("c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d"),
             15,30, byrow = TRUE)
  
)

d.models <- list(
  
  d1 = covariates_all_years_steelhead_new[c(121:135,61:75),], # hatchery and flow
  d2 = covariates_all_years_steelhead_new[106:135,]  #resid and hatchery
  
)

names(d.models) <- c(
  "Hatchery and flow",
  "Residuals and hatchery"
)


out.tab_steelhead3 <- NULL
fits_steelhead3 <- list()
for(i in 1:length(d.models)){
  
  print(names(d.models)[i])
  
  fit.model = c(list(d= d.models[[i]]), mod.list_2_2)
  fit <- MARSS(subset_perhour_steelhead, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(d.models)[i],
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_steelhead3=rbind(out.tab_steelhead3,out)
  fits_steelhead3=c(fits_steelhead3,list(fit))
  
  
}




min.AICc <- order(out.tab_steelhead3$AICc)
out.tab.steelhead3 <- out.tab_steelhead3[min.AICc, ]
out.tab.steelhead3


# final model has hatchery in d and photoperiod sqaure

ci_final_steelhead <- tidy(fits_steelhead1[[7]])

# Detect the number of available cores and create cluster
# cl <- parallel::makeCluster(detectCores())
# # Run parallel computation
# time_parallel <- system.time(
#   parallel::parLapply(cl,
#                       data_list,
#                       mean)
# )
# # Close cluster
# parallel::stopCluster(cl)

ci_steelhead_final <- MARSSparamCIs(fits_steelhead1[[7]], method = "parametric", alpha = 0.05, 
                                    nboot = 500, silent = FALSE) #, hessian.fun = "optim")


ci_steelhead_final_tab <- tidy(ci_steelhead_final)


ggplot(ci_steelhead_final_tab[18:32,], aes(x = c(2005:2014, 2016:2020), y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Estimate of Hatchery effect") +
  ggtitle("Effect of hatchery Steelhead smolts") + 
  theme(plot.title = element_text(size = 20))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

ggsave(here("output","steelheadsmolts_hatchery_effect.jpeg"), width = 10, height = 8)


ggplot(ci_steelhead_final_tab[33:34,], aes(x = c("Photoperiod", "Photoperiod squared"), y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of environmental effect") +
  ggtitle("Environmental effects on Steelhead smolts") + 
  theme(plot.title = element_text(size = 20))+
  theme(axis.text.x=element_text(size=14),axis.title.y=element_text(size=14))

ggsave(here("output","steelheadsmolts_environmental_effect.jpeg"), width = 10, height = 8)

autoplot(fits_steelhead1[[7]])


autoplot(fits_steelhead1[[7]], plot.type = "fitted.ytT")
ggsave(here("output","fitted_y_steelhead.jpeg"), width = 10, height = 8)


#since it is only photoperiod, I am adding 2015 back in



data <- read.csv(here("data","dungeness_subset_covariates.csv"), header = TRUE)

data$steelheadsmolt_hatchery_perhour_interpolate <- na.approx(data$steelheadsmolt_hatchery_perhour)

lm_temp <- lm(temp ~ 1+photoperiod, data)


data <- data %>% add_residuals(lm_temp)

plot(data$doy[data$doy >120 & data$doy < 160], 
     data$steelheadsmolt_wild_perhour[data$doy >120 & data$doy < 160])

plot(data$doy[data$doy >120 & data$doy < 160], 
     data$steelheadsmolt_hatchery_perhour[data$doy >120 & data$doy < 160])


subset_perhour_steelhead <- arrange(data,doy) %>%
  filter(doy > 120 & doy <= 160) %>%
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
  filter(doy >120 & doy <= 160) %>%
  dplyr::select(year,doy, temp, atu_april, photoperiod, lunar_phase, resid, flow, steelheadsmolt_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = year, values_from = c(
    temp, atu_april, photoperiod, lunar_phase, resid, flow, steelheadsmolt_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


for(i in 1:(dim(covariates_all_years_steelhead)[1]-48)){
  covariates_all_years[i,] = scale(covariates_all_years[i,])[,1]
}


for(i in 81:96){
  covariates_all_years_steelhead[i,] = scale(covariates_all_years_steelhead[i,])[,1]
}

for(i in 97:109){
  covariates_all_years_steelhead[i,] = scale(covariates_all_years_steelhead[i,], center = FALSE, scale= TRUE)[,1]
}


for(i in 111:112){
  covariates_all_years_steelhead[i,] = scale(covariates_all_years_steelhead[i,], center = FALSE, scale= TRUE)[,1]
}

#change numbers in this
make_poly_columns <- function(covariates_all_years){
  
  covariates_all_years_new <- matrix(NA,160,40)
  
  for(i in 1:16){
    if(i==11){
      covariates_all_years_new[i,] <- covariates_all_years[11,]
      covariates_all_years_new[16+i,] <- covariates_all_years[11,]
      list1 <- c(list1, rownames(covariates_all_years)[i])
      list2 <- c(list2, paste0(rownames(covariates_all_years)[i], "_square"))
    }
    else{
      covariates_all_years_new[i,] <- t(poly(covariates_all_years[i,],2)[,1])
      
      covariates_all_years_new[16+i,] <- t(poly(covariates_all_years[i,],2)[,2])
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
    }
    
    
    
  
  
  for(i in 33:48){
    covariates_all_years_new[i,] <- t(poly(covariates_all_years[i,],2)[,1])
    
    covariates_all_years_new[16+i,] <- t(poly(covariates_all_years[i,],2)[,2])
    print(rownames(covariates_all_years)[i])
    
    if(i==33){
      list3 <- rownames(covariates_all_years)[i]
      list4 <- paste0(rownames(covariates_all_years)[i], "_square")
    }
    else{
      list3 <- c(list3, rownames(covariates_all_years)[i])
      list4 <- c(list4, paste0(rownames(covariates_all_years)[i], "_square"))
    }
  }
  
  for(i in 81:96){
    covariates_all_years_new[i-15,] <- t(poly(covariates_all_years[i,],2)[,1])
    
    covariates_all_years_new[i,] <- t(poly(covariates_all_years[i,],2)[,2])
    print(rownames(covariates_all_years)[i])
    
    if(i==81){
      list5 <- rownames(covariates_all_years)[i]
      list6 <- paste0(rownames(covariates_all_years)[i], "_square")
    }
    else{
      list5 <- c(list5, rownames(covariates_all_years)[i])
      list6 <- c(list6, paste0(rownames(covariates_all_years)[i], "_square"))
    }
  }
  
  
  
  colnames(covariates_all_years_new) <- colnames(covariates_all_years)
  
  
  
  covariates_all_years_new[97:112,]<- covariates_all_years[17:32,]
  covariates_all_years_new[145:160,]<- covariates_all_years[49:64,]
  covariates_all_years_new[113:128,]<- covariates_all_years[65:80,]
  covariates_all_years_new[129:144,]<- covariates_all_years[97:112,]
  
  
  rownames(covariates_all_years_new) <- c(list1,list2, list3, list4, list5, list6,
                                          rownames(covariates_all_years)[17:32],
                                          rownames(covariates_all_years)[65:80],
                                          rownames(covariates_all_years)[97:112],
                                          rownames(covariates_all_years)[49:64])
  
  
  return(covariates_all_years_new)
  
}


covariates_all_years_steelhead_new <- make_poly_columns(covariates_all_years_steelhead)


#final steelhead


mod.list_2_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"photoperiod^2"),
             16,32, byrow = TRUE),
  c = covariates_all_years_steelhead_new[33:64,],
  D = "diagonal and unequal",
  d = covariates_all_years_steelhead_new[129:144,] 
)

fit_final_steelhead <- MARSS(subset_perhour_steelhead, model=mod.list_2_1, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
autoplot(fit_final_steelhead)

autoplot(fit_final_steelhead, plot.type = "fitted.ytT")
ggsave(here("output","fitted_y_steelhead.jpeg"), width = 10, height = 8)


ci_steelhead_final <- MARSSparamCIs(fit_final_steelhead, method = "parametric", alpha = 0.05, 
                                    nboot = 500, silent = FALSE) #, hessian.fun = "optim")


ci_steelhead_final_tab <- tidy(ci_steelhead_final)
