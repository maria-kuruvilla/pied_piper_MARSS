# goal to try MARSS for all years of chinook0
#first try for both spring and summer runs together - we do not have most covariates for spring
#if that does not work, try separately

#packages
library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)


#readng the data

data <- read.csv(here("data","dungeness_subset_covariates.csv"), header = TRUE)

data$chinook0_hatchery_num_interpolate <- na.approx(data$chinook0_hatchery_num)

plot(data$doy, data$chinook0_wild_num)

plot(data$doy[data$doy >25 & data$doy <100], data$chinook0_wild_num[data$doy >25 & data$doy <100])
plot(data$doy[data$doy >120 & data$doy <225], data$chinook0_wild_num[data$doy >120 & data$doy <225])

#subset
subset_chinook <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 25 & doy <= 225) %>%
  mutate(log.value = ifelse(chinook0_wild_num == 0, NA, log(chinook0_wild_num))) %>%
  select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scale the values

for(i in 1:dim(subset_chinook)[1]){
  subset_chinook[i,] = scale(subset_chinook[i,])[,1]
}



covariates_all_years <- arrange(data,doy) %>%
  filter(year != 2015 & doy >25 & doy <= 225) %>%
  select(year,doy, temp, atu_april, photoperiod, flow) %>%
  pivot_wider(names_from = year, values_from = c(temp, atu_april, photoperiod, flow)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


mod.list <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = "diagonal and unequal"
)

mod.list0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal"
)


fit_chinook0 <- MARSS(subset_chinook, model = mod.list, method = "BFGS")

autoplot(fit_chinook0)

fit_chinook0$AICc #4591





#summer

#subset
subset_chinook_summer <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 120 & doy <= 209) %>%
  mutate(log.value = ifelse(chinook0_wild_num == 0, NA, log(chinook0_wild_num))) %>%
  select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


covariates_all_years_summer <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 120 & doy <= 209) %>%
  select(year,doy, temp, atu_april, photoperiod, flow) %>%
  pivot_wider(names_from = year, values_from = c(temp, atu_april, photoperiod, flow)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()



#putting it into a loop

c.models <- list( 
  c1 = "zero", # No covariates
  c2 = zscore(covariates_all_years_summer[1:15,]), # Temperature
  c3 = zscore(covariates_all_years_summer[16:30,]), # ATU April
  c4 = zscore(covariates_all_years_summer[31:45,]),  #Photoperiod
  c5 = zscore(covariates_all_years_summer[46:60,])  #Flow
)
names(c.models) <- c("No covariates", "Temperature", "ATU", "Photoperiod", "Flow")


out.tab <- NULL
fits <- list()
for(i in 1:length(c.models)){
  if(i == 1){
    fit.model = c(list(c= c.models[[i]]), mod.list0)
  }
  else{
    fit.model = c(list(c= c.models[[i]]), mod.list)
  }
  
  fit <- MARSS(subset_chinook_summer, model=fit.model, silent = TRUE,
                   control=list(maxit=2000, allow.degen=TRUE))

      
      out=data.frame(c=names(c.models)[i], 
                     logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                     num.iter=fit$numIter, converged=!fit$convergence,
                     stringsAsFactors = FALSE)
      out.tab=rbind(out.tab,out)
      fits=c(fits,list(fit))
      
}

 

min.AICc <- order(out.tab$AICc)
out.tab.1 <- out.tab[min.AICc, ]
out.tab.1


#model with flow is best

#now trying flow with temp, atu and photoperiod
c.models2 <- list( 
  c1 = zscore(covariates_all_years_summer[c(1:15,46:60),]), # Temperature and flow
  c2 = zscore(covariates_all_years_summer[c(16:30,46:60),]), # ATU April and flow
  c3 = zscore(covariates_all_years_summer[c(31:60),])  #Photoperiod and flow
  
)
names(c.models2) <- c("Temperature and Flow", "ATU and Flow", "Photoperiod and Flow")

mod.list2 <- mod.list.photo_flow_all_years <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = matrix(list("p2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"p2006",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"p2007",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"p2008",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"p2009",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"p2010",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"p2011",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"p2012",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"p2013",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"p2014",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"p2016",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"p2017",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"p2018",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"p2019",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p2020",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2020"),
             15,30, byrow = TRUE)
)


out.tab2 <- NULL
fits2 <- list()
for(i in 1:length(c.models2)){
  
  fit.model2 = c(list(c= c.models2[[i]]), mod.list2)
  
  
  fit2 <- MARSS(subset_chinook_summer, model=fit.model2, silent = TRUE,
               control=list(maxit=2000, allow.degen=TRUE))
  
  
  out=data.frame(c=names(c.models2)[i], 
                 logLik=fit2$logLik, AICc=fit2$AICc, num.param=fit2$num.params,
                 num.iter=fit2$numIter, converged=!fit2$convergence,
                 stringsAsFactors = FALSE)
  out.tab2=rbind(out.tab2,out)
  fits2=c(fits2,list(fit2))
  
}



min.AICc <- order(out.tab2$AICc)
out.tab.2 <- out.tab2[min.AICc, ]
out.tab.2



#trying flow and hatchery

covariates_hatchery_all_years_summer <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 120 & doy <= 209) %>%
  select(year,doy, chinook0_hatchery_num_interpolate) %>%
  pivot_wider(names_from = year, values_from = chinook0_hatchery_num_interpolate) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

c.models3 <- list( 
  
  c1 = rbind(covariates_hatchery_all_years_summer,zscore(covariates_all_years_summer[c(46:60),]))  #hatchery and flow
  
)
names(c.models3) <- c("Hatchery and Flow")

mod.list2 <- mod.list.photo_flow_all_years <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = matrix(list("p2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"p2006",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"p2007",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"p2008",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"p2009",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"p2010",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"p2011",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"p2012",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"p2013",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"p2014",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"p2016",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"p2017",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"p2018",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"p2019",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p2020",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2020"),
             15,30, byrow = TRUE),
  c = c.models3[[1]]
)


fit_chinook0_hatchery <- MARSS(subset_chinook_summer, model = mod.list2, method = "BFGS")
