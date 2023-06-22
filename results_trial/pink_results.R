# repeat same MARSS analysis for pink


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

data$pink0_hatchery_perhour_interpolate <- na.approx(data$pink0_hatchery_perhour)

plot(data$doy[data$doy > 60 & data$doy < 140],data$pink0_wild_perhour[data$doy > 60 & data$doy < 140])

year = 2007
plot(data$doy[data$doy > 60 & data$doy < 140 & data$year == year],
     data$pink0_hatchery_perhour_interpolate[data$doy > 60 & data$doy < 140 & data$year == year])
#only 2008 and 2010 have hatchery

plot(data$doy[data$doy > 60 & data$doy < 140 & data$year == year],
     data$pink0_wild_perhour[data$doy > 60 & data$doy < 140 & data$year == year])


lm_temp <- lm(temp ~ 1+photoperiod, data)


data <- data %>% add_residuals(lm_temp)

#need to leave out 2005, 2015, 2016 and atu

subset_perhour_pink <- arrange(data,doy) %>%
  filter(year != 2015 & year != 2005 & year != 2016 & doy > 60 & doy <= 140) %>%
  mutate(log.value = log(pink0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

#scale the values

for(i in 1:dim(subset_perhour_pink)[1]){
  subset_perhour_pink[i,] = scale(subset_perhour_pink[i,])[,1]
}


covariates_all_years_pink <- arrange(data,doy) %>%
  filter(year != 2015 & year != 2005 & year != 2016 & doy > 60 & doy <= 140) %>%
  dplyr::select(year,doy, temp, photoperiod, lunar_phase, resid, flow, pink0_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = year, values_from = c(
    temp, photoperiod, lunar_phase, resid, flow, pink0_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


for(i in 1:(dim(covariates_all_years_pink)[1]-39)){
  covariates_all_years_pink[i,] = scale(covariates_all_years_pink[i,])[,1]
}


for(i in 53:65){
  covariates_all_years_pink[i,] = scale(covariates_all_years_pink[i,])[,1]
}

for(i in c(68,70)){
  covariates_all_years_pink[i,] = scale(covariates_all_years_pink[i,], center = FALSE, scale= TRUE)[,1]
}


make_poly_columns <- function(covariates_all_years){
  
  covariates_all_years_new <- matrix(NA,117,80)
  
  for(i in 1:13){
    covariates_all_years_new[i,] <- t(poly(covariates_all_years[i,],2)[,1])
    
    covariates_all_years_new[13+i,] <- t(poly(covariates_all_years[i,],2)[,2])
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
  
  
  for(i in 14:26){
    covariates_all_years_new[i+13,] <- t(poly(covariates_all_years[i,],2)[,1])
    
    covariates_all_years_new[26+i,] <- t(poly(covariates_all_years[i,],2)[,2])
    print(rownames(covariates_all_years)[i])
    
    if(i==14){
      list3 <- rownames(covariates_all_years)[i]
      list4 <- paste0(rownames(covariates_all_years)[i], "_square")
    }
    else{
      list3 <- c(list3, rownames(covariates_all_years)[i])
      list4 <- c(list4, paste0(rownames(covariates_all_years)[i], "_square"))
    }
  }
  
  for(i in 53:65){
    covariates_all_years_new[i,] <- t(poly(covariates_all_years[i,],2)[,1])
    
    covariates_all_years_new[i+13,] <- t(poly(covariates_all_years[i,],2)[,2])
    print(rownames(covariates_all_years)[i])
    
    if(i==53){
      list5 <- rownames(covariates_all_years)[i]
      list6 <- paste0(rownames(covariates_all_years)[i], "_square")
    }
    else{
      list5 <- c(list5, rownames(covariates_all_years)[i])
      list6 <- c(list6, paste0(rownames(covariates_all_years)[i], "_square"))
    }
  }
  
  
  
  colnames(covariates_all_years_new) <- colnames(covariates_all_years)
  
  
  
  covariates_all_years_new[79:91,]<- covariates_all_years[27:39,]
  covariates_all_years_new[92:104,]<- covariates_all_years[40:52,]
  covariates_all_years_new[105:117,]<- covariates_all_years[66:78,]
  

  
  
  
  
  rownames(covariates_all_years_new) <- c(list1,list2, list3, list4, list5, list6,
                                          rownames(covariates_all_years)[27:39],
                                          rownames(covariates_all_years)[40:52],
                                          rownames(covariates_all_years)[66:78])
  
  return(covariates_all_years_new)
  
}


covariates_all_years_pink_new <- make_poly_columns(covariates_all_years_pink)




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
  d = covariates_all_years_pink_new[105:117,] 
)

mod.list_1_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = "diagonal and equal",
  D = "diagonal and unequal",
  d = covariates_all_years_pink_new[105:117,] 
)


# photo square, temp square, flow square in d


mod.list_2_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,"p"),
             13,26, byrow = TRUE),
  D = "diagonal and unequal",
  d = covariates_all_years_pink_new[105:117,] 
)


c.models <- list(
  c1 = "zero",
  c2 = "zero",
  c3 = covariates_all_years_pink_new[1:13,], # Temperature 
  c4 = covariates_all_years_pink_new[27:39,], # photoperiod 
  c6 = covariates_all_years_pink_new[1:26,], # Temperature sqaure
  c7 = covariates_all_years_pink_new[27:52,], # photoperiod sqaure
  c8 = covariates_all_years_pink_new[53:78,]  #flow sqaure 
  
)
names(c.models) <- c("No covariates", 
                     "No covariates in c, hatchery in d", 
                     "Temperature, hatchery in d", 
                     "Photoperiod, hatchery in d", 
                     "Temperature sqaure, hatchery in d", 
                     "Photoperiod square, hatchery in d", 
                     "Flow square, hatchery in d")

c_num = c(0,0,1,1,2,2,2)
d_num = c(0,1,1,1,1,1,1)

out.tab_pink1 <- NULL
fits_pink1 <- list()
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
  
  
  
  
  
  
  fit <- MARSS(subset_perhour_pink, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i],
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_pink1=rbind(out.tab_pink1,out)
  fits_pink1=c(fits_pink1,list(fit))
  
  
}




min.AICc <- order(out.tab_pink1$AICc)
out.tab.pink1 <- out.tab_pink1[min.AICc, ]
out.tab.pink1



# then add lunar phase, resid, flow in c (along with photoperiod) with hatchery in d



c.models <- list(
  
  c9 = covariates_all_years_pink_new[c(27:39,79:91),], # photo and lunar phase
  c10 = covariates_all_years_pink_new[c(27:39,92:104),],  #photo sqaure and resid
  c11 = covariates_all_years_pink_new[c(27:39,53:65),]  # photo square and flow
  
)

names(c.models) <- c(
  "Photoperiod, lunar phase, hatchery in d",
  "Photoperiod, residuals, hatchery in d",
  "Photoperiod, flow, hatchery in d"
)



out.tab_pink2 <- NULL
fits_pink2 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  
  fit.model = c(list(c= c.models[[i]]), mod.list_2_1)
  fit <- MARSS(subset_perhour_pink, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i],
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_pink2=rbind(out.tab_pink2,out)
  fits_pink2=c(fits_pink2,list(fit))
  
  
}




min.AICc <- order(out.tab_pink2$AICc)
out.tab.pink2 <- out.tab_pink2[min.AICc, ]
out.tab.pink2

autoplot(fits_pink2[[1]])

tidy(fits_pink2[[1]])
