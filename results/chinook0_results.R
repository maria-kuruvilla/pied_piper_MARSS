# goal 

# add new covariates - derivative of flow, temp, and photoperiod

# redo analysis with hatchery in c covariate


library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)
library(MASS)
library(modelr)
library(broom)
library(ggfortify)
library(qpcR)

#wrangling data

#########

data <- read.csv(here("data","dungeness_subset_covariates.csv"), header = TRUE)

data$chinook0_hatchery_perhour_interpolate <- na.approx(data$chinook0_hatchery_perhour)

lm_temp <- lm(temp ~ 1+photoperiod, data)

data <- data %>% add_residuals(lm_temp)

data$temp_diff <- c(NA, diff(data$temp))

data$flow_diff <- c(NA, diff(data$flow))

data$photo_diff <- c(NA, diff(data$photoperiod))

covariates_all_years_chinook0 <- arrange(data,doy) %>%
  filter(year != 2015 & doy >130 & doy <= 200) %>%
  dplyr::select(year,doy, temp, flow, photoperiod, atu_april, 
                lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
                chinook0_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = year, values_from = c(
    temp, flow, photoperiod, atu_april, lunar_phase, 
    resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


for(i in 1:(dim(covariates_all_years_chinook0)[1]-75)){ # everything except resid, diffs and hatchery
  covariates_all_years_chinook0[i,] = scale(covariates_all_years_chinook0[i,])[,1]
}


for(i in 76:135){
  covariates_all_years_chinook0[i,] = scale(covariates_all_years_chinook0[i,], center = FALSE, scale= TRUE)[,1]
}

for(i in 137:150){ # not starting at 136 beacuse 2005 does not have hatchery
  covariates_all_years_chinook0[i,] = scale(covariates_all_years_chinook0[i,], center = FALSE, scale= TRUE)[,1]
}

#subset
subset_chinook_summer_perhour <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 130 & doy <= 200) %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


for(i in 1:dim(subset_chinook_summer_perhour)[1]){
  subset_chinook_summer_perhour[i,] = scale(subset_chinook_summer_perhour[i,])[,1]
}



##########

#making poly

######

make_poly <- function(covariates_all_years_chinook0){
  
  covariates_all_years_new_chinook <- matrix(NA,195,70)
  
  #we want 1-15 to be temp. 16-30 to be temp^2
  # 31-45 to be flow. 46-60 to be flow^2
  # 61-75 to be photoperiod. 76-90 to be photoperiod^2
  
  
  for(i in 1:15){
    covariates_all_years_new_chinook[i,] <- t(poly(covariates_all_years_chinook0[i,],2)[,1])
    
    covariates_all_years_new_chinook[15+i,] <- t(poly(covariates_all_years_chinook0[i,],2)[,2])
    print(rownames(covariates_all_years_chinook0)[i])
    
    if(i==1){
      list1 <- rownames(covariates_all_years_chinook0)[i]
      list2 <- paste0(rownames(covariates_all_years_chinook0)[i], "_square")
    }
    else{
      list1 <- c(list1, rownames(covariates_all_years_chinook0)[i])
      list2 <- c(list2, paste0(rownames(covariates_all_years_chinook0)[i], "_square"))
    }
  }
  
  
  for(i in 16:30){
    covariates_all_years_new_chinook[i+15,] <- t(poly(covariates_all_years_chinook0[i,],2)[,1])
    
    covariates_all_years_new_chinook[30+i,] <- t(poly(covariates_all_years_chinook0[i,],2)[,2])
    print(rownames(covariates_all_years_chinook0)[i])
    
    if(i==16){
      list3 <- rownames(covariates_all_years_chinook0)[i]
      list4 <- paste0(rownames(covariates_all_years_chinook0)[i], "_square")
    }
    else{
      list3 <- c(list3, rownames(covariates_all_years_chinook0)[i])
      list4 <- c(list4, paste0(rownames(covariates_all_years_chinook0)[i], "_square"))
    }
  }
  
  
  for(i in 31:45){
    covariates_all_years_new_chinook[i+30,] <- t(poly(covariates_all_years_chinook0[i,],2)[,1])
    
    covariates_all_years_new_chinook[i+45,] <- t(poly(covariates_all_years_chinook0[i,],2)[,2])
    print(rownames(covariates_all_years_chinook0)[i])
    
    if(i==31){
      list5 <- rownames(covariates_all_years_chinook0)[i]
      list6 <- paste0(rownames(covariates_all_years_chinook0)[i], "_square")
    }
    else{
      list5 <- c(list5, rownames(covariates_all_years_chinook0)[i])
      list6 <- c(list6, paste0(rownames(covariates_all_years_chinook0)[i], "_square"))
    }
  }
  
  
  
  covariates_all_years_new_chinook[91:195,]<- covariates_all_years_chinook0[46:150,]
  
  
  rownames(covariates_all_years_new_chinook) <- c(list1,list2, list3, list4, list5, list6,
                                                  rownames(covariates_all_years_chinook0)[46:150])
  
  
  colnames(covariates_all_years_new_chinook) <- colnames(covariates_all_years_chinook0)
  return(covariates_all_years_new_chinook)
  
}


covariates_all_years_new_chinook <- make_poly(covariates_all_years_chinook0)


#list of models to try

# no covariates in c, no covariates in d 

# temp, photoperiod, atu, photo_diff, 



mod.list_0_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal"
)

# no covariates in c but flow in d

mod.list_0_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  D = "diagonal and equal",
  d = covariates_all_years_new_chinook[31:45,] 
)



mod.list_1_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = "diagonal and equal"
)


mod.list_1_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = "diagonal and equal",
  D = "diagonal and equal",
  d = covariates_all_years_new_chinook[31:45,] 
)


mod.list_2_0 <- list(
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
             15,30, byrow = TRUE)
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
  D = "diagonal and equal",
  d = covariates_all_years_new_chinook[31:45,] 
)


mod.list_3_0 <- list(
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
             15,45, byrow = TRUE)
  
)


mod.list_3_1 <- list(
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
  D = "diagonal and equal",
  d = covariates_all_years_new_chinook[31:45,] 
)

c.models <- list(
  c1 = "zero",
  c2 = "zero",
  c3 = covariates_all_years_new_chinook[1:15,], # Temperature 
  c4 = covariates_all_years_new_chinook[31:45,], # flow
  c5 = covariates_all_years_new_chinook[61:75,], # photoperiod
  c6 = covariates_all_years_new_chinook[91:105,], # atu
  c7 = covariates_all_years_new_chinook[106:120,],  # lunar phase
  c8 = covariates_all_years_new_chinook[121:135,],  # resid
  c9 = covariates_all_years_new_chinook[136:150,],  # temp diff
  c10 = covariates_all_years_new_chinook[151:165,],  # flowdiff
  c11 = covariates_all_years_new_chinook[166:180,],  # photo diff
  c12 = covariates_all_years_new_chinook[181:195,]  # hatchery
)
names(c.models) <- c("No covariates", 
                     "No covariates in c, flow in d", 
                     "Temperature", 
                     "Flow",
                     "Photoperiod",
                     "ATU", 
                     "Lunar phase", 
                     "Temperature residuals", 
                     "Temperature difference",
                     "Flow difference",
                     "Photoperiod difference",
                     "Hatchery")

d.models <- list(
  d1 = "zero",
  d2 = covariates_all_years_new_chinook[31:45,], # flow
  d3 = "zero",
  d4 = "zero",
  d5 = "zero",
  d6 = "zero",
  d7 = "zero",
  d8 = "zero",
  d9 = "zero",
  d10 = "zero",
  d11 = "zero",
  d12 = "zero"
)

names(d.models) <- c("No covariates", 
                     "Flow", 
                     "No covariates",
                     "No covariates",
                     "No covariates",
                     "No covariates",
                     "No covariates",
                     "No covariates",
                     "No covariates",
                     "No covariates",
                     "No covariates",
                     "No covariates")

c_num = c(0,0,1,1,1,1,1,1,1,1,1,1)
d_num = c(0,1,0,0,0,0,0,0,0,0,0,0)

out.tab_chinook1 <- NULL
fits_chinook1 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 0 & d_num[i] == 0){
    
    fit.model = mod.list_0_0
  }
  
  else if(c_num[i] == 0 & d_num[i] == 1){
    
    fit.model = mod.list_0_1
  }
  
  else if(c_num[i] == 1 & d_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_1_0)
  }
  
  else if(c_num[i] == 1 & d_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_1_1)
  }
  
  else if(c_num[i] == 2 & d_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_2_1)
  }
  
  
  
  
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d=names(d.models)[i],
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook1=rbind(out.tab_chinook1,out)
  fits_chinook1=c(fits_chinook1,list(fit))
  
  
}




min.AICc <- order(out.tab_chinook1$AICc)
out.tab.chinook1 <- out.tab_chinook1[min.AICc, ]
out.tab.chinook1



c.models <- list(
  c3 = covariates_all_years_new_chinook[1:15,], # Temperature 
  c4 = covariates_all_years_new_chinook[31:45,], # flow
  c5 = covariates_all_years_new_chinook[61:75,], # photoperiod
  c6 = covariates_all_years_new_chinook[91:105,], # atu
  c7 = covariates_all_years_new_chinook[106:120,],  # lunar phase
  c8 = covariates_all_years_new_chinook[121:135,],  # resid
  c9 = covariates_all_years_new_chinook[136:150,],  # temp diff
  c10 = covariates_all_years_new_chinook[151:165,],  # flowdiff
  c11 = covariates_all_years_new_chinook[166:180,] # photo diff
  
)
names(c.models) <- c("Temperature", 
                     "Flow",
                     "Photoperiod",
                     "ATU", 
                     "Lunar phase", 
                     "Temperature residuals", 
                     "Temperature difference",
                     "Flow difference",
                     "Photoperiod difference")


d.models <- list(
  c1 = covariates_all_years_new_chinook[31:45,], # flow
  c2 = covariates_all_years_new_chinook[31:45,], # flow
  c3 = covariates_all_years_new_chinook[31:45,], # flow
  c4 = covariates_all_years_new_chinook[31:45,], # flow
  c5 = covariates_all_years_new_chinook[31:45,], # flow
  c6 = covariates_all_years_new_chinook[31:45,], # flow
  c7 = covariates_all_years_new_chinook[31:45,], # flow
  c8 = covariates_all_years_new_chinook[31:45,], # flow
  c9 = covariates_all_years_new_chinook[31:45,] # flow
  
)
names(d.models) <- c("Flow", 
                     "Flow",
                     "Flow",
                     "Flow",
                     "Flow",
                     "Flow",
                     "Flow",
                     "Flow",
                     "Flow")

c_num = c(1,1,1,1,1,1,1,1,1)
d_num = c(1,1,1,1,1,1,1,1,1)

out.tab_chinook2 <- NULL
fits_chinook2 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 0 & d_num[i] == 0){
    
    fit.model = mod.list_0_0
  }
  
  else if(c_num[i] == 0 & d_num[i] == 1){
    
    fit.model = mod.list_0_1
  }
  
  else if(c_num[i] == 1 & d_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_1_0)
  }
  
  else if(c_num[i] == 1 & d_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_1_1)
  }
  
  else if(c_num[i] == 2 & d_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_2_1)
  }
  
  
  
  
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d=names(d.models)[i],
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook2=rbind(out.tab_chinook2,out)
  fits_chinook2=c(fits_chinook2,list(fit))
  
  
}




min.AICc <- order(out.tab_chinook2$AICc)
out.tab.chinook2 <- out.tab_chinook2[min.AICc, ]
out.tab.chinook2

# let us have three global models with flow in observation covariate

# 1 - temp, flow, temp diff, flow diff, photo diff, resid, lunar phase
# 2 - atu, flow, temp diff, flow diff, photo diff, resid, lunar phase
# 3 - photoperiod, flow, temp diff, flow diff, photo diff, resid, lunar phase



mod.list_7_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"g"),
             15,105, byrow = TRUE),
  D = "diagonal and equal",
  d = covariates_all_years_new_chinook[31:45,] 
)


mod.list_7_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"g"),
             15,105, byrow = TRUE)
)


c.models <- list(
  c1 = covariates_all_years_new_chinook[c(1:15,31:45,106:180),], # Temperature 
  c2 = covariates_all_years_new_chinook[c(61:75,31:45,106:180),], # photoperiod
  c3 = covariates_all_years_new_chinook[c(91:105,31:45,106:180),], # atu
  c4 = covariates_all_years_new_chinook[c(1:15,31:45,106:180),], # Temperature 
  c5 = covariates_all_years_new_chinook[c(61:75,31:45,106:180),], # photoperiod
  c6 = covariates_all_years_new_chinook[c(91:105,31:45,106:180),] # atu
  
  
  
)
names(c.models) <- c("Temperature, Flow, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference",
                     "Photoperiod, Flow, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference",
                     "ATU, Flow, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference",
                     "Temperature, Flow, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference",
                     "Photoperiod, Flow, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference",
                     "ATU, Flow, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference"
                     )


d.models <- list(
  d1 = covariates_all_years_new_chinook[31:45,], # flow
  d2 = covariates_all_years_new_chinook[31:45,], # flow
  d3 = covariates_all_years_new_chinook[31:45,],
  d4 = "zero",
  d5 = "zero",
  d6 = "zero"
  
)
names(d.models) <- c("Flow", 
                     "Flow",
                     "Flow",
                     "No covariates",
                     "No covariates",
                     "No covariates")

c_num = c(7,7,7,7,7,7)
d_num = c(1,1,1,0,0,0)


out.tab_chinook3 <- NULL
fits_chinook3 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 7 & d_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_7_1)
  }
  
  else if(c_num[i] == 7 & d_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_7_0)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d=names(d.models)[i],
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook3=rbind(out.tab_chinook3,out)
  fits_chinook3=c(fits_chinook3,list(fit))
  
  
}




min.AICc <- order(out.tab_chinook3$AICc)
out.tab.chinook3 <- out.tab_chinook3[min.AICc, ]
out.tab.chinook3



mod.list_9_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i"),
             15,135, byrow = TRUE)
)


mod.list_9_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i"),
             15,135, byrow = TRUE),
  D = "diagonal and equal",
  d = covariates_all_years_new_chinook[31:45,] 
)

c.models <- list(
  c1 = covariates_all_years_new_chinook[c(1:60,106:180),], # Temperature 
  c2 = covariates_all_years_new_chinook[c(61:90,31:60,106:180),], # photoperiod
  c3 = covariates_all_years_new_chinook[c(91:105,31:60,106:180),], # atu
  c4 = covariates_all_years_new_chinook[c(1:60,106:180),], # Temperature 
  c5 = covariates_all_years_new_chinook[c(61:90,31:60,106:180),], # photoperiod
  c6 = covariates_all_years_new_chinook[c(91:105,31:60,106:180),] # atu
  
  
  
)
names(c.models) <- c("Temperature, Temperature^2, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference",
                     "Photoperiod, Photoperiod^2, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference",
                     "ATU, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference",
                     "Temperature, Temperature^2, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference",
                     "Photoperiod, Photoperiod^2, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference",
                     "ATU, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference"
)


d.models <- list(
  d1 = covariates_all_years_new_chinook[31:45,], # flow
  d2 = covariates_all_years_new_chinook[31:45,], # flow
  d3 = covariates_all_years_new_chinook[31:45,],
  d4 = "zero",
  d5 = "zero",
  d6 = "zero"
  
)
names(d.models) <- c("Flow", 
                     "Flow",
                     "Flow",
                     "No covariates",
                     "No covariates",
                     "No covariates")

c_num = c(9,9,8,9,9,8)
d_num = c(1,1,1,0,0,0)



out.tab_chinook4 <- NULL
fits_chinook4 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 9 & d_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_9_1)
  }
  
  else if(c_num[i] == 9 & d_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_9_0)
  }
  else if(c_num[i] == 8 & d_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_8_0)
  }
  else if(c_num[i] == 8 & d_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_8_1)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d=names(d.models)[i],
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook4=rbind(out.tab_chinook4,out)
  fits_chinook4=c(fits_chinook4,list(fit))
  
  
}




min.AICc <- order(out.tab_chinook4$AICc)
out.tab.chinook4 <- out.tab_chinook4[min.AICc, ]
out.tab.chinook4






mod.list_8_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h"),
             15,120, byrow = TRUE)
)

mod.list_8_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h"),
             15,120, byrow = TRUE),
  D = "diagonal and equal",
  d = covariates_all_years_new_chinook[31:45,] 
)



mod.list_9_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2020"),
             15,135, byrow = TRUE)
)


mod.list_9_1_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2020"),
             15,135, byrow = TRUE),
  D = "diagonal and equal",
  d = covariates_all_years_new_chinook[31:45,] 
)


mod.list_10_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "h2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"h2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"h2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"h2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"h2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"h2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"h2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"h2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"h2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"h2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"h2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"h2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"h2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"h2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2020"),
             15,150, byrow = TRUE)
)


mod.list_10_1 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "h2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"h2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"h2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"h2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"h2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"h2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"h2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"h2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"h2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"h2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"h2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"h2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"h2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"h2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"g",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"i",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2020"),
             15,150, byrow = TRUE),
  D = "diagonal and equal",
  d = covariates_all_years_new_chinook[31:45,] 
)


c.models <- list(
  c1 = covariates_all_years_new_chinook[c(1:60,106:195),], # Temperature 
  c2 = covariates_all_years_new_chinook[c(61:90,31:60,106:195),], # photoperiod
  c3 = covariates_all_years_new_chinook[c(91:105,31:60,106:195),], # atu
  c4 = covariates_all_years_new_chinook[c(1:60,106:195),], # Temperature 
  c5 = covariates_all_years_new_chinook[c(61:90,31:60,106:195),], # photoperiod
  c6 = covariates_all_years_new_chinook[c(91:105,31:60,106:195),] # atu
  
  
  
)
names(c.models) <- c("Temperature, Temperature^2, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference, Hatchery",
                     "Photoperiod, Photoperiod^2, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference, Hatchery",
                     "ATU, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference, Hatchery",
                     "Temperature, Temperature^2, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference, Hatchery",
                     "Photoperiod, Photoperiod^2, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference, Hatchery",
                     "ATU, Flow, Flow^2, Lunar phase, Temperature residuals, 
                     Temperature difference, Flow difference, Photoperiod difference, Hatchery"
)


d.models <- list(
  d1 = covariates_all_years_new_chinook[31:45,], # flow
  d2 = covariates_all_years_new_chinook[31:45,], # flow
  d3 = covariates_all_years_new_chinook[31:45,],
  d4 = "zero",
  d5 = "zero",
  d6 = "zero"
  
)
names(d.models) <- c("Flow", 
                     "Flow",
                     "Flow",
                     "No covariates",
                     "No covariates",
                     "No covariates")

c_num = c(10,10,9,10,10,9)
d_num = c(1,1,1,0,0,0)

out.tab_chinook5 <- NULL
fits_chinook5 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 9 & d_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_9_1_h)
  }
  
  else if(c_num[i] == 9 & d_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_9_0_h)
  }
  else if(c_num[i] == 10 & d_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_10_0)
  }
  else if(c_num[i] == 10 & d_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_10_1)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d=names(d.models)[i],
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook5=rbind(out.tab_chinook5,out)
  fits_chinook5=c(fits_chinook5,list(fit))
  
  
}




min.AICc <- order(out.tab_chinook5$AICc)
out.tab.chinook5 <- out.tab_chinook5[min.AICc, ]
out.tab.chinook5


################################
#note for later - standard deviation after using poly is not 1



#################################


pairs(data[,29:39])

#temp, ATU, photo, resid, photo_diff are correlated
#flow and flow_diff and correleated

# let's leave out temp square, flow_square and photoperiod square for now

# use covariates all years chinook0

#make different covariate matrices from 5 to 0. Let's leave out observation covariates for now


#matrix with 5 cov (w hatchery)

#############


mod.list_5_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2020"),
             15,75, byrow = TRUE)
)


mod.list_5_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"e"),
             15,75, byrow = TRUE)
)




##########

# matrix with 4 cov (w and wo hatchery)

#####


mod.list_4_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2020"),
             15,60, byrow = TRUE)
)


mod.list_4_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d"),
             15,60, byrow = TRUE)
)


######

# matrix with 3 cov (w and wo hatchery)


#####


mod.list_3_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  "d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"d",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2020"),
             15,45, byrow = TRUE)
)


mod.list_3_0 <- list(
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
                  
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c"),
             15,45, byrow = TRUE)
)


######

# matrix with 2 cov (w and wo hatchery)


#####


mod.list_2_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = matrix(list("a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2020"),
             15,30, byrow = TRUE)
)


mod.list_2_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
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
             15,30, byrow = TRUE)
)


######


#matrix with 1 cov (w and wo hatchery)

mod.list_1_0_h <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = "diagonal and unequal"
)

mod.list_1_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal",
  C = "diagonal and equal"
)

mod.list_0_0 <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and equal"
)


c.models <- list(
  c1 = covariates_all_years_chinook0[c(1:30,61:75, 91:105,136:150),], #temp, flow, lp, temp diff,  hatchery 
  c2 = covariates_all_years_chinook0[c(31:45, 16:30, 61:75, 91:105,136:150),], #photoperiod, flow, lp, temp diff,  hatchery
  c3 = covariates_all_years_chinook0[c(46:60,16:30, 61:75, 91:105,136:150),], # atu, flow, lp, temp diff,  hatchery
  c4 = covariates_all_years_chinook0[c(76:90,16:30, 61:75, 91:105,136:150),], #resid,  flow, lp, temp diff,  hatchery
  c5 = covariates_all_years_chinook0[c(121:135,16:30, 61:75, 91:105,136:150),] # photodiff, flow, lp, temp diff,  hatchery
  
  
  
  
)
names(c.models) <- c("Temperature, Flow, Lunar phase, Temperature difference, Hatchery",
                     "Photoperiod, Flow, Lunar phase, Temperature difference, Hatchery",
                     "ATU, Flow, Lunar phase, Temperature difference, Hatchery",
                     "Temperature residuals, Flow, Lunar phase, Temperature difference, Hatchery",
                     "Photoperiod difference, Flow, Lunar phase, Temperature difference, Hatchery"
)


out.tab_chinook_new_0 <- NULL
fits_chinook_new_0 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  
    
  fit.model = c(list(c= c.models[[i]]), mod.list_5_0_h)
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_0=rbind(out.tab_chinook_new_0,out)
  fits_chinook_new_0=c(fits_chinook_new_0,list(fit))
  
  
}




min.AICc <- order(out.tab_chinook_new_0$AICc)
out.tab.chinook_new_0 <- out.tab_chinook_new_0[min.AICc, ]
out.tab.chinook_new_0

autoplot(fits_chinook_new_0[[5]])

ci <- tidy(fits_chinook_new_0[[5]]) # does not work

# leaving out temp diff

c.models <- list(
  c1 = covariates_all_years_chinook0[c(1:30,61:75, 91:105),], #temp, flow, lp, temp diff, 
  c2 = covariates_all_years_chinook0[c(31:45, 16:30, 61:75, 91:105),], #photoperiod, flow, lp, temp diff,  
  c3 = covariates_all_years_chinook0[c(46:60,16:30, 61:75, 91:105),], # atu, flow, lp, temp diff,  
  c4 = covariates_all_years_chinook0[c(76:90,16:30, 61:75, 91:105),], #resid,  flow, lp, temp diff,  
  c5 = covariates_all_years_chinook0[c(121:135,16:30, 61:75, 91:105),], # photodiff, flow, lp, temp diff,
  c6 = covariates_all_years_chinook0[c(1:30,61:75,  136:150),], #temp, flow, lp,    hatchery 
  c7 = covariates_all_years_chinook0[c(31:45, 16:30, 61:75,  136:150),], #photoperiod, flow, lp,    hatchery
  c8 = covariates_all_years_chinook0[c(46:60,16:30, 61:75,  136:150),], # atu, flow, lp,    hatchery
  c9 = covariates_all_years_chinook0[c(76:90,16:30, 61:75,  136:150),], #resid,  flow, lp,    hatchery
  c10 = covariates_all_years_chinook0[c(121:135,16:30, 61:75,  136:150),], # photodiff, flow, lp,    hatchery
  c11 = covariates_all_years_chinook0[c(1:30,61:75 ),], #temp, flow, lp,   
  c12 = covariates_all_years_chinook0[c(31:45, 16:30, 61:75 ),], #photoperiod, flow, lp,    
  c13 = covariates_all_years_chinook0[c(46:60,16:30, 61:75),], # atu, flow, lp,    
  c14 = covariates_all_years_chinook0[c(76:90,16:30, 61:75),], #resid,  flow, lp,    
  c15 = covariates_all_years_chinook0[c(121:135,16:30, 61:75),] # photodiff, flow, lp,  
  
  
  
  
)
names(c.models) <- c("Temperature, Flow, Lunar phase, Temperature difference ",
                     "Photoperiod, Flow, Lunar phase, Temperature difference ",
                     "ATU, Flow, Lunar phase, Temperature difference ",
                     "Temperature residuals, Flow, Lunar phase, Temperature difference ",
                     "Photoperiod difference, Flow, Lunar phase, Temperature difference ",
                     
                     "Temperature, Flow, Lunar phase,  Hatchery",
                     "Photoperiod, Flow, Lunar phase,  Hatchery",
                     "ATU, Flow, Lunar phase,  Hatchery",
                     "Temperature residuals, Flow, Lunar phase,  Hatchery",
                     "Photoperiod difference, Flow, Lunar phase,  Hatchery",
                     
                     "Temperature, Flow, Lunar phase",
                     "Photoperiod, Flow, Lunar phase",
                     "ATU, Flow, Lunar phase",
                     "Temperature residuals, Flow, Lunar phase",
                     "Photoperiod difference, Flow, Lunar phase"
)

c_num = c(4,4,4,4,4,4,4,4,4,4,3,3,3,3,3)
h_num = c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0)


out.tab_chinook_new_1 <- NULL
fits_chinook_new_1 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 4 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_4_0)
  }
  
  else if(c_num[i] == 4 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_4_0_h)
  }
  else if(c_num[i] == 3 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_3_0)
  }
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_1=rbind(out.tab_chinook_new_1,out)
  fits_chinook_new_1=c(fits_chinook_new_1,list(fit))
  
  
}



out.tab_chinook_new_1 <- rbind(out.tab_chinook_new_1,out.tab_chinook_new_0)
min.AICc <- order(out.tab_chinook_new_1$AICc)
out.tab.chinook_new_1 <- out.tab_chinook_new_1[min.AICc, ]
out.tab.chinook_new_1

#leaving out lunar phase

c.models <- list(
  c1 = covariates_all_years_chinook0[c(1:30,  91:105),], #temp, flow, temp diff, 
  c2 = covariates_all_years_chinook0[c(31:45, 16:30,   91:105),], #photoperiod, flow, temp diff,  
  c3 = covariates_all_years_chinook0[c(46:60,16:30,   91:105),], # atu, flow, temp diff,  
  c4 = covariates_all_years_chinook0[c(76:90,16:30,   91:105),], #resid,  flow, temp diff,  
  c5 = covariates_all_years_chinook0[c(121:135,16:30,   91:105),], # photodiff, flow, temp diff,
  c6 = covariates_all_years_chinook0[c(1:30,  91:105,   136:150),], #temp, flow,    hatchery 
  c7 = covariates_all_years_chinook0[c(31:45, 16:30,  91:105,    136:150),], #photoperiod, flow,    hatchery
  c8 = covariates_all_years_chinook0[c(46:60,16:30,  91:105,    136:150),], # atu, flow,    hatchery
  c9 = covariates_all_years_chinook0[c(76:90,16:30,  91:105,    136:150),], #resid,  flow,    hatchery
  c10 = covariates_all_years_chinook0[c(121:135,16:30,  91:105,    136:150),] # photodiff, flow,    hatchery
  
  
  
  
)

names(c.models) <- c("Temperature, Flow,  Temperature difference ",
                     "Photoperiod, Flow,  Temperature difference ",
                     "ATU, Flow,  Temperature difference ",
                     "Temperature residuals, Flow,  Temperature difference ",
                     "Photoperiod difference, Flow,  Temperature difference ",
                     
                     "Temperature, Flow,  Temperature difference,   Hatchery",
                     "Photoperiod, Flow,  Temperature difference,   Hatchery",
                     "ATU, Flow,  Temperature difference,   Hatchery",
                     "Temperature residuals, Flow,  Temperature difference,   Hatchery",
                     "Photoperiod difference, Flow,  Temperature difference,   Hatchery"
                     
)


c_num = c(3,3,3,3,3,4,4,4,4,4)
h_num = c(0,0,0,0,0,1,1,1,1,1)

out.tab_chinook_new_2 <- NULL
fits_chinook_new_2 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 3 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_3_0)
  }
  
  else if(c_num[i] == 4 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_4_0_h)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_2=rbind(out.tab_chinook_new_2,out)
  fits_chinook_new_2=c(fits_chinook_new_2,list(fit))
  
  
}



out.tab_chinook_new_2 <- rbind(out.tab_chinook_new_2,out.tab_chinook_new_1)
min.AICc <- order(out.tab_chinook_new_2$AICc)
out.tab.chinook_new_2 <- out.tab_chinook_new_2[min.AICc, ]
out.tab.chinook_new_2




#leaving out flow


c.models <- list(
  c1 = covariates_all_years_chinook0[c(1:15,61:75, 91:105,136:150),], #temp, lp, temp diff,  hatchery 
  c2 = covariates_all_years_chinook0[c(31:45, 61:75, 91:105,136:150),], #photoperiod, lp, temp diff,  hatchery
  c3 = covariates_all_years_chinook0[c(46:60,61:75, 91:105,136:150),], # atu, lp, temp diff,  hatchery
  c4 = covariates_all_years_chinook0[c(76:90,61:75, 91:105,136:150),], #resid,  lp, temp diff,  hatchery
  c5 = covariates_all_years_chinook0[c(121:135,61:75, 91:105,136:150),], # photodiff, lp, temp diff,  hatchery
  c6 = covariates_all_years_chinook0[c(1:15,61:75, 91:105),], #temp, lp, temp diff,  hatchery 
  c7 = covariates_all_years_chinook0[c(31:45, 61:75, 91:105),], #photoperiod, lp, temp diff, 
  c8 = covariates_all_years_chinook0[c(46:60,61:75, 91:105),], # atu, lp, temp diff,  
  c9 = covariates_all_years_chinook0[c(76:90,61:75, 91:105),], #resid,  lp, temp diff,  
  c10 = covariates_all_years_chinook0[c(121:135,61:75, 91:105),] # photodiff, lp, temp diff, 
  
  
  
  
  
)
names(c.models) <- c("Temperature,  Lunar phase, Temperature difference,   Hatchery",
                     "Photoperiod,  Lunar phase, Temperature difference,   Hatchery",
                     "ATU,  Lunar phase, Temperature difference,   Hatchery",
                     "Temperature residuals,  Lunar phase, Temperature difference,   Hatchery",
                     "Photoperiod difference,  Lunar phase, Temperature difference,   Hatchery",
                     
                     "Temperature,  Lunar phase, Temperature difference ",
                     "Photoperiod,  Lunar phase, Temperature difference ",
                     "ATU,  Lunar phase, Temperature difference ",
                     "Temperature residuals,  Lunar phase, Temperature difference ",
                     "Photoperiod difference,  Lunar phase, Temperature difference "
                     
                     
                     
)


c_num = c(4,4,4,4,4,3,3,3,3,3)
h_num = c(1,1,1,1,1,0,0,0,0,0)


out.tab_chinook_new_3 <- NULL
fits_chinook_new_3 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 3 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_3_0)
  }
  
  else if(c_num[i] == 4 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_4_0_h)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_3=rbind(out.tab_chinook_new_3,out)
  fits_chinook_new_3=c(fits_chinook_new_3,list(fit))
  
  
}



out.tab_chinook_new_3 <- rbind(out.tab_chinook_new_3,out.tab_chinook_new_2)
min.AICc <- order(out.tab_chinook_new_3$AICc)
out.tab.chinook_new_3 <- out.tab_chinook_new_3[min.AICc, ]
out.tab.chinook_new_3


#leaving out temp/atu/photoperiod/resid/photodiff


c.models <- list(
  c1 = covariates_all_years_chinook0[c(16:30,61:75, 91:105,136:150),], #flow, lp, temp diff,  hatchery 
  
  c2 = covariates_all_years_chinook0[c(16:30,61:75, 91:105),] # flow, lp, temp diff, 
  
)

names(c.models) <- c("Flow,  Lunar phase, Temperature difference,  Hatchery",
                     
                     
                     "Flow,  Lunar phase, Temperature difference"              
                     
)


c_num = c(4,3)
h_num = c(1,0)



out.tab_chinook_new_4 <- NULL
fits_chinook_new_4 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 3 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_3_0)
  }
  
  else if(c_num[i] == 4 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_4_0_h)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_4=rbind(out.tab_chinook_new_4,out)
  fits_chinook_new_4=c(fits_chinook_new_4,list(fit))
  
  
}



out.tab_chinook_new_4 <- rbind(out.tab_chinook_new_4,out.tab_chinook_new_3)
min.AICc <- order(out.tab_chinook_new_4$AICc)
out.tab.chinook_new_4 <- out.tab_chinook_new_4[min.AICc, ]
out.tab.chinook_new_4



# leaving out flow and lunar phase




c.models <- list(
  c1 = covariates_all_years_chinook0[c(1:15, 91:105,136:150),], #temp,  temp diff,  hatchery 
  c2 = covariates_all_years_chinook0[c(31:45,  91:105,136:150),], #photoperiod,  temp diff,  hatchery
  c3 = covariates_all_years_chinook0[c(46:60, 91:105,136:150),], # atu,  temp diff,  hatchery
  c4 = covariates_all_years_chinook0[c(76:90, 91:105,136:150),], #resid,   temp diff,  hatchery
  c5 = covariates_all_years_chinook0[c(121:135, 91:105,136:150),], # photodiff,  temp diff,  hatchery
  c6 = covariates_all_years_chinook0[c(1:15, 91:105),], #temp,  temp diff,  hatchery 
  c7 = covariates_all_years_chinook0[c(31:45,  91:105),], #photoperiod,  temp diff, 
  c8 = covariates_all_years_chinook0[c(46:60, 91:105),], # atu,  temp diff,  
  c9 = covariates_all_years_chinook0[c(76:90, 91:105),], #resid,   temp diff,  
  c10 = covariates_all_years_chinook0[c(121:135, 91:105),] # photodiff,  temp diff, 
  

)
names(c.models) <- c("Temperature, Temperature difference,   Hatchery",
                     "Photoperiod, Temperature difference,   Hatchery",
                     "ATU, Temperature difference,   Hatchery",
                     "Temperature residuals, Temperature difference,   Hatchery",
                     "Photoperiod difference, Temperature difference,   Hatchery",
                     
                     "Temperature, Temperature difference ",
                     "Photoperiod, Temperature difference ",
                     "ATU, Temperature difference ",
                     "Temperature residuals, Temperature difference ",
                     "Photoperiod difference, Temperature difference "
                     
                     
                     
)


c_num = c(3,3,3,3,3,2,2,2,2,2)
h_num = c(1,1,1,1,1,0,0,0,0,0)


out.tab_chinook_new_5 <- NULL
fits_chinook_new_5 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 3 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_3_0_h)
  }
  
  else if(c_num[i] == 2 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_2_0)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_5=rbind(out.tab_chinook_new_5,out)
  fits_chinook_new_5=c(fits_chinook_new_5,list(fit))
  
  
}



out.tab_chinook_new_5 <- rbind(out.tab_chinook_new_5,out.tab_chinook_new_4)
min.AICc <- order(out.tab_chinook_new_5$AICc)
out.tab.chinook_new_5 <- out.tab_chinook_new_5[min.AICc, ]
out.tab.chinook_new_5



#leaving out temp diff and flow


c.models <- list(
  c1 = covariates_all_years_chinook0[c(1:15,61:75, 136:150),], #temp, lp,   hatchery 
  c2 = covariates_all_years_chinook0[c(31:45, 61:75, 136:150),], #photoperiod, lp,   hatchery
  c3 = covariates_all_years_chinook0[c(46:60,61:75, 136:150),], # atu, lp,   hatchery
  c4 = covariates_all_years_chinook0[c(76:90,61:75, 136:150),], #resid,  lp,   hatchery
  c5 = covariates_all_years_chinook0[c(121:135,61:75, 136:150),], # photodiff, lp,   hatchery
  c6 = covariates_all_years_chinook0[c(1:15,61:75),], #temp, lp,   hatchery 
  c7 = covariates_all_years_chinook0[c(31:45, 61:75),], #photoperiod, lp,  
  c8 = covariates_all_years_chinook0[c(46:60,61:75),], # atu, lp,   
  c9 = covariates_all_years_chinook0[c(76:90,61:75),], #resid,  lp,   
  c10 = covariates_all_years_chinook0[c(121:135,61:75),] # photodiff, lp,  
  
  
)
names(c.models) <- c("Temperature,  Lunar phase,   Hatchery",
                     "Photoperiod,  Lunar phase,   Hatchery",
                     "ATU,  Lunar phase,   Hatchery",
                     "Temperature residuals,  Lunar phase,   Hatchery",
                     "Photoperiod difference,  Lunar phase,   Hatchery",
                     
                     "Temperature,  Lunar phase ",
                     "Photoperiod,  Lunar phase ",
                     "ATU,  Lunar phase ",
                     "Temperature residuals,  Lunar phase ",
                     "Photoperiod difference,  Lunar phase "
                     
                     
                     
)



c_num = c(3,3,3,3,3,2,2,2,2,2)
h_num = c(1,1,1,1,1,0,0,0,0,0)


out.tab_chinook_new_6 <- NULL
fits_chinook_new_6 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 3 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_3_0_h)
  }
  
  else if(c_num[i] == 2 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_2_0)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_6=rbind(out.tab_chinook_new_6,out)
  fits_chinook_new_6=c(fits_chinook_new_6,list(fit))
  
  
}



out.tab_chinook_new_6 <- rbind(out.tab_chinook_new_6,out.tab_chinook_new_5)
min.AICc <- order(out.tab_chinook_new_6$AICc)
out.tab.chinook_new_6 <- out.tab_chinook_new_6[min.AICc, ]
out.tab.chinook_new_6



#leaving out temp diff and temp/atu/...



c.models <- list(
  c1 = covariates_all_years_chinook0[c(16:30,61:75, 136:150),], #flow, lp,   hatchery 
  c2 = covariates_all_years_chinook0[c(16:30,61:75),] #flow, lp

  
)

names(c.models) <- c("Flow,  Lunar phase,   Hatchery",
                     "Flow,  Lunar phase"
                     
)




c_num = c(3,2)
h_num = c(1,0)

out.tab_chinook_new_7 <- NULL
fits_chinook_new_7 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 3 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_3_0_h)
  }
  
  else if(c_num[i] == 2 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_2_0)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_7=rbind(out.tab_chinook_new_7,out)
  fits_chinook_new_7=c(fits_chinook_new_7,list(fit))
  
  
}



out.tab_chinook_new_7 <- rbind(out.tab_chinook_new_7,out.tab_chinook_new_6)
min.AICc <- order(out.tab_chinook_new_7$AICc)
out.tab.chinook_new_7 <- out.tab_chinook_new_7[min.AICc, ]
out.tab.chinook_new_7


#leaving out temp/atu ... and lp
#leaving out temp/atu... and flow




c.models <- list(
  c1 = covariates_all_years_chinook0[c(16:30,91:105, 136:150),], #flow, temp diff,   hatchery 
  c2 = covariates_all_years_chinook0[c(16:30,91:105),], #flow, temp diff
  c3 = covariates_all_years_chinook0[c(61:75,91:105, 136:150),], #lp,  temp diff. hatchery 
  c4 = covariates_all_years_chinook0[c(61:75,91:105),] #lp, temp diff
  
  
)

names(c.models) <- c("Flow,  Temperature difference,   Hatchery",
                     "Flow,  Temperature difference",
                     "Lunar phase,  Temperature difference,   Hatchery",
                     "Lunar phase,  Temperature difference"
                     
)




c_num = c(3,2,3,2)
h_num = c(1,0,1,0)

out.tab_chinook_new_8 <- NULL
fits_chinook_new_8 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 3 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_3_0_h)
  }
  
  else if(c_num[i] == 2 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_2_0)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_8=rbind(out.tab_chinook_new_8,out)
  fits_chinook_new_8=c(fits_chinook_new_8,list(fit))
  
  
}



out.tab_chinook_new_8 <- rbind(out.tab_chinook_new_8,out.tab_chinook_new_7)
min.AICc <- order(out.tab_chinook_new_8$AICc)
out.tab.chinook_new_8 <- out.tab_chinook_new_8[min.AICc, ]
out.tab.chinook_new_8


#leaving out temp diff and lp


c.models <- list(
  c1 = covariates_all_years_chinook0[c(1:15,16:30, 136:150),], #temp, flow,   hatchery 
  c2 = covariates_all_years_chinook0[c(31:45,16:30, 136:150),], #photoperiod, flow,   hatchery
  c3 = covariates_all_years_chinook0[c(46:60,16:30, 136:150),], # atu, flow,   hatchery
  c4 = covariates_all_years_chinook0[c(76:90,16:30, 136:150),], #resid, flow,   hatchery
  c5 = covariates_all_years_chinook0[c(121:135,16:30, 136:150),], # photodiff, flow,   hatchery
  c6 = covariates_all_years_chinook0[c(1:15,16:30),], #temp, flow,   
  c7 = covariates_all_years_chinook0[c(31:45,16:30),], #photoperiod, flow,  
  c8 = covariates_all_years_chinook0[c(46:60,16:30),], # atu, flow,   
  c9 = covariates_all_years_chinook0[c(76:90,16:30),], #resid, flow,   
  c10 = covariates_all_years_chinook0[c(121:135,16:30),] # photodiff, flow,  
  
  
)

names(c.models) <- c("Temperature, Flow,   Hatchery",
                     "Photoperiod, Flow,   Hatchery",
                     "ATU, Flow,   Hatchery",
                     "Temperature residuals, Flow,   Hatchery",
                     "Photoperiod difference, Flow,   Hatchery",
                     
                     "Temperature, Flow ",
                     "Photoperiod, Flow ",
                     "ATU, Flow ",
                     "Temperature residuals, Flow ",
                     "Photoperiod difference, Flow "
                     
)



c_num = c(3,3,3,3,3,2,2,2,2,2)
h_num = c(1,1,1,1,1,0,0,0,0,0)


out.tab_chinook_new_9 <- NULL
fits_chinook_new_9 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 3 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_3_0_h)
  }
  
  else if(c_num[i] == 2 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_2_0)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_9=rbind(out.tab_chinook_new_9,out)
  fits_chinook_new_9=c(fits_chinook_new_9,list(fit))
  
  
}



out.tab_chinook_new_9 <- rbind(out.tab_chinook_new_9,out.tab_chinook_new_8)
min.AICc <- order(out.tab_chinook_new_9$AICc)
out.tab.chinook_new_9 <- out.tab_chinook_new_9[min.AICc, ]
out.tab.chinook_new_9


#only having temp/atu/.. 

c.models <- list(
  c1 = covariates_all_years_chinook0[c(1:15, 136:150),], #temp   hatchery 
  c2 = covariates_all_years_chinook0[c(31:45, 136:150),], #photoperiod   hatchery
  c3 = covariates_all_years_chinook0[c(46:60, 136:150),], # atu   hatchery
  c4 = covariates_all_years_chinook0[c(76:90, 136:150),], #resid   hatchery
  c5 = covariates_all_years_chinook0[c(121:135, 136:150),], # photodiff   hatchery
  c6 = covariates_all_years_chinook0[c(1:15),], #temp   
  c7 = covariates_all_years_chinook0[c(31:45),], #photoperiod  
  c8 = covariates_all_years_chinook0[c(46:60),], # atu   
  c9 = covariates_all_years_chinook0[c(76:90),], #resid   
  c10 = covariates_all_years_chinook0[c(121:135),] # photodiff  
  
  
)

names(c.models) <- c("Temperature,   Hatchery",
                     "Photoperiod,   Hatchery",
                     "ATU,   Hatchery",
                     "Temperature residuals,   Hatchery",
                     "Photoperiod difference,   Hatchery",
                     
                     "Temperature ",
                     "Photoperiod ",
                     "ATU ",
                     "Temperature residuals ",
                     "Photoperiod difference "
                     
)


c_num = c(2,2,2,2,2,1,1,1,1,1)
h_num = c(1,1,1,1,1,0,0,0,0,0)


out.tab_chinook_new_10 <- NULL
fits_chinook_new_10 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 2 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_2_0_h)
  }
  
  else if(c_num[i] == 1 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_1_0)
  }
  
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_10=rbind(out.tab_chinook_new_10,out)
  fits_chinook_new_10=c(fits_chinook_new_10,list(fit))
  
  
}



out.tab_chinook_new_10 <- rbind(out.tab_chinook_new_10,out.tab_chinook_new_9)
min.AICc <- order(out.tab_chinook_new_10$AICc)
out.tab.chinook_new_10 <- out.tab_chinook_new_10[min.AICc, ]
out.tab.chinook_new_10


#only have temp diff
#only have flow
#only have lp

c.models <- list(
  c1 = covariates_all_years_chinook0[c(91:105, 136:150),], #temp diff,   hatchery 
  c2 = covariates_all_years_chinook0[c(16:30, 136:150),], #flow,  hatchery
  c3 = covariates_all_years_chinook0[c(61:75, 136:150),], #lp  hatchery
  c4 = covariates_all_years_chinook0[c(136:150),], #hatchery
  c5 = covariates_all_years_chinook0[c(91:105),], #temp diff 
  c6 = covariates_all_years_chinook0[c(16:30),], #flow
  c7 = covariates_all_years_chinook0[c(61:75),], #lp
  c8 = "zero"
  
  
)

names(c.models) <- c("Temperature difference,   Hatchery",
                     "Flow,   Hatchery",
                     "Lunar phase,   Hatchery",
                     "Hatchery",
                     "Temperature difference",
                     "Flow",
                     "Lunar phase",
                     "None"
                     
)


c_num = c(2,2,2,1,1,1,1,0)
h_num = c(1,1,1,1,0,0,0,0)



out.tab_chinook_new_11 <- NULL
fits_chinook_new_11 <- list()
for(i in 1:length(c.models)){
  
  print(names(c.models)[i])
  if(c_num[i] == 2 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_2_0_h)
  }
  
  else if(c_num[i] == 1 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_1_0)
  }
  
  else if(c_num[i] == 1 & h_num[i] == 1){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_1_0_h)
  }
  else if(c_num[i] == 0 & h_num[i] == 0){
    
    fit.model = c(list(c= c.models[[i]]), mod.list_0_0)
  }
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=names(c.models)[i], d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_chinook_new_11=rbind(out.tab_chinook_new_11,out)
  fits_chinook_new_11=c(fits_chinook_new_11,list(fit))
  
  
}



out.tab_chinook_new_11 <- rbind(out.tab_chinook_new_11,out.tab_chinook_new_10)
min.AICc <- order(out.tab_chinook_new_11$AICc)
out.tab.chinook_new_11 <- out.tab_chinook_new_11[min.AICc, ]
out.tab.chinook_new_11

weights <- akaike.weights(out.tab.chinook_new_11$AICc)

out.tab.chinook_new_11$deltaAICc <- weights$deltaAIC
out.tab.chinook_new_11$rel.LL <- weights$rel.LL
out.tab.chinook_new_11$weights <- weights$weights

searchString <- function(search_string) {
  matching_rows <- grep(search_string, out.tab.chinook_new_11$c, ignore.case = FALSE)
  return(out.tab.chinook_new_11[matching_rows, ])
}

rows <- searchString("ATU")
dim(rows)[1]

relative_importance <- data.frame(variable = c("Temperature difference", "Flow", "Lunar phase", "Hatchery"), 
                                               relative_importance_weight = c(
                                                 sum(searchString("Temperature difference")$weights),
                                  sum(searchString("Flow")$weights),
                                  sum(searchString("Lunar phase")$weights),
                                  sum(searchString("Hatchery")$weights))
                                  )

write.csv(out.tab.chinook_new_11, here("output","chinook_model_selection.csv"))


# trying diagonal and equal for hatchery


mod.list_3_0_trial <- list(
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
                  
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"a",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"b",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"c"),
             15,45, byrow = TRUE),
  c = covariates_all_years_chinook0[c(121:135, 91:105,136:150),]# photodiff,  temp diff,  hatchery
)

fit_trial <- MARSS(subset_chinook_summer_perhour, model=mod.list_3_0_trial, 
                   silent = TRUE, method = "BFGS",
                   control=list(maxit=2000))


#maybe we should try all the models but with same effect for hatchery for all years


get_covariate_combinations <- function(covariates) {
  n <- length(covariates)
  combinations <- lapply(1:n, function(x) combn(covariates, x, simplify = FALSE))
  unlist(combinations, recursive = FALSE)
}


list_combinations <- get_covariate_combinations(1:5)
out.tab<- NULL
fits <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  for(j in covariates){
    if(j == 1){
      k = 1
    }
    else if(j==2){
      k = 7
    }
    else if(j==3){
      k = 2
    }
    else if(j==4){
      k = 5
    }
    else if(j==5){
      k = 10
    }
    
    c = rbind(c,covariates_all_years_chinook0[((1+(k-1)*15):(k*15)),])
    name_long = rownames(covariates_all_years_chinook0)[1+(k-1)*15]
    name = paste(name, substr(name_long,1,nchar(name_long)-5))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(c_num == 1){
    
    fit.model = c(list(c= c), mod.list_1_0)
  }
  
  else if(c_num == 2){
    
    fit.model = c(list(c= c), mod.list_2_0)
  }
  else if(c_num == 3){
    
    fit.model = c(list(c= c), mod.list_3_0)
  }
  else if(c_num == 4){
    
    fit.model = c(list(c= c), mod.list_4_0)
  }
  else if(c_num == 5){
    
    fit.model = c(list(c= c), mod.list_5_0)
  }
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, d = "None",
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab=rbind(out.tab,out)
  fits=c(fits,list(fit))

}
out.tab

#trying with multiple temperature type variables


#######


list_combinations <- get_covariate_combinations(1:5)
out.tab<- NULL
fits <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  if(covariates[1] == 1){
    
    for(kk in c(1,3,4,6,9)){
      c = NULL
      name = NULL
      for(j in covariates){
        if(j == 1){
          k = kk
        }
        else if(j==2){
          k = 7
        }
        else if(j==3){
          k = 2
        }
        else if(j==4){
          k = 5
        }
        else if(j==5){
          k = 10
        }
        
        c = rbind(c,covariates_all_years_chinook0[((1+(k-1)*15):(k*15)),])
        name_long = rownames(covariates_all_years_chinook0)[1+(k-1)*15]
        name = paste(name, substr(name_long,1,nchar(name_long)-5))
        
      }
      # print(c)
      print(name)
      c_num <- length(covariates)
      if(c_num == 1){
        
        fit.model = c(list(c= c), mod.list_1_0)
      }
      
      else if(c_num == 2){
        
        fit.model = c(list(c= c), mod.list_2_0)
      }
      else if(c_num == 3){
        
        fit.model = c(list(c= c), mod.list_3_0)
      }
      else if(c_num == 4){
        
        fit.model = c(list(c= c), mod.list_4_0)
      }
      else if(c_num == 5){
        
        fit.model = c(list(c= c), mod.list_5_0)
      }
      
      fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                   control=list(maxit=2000))
      
      
      out=data.frame(c=name, d = "None",
                     logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                     num.iter=fit$numIter, converged=!fit$convergence,
                     stringsAsFactors = FALSE)
      out.tab=rbind(out.tab,out)
      fits=c(fits,list(fit))
    }
    
  }
  else{
    for(j in covariates){
      if(j == 1){
        k = 1
      }
      else if(j==2){
        k = 7
      }
      else if(j==3){
        k = 2
      }
      else if(j==4){
        k = 5
      }
      else if(j==5){
        k = 10
      }
      
      c = rbind(c,covariates_all_years_chinook0[((1+(k-1)*15):(k*15)),])
      name_long = rownames(covariates_all_years_chinook0)[1+(k-1)*15]
      name = paste(name, substr(name_long,1,nchar(name_long)-5))
      
    }
    # print(c)
    print(name)
    c_num <- length(covariates)
    if(c_num == 1){
      
      fit.model = c(list(c= c), mod.list_1_0)
    }
    
    else if(c_num == 2){
      
      fit.model = c(list(c= c), mod.list_2_0)
    }
    else if(c_num == 3){
      
      fit.model = c(list(c= c), mod.list_3_0)
    }
    else if(c_num == 4){
      
      fit.model = c(list(c= c), mod.list_4_0)
    }
    else if(c_num == 5){
      
      fit.model = c(list(c= c), mod.list_5_0)
    }
    
    fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                 control=list(maxit=2000))
    
    
    out=data.frame(c=name, d = "None",
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence,
                   stringsAsFactors = FALSE)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
    
  }
  
 
}
out.tab$deltaAICc <- out.tab$AICc - min(out.tab$AICc)
min.AICc <- order(out.tab$AICc)
out.tab.ordered <- out.tab[min.AICc, ]
out.tab.ordered
write.csv(out.tab.ordered, here("output","chinook0_model_selection.csv"))


#####



searchString <- function(search_string) {
  matching_rows <- grep(search_string, out.tab.ordered$c, ignore.case = FALSE)
  return(out.tab.ordered[matching_rows, ])
}

rows <- searchString("temp ")
dim(rows)[1]

# get balanced set with just photoperiod difference and with dataframe of indicator variables

#######




get_covariate_combinations <- function(covariates) {
  n <- length(covariates)
  combinations <- lapply(1:n, function(x) combn(covariates, x, simplify = FALSE))
  unlist(combinations, recursive = FALSE)
}


list_combinations <- get_covariate_combinations(1:5)
out.tab_photo_diff<- NULL
fits_photo_diff <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  photoperiod_difference = 0
  temperature_difference = 0
  flow = 0
  lunar_phase = 0
  hatchery = 0
  for(j in covariates){
    if(j == 1){
      k = 9
      photoperiod_difference =1
    }
    else if(j==2){
      k = 7
      temperature_difference = 1
    }
    else if(j==3){
      k = 2
      flow = 1
    }
    else if(j==4){
      k = 5
      lunar_phase = 1
    }
    else if(j==5){
      k = 10
      hatchery = 1
    }
    
    c = rbind(c,covariates_all_years_chinook0[((1+(k-1)*15):(k*15)),])
    name_long = rownames(covariates_all_years_chinook0)[1+(k-1)*15]
    name = paste(name, substr(name_long,1,nchar(name_long)-5))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(c_num == 1){
    
    fit.model = c(list(c= c), mod.list_1_0)
  }
  
  else if(c_num == 2){
    
    fit.model = c(list(c= c), mod.list_2_0)
  }
  else if(c_num == 3){
    
    fit.model = c(list(c= c), mod.list_3_0)
  }
  else if(c_num == 4){
    
    fit.model = c(list(c= c), mod.list_4_0)
  }
  else if(c_num == 5){
    
    fit.model = c(list(c= c), mod.list_5_0)
  }
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, photoperiod_difference = photoperiod_difference,
                 temperature_difference = temperature_difference, flow = flow,
                 lunar_phase = lunar_phase, hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_photo_diff=rbind(out.tab_photo_diff,out)
  fits_photo_diff=c(fits_photo_diff,list(fit))
  
}
out.tab_photo_diff


#######

#temp

######

list_combinations <- get_covariate_combinations(1:5)
out.tab_temp<- NULL
fits_temp <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  temperature = 0
  temperature_difference = 0
  flow = 0
  lunar_phase = 0
  hatchery = 0
  for(j in covariates){
    if(j == 1){
      k = 1
      temperature =1
    }
    else if(j==2){
      k = 7
      temperature_difference = 1
    }
    else if(j==3){
      k = 2
      flow = 1
    }
    else if(j==4){
      k = 5
      lunar_phase = 1
    }
    else if(j==5){
      k = 10
      hatchery = 1
    }
    
    c = rbind(c,covariates_all_years_chinook0[((1+(k-1)*15):(k*15)),])
    name_long = rownames(covariates_all_years_chinook0)[1+(k-1)*15]
    name = paste(name, substr(name_long,1,nchar(name_long)-5))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(c_num == 1){
    
    fit.model = c(list(c= c), mod.list_1_0)
  }
  
  else if(c_num == 2){
    
    fit.model = c(list(c= c), mod.list_2_0)
  }
  else if(c_num == 3){
    
    fit.model = c(list(c= c), mod.list_3_0)
  }
  else if(c_num == 4){
    
    fit.model = c(list(c= c), mod.list_4_0)
  }
  else if(c_num == 5){
    
    fit.model = c(list(c= c), mod.list_5_0)
  }
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, temperature = temperature,
                 temperature_difference = temperature_difference, flow = flow,
                 lunar_phase = lunar_phase, hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_temp=rbind(out.tab_temp,out)
  fits_temp=c(fits_temp,list(fit))
  
}
out.tab_temp

#######



#photo

######

list_combinations <- get_covariate_combinations(1:5)
out.tab_photo<- NULL
fits_photo <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  photoperiod= 0
  temperature_difference = 0
  flow = 0
  lunar_phase = 0
  hatchery = 0
  for(j in covariates){
    if(j == 1){
      k = 3
      photoperiod=1
    }
    else if(j==2){
      k = 7
      temperature_difference = 1
    }
    else if(j==3){
      k = 2
      flow = 1
    }
    else if(j==4){
      k = 5
      lunar_phase = 1
    }
    else if(j==5){
      k = 10
      hatchery = 1
    }
    
    c = rbind(c,covariates_all_years_chinook0[((1+(k-1)*15):(k*15)),])
    name_long = rownames(covariates_all_years_chinook0)[1+(k-1)*15]
    name = paste(name, substr(name_long,1,nchar(name_long)-5))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(c_num == 1){
    
    fit.model = c(list(c= c), mod.list_1_0)
  }
  
  else if(c_num == 2){
    
    fit.model = c(list(c= c), mod.list_2_0)
  }
  else if(c_num == 3){
    
    fit.model = c(list(c= c), mod.list_3_0)
  }
  else if(c_num == 4){
    
    fit.model = c(list(c= c), mod.list_4_0)
  }
  else if(c_num == 5){
    
    fit.model = c(list(c= c), mod.list_5_0)
  }
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, photoperiod= photoperiod,
                 temperature_difference = temperature_difference, flow = flow,
                 lunar_phase = lunar_phase, hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_photo=rbind(out.tab_photo,out)
  fits_photo=c(fits_photo,list(fit))
  
}
out.tab_photo

#######


#atu

######

list_combinations <- get_covariate_combinations(1:5)
out.tab_atu<- NULL
fits_atu <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  atu= 0
  temperature_difference = 0
  flow = 0
  lunar_phase = 0
  hatchery = 0
  for(j in covariates){
    if(j == 1){
      k = 4
      atu=1
    }
    else if(j==2){
      k = 7
      temperature_difference = 1
    }
    else if(j==3){
      k = 2
      flow = 1
    }
    else if(j==4){
      k = 5
      lunar_phase = 1
    }
    else if(j==5){
      k = 10
      hatchery = 1
    }
    
    c = rbind(c,covariates_all_years_chinook0[((1+(k-1)*15):(k*15)),])
    name_long = rownames(covariates_all_years_chinook0)[1+(k-1)*15]
    name = paste(name, substr(name_long,1,nchar(name_long)-5))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(c_num == 1){
    
    fit.model = c(list(c= c), mod.list_1_0)
  }
  
  else if(c_num == 2){
    
    fit.model = c(list(c= c), mod.list_2_0)
  }
  else if(c_num == 3){
    
    fit.model = c(list(c= c), mod.list_3_0)
  }
  else if(c_num == 4){
    
    fit.model = c(list(c= c), mod.list_4_0)
  }
  else if(c_num == 5){
    
    fit.model = c(list(c= c), mod.list_5_0)
  }
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, atu= atu,
                 temperature_difference = temperature_difference, flow = flow,
                 lunar_phase = lunar_phase, hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_atu=rbind(out.tab_atu,out)
  fits_atu=c(fits_atu,list(fit))
  
}
out.tab_atu

#######




#resid

######

list_combinations <- get_covariate_combinations(1:5)
out.tab_resid<- NULL
fits_resid <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  resid= 0
  temperature_difference = 0
  flow = 0
  lunar_phase = 0
  hatchery = 0
  for(j in covariates){
    if(j == 1){
      k = 6
      resid=1
    }
    else if(j==2){
      k = 7
      temperature_difference = 1
    }
    else if(j==3){
      k = 2
      flow = 1
    }
    else if(j==4){
      k = 5
      lunar_phase = 1
    }
    else if(j==5){
      k = 10
      hatchery = 1
    }
    
    c = rbind(c,covariates_all_years_chinook0[((1+(k-1)*15):(k*15)),])
    name_long = rownames(covariates_all_years_chinook0)[1+(k-1)*15]
    name = paste(name, substr(name_long,1,nchar(name_long)-5))
    
  }
  # print(c)
  print(name)
  c_num <- length(covariates)
  if(c_num == 1){
    
    fit.model = c(list(c= c), mod.list_1_0)
  }
  
  else if(c_num == 2){
    
    fit.model = c(list(c= c), mod.list_2_0)
  }
  else if(c_num == 3){
    
    fit.model = c(list(c= c), mod.list_3_0)
  }
  else if(c_num == 4){
    
    fit.model = c(list(c= c), mod.list_4_0)
  }
  else if(c_num == 5){
    
    fit.model = c(list(c= c), mod.list_5_0)
  }
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, resid= resid,
                 temperature_difference = temperature_difference, flow = flow,
                 lunar_phase = lunar_phase, hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  out.tab_resid=rbind(out.tab_resid,out)
  fits_resid=c(fits_resid,list(fit))
  
}
out.tab_resid

#######

#resid


#####

weights <- akaike.weights(out.tab_resid$AICc)

out.tab_resid$deltaAICc <- weights$deltaAIC
out.tab_resid$rel.LL <- weights$rel.LL
out.tab_resid$weights <- weights$weights


min.AICc <- order(out.tab_resid$AICc)
out.tab_resid.ordered <- out.tab_resid[min.AICc, ]
out.tab_resid.ordered

out.tab_resid.ordered$cumulative_weights <- cumsum(out.tab_resid.ordered$weights)

relative_importance_resid <- sum(out.tab_resid$weights[out.tab_resid$resid==1])
relative_importance_temperature_difference <- sum(out.tab_resid$weights[out.tab_resid$temperature_difference==1])
relative_importance_flow <- sum(out.tab_resid$weights[out.tab_resid$flow==1])
relative_importance_lunar_phase <- sum(out.tab_resid$weights[out.tab_resid$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab_resid$weights[out.tab_resid$hatchery==1])

riv_residuals <- data.frame(variable = c("temperature residuals",
                                         "temperature difference",
                                         "flow",
                                         "lunar phase",
                                         "hatchery"),
                            relative_importance = c(relative_importance_resid,
                                                    relative_importance_temperature_difference,
                                                    relative_importance_flow,
                                                    relative_importance_lunar_phase,
                                                    relative_importance_hatchery))

######


#temp


#####
weights <- akaike.weights(out.tab_temp$AICc)

out.tab_temp$deltaAICc <- weights$deltaAIC
out.tab_temp$rel.LL <- weights$rel.LL
out.tab_temp$weights <- weights$weights


min.AICc <- order(out.tab_temp$AICc)
out.tab_temp.ordered <- out.tab_temp[min.AICc, ]
out.tab_temp.ordered

out.tab_temp.ordered$cumulative_weights <- cumsum(out.tab_temp.ordered$weights)

relative_importance_temp <- sum(out.tab_temp$weights[out.tab_temp$temp==1])
relative_importance_temperature_difference <- sum(out.tab_temp$weights[out.tab_temp$temperature_difference==1])
relative_importance_flow <- sum(out.tab_temp$weights[out.tab_temp$flow==1])
relative_importance_lunar_phase <- sum(out.tab_temp$weights[out.tab_temp$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab_temp$weights[out.tab_temp$hatchery==1])

riv_temp <- data.frame(variable = c("temperature",
                                         "temperature difference",
                                         "flow",
                                         "lunar phase",
                                         "hatchery"),
                            relative_importance = c(relative_importance_temp,
                                                    relative_importance_temperature_difference,
                                                    relative_importance_flow,
                                                    relative_importance_lunar_phase,
                                                    relative_importance_hatchery))



######

#photoperiod


########


weights <- akaike.weights(out.tab_photo$AICc)

out.tab_photo$deltaAICc <- weights$deltaAIC
out.tab_photo$rel.LL <- weights$rel.LL
out.tab_photo$weights <- weights$weights


min.AICc <- order(out.tab_photo$AICc)
out.tab_photo.ordered <- out.tab_photo[min.AICc, ]
out.tab_photo.ordered

out.tab_photo.ordered$cumulative_weights <- cumsum(out.tab_photo.ordered$weights)

relative_importance_photo <- sum(out.tab_photo$weights[out.tab_photo$photo==1])
relative_importance_temperature_difference <- sum(out.tab_photo$weights[out.tab_photo$temperature_difference==1])
relative_importance_flow <- sum(out.tab_photo$weights[out.tab_photo$flow==1])
relative_importance_lunar_phase <- sum(out.tab_photo$weights[out.tab_photo$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab_photo$weights[out.tab_photo$hatchery==1])

riv_photo <- data.frame(variable = c("photoperiod",
                                    "temperature difference",
                                    "flow",
                                    "lunar phase",
                                    "hatchery"),
                       relative_importance = c(relative_importance_photo,
                                               relative_importance_temperature_difference,
                                               relative_importance_flow,
                                               relative_importance_lunar_phase,
                                               relative_importance_hatchery))



#############



#atu


#####



weights <- akaike.weights(out.tab_atu$AICc)

out.tab_atu$deltaAICc <- weights$deltaAIC
out.tab_atu$rel.LL <- weights$rel.LL
out.tab_atu$weights <- weights$weights


min.AICc <- order(out.tab_atu$AICc)
out.tab_atu.ordered <- out.tab_atu[min.AICc, ]
out.tab_atu.ordered

out.tab_atu.ordered$cumulative_weights <- cumsum(out.tab_atu.ordered$weights)

relative_importance_atu <- sum(out.tab_atu$weights[out.tab_atu$atu==1])
relative_importance_temperature_difference <- sum(out.tab_atu$weights[out.tab_atu$temperature_difference==1])
relative_importance_flow <- sum(out.tab_atu$weights[out.tab_atu$flow==1])
relative_importance_lunar_phase <- sum(out.tab_atu$weights[out.tab_atu$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab_atu$weights[out.tab_atu$hatchery==1])

riv_atu <- data.frame(variable = c("atu",
                                     "temperature difference",
                                     "flow",
                                     "lunar phase",
                                     "hatchery"),
                        relative_importance = c(relative_importance_atu,
                                                relative_importance_temperature_difference,
                                                relative_importance_flow,
                                                relative_importance_lunar_phase,
                                                relative_importance_hatchery))



#######

#photoperiod difference



######


weights <- akaike.weights(out.tab_photo_diff$AICc)

out.tab_photo_diff$deltaAICc <- weights$deltaAIC
out.tab_photo_diff$rel.LL <- weights$rel.LL
out.tab_photo_diff$weights <- weights$weights


min.AICc <- order(out.tab_photo_diff$AICc)
out.tab_photo_diff.ordered <- out.tab_photo_diff[min.AICc, ]
out.tab_photo_diff.ordered

out.tab_photo_diff.ordered$cumulative_weights <- cumsum(out.tab_photo_diff.ordered$weights)

relative_importance_photo_diff <- sum(out.tab_photo_diff$weights[out.tab_photo_diff$photoperiod_difference==1])
relative_importance_temperature_difference <- sum(out.tab_photo_diff$weights[out.tab_photo_diff$temperature_difference==1])
relative_importance_flow <- sum(out.tab_photo_diff$weights[out.tab_photo_diff$flow==1])
relative_importance_lunar_phase <- sum(out.tab_photo_diff$weights[out.tab_photo_diff$lunar_phase==1])
relative_importance_hatchery <- sum(out.tab_photo_diff$weights[out.tab_photo_diff$hatchery==1])

riv_photo_diff <- data.frame(variable = c("photoperiod difference",
                                   "temperature difference",
                                   "flow",
                                   "lunar phase",
                                   "hatchery"),
                      relative_importance = c(relative_importance_photo_diff,
                                              relative_importance_temperature_difference,
                                              relative_importance_flow,
                                              relative_importance_lunar_phase,
                                              relative_importance_hatchery))







######

#adding model with no variables

name = "None"
resid= 0
temperature_difference = 0
flow = 0
lunar_phase = 0
hatchery = 0


fit <- MARSS(subset_chinook_summer_perhour, model=mod.list_0_0, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))


out=data.frame(c=name, resid= resid,
               temperature_difference = temperature_difference, flow = flow,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)

out.tab_resid$deltaAICc <- NULL
out.tab_resid$rel.LL <- NULL
out.tab_resid$weights <- NULL

out.tab_resid=rbind(out.tab_resid,out)
fits_resid=c(fits_resid,list(fit))


name = "None"
resid= 0
temperature_difference = 0
flow = 0
lunar_phase = 0
hatchery = 0
temperature = 0

out=data.frame(c=name, temperature= temperature,
               temperature_difference = temperature_difference, flow = flow,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)

out.tab_temp$deltaAICc <- NULL
out.tab_temp$rel.LL <- NULL
out.tab_temp$weights <- NULL


out.tab_temp=rbind(out.tab_temp,out)
fits_temp=c(fits_temp,list(fit))


atu = 0
out=data.frame(c=name, atu=atu,
               temperature_difference = temperature_difference, flow = flow,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)

out.tab_atu$deltaAICc <- NULL
out.tab_atu$rel.LL <- NULL
out.tab_atu$weights <- NULL


out.tab_atu=rbind(out.tab_atu,out)
fits_atu=c(fits_atu,list(fit))



photoperiod = 0
out=data.frame(c=name, photoperiod = photoperiod,
               temperature_difference = temperature_difference, flow = flow,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)

out.tab_photo$deltaAICc <- NULL
out.tab_photo$rel.LL <- NULL
out.tab_photo$weights <- NULL




out.tab_photo=rbind(out.tab_photo,out)
fits_photo=c(fits_photo,list(fit))



photoperiod_difference = 0
out=data.frame(c=name, photoperiod_difference = photoperiod_difference,
               temperature_difference = temperature_difference, flow = flow,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)

out.tab_photo_diff$deltaAICc <- NULL
out.tab_photo_diff$rel.LL <- NULL
out.tab_photo_diff$weights <- NULL

out.tab_photo_diff=rbind(out.tab_photo_diff,out)
fits_photo_diff=c(fits_photo_diff,list(fit))


write.csv(out.tab_temp.ordered, here("output","chinook0_model_selection_temp.csv"))
write.csv(out.tab_resid.ordered, here("output","chinook0_model_selection_resid.csv"))
write.csv(out.tab_atu.ordered, here("output","chinook0_model_selection_atu.csv"))
write.csv(out.tab_photo.ordered, here("output","chinook0_model_selection_photo.csv"))
write.csv(out.tab_photo_diff.ordered, here("output","chinook0_model_selection_photo_diff.csv"))


riv_residuals = riv_residuals[order(riv_residuals$relative_importance, decreasing = TRUE), ]

write.csv(riv_residuals, here("output","riv_resid.csv"))



riv_temp = riv_temp[order(riv_temp$relative_importance, decreasing = TRUE), ]

write.csv(riv_temp, here("output","riv_temp.csv"))



riv_atu = riv_atu[order(riv_atu$relative_importance, decreasing = TRUE), ]

write.csv(riv_atu, here("output","riv_atu.csv"))


riv_photo = riv_photo[order(riv_photo$relative_importance, decreasing = TRUE), ]

write.csv(riv_photo, here("output","riv_photo.csv"))



riv_photo_diff = riv_photo_diff[order(riv_photo_diff$relative_importance, decreasing = TRUE), ]

write.csv(riv_photo_diff, here("output","riv_photo_diff.csv"))


out.tab_photo_diff.ordered
ci_best_photo_diff <- tidy(fits_photo_diff[[18]])

ggplot(ci_best_photo_diff[c(18,19,20),], 
       aes(x = c("Photoperiod\n difference", "Temperature\n difference", "Hatchery"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  ggtitle("Chinook subyearlings") + 
  theme(plot.title = element_text(size = 24))+
  theme(axis.text.x=element_text(size=24),axis.title.y=element_text(size=24))

ggsave(here("output","chinook0_effect_final_photo_diff.jpeg"), width = 8, height = 8)
autoplot(fits_photo_diff[[18]], plot.type = "fitted.ytT")

ggsave(here("output","fitted_y_chinook_best_model.jpeg"), width = 8, height = 8)


