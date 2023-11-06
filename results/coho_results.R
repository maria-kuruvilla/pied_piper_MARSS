


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

data <- read.csv(here("data","dungeness_subset_covariates_diff.csv"), header = TRUE)
data$coho1_hatchery_perhour_interpolate <- na.approx(data$coho1_hatchery_perhour)

#subset
subset_coho_perhour <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 120 & doy <= 160) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


for(i in 1:dim(subset_coho_perhour)[1]){
  subset_coho_perhour[i,] = scale(subset_coho_perhour[i,])[,1]
}


covariates_all_years_coho <- arrange(data,doy) %>%
  filter(year != 2015 & doy >120 & doy <= 160) %>%
  dplyr::select(year,doy, temp, flow, photoperiod, atu_april, 
                lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
                coho1_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = year, values_from = c(
    temp, flow, photoperiod, atu_april, lunar_phase, 
    resid, temp_diff, flow_diff, photo_diff, coho1_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


for(i in 1:(dim(covariates_all_years_coho)[1]-75)){ # everything except resid, diffs and hatchery
  covariates_all_years_coho[i,] = scale(covariates_all_years_coho[i,])[,1]
}

#scaling but not centering the rest

for(i in 76:135){
  covariates_all_years_coho[i,] = scale(covariates_all_years_coho[i,], center = FALSE, scale= TRUE)[,1]
}

for(i in 136:150){
  covariates_all_years_coho[i,] = scale(covariates_all_years_coho[i,], center = FALSE, scale= TRUE)[,1]
}




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




list_combinations <- get_covariate_combinations(1:5)
out.tab_coho <- NULL
fits_coho <- list()
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
        
        c = rbind(c,covariates_all_years_coho[((1+(k-1)*15):(k*15)),])
        name_long = rownames(covariates_all_years_coho)[1+(k-1)*15]
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
      
      fit <- MARSS(subset_coho_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                   control=list(maxit=2000))
      
      
      out=data.frame(c=name, d = "None",
                     logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                     num.iter=fit$numIter, converged=!fit$convergence,
                     stringsAsFactors = FALSE)
      out.tab_coho=rbind(out.tab_coho,out)
      fits_coho=c(fits_coho,list(fit))
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
      
      c = rbind(c,covariates_all_years_coho[((1+(k-1)*15):(k*15)),])
      name_long = rownames(covariates_all_years_coho)[1+(k-1)*15]
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
    
    fit <- MARSS(subset_coho_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                 control=list(maxit=2000))
    
    
    out=data.frame(c=name, d = "None",
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence,
                   stringsAsFactors = FALSE)
    out.tab_coho=rbind(out.tab_coho,out)
    fits_coho=c(fits_coho,list(fit))
    
  }
  
  
}
out.tab_coho$deltaAICc <- out.tab_coho$AICc - min(out.tab_coho$AICc)
min.AICc <- order(out.tab_coho$AICc)
out.tab_coho.ordered <- out.tab_coho[min.AICc, ]
out.tab_coho.ordered
write.csv(out.tab_coho.ordered, here("output","coho_model_selection.csv"))


#separate set of models to run model selection on 


#photoperiod difference

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
    
    c = rbind(c,covariates_all_years_coho[((1+(k-1)*15):(k*15)),])
    name_long = rownames(covariates_all_years_coho)[1+(k-1)*15]
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
  
  fit <- MARSS(subset_coho_perhour, model=fit.model, silent = TRUE, method = "BFGS",
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


photoperiod_difference = 0
name = "None"
resid= 0
temperature_difference = 0
flow = 0
lunar_phase = 0
hatchery = 0


fit <- MARSS(subset_coho_perhour, model=mod.list_0_0, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

out=data.frame(c=name, photoperiod_difference = photoperiod_difference,
               temperature_difference = temperature_difference, flow = flow,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)

out.tab_photo_diff=rbind(out.tab_photo_diff,out)
fits_photo_diff=c(fits_photo_diff,list(fit))


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


riv_photo_diff = riv_photo_diff[order(riv_photo_diff$relative_importance, decreasing = TRUE), ]




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
    
    c = rbind(c,covariates_all_years_coho[((1+(k-1)*15):(k*15)),])
    name_long = rownames(covariates_all_years_coho)[1+(k-1)*15]
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
  
  fit <- MARSS(subset_coho_perhour, model=fit.model, silent = TRUE, method = "BFGS",
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
    
    c = rbind(c,covariates_all_years_coho[((1+(k-1)*15):(k*15)),])
    name_long = rownames(covariates_all_years_coho)[1+(k-1)*15]
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
  
  fit <- MARSS(subset_coho_perhour, model=fit.model, silent = TRUE, method = "BFGS",
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
    
    c = rbind(c,covariates_all_years_coho[((1+(k-1)*15):(k*15)),])
    name_long = rownames(covariates_all_years_coho)[1+(k-1)*15]
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
  
  fit <- MARSS(subset_coho_perhour, model=fit.model, silent = TRUE, method = "BFGS",
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
    
    c = rbind(c,covariates_all_years_coho[((1+(k-1)*15):(k*15)),])
    name_long = rownames(covariates_all_years_coho)[1+(k-1)*15]
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
  
  fit <- MARSS(subset_coho_perhour, model=fit.model, silent = TRUE, method = "BFGS",
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


######

#adding model with no variables


fit <- MARSS(subset_chinook_summer_perhour, model=mod.list_0_0, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))


name = "None"
resid= 0
temperature_difference = 0
flow = 0
lunar_phase = 0
hatchery = 0



out=data.frame(c=name, resid= resid,
               temperature_difference = temperature_difference, flow = flow,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)


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



out.tab_temp=rbind(out.tab_temp,out)
fits_temp=c(fits_temp,list(fit))


atu = 0
out=data.frame(c=name, atu=atu,
               temperature_difference = temperature_difference, flow = flow,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)



out.tab_atu=rbind(out.tab_atu,out)
fits_atu=c(fits_atu,list(fit))



photoperiod = 0
out=data.frame(c=name, photoperiod = photoperiod,
               temperature_difference = temperature_difference, flow = flow,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)





out.tab_photo=rbind(out.tab_photo,out)
fits_photo=c(fits_photo,list(fit))


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


riv_residuals = riv_residuals[order(riv_residuals$relative_importance, decreasing = TRUE), ]


#temperature
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

riv_temp = riv_temp[order(riv_temp$relative_importance, decreasing = TRUE), ]

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
riv_photo = riv_photo[order(riv_photo$relative_importance, decreasing = TRUE), ]



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

riv_atu = riv_atu[order(riv_atu$relative_importance, decreasing = TRUE), ]





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


riv_photo_diff = riv_photo_diff[order(riv_photo_diff$relative_importance, decreasing = TRUE), ]


write.csv(riv_temp, here("output","riv_temp_coho.csv"))
write.csv(riv_photo, here("output","riv_photo_coho.csv"))
write.csv(riv_atu, here("output","riv_atu_coho.csv"))
write.csv(riv_residuals, here("output","riv_resid_coho.csv"))
write.csv(riv_photo_diff, here("output","riv_photo_diff_coho.csv"))

write.csv(out.tab_temp.ordered, here("output","coho_model_selection_temp.csv"))
write.csv(out.tab_resid.ordered, here("output","coho_model_selection_resid.csv"))
write.csv(out.tab_atu.ordered, here("output","coho_model_selection_atu.csv"))
write.csv(out.tab_photo.ordered, here("output","coho_model_selection_photo.csv"))
write.csv(out.tab_photo_diff.ordered, here("output","coho_model_selection_photo_diff.csv"))


out_together <- bind_rows(out.tab_temp, out.tab_resid, out.tab_atu, out.tab_photo, out.tab_photo_diff)

out_together$total_delta_AICc <- out_together$AICc - min(out_together$AICc)

out_together.ordered <- out_together[order(out_together$total_delta_AICc),]

ci_best_photo <- tidy(fits_photo[[7]])

ci_w_hatchery <- tidy(fits_photo[[20]])

ggplot(ci_best_photo[c(18,19),], 
       aes(x = c("Photoperiod", "Flow"),
           y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  ggtitle("Coho yearlings") + 
  theme(plot.title = element_text(size = 24))+
  theme(axis.text.x=element_text(size=24),axis.title.y=element_text(size=24))

ggsave(here("output","coho_effect_final_photo.jpeg"), width = 8, height = 8)
autoplot(fits_photo[[7]], plot.type = "fitted.ytT")

ggsave(here("output","fitted_y_coho_best_model.jpeg"), width = 8, height = 8)




