# try marss model with covariates

#packages
library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)

#readng the data

data <- read.csv(here("data","dungeness_subset_covariates.csv"), header = TRUE)

#subset
subset <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2010 & doy > 100 & doy <=200) %>%
  mutate(log.value = ifelse(coho1_wild_num == 0, NA, log(coho1_wild_num))) %>%
  select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


subset_covariates <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2010 & doy > 100 & doy <=200) %>%
  select(temp,year,doy) %>%
  pivot_wider(names_from = year, values_from = temp) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


mod.list.cov <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = "diagonal and unequal",
  c = subset_covariates
)



fit_covariates <- MARSS(subset, model = mod.list.cov, method = "BFGS")

autoplot(fit_covariates)

AIC(fit_covariates)
#AIC with temp  = 1025


#without temp

mod.list <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "diagonal and unequal"
)


fit <- MARSS(subset, model = mod.list, method = "BFGS")
AIC(fit) #1028


#hatchery

data$coho1_hatchery_num_interpolate <- na.approx(data$coho1_hatchery_num)

subset_covariates_hatchery <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2010 & doy > 100 & doy <=200) %>%
  select(coho1_hatchery_num_interpolate,year,doy) %>%
  pivot_wider(names_from = year, values_from = coho1_hatchery_num_interpolate) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

mod.list.hatchery <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = "diagonal and unequal",
  c = subset_covariates_hatchery
)



fit_covariates_hatchery <- MARSS(subset, model = mod.list.hatchery, method = "BFGS")

autoplot(fit_covariates_hatchery)

AIC(fit_covariates_hatchery) # 1032


#temp and hatchery as covariates

covariates_hatchery_temp <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2010 & doy > 100 & doy <=200) %>%
  select(coho1_hatchery_num_interpolate,year,doy, temp) %>%
  pivot_wider(names_from = year, values_from = c(coho1_hatchery_num_interpolate,temp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


mod.list.hatchery_temp <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = matrix(list("h2005",0,0,0,0,"t2005",0,0,0,0,
                  0,"h2006",0,0,0,0,"t2006",0,0,0,
                  0,0,"h2007",0,0,0,0,"t2007",0,0,
                  0,0,0,"h2008",0,0,0,0,"t2008",0,
                  0,0,0,0,"h2009",0,0,0,0,"t2009"),5,10, byrow = TRUE),
  c = covariates_hatchery_temp
)

fit_temp_hatchery <- MARSS(subset, model = mod.list.hatchery_temp, method = "BFGS")

autoplot(fit_temp_hatchery)

AIC(fit_temp_hatchery) #1029



#try for all years
coho <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 100 & doy <=200) %>%
  mutate(log.value = ifelse(coho1_wild_num == 0, NA, log(coho1_wild_num))) %>%
  select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


covariates_coho_hatchery_temp_all_years <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 100 & doy <=200) %>%
  select(coho1_hatchery_num_interpolate,year,doy, temp) %>%
  pivot_wider(names_from = year, values_from = c(coho1_hatchery_num_interpolate,temp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

mod.list.hatchery_temp_all_years <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = matrix(list("h2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"h2006",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"h2007",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"h2008",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"h2009",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"h2010",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"h2011",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"h2012",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"h2013",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"h2014",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"h2016",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"h2017",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"h2018",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"h2019",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2020",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2020"),
             15,30, byrow = TRUE),
  c = covariates_coho_hatchery_temp_all_years
)

fit_temp_hatchery_all_years <- MARSS(coho, model = mod.list.hatchery_temp_all_years, method = "BFGS")

autoplot(fit_temp_hatchery_all_years)

AIC(fit_temp_hatchery_all_years) #2623


#just temp

covariates_temp_all_years <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 100 & doy <=200) %>%
  select(year,doy, temp) %>%
  pivot_wider(names_from = year, values_from = temp) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

mod.list.temp_all_years <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = matrix(list("t2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"t2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"t2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"t2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"t2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"t2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"t2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"t2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"t2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"t2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"t2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"t2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"t2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"t2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"t2020"),
             15,15, byrow = TRUE),
  c = covariates_temp_all_years
)

fit_temp_all_years <- MARSS(coho, model = mod.list.temp_all_years, method = "BFGS")

autoplot(fit_temp_all_years)

AIC(fit_temp_all_years) #2608



#############


#modlist = list()
#fit = MARSS(,,fit = FALSE)
#summary(fit$model)
#toLatex(fit$model)

#if we demean all the data (y,c, and d), we do not need any intercept
#first log then demean

#do not usually standardize the means
# the q matrix is usually accounting for the variance

##########

#subset demean
subset <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2010 & doy > 100 & doy <=200) %>%
  mutate(log.value = ifelse(coho1_wild_num == 0, NA, log(coho1_wild_num))) %>%
  select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset)[1]){
  subset[i,] = subset[i,] - mean(subset[i,], na.rm = TRUE)
}

covariates_temp_zscore <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2010 & doy > 100 & doy <=200) %>%
  mutate(temp_z = zscore(temp)) %>%
  select(year,doy, temp_z) %>%
  pivot_wider(names_from = year, values_from = temp_z) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

mod.list <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal"
)

fit <- MARSS(subset, model = mod.list, method = "BFGS")

autoplot(fit)

AIC(fit)

mod.list_temp <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = "diagonal and unequal",
  c = covariates_temp_zscore
)

fit_temp <- MARSS(subset, model = mod.list_temp, method = "BFGS")

autoplot(fit_temp)

AIC(fit_temp)


covariates_temp_zscore_hatchery <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2010 & doy > 100 & doy <=200) %>%
  mutate(temp_z = zscore(temp)) %>%
  select(year,doy, temp_z, coho1_hatchery_num_interpolate) %>%
  pivot_wider(names_from = year, values_from = c(temp_z,coho1_hatchery_num_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

mod.list_temp_hatchery <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C =  matrix(list("h2005",0,0,0,0,"t2005",0,0,0,0,
                      0,"h2006",0,0,0,0,"t2006",0,0,0,
                      0,0,"h2007",0,0,0,0,"t2007",0,0,
                      0,0,0,"h2008",0,0,0,0,"t2008",0,
                      0,0,0,0,"h2009",0,0,0,0,"t2009"
                      ),
                 5,10, byrow = TRUE),
  c = covariates_temp_zscore_hatchery
)
fit_temp_hatchery <- MARSS(subset, model = mod.list_temp_hatchery, method = "BFGS")

autoplot(fit_temp_hatchery)

AIC(fit_temp_hatchery)


