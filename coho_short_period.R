# Try MARSS for coho with a smaller restriction on period to see if 2012 and 2017 has bettter results

#packages
library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)


#readng the data

data <- read.csv(here("data","dungeness_subset_covariates.csv"), header = TRUE)

data$coho1_hatchery_num_interpolate <- na.approx(data$coho1_hatchery_num)


ggplot() + geom_line(aes(data$doy, data$coho1_wild_num), data = data) + facet_wrap(data$year,nrow = 4, ncol = 4)

#subset
subset <- arrange(data,doy) %>%
  filter(year != 2015 & doy > 100 & doy <= 185) %>%
  mutate(log.value = ifelse(coho1_wild_num == 0, NA, log(coho1_wild_num))) %>%
  select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

#scale the values

for(i in 1:dim(subset)[1]){
  subset[i,] = scale(subset[i,])[,1]
}


mod.list <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal"
)

fit <- MARSS(subset, model = mod.list, method = "BFGS")

autoplot(fit)

AIC(fit) #2084

covariates_temp_all_years <- arrange(data,doy) %>%
  filter(year != 2015 & doy >100 & doy <= 185) %>%
  select(year,doy, temp) %>%
  pivot_wider(names_from = year, values_from = temp) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(covariates_temp_all_years)[1]){
  covariates_temp_all_years[i,] = scale(covariates_temp_all_years[i,])[,1]
}


mod.list.temp_all_years <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = "diagonal and unequal",
  c = covariates_temp_all_years
)

fit_temp_all_years <- MARSS(subset, model = mod.list.temp_all_years, method = "BFGS")

autoplot(fit_temp_all_years)

AIC(fit_temp_all_years) #2086
#model is worse with temp


#let's try atu

covariates_atu_all_years <- arrange(data,doy) %>%
  filter(year != 2015 & doy >100 & doy <= 185) %>%
  select(year,doy, atu_april) %>%
  pivot_wider(names_from = year, values_from = atu_april) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(covariates_atu_all_years)[1]){
  covariates_atu_all_years[i,] = scale(covariates_atu_all_years[i,])[,1]
}

mod.list.atu_all_years <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = "diagonal and unequal",
  c = covariates_atu_all_years
)

fit_atu_all_years <- MARSS(subset, model = mod.list.atu_all_years, method = "BFGS")

autoplot(fit_atu_all_years)

AIC(fit_atu_all_years) #2067
#model is much better with atu

#let's try photoperiod

covariates_photo_all_years <- arrange(data,doy) %>%
  filter(year != 2015 & doy >100 & doy <= 185) %>%
  select(year,doy, photoperiod) %>%
  pivot_wider(names_from = year, values_from = photoperiod) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(covariates_photo_all_years)[1]){
  covariates_photo_all_years[i,] = scale(covariates_photo_all_years[i,])[,1]
}

mod.list.photo_all_years <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = "diagonal and unequal",
  c = covariates_photo_all_years
)

fit_photo_all_years <- MARSS(subset, model = mod.list.photo_all_years, method = "BFGS")

autoplot(fit_photo_all_years)

AIC(fit_photo_all_years) #2057
#model is much better with photoperiod than with atu and temp

lm_photo_temp <- lm(temp~photoperiod, data)
summary(lm_photo_temp)
plot(residuals(lm_photo_temp), fitted(lm_photo_temp))
#residuals look really bad so maybe try MARSS model 



covariates_photo_flow_all_years <- arrange(data,doy) %>%
  filter(year != 2015 & doy >100 & doy <= 185) %>%
  select(year,doy, photoperiod, flow) %>%
  pivot_wider(names_from = year, values_from = c(photoperiod,flow)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(covariates_photo_flow_all_years)[1]){
  covariates_photo_flow_all_years[i,] = scale(covariates_photo_flow_all_years[i,])[,1]
}

mod.list.photo_flow_all_years <- list(
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
  c = covariates_photo_flow_all_years
)

fit_photo_flow_all_years <- MARSS(subset, model = mod.list.photo_flow_all_years, method = "BFGS")

autoplot(fit_photo_flow_all_years)

AIC(fit_photo_flow_all_years) #2057
#adding flow does not make the model better


#let's add flow to d covariate


covariates_flow_all_years <- arrange(data,doy) %>%
  filter(year != 2015 & doy >100 & doy <= 185) %>%
  select(year,doy,  flow) %>%
  pivot_wider(names_from = year, values_from = flow) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(covariates_flow_all_years)[1]){
  covariates_flow_all_years[i,] = scale(covariates_flow_all_years[i,])[,1]
}

mod.list.photo_flow_all_years_d <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = "diagonal and unequal",
  c = covariates_photo_all_years,
  D = "diagonal and unequal",
  d = covariates_flow_all_years
)

fit_photo_flow_all_years_d <- MARSS(subset, model = mod.list.photo_flow_all_years_d, method = "BFGS")

autoplot(fit_photo_flow_all_years_d)

AIC(fit_photo_flow_all_years_d) #2040
#adding flow to observation covariates makes it better


#adding flow to both d and c

mod.list.photo_flow_all_years_c_d <- list(
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
  c = covariates_photo_flow_all_years,
  D = "diagonal and unequal",
  d = covariates_flow_all_years
)
fit_photo_flow_all_years_c_d <- MARSS(subset, model = mod.list.photo_flow_all_years_c_d, method = "BFGS")

autoplot(fit_photo_flow_all_years_c_d)

AIC(fit_photo_flow_all_years_c_d) #2035
#adding flow to observation and process covariates makes it better but some x residuals look weird




#finally adding hatchery


covariates_photo_flow_hatchery_all_years <- arrange(data,doy) %>%
  filter(year != 2015 & doy >100 & doy <= 185) %>%
  select(year,doy,  photoperiod, flow, coho1_hatchery_num_interpolate) %>%
  pivot_wider(names_from = year, values_from = c(photoperiod, flow, coho1_hatchery_num_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:30){ #only for photoperiod and flow
  covariates_photo_flow_hatchery_all_years[i,] = scale(covariates_photo_flow_hatchery_all_years[i,])[,1]
}


mod.list.photo_flow_hatchery_all_years_c_d <- list(
  U = "zero",
  R = "diagonal and equal",
  Q = "diagonal and unequal",
  C = matrix(list("p2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2005",0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,"p2006",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2006",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2006",0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,"p2007",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2007",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2007",0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,"p2008",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2008",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2008",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,"p2009",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2009",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2009",0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"p2010",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2010",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2010",0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,"p2011",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2011",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2011",0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,"p2012",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2012",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2012",0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,"p2013",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2013",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2013",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,"p2014",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2014",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2014",0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"p2016",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2016",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2016",0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,"p2017",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2017",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2017",0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,"p2018",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2018",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2018",0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,"p2019",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2019",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2019",0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,"p2020",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"f2020",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"h2020"),
             15,45, byrow = TRUE),
  c = covariates_photo_flow_hatchery_all_years,
  D = "diagonal and unequal",
  d = covariates_flow_all_years
)
fit_photo_flow_hathchery_all_years_c_d <- MARSS(subset, model = mod.list.photo_flow_hatchery_all_years_c_d, method = "BFGS")

fit_photo_flow_hatchery_all_years_c_d <- fit_photo_flow_hathchery_all_years_c_d 
autoplot(fit_photo_flow_hatchery_all_years_c_d)

AIC(fit_photo_flow_hatchery_all_years_c_d) #2036
#adding hatchery did not make it better


#put it all into table

results <- matrix(NA, nrow = 9, ncol = 3)

results[1,1] <- "C"
results[1,2] <- "D"
results[1,3] <- "AIC"

results[2,1] <- "temperature"
results[3,1] <- "atu"
results[4,1] <- "photoperiod"
results[5,1] <- "photoperiod, flow"
results[6,1] <- "photoperiod"
results[6,2] <- "flow"
results[7,1] <- "photoperiod, flow"
results[7,2] <- "flow"
results[8,1] <- "photoperiod, flow, hatchery"
results[8,2] <- "flow"


results[2,3] <- AIC(fit_temp_all_years)
results[3,3] <- AIC(fit_atu_all_years)
results[4,3] <- AIC(fit_photo_all_years)
results[5,3] <- AIC(fit_photo_flow_all_years)
results[6,3] <- AIC(fit_photo_flow_all_years_d)
results[7,3] <- AIC(fit_photo_flow_all_years_c_d)
results[8,3] <- AIC(fit_photo_flow_hatchery_all_years_c_d)
results[9,3] <- AIC(fit)

results
