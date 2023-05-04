#packages
library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)


#readng the data

data <- read.csv(here("data","dungeness_all_days.csv"), header = TRUE)



#filtering only 2005-2009 data for coho wild

subset <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2010) %>%
  select(coho1_wild_num,year,doy) %>%
  pivot_wider(names_from = year, values_from = coho1_wild_num) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

mod.list <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "diagonal and unequal"
  )

fit <- MARSS(subset, model = mod.list, method = "BFGS")

plot(fit, plot.type = "model.resids.ytt1") # not good - try log
plot(fit, plot.type = "qqplot.std.model.resids.ytt1") # not good 
plot(fit, plot.type = "acf.std.model.resids.ytt1") # not good


#filtering on doy and removing zeros and logging the data
subset2 <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2010 & doy > 100 & doy <=200) %>%
  select(coho1_wild_num,year,doy) %>%
  pivot_wider(names_from = year, values_from = coho1_wild_num) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


subset3 <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2010 & doy > 100 & doy <=200) %>%
  mutate(log.value = ifelse(coho1_wild_num == 0, NA, log(coho1_wild_num))) %>%
  select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()



mod.list <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "diagonal and unequal"
)

fit3 <- MARSS(subset3, model = mod.list, method = "BFGS")

autoplot(fit3)


#adding covariates
#interplating temp values - need to change


subset_covariates <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2010 & doy > 100 & doy <=200) %>%
  mutate(temp =na.approx(temp, na.rm = FALSE, maxgap = 5)) %>%
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


fit_covariates <- MARSS(subset3, model = mod.list.cov, method = "BFGS")

autoplot(fit_covariates)