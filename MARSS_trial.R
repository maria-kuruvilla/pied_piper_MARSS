#packages
library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)


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
