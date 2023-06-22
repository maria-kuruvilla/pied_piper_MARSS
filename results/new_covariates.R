#goal - make new covariates

# take derivative of flow, temperature and day of year


library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)
library(MASS)
library(modelr)

data <- read.csv(here("data","dungeness_subset_covariates.csv"), header = TRUE)

data$chinook0_hatchery_perhour_interpolate <- na.approx(data$chinook0_hatchery_perhour)

lm_temp <- lm(temp ~ 1+photoperiod, data)

data <- data %>% add_residuals(lm_temp)

data$temp_diff <- c(NA, diff(data$temp))

data$flow_diff <- c(NA, diff(data$flow))

data$photo_diff <- c(NA, diff(data$photoperiod))




#save csv file

write.csv(data,here("data","dungeness_subset_covariates_diff.csv"))

write.csv(data,here("..","..","..","OneDrive","Documents","data","pied_piper",
                    "dungeness","dungeness_subset_covariates_diff.csv"))

