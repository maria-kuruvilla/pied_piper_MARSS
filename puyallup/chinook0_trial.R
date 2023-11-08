# Goal 

# try MARSS analysis with Puyallup data fo chinook0

#load packages

library(here)
library(MARSS)
library(ggplot2)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)
library(tidyverse)


# load data

data <- read.csv(here("puyallup", "data","puyallup_2004-2018_all_days_w_covariates.csv"))
