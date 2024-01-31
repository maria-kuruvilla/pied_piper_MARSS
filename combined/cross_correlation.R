# Goal - to plot the cross correlation between hatchery fish and wild fish
# for all rivers and all species

# Load libraries
library(ggplot2)
library(tidyverse)
library(forecast)
library(here)
library(zoo)
library(cowplot)

# Load data
#load data
data <- read.csv(here("puyallup", "data","puyallup_2004-2021_all_days_w_covariates_w_temp.csv"),
                 na.strings = c("NA",""))

data$Date = as.Date(data$Date, format = "%Y-%m-%d")
data$year = year(data$Date)

data$chinook0_hatchery_perhour[dim(data)[1]] <- 0
data$chinook0_wild_perhour[dim(data)[1]] <- 0

ccf_puyallup_data <- data %>% 
  dplyr::select(doy,year, chinook0_wild_perhour, chinook0_hatchery_perhour) %>%
  mutate(chinook0_wild_perhour_inp = na.approx(chinook0_wild_perhour, na.rm = FALSE),
         chinook0_hatchery_perhour_inp = na.approx(chinook0_hatchery_perhour, na.rm = FALSE))

puyallup_chinook <- ccf(ccf_data$chinook0_wild_perhour_inp, ccf_data$chinook0_hatchery_perhour_inp,
    lag.max = 10, plot = TRUE, ylab = "cross-correlation",
    main = "Chinook, Puyallup")


ccf_puyallup_data_hatchery <- arrange(data,doy) %>% 
  dplyr::select(doy,year,chinook0_hatchery_perhour) %>%
  filter(doy>100 & doy<=200) %>%
  mutate(chinook0_hatchery_perhour_inp = na.approx(chinook0_hatchery_perhour, na.rm = FALSE), chinook0_hatchery_perhour = NULL) %>%
  pivot_wider(names_from = c(year), values_from = c(chinook0_hatchery_perhour_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

ccf_puyallup_data_wild<- arrange(data,doy) %>% 
  dplyr::select(doy,year,chinook0_wild_perhour) %>%
  filter(doy>100 & doy<=200) %>%
  mutate(chinook0_wild_perhour_inp = na.approx(chinook0_wild_perhour, na.rm = FALSE), chinook0_wild_perhour = NULL) %>%
  pivot_wider(names_from = c(year), values_from = c(chinook0_wild_perhour_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()
ccf_puyallup_data_hatchery <- ts(ccf_puyallup_data_hatchery, start = 100, end = 200)
ccf_puyallup_data_wild <- ts(ccf_puyallup_data_wild, frequency = 100)
ggCcf(ccf_puyallup_data_hatchery,
      ccf_puyallup_data_wild)
      
ggCcf(ccf_puyallup_data_hatchery,
        ccf_puyallup_data_wild, lag.max = 10, plot = TRUE)+
  theme_classic()

#dungeness

# Load data
data_dungeness <- read.csv(here("data","dungeness_aggregated_w_trap_efficiency.csv"),
                 na.strings = c("NA",""))

data_dungeness$Date = as.Date(data_dungeness$Date, format = "%Y-%m-%d")
unique(data_dungeness$year)

ccf_dungeness_data <- data_dungeness %>% 
  dplyr::select(doy,year, chinook0_wild_perhour, chinook0_hatchery_perhour) %>%
  mutate(chinook0_wild_perhour_inp = na.approx(chinook0_wild_perhour, na.rm = FALSE),
         chinook0_hatchery_perhour_inp = na.approx(chinook0_hatchery_perhour, na.rm = FALSE))

dungeness_chinook <- ccf(ccf_dungeness_data$chinook0_wild_perhour_inp,
    ccf_dungeness_data$chinook0_hatchery_perhour_inp,
    lag.max = 10, plot = TRUE, 
    ylab = "cross-correlation", main = "Chinook, Dungeness")

ggplot(ccf_dungeness_data, aes(x = chinook0_wild_perhour_inp, y = chinook0_hatchery_perhour_inp))+
  geom_point()+
  theme_classic()+
  labs(x = "Wild", y = "Hatchery", title = "Chinook, Dungeness")

ccf_dungeness_coho_data <- data_dungeness %>% 
  dplyr::select(doy,year, coho1_wild_perhour, coho1_hatchery_perhour) %>%
  mutate(coho1_wild_perhour_inp = na.approx(coho1_wild_perhour, na.rm = FALSE),
         coho1_hatchery_perhour_inp = na.approx(coho1_hatchery_perhour, na.rm = FALSE))

dungeness_coho <- ccf(ccf_dungeness_coho_data$coho1_wild_perhour_inp,
    ccf_dungeness_coho_data$coho1_hatchery_perhour_inp,
    lag.max = 10, plot = TRUE, 
    ylab = "cross-correlation", main = "Coho, Dungeness")

#skagit
ccf_skagit_chinook_data <- data_day_night %>% 
  filter(daytime_category == 'night') %>%
  mutate(chinook0_wild_perhour = chinook0_wild_num/In, chinook0_hatchery_perhour = chinook0_hatchery_num/In) %>%
  dplyr::select(doy,year, chinook0_wild_perhour, chinook0_hatchery_perhour) %>%
  mutate(chinook0_wild_perhour_inp = na.approx(chinook0_wild_perhour, na.rm = FALSE),
         chinook0_hatchery_perhour_inp = na.approx(chinook0_hatchery_perhour, na.rm = FALSE)) %>% 
  filter(doy > 25 & doy < 189)

skagit_chinook <- ccf(ccf_skagit_chinook_data$chinook0_wild_perhour_inp,
    ccf_skagit_chinook_data$chinook0_hatchery_perhour_inp,
    lag.max = 10, plot = TRUE, 
    ylab = "cross-correlation", main = "Chinook, Skagit")

ccf_skagit_coho_data <- data_day_night %>% 
  filter(daytime_category == 'night') %>%
  mutate(coho1_wild_perhour = coho1_wild_num/In, coho1_hatchery_perhour = coho1_hatchery_num/In) %>%
  dplyr::select(doy,year, coho1_wild_perhour, coho1_hatchery_perhour) %>%
  mutate(coho1_wild_perhour_inp = na.approx(coho1_wild_perhour, na.rm = FALSE),
         coho1_hatchery_perhour_inp = na.approx(coho1_hatchery_perhour, na.rm = FALSE)) %>% 
  filter(doy > 25 & doy < 189)

skagit_coho <- ccf(ccf_skagit_coho_data$coho1_wild_perhour_inp,
    ccf_skagit_coho_data$coho1_hatchery_perhour_inp,
    lag.max = 10, plot = TRUE, 
    ylab = "cross-correlation", main = "Coho, Skagit")

attach(ccf_skagit_coho_data)
layout(matrix(c(1), 1, 1, byrow = TRUE))
ccf(coho1_wild_perhour_inp, coho1_hatchery_perhour_inp, lag.max = 10, plot = TRUE, 
    ylab = "cross-correlation", main = "Coho, Skagit")
hist(mpg)
hist(disp)


