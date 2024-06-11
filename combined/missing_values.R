#goal - calculate number of missing values in the dataset

library(tidyverse)
library(here)


# load the dataset

data_day <- read.csv(here("data","dungeness_unaggregated_day_w_covariates.csv"))
data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"

data_day_night <- rbind(data_day,data_night)

glimpse(data_day_night)

# calculate number of missing values in the dataset



missing_values_chinook_dungeness <- data_day_night %>%
  filter(doy > 130, doy < 200, year != 2015) %>%
  select(chinook0_wild_perhour) %>% 
  summarise_all(funs(sum(is.na(.)))) %>%
  gather(key = "variable", value = "missing_values") %>%
  arrange(desc(missing_values))

nonmissing_values_chinook_dungeness <- data_day_night %>%
  filter(doy > 130, doy < 200, year != 2015) %>%
  select(chinook0_wild_perhour) %>% 
  summarise_all(funs(sum(!is.na(.)))) %>%
  gather(key = "variable", value = "missing_values") %>%
  arrange(desc(missing_values))

nonmissing_values_chinook_dungeness[1,2] + missing_values_chinook_dungeness[1,2]

#335/2070


missing_values_coho_dungeness <- data_day_night %>%
  filter(doy > 120, doy < 160, year != 2015) %>%
  select(coho1_wild_perhour) %>% 
  summarise_all(funs(sum(is.na(.)))) %>%
  gather(key = "variable", value = "missing_values") %>%
  arrange(desc(missing_values))

nonmissing_values_coho_dungeness <- data_day_night %>%
  filter(doy > 120, doy < 160, year != 2015) %>%
  select(coho1_wild_perhour) %>% 
  summarise_all(funs(sum(!is.na(.)))) %>%
  gather(key = "variable", value = "missing_values") %>%
  arrange(desc(missing_values))


missing_values_coho_dungeness[1,2] + nonmissing_values_coho_dungeness[1,2]

#191/1170


#load puyallup dataset

puyallup <- read.csv(here("puyallup","data","puyallup_final.csv"))

missing_values_chinook_puyallup <- puyallup %>% 
  filter(doy >130, doy < 218) %>% 
  select(chinook0_wild_perhour_day) %>% 
  summarise_all(funs(sum(is.na(.)))) %>%
  gather(key = "variable", value = "missing_values")

missing_values_chinook_puyallup_night <- puyallup %>% 
  filter(doy >130, doy < 218) %>% 
  select(chinook0_wild_perhour_night) %>% 
  summarise_all(funs(sum(is.na(.)))) %>%
  gather(key = "variable", value = "missing_values")



nonmissing_values_chinook_puyallup <- puyallup %>% 
  filter(doy >130, doy < 218) %>% 
  select(chinook0_wild_perhour_day) %>% 
  summarise_all(funs(sum(!is.na(.)))) %>%
  gather(key = "variable", value = "missing_values")



nonmissing_values_chinook_puyallup_night <- puyallup %>% 
  filter(doy >130, doy < 218) %>% 
  select(chinook0_wild_perhour_night) %>% 
  summarise_all(funs(sum(!is.na(.)))) %>%
  gather(key = "variable", value = "missing_values")
  
missing_values_chinook_puyallup[1,2]+ nonmissing_values_chinook_puyallup[1,2] 


missing_values_chinook_puyallup_night[1,2] + nonmissing_values_chinook_puyallup_night[1,2]


#552/1566 

#596/1566

#1148/3132

#skagit data

skagit <- read.csv(here("skagit","data","skagit_2010-2022_w_covariates.csv"))

missing_values_chinook_skagit <- skagit %>% 
  filter(doy >150, doy < 189, daytime_category == "night") %>% 
  select(chinook0_wild_num) %>% 
  summarise_all(funs(sum(is.na(.)))) %>%
  gather(key = "variable", value = "missing_values")


nonmissing_values_chinook_skagit <- skagit %>%
  filter(doy >150, doy < 189, daytime_category == "night") %>% 
  select(chinook0_wild_num) %>% 
  summarise_all(funs(sum(!is.na(.)))) %>%
  gather(key = "variable", value = "missing_values")

missing_values_chinook_skagit[1,2] + nonmissing_values_chinook_skagit[1,2]

#159/912

missing_values_coho_skagit <- skagit %>% 
  filter(doy >100, doy < 150, daytime_category == "night") %>% 
  select(coho1_wild_num) %>% 
  summarise_all(funs(sum(is.na(.)))) %>%
  gather(key = "variable", value = "missing_values")

nonmissing_values_coho_skagit <- skagit %>%
  filter(doy >100, doy < 150, daytime_category == "night") %>% 
  select(coho1_wild_num) %>% 
  summarise_all(funs(sum(!is.na(.)))) %>%
  gather(key = "variable", value = "missing_values")


missing_values_coho_skagit[1,2] + nonmissing_values_coho_skagit[1,2]


# 108/1176




