#Goal -to correct the puyallup coho numbers in the summary table

library(ggplot2)
library(tidyverse)
library(here)

#load data
dungeness <- read.csv(here("data", "dungeness_subset_covariates_diff.csv"))

dungeness_day <- read.csv(here("data", "dungeness_unaggregated_day_w_covariates.csv"))

dungeness_night <- read.csv(here("data", "dungeness_unaggregated_night_w_covariates.csv"))

puyallup <- read.csv(here("puyallup","data", "puyallup_final.csv"))

skagit <- read.csv(here("skagit","data", "skagit_2010-2022_w_covariates.csv"))


#make column for river in all data and then rbind all data

dungeness$river <- "Dungeness"
dungeness_day$river <- "Dungeness"
dungeness_night$river <- "Dungeness"
puyallup$river <- "Puyallup"
skagit$river <- "Skagit"


dungeness_day$daytime_category <- "day"
dungeness_night$daytime_category <- "night"

dungeness_day_night <- rbind(dungeness_day, dungeness_night)



dungeness_subset <- dungeness_day_night %>% 
  select(c("Date","doy","year","river",
           "chinook0_wild_num","chinook0_hatchery_num",
           "coho1_wild_num","coho1_hatchery_num","daytime_category")) %>% 
  pivot_longer(cols = c(chinook0_wild_num,chinook0_hatchery_num,
                        coho1_wild_num,coho1_hatchery_num),
               names_to = c("species","origin"),
               names_pattern = "(.*)_(.*)_num",
               values_to = "num")

puyallup_subset <- puyallup %>% 
  select(c("Date","doy","year","river",
           "chinook0_wild_num_day","chinook0_hatchery_num_day",
           "chinook0_wild_num_night","chinook0_hatchery_num_night",
           "coho1_wild_num_day","coho1_hatchery_num_day",
           "coho1_wild_num_night","coho1_hatchery_num_night"
  ))%>% 
  pivot_longer(cols = c(chinook0_wild_num_day,chinook0_hatchery_num_day,
                        chinook0_wild_num_night,chinook0_hatchery_num_night,
                        coho1_wild_num_day,coho1_hatchery_num_day,
                        coho1_wild_num_night,coho1_hatchery_num_night
  ),
  names_to = c("species","origin","daytime_category"),
  names_pattern = "(.*)_(.*)_num_(.*)",
  values_to = "num")




skagit_subset <- skagit %>%
  select(c("Date","doy","year","river",
           "chinook0_wild_num","chinook0_hatchery_num",
           "coho1_wild_num","coho1_hatchery_num",
           "trap", "daytime_category")) %>% 
  pivot_longer(cols = c(chinook0_wild_num,chinook0_hatchery_num,
                        coho1_wild_num,coho1_hatchery_num),
               names_to = c("species","origin"),
               names_pattern = "(.*)_(.*)_num",
               values_to = "num") %>%
  group_by(Date, doy, year, river, species, origin, daytime_category) %>% 
  summarize(num = sum(num)) %>%
  ungroup()


all_rivers <- rbind(dungeness_subset, puyallup_subset, skagit_subset)

#make a column for proportion of fish caught during the day and night

all_rivers <- all_rivers %>% 
  group_by(Date, doy, year, river, species, origin) %>% 
  mutate(total = sum(num)) %>%
  ungroup() %>%
  group_by(Date, doy, year, river, species, origin, daytime_category) %>% 
  mutate(proportion = num/total) %>%
  ungroup() %>%
  select(-total)



dungeness_chinook_summary <- all_rivers %>% 
  filter(doy >= 130 & doy <= 200 & year != 2015 & river == "Dungeness" & species == "chinook0") %>%
  group_by(river, species, origin) %>% 
  summarize(total = sum(num, na.rm = TRUE),
            max = max(num, na.rm = TRUE))


dungeness_coho_summary <- all_rivers %>% 
  filter(doy >= 120 & doy <= 160 & year != 2015 & river == "Dungeness" & species == "coho1") %>%
  group_by(river, species, origin) %>% 
  summarize(total = sum(num, na.rm = TRUE),
            max = max(num, na.rm = TRUE))

puyallup_chinook_summary <- all_rivers %>% 
  filter(doy >= 130 & doy <= 218 & river == "Puyallup" & species == "chinook0") %>%
  group_by(river, species, origin) %>% 
  summarize(total = sum(num, na.rm = TRUE),
            max = max(num, na.rm = TRUE))

skagit_chinook_summary <- all_rivers %>% 
  filter(doy >= 150 & doy <= 189 & river == "Skagit" & species == "chinook0") %>%
  group_by(river, species, origin) %>% 
  summarize(total = sum(num, na.rm = TRUE),
            max = max(num, na.rm = TRUE))


skagit_coho_summary <- all_rivers %>% 
  filter(doy >= 100 & doy <= 150 & river == "Skagit" & species == "coho1") %>%
  group_by(river, species, origin) %>% 
  summarize(total = sum(num, na.rm = TRUE),
            max = max(num, na.rm = TRUE))


puyallup_coho_summary <- all_rivers %>% 
  filter(doy >= 90 & doy <= 160 & river == "Puyallup" & species == "coho1") %>%
  group_by(river, species, origin) %>% 
  summarize(total = sum(num, na.rm = TRUE),
            max = max(num, na.rm = TRUE)) %>% 
  mutate(total = round(total, 0),
         max = round(max, 0))

table_summary <- rbind(dungeness_chinook_summary, dungeness_coho_summary,
                       puyallup_chinook_summary, puyallup_coho_summary,
                       skagit_chinook_summary, skagit_coho_summary)




puyallup_coho_check <- all_rivers %>% 
  filter(doy >= 100 & doy <= 120 & river == "Puyallup" & species == "coho1" & origin == "wild")



table_summary$species <- factor(table_summary$species, levels = c("chinook0","coho1"))

#change chinook0 to Chinook and coho1 to Coho

table_summary$species <- fct_recode(table_summary$species, 
                                    "Chinook" = "chinook0",
                                    "Coho" = "coho1")
table_summary

#save table

write.csv(table_summary, here("output",
                              "summary_all_rivers_numbers_new_corrected.csv"), 
          row.names = FALSE)
