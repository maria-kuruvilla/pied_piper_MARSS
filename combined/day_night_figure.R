
#Goal - to make one figure that shows the proportion 
#of fish caught during the day and night for each species and origin
#for all 3 rivers

#load data
dungeness <- read.csv(here("data", "dungeness_subset_covariates_diff.csv"))

dungeness_day <- read.csv(here("data", "dungeness_unaggregated_day_w_covariates.csv"))

dungeness_night <- read.csv(here("data", "dungeness_unaggregated_night_w_covariates.csv"))

puyallup <- read.csv(here("puyallup","data", "puyallup_2004-2021_all_days_w_covariates_edited.csv"))

skagit <- read.csv(here("skagit","data", "skagit_2010-2022_w_covariates.csv"))


#make column for river in all data and then rbind all data

dungeness$river <- "Dungeness"
dungeness_day$river <- "Dungeness"
dungeness_night$river <- "Dungeness"
puyallup$river <- "Puyallup"
skagit$river <- "Skagit"


dungeness_day$daytime_category <- "day"
dungeness_night$daytime_category <- "night"

#check unique secchi_depth values in puyallup
unique(puyallup$secchi_depth_night)

#In puyallup when secchi_depth is NA, average secchi_depth_day and secchi_depth_night
#unless both are NA, then secchi_depth is NA
#first comvert all '-' to NA
puyallup <- puyallup %>%
  mutate(secchi_depth_day = ifelse(secchi_depth_day == " - ", NA, as.numeric(secchi_depth_day)),
         secchi_depth_night = ifelse(secchi_depth_night == " - ", NA, as.numeric(secchi_depth_night))) %>%
  mutate(flag = ifelse(is.na(secchi_depth_day) & is.na(secchi_depth_night), 1, 0)) %>%
  mutate(secchi_depth = ifelse(is.na(secchi_depth) & flag == 0,
                               (secchi_depth_day + secchi_depth_night)/2, secchi_depth))


#combine dungeness_day and dungeness_night

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


#the last row of chinook0_hatchery_perhour in puyallup should be set to 0


puyallup_subset <- puyallup %>% 
  select(c("Date","doy","year","river",
           "chinook0_wild_num_day","chinook0_hatchery_num_day",
           "chinook0_wild_num_night","chinook0_hatchery_num_night"
           ))%>% 
  pivot_longer(cols = c(chinook0_wild_num_day,chinook0_hatchery_num_day,
                        chinook0_wild_num_night,chinook0_hatchery_num_night),
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


#make columns to find the total number, the average and the maximum 
#number of fish caught during doy 100 to 200 for each species and origin
#and river

all_rivers_summary <- all_rivers %>% 
  filter(doy >= 100 & doy <= 200) %>%
  group_by(river, species, origin) %>% 
  summarize(total = sum(num, na.rm = TRUE),
            max = max(num, na.rm = TRUE))

dungeness_chinook_summary <- all_rivers %>% 
  filter(doy >= 130 & doy <= 200 & year != 2015 & river == "Dungeness" & species == "chinook0") %>%
  group_by(river, species, origin) %>% 
  summarize(total = sum(num, na.rm = TRUE),
            max = max(num, na.rm = TRUE))

dungeness_coho_summary <- all_rivers %>% 
  filter(doy >= 120 & doy <= 1600 & year != 2015 & river == "Dungeness" & species == "coho1") %>%
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


table_summary <- rbind(dungeness_chinook_summary, dungeness_coho_summary,
      puyallup_chinook_summary, skagit_chinook_summary, skagit_coho_summary)

#change chinook0 to Chinook and coho1 to Coho

table_summary$species <- factor(table_summary$species, levels = c("chinook0","coho1"))

#change chinook0 to Chinook and coho1 to Coho

table_summary$species <- fct_recode(table_summary$species, 
                                    "Chinook" = "chinook0",
                                    "Coho" = "coho1")

#save table

write.csv(table_summary, here("output",
                              "summary_all_rivers_numbers.csv"), 
          row.names = FALSE)


#take the average of the proportions by species, origin, daytime_category, 
#and river 
#then make plot

all_rivers_avg <- all_rivers %>% 
  group_by(river, species, origin, daytime_category) %>% 
  summarize(mean_proportion = mean(proportion, na.rm = TRUE), 
            n = n(), 
            se = sd(proportion, na.rm = T)/sqrt(n)) %>%
  ungroup() %>%
  mutate(daytime_category = factor(daytime_category, levels = c("day","night")))

#convert species to factor

all_rivers_avg$species <- factor(all_rivers_avg$species, levels = c("chinook0","coho1"))

#change chinook0 to Chinook and coho1 to Coho

all_rivers_avg$species <- fct_recode(all_rivers_avg$species, 
                                     "Chinook" = "chinook0",
                                     "Coho" = "coho1")

#plot
#change name of species to Chinook and Coho
ggplot(all_rivers_avg, aes(x = daytime_category, y = mean_proportion, fill = origin)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_proportion - se, ymax = mean_proportion + se),
                position = position_dodge(width = 0.9), width = 0.25) +
  facet_grid(species~river) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0)) +
  labs(x = "", y = "Proportion of fish caught", fill = "Origin") +
  scale_fill_manual(values = c("cadetblue","salmon")) +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.background = element_blank(),
        )

ggsave(here("output",
            "all_rivers_proportion_day_night.png"),
       width = 13, height = 8, units = "in", dpi = 300)
########




data <- read.csv(here("puyallup", "data","puyallup_2004-2021_all_days_w_covariates.csv"),
                 na.strings = c("NA",""))

data_subset <- data %>% 
  select(chinook0_wild_day_fraction,chinook0_hatchery_day_fraction)

data_subset_long <- data_subset %>% 
  pivot_longer(cols = c(chinook0_wild_day_fraction,chinook0_hatchery_day_fraction),
               names_to = c("species","origin"),
               names_pattern = "(.*)_(.*)_day_fraction",
               values_to = "proportion_day") %>%
  mutate(proportion_night = 1 - proportion_day) %>%
  pivot_longer(cols = c(proportion_day,proportion_night),
               names_to = "daytime_category",
               names_pattern = "proportion_(.*)",
               values_to = "proportion")

data_long_subset_avg <-   data_subset_long %>% 
  filter(!is.na(proportion)) %>%
  group_by(species,daytime_category, origin) %>%
  summarize(average_prop = mean(proportion, na.rm = T), n = n(), se = sd(proportion, na.rm = T)/sqrt(n))

ggplot(data = data_long_subset_avg, aes(y = average_prop, x = daytime_category, fill = origin)) +
  geom_col(position = "dodge")+
  facet_wrap(~species)+
  geom_errorbar(aes(ymin = average_prop-se, ymax = average_prop+se), 
                position = "dodge", alpha = 0.5)+
  theme(axis.text.x=element_text(size=24),axis.title.y=element_text(size=24))+
  labs(x = "",
       y = "Average proportion of \nfish caught per hour", 
       fill = "Origin")+
  scale_fill_manual(values = c("cadetblue","salmon"))

ggsave(here("puyallup","output","chinook_all_years_day_night_proportion.png"), width = 10, height = 8, units = "in", dpi = 300)
