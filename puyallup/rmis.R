

rmis_puyallup <- read_csv("https://www.rmis.org/reports/CSV27588.txt")

#read first_release_date and make the first four digits the year
#next two digits as the month, and the last two digits as the day
rmis_puyallup <- rmis_puyallup %>% 
  mutate(year = as.numeric(substr(first_release_date, 1, 4)),
         month = substr(first_release_date, 5, 6),
         day = substr(first_release_date, 7, 8),
         date = as.Date(paste(year, month, day, sep = "-")),
         species_name = ifelse(species == 1, "chinook", "coho"),
         age = year - brood_year,
         total = untagged_unclipped + untagged_adclipped + tagged_unclipped + tagged_adclipped + untagged_unknown
         )

chinook_subset <- rmis_puyallup %>% 
  filter(species == 1, year >= 2004, year <= 2021) %>% 
  select(date, age, total, untagged_unclipped) %>% 
  group_by(date, age) %>% 
  summarise(total = sum(total), untagged_unclipped = sum(untagged_unclipped)) %>% 
  mutate(prop_unmarked = untagged_unclipped / total)

coho_subset <- rmis_puyallup %>% 
  filter(species == 2, year >= 2004, year <= 2021) %>% 
  select(date, age, total, untagged_unclipped) %>% 
  group_by(date, age) %>% 
  summarise(total = sum(total), untagged_unclipped = sum(untagged_unclipped)) %>% 
  mutate(prop_unmarked = untagged_unclipped / total)


print(chinook_subset)
