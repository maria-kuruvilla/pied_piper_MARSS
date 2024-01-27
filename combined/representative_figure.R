
#load libraries


library(ggplot2)
library(here)
library(zoo)
library(MASS)
library(GGally)
library(tidyverse)
library(broom)
library(ggpubr)

#load data

#dungeness



if(.Platform$OS.type == "unix") {
  data_string = here("..","..","data","pied_piper","dungeness")
} else {
  data_string = here("..","..","..","OneDrive","Documents","data","pied_piper","dungeness")
}



data_day <- read.csv(here("data","dungeness_unaggregated_day_w_covariates.csv"))
data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"

data_day_night <- rbind(data_day,data_night)


covariates_chinook0_dungeness <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy >130 & doy <= 200) %>%
  dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_april, 
                lunar_phase, resid, temp_diff, flow_diff, photo_diff, 
                chinook0_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = c(
    temp, flow, photoperiod, atu_april, lunar_phase, 
    resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

num_years = 2020-2005
num_rows = num_years*2


dim(covariates_chinook0_dungeness)

#scale the covariates

for(i in 1:(dim(covariates_chinook0_dungeness)[1]/2)){ # everything except resid, diffs and hatchery
  covariates_chinook0_dungeness[i,] = scale(covariates_chinook0_dungeness[i,])[,1]
}

#skip 2005 day
for(i in (dim(covariates_chinook0_dungeness)[1]/2 +1):(dim(covariates_chinook0_dungeness)[1] - num_rows)){
  covariates_chinook0_dungeness[i,] = scale(covariates_chinook0_dungeness[i,], center = FALSE, scale= TRUE)[,1]
}

#skip 2005 night
for(i in (dim(covariates_chinook0_dungeness)[1] - num_rows + 2):(dim(covariates_chinook0_dungeness)[1] - num_rows/2)){
  covariates_chinook0_dungeness[i,] = scale(covariates_chinook0_dungeness[i,], center = FALSE, scale= TRUE)[,1]
}


for(i in (dim(covariates_chinook0_dungeness)[1] - num_rows/2 + 2):(dim(covariates_chinook0_dungeness)[1])){
  covariates_chinook0_dungeness[i,] = scale(covariates_chinook0_dungeness[i,], center = FALSE, scale= TRUE)[,1]
}




#subset response variable
chinook0_dungeness <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy > 130 & doy <= 200) %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

chinook0_dungeness_wild_hatchery <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy > 130 & doy <= 200) %>%
  mutate(log.wild = log(chinook0_wild_perhour + 1), log.hatchery = log(chinook0_hatchery_perhour +1),
         log.hatchery.interpolate = log(chinook0_hatchery_perhour_interpolate +1)) %>%
         
  dplyr::select(log.wild, log.hatchery, log.hatchery.interpolate, year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = c(log.wild,log.hatchery,log.hatchery.interpolate))



colnames(chinook0_dungeness_wild_hatchery)

for(i in 2:(dim(chinook0_dungeness_wild_hatchery)[2]/3+1)){
  chinook0_dungeness_wild_hatchery[,i] = scale(chinook0_dungeness_wild_hatchery[,i])[,1]
  print(i)
  print(colnames(chinook0_dungeness_wild_hatchery)[i])
}

for(i in 31:dim(chinook0_dungeness_wild_hatchery)[2]){
  chinook0_dungeness_wild_hatchery[,i] = scale(chinook0_dungeness_wild_hatchery[,i],center = FALSE)[,1]
  print(i)
  print(colnames(chinook0_dungeness_wild_hatchery)[i])
}

chinook0_dungeness_wild_hatchery <- chinook0_dungeness_wild_hatchery %>% 
  #if hatchery is NA, then interpolate keep the value, otherwise replace with NA
  mutate(log.hatchery.interpolate_2005_Night = ifelse(is.na(log.hatchery_2005_Night), log.hatchery.interpolate_2005_Night, NA),
         log.hatchery.interpolate_2006_Night = ifelse(is.na(log.hatchery_2006_Night), log.hatchery.interpolate_2006_Night, NA),
         log.hatchery.interpolate_2007_Night = ifelse(is.na(log.hatchery_2007_Night), log.hatchery.interpolate_2007_Night, NA),
         log.hatchery.interpolate_2008_Night = ifelse(is.na(log.hatchery_2008_Night), log.hatchery.interpolate_2008_Night, NA),
         log.hatchery.interpolate_2009_Night = ifelse(is.na(log.hatchery_2009_Night), log.hatchery.interpolate_2009_Night, NA),
         log.hatchery.interpolate_2010_Night = ifelse(is.na(log.hatchery_2010_Night), log.hatchery.interpolate_2010_Night, NA),
         log.hatchery.interpolate_2011_Night = ifelse(is.na(log.hatchery_2011_Night), log.hatchery.interpolate_2011_Night, NA),
         log.hatchery.interpolate_2012_Night = ifelse(is.na(log.hatchery_2012_Night), log.hatchery.interpolate_2012_Night, NA),
         log.hatchery.interpolate_2013_Night = ifelse(is.na(log.hatchery_2013_Night), log.hatchery.interpolate_2013_Night, NA),
         log.hatchery.interpolate_2014_Night = ifelse(is.na(log.hatchery_2014_Night), log.hatchery.interpolate_2014_Night, NA),
         log.hatchery.interpolate_2016_Night = ifelse(is.na(log.hatchery_2016_Night), log.hatchery.interpolate_2016_Night, NA),
         log.hatchery.interpolate_2017_Night = ifelse(is.na(log.hatchery_2017_Night), log.hatchery.interpolate_2017_Night, NA),
         log.hatchery.interpolate_2018_Night = ifelse(is.na(log.hatchery_2018_Night), log.hatchery.interpolate_2018_Night, NA),
         log.hatchery.interpolate_2019_Night = ifelse(is.na(log.hatchery_2019_Night), log.hatchery.interpolate_2019_Night, NA),
         log.hatchery.interpolate_2020_Night = ifelse(is.na(log.hatchery_2020_Night), log.hatchery.interpolate_2020_Night, NA),
         
         
  )
  


ggplot(chinook0_dungeness_wild_hatchery) +
  geom_line(aes(x = doy, y = log.wild_2016_Night, color = "wild"), linewidth = 1, alpha = 0.8)+
  geom_line(aes(x = doy, y = log.hatchery_2016_Night, color = "hatchery"), linewidth = 1, alpha = 0.8)+
  geom_point(aes(x = doy, y = log.hatchery.interpolate_2016_Night, color = "hatchery"), size = 0.2, alpha = 0.5)+
  scale_color_manual(values = c("wild" = "salmon", "hatchery" = "cadetblue"))+
  theme_classic()+
  labs(x = "Day of Year", y = "Log Chinook per Hour", color = "Origin")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
        )

ggplot(chinook0_dungeness_wild_hatchery) +
  geom_line(aes(x = doy, y = log.wild_2008_Night, color = "wild"), linewidth = 1, alpha = 0.8)+
  geom_line(aes(x = doy, y = log.hatchery_2008_Night, color = "hatchery"), linewidth = 1, alpha = 0.8)+
  geom_point(aes(x = doy, y = log.hatchery.interpolate_2008_Night, color = "hatchery"), size = 0.2, alpha = 0.5)+
  scale_color_manual(values = c("wild" = "salmon", "hatchery" = "cadetblue"))+
  theme_classic()+
  labs(x = "Day of Year", y = "Log Chinook per Hour", color = "Origin")+
  theme(legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
)

#2006
ggplot(chinook0_dungeness_wild_hatchery) +
  geom_line(aes(x = doy, y = log.wild_2006_Night, color = "wild"), linewidth = 1, alpha = 0.8)+
  geom_line(aes(x = doy, y = log.hatchery_2006_Night, color = "hatchery"), linewidth = 1, alpha = 0.8)+
  geom_point(aes(x = doy, y = log.hatchery.interpolate_2006_Night, color = "hatchery"), size = 0.2, alpha = 0.5)+
  scale_color_manual(values = c("wild" = "salmon", "hatchery" = "cadetblue"))+
  theme_classic()+
  labs(x = "Day of Year", y = "Log Chinook per Hour", color = "Origin")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
  )

#2007
ggplot(chinook0_dungeness_wild_hatchery) +
  geom_line(aes(x = doy, y = log.wild_2007_Night, color = "wild"), linewidth = 1, alpha = 0.8)+
  geom_line(aes(x = doy, y = log.hatchery_2007_Night, color = "hatchery"), linewidth = 1, alpha = 0.8)+
  geom_point(aes(x = doy, y = log.hatchery.interpolate_2007_Night, color = "hatchery"), size = 0.2, alpha = 0.5)+
  scale_color_manual(values = c("wild" = "salmon", "hatchery" = "cadetblue"))+
  theme_classic()+
  labs(x = "Day of Year", y = "Log Chinook per Hour", color = "Origin")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
  )

#2009

ggplot(chinook0_dungeness_wild_hatchery) +
  geom_line(aes(x = doy, y = log.wild_2009_Night, color = "wild"), linewidth = 1, alpha = 0.8)+
  geom_line(aes(x = doy, y = log.hatchery_2009_Night, color = "hatchery"), linewidth = 1, alpha = 0.8)+
  scale_color_manual(values = c("wild" = "salmon", "hatchery" = "cadetblue"))+
  theme_classic()+
  labs(x = "Day of Year", y = "Log Chinook per Hour", color = "Origin")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
  )

#2010

ggplot(chinook0_dungeness_wild_hatchery) +
  geom_line(aes(x = doy, y = log.wild_2010_Night, color = "wild"), linewidth = 1, alpha = 0.8)+
  geom_line(aes(x = doy, y = log.hatchery_2010_Night, color = "hatchery"), linewidth = 1, alpha = 0.8)+
  geom_point(aes(x = doy, y = log.hatchery.interpolate_2010_Night, color = "hatchery"), size = 0.2, alpha = 0.5)+
  scale_color_manual(values = c("wild" = "salmon", "hatchery" = "cadetblue"))+
  theme_classic()+
  labs(x = "Day of Year", y = "Log Chinook per Hour", color = "Origin")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
  )

#make the chinook0_dungeness_wild_hatchery dataframe long
#select log.wild_2006_Night, log.wild_2007_Night, log.wild_2008_Night, log.wild_2009_Night, log.wild_2010_Night, etc
#pivot_longer() to make a new dataframe with
#with doy, year, daytime_category, origin, wild, hatchery, interpolate columns

chinook0_dungeness_wild_hatchery_long <- chinook0_dungeness_wild_hatchery %>%
  pivot_longer(cols = -doy,
               names_to = c("origin", "year", "daytime_category"),
               names_sep = "_",
               values_to = "log.chinook_per_hour")

ggplot(chinook0_dungeness_wild_hatchery_long %>% 
         filter(daytime_category == 'Night')) +
  geom_line(aes(x = doy, y = log.chinook_per_hour, color = origin), alpha = 0.5, linewidth = 0.8)+
  facet_wrap(~year, ncol = 5)+
  scale_color_manual(values = c("log.wild" = "salmon", "log.hatchery" = "cadetblue"))+
  theme_classic()+
  labs(x = "Day of Year", y = "Log Chinook per Hour", color = "Origin")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white", color = "white")
  )    

ggsave(here("output","dungeness_chinook_wild_hatchery_all_years.png"), width = 12, height = 8, units = "in")    


chinook0_dungeness_wild_hatchery_linear <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy > 130 & doy <= 200) %>%
  mutate(wild = chinook0_wild_perhour, hatchery = chinook0_hatchery_perhour,
         hatchery.interpolate = chinook0_hatchery_perhour_interpolate) %>%
  
  dplyr::select(wild, hatchery, hatchery.interpolate, year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = c(wild,hatchery,hatchery.interpolate))


for(i in 2:(dim(chinook0_dungeness_wild_hatchery_linear)[2]/3+1)){
  chinook0_dungeness_wild_hatchery_linear[,i] = scale(chinook0_dungeness_wild_hatchery_linear[,i])[,1]
  print(i)
  print(colnames(chinook0_dungeness_wild_hatchery_linear)[i])
}

for(i in 31:dim(chinook0_dungeness_wild_hatchery_linear)[2]){
  chinook0_dungeness_wild_hatchery_linear[,i] = scale(chinook0_dungeness_wild_hatchery_linear[,i],center = FALSE)[,1]
  print(i)
  print(colnames(chinook0_dungeness_wild_hatchery_linear)[i])
}

chinook0_dungeness_wild_hatchery_linear_long <- chinook0_dungeness_wild_hatchery_linear %>%
  pivot_longer(cols = -doy,
               names_to = c("origin", "year", "daytime_category"),
               names_sep = "_",
               values_to = "chinook_per_hour")

ggplot(chinook0_dungeness_wild_hatchery_linear_long %>% 
         filter(daytime_category == 'Night', origin != 'hatchery.interpolate')) +
  geom_line(aes(x = doy, y = chinook_per_hour, color = origin), alpha = 0.5, linewidth = 0.8)+
  facet_wrap(~year, ncol = 5)+
  scale_color_manual(values = c("wild" = "salmon", "hatchery" = "cadetblue"))+
  theme_classic()+
  labs(x = "Day of Year", y = "Scaled Chinook per Hour", color = "Origin")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white", color = "white")
  )    

ggsave(here("output","dungeness_chinook_wild_hatchery_all_years_linear.png"), width = 12, height = 8, units = "in")

ggplot(chinook0_dungeness_wild_hatchery_linear_long %>% 
         arrange(origin, doy ,year) %>%
         filter(daytime_category == 'Night', origin != 'hatchery.interpolate', year == 2016 | year == 2020)) +
  geom_line(aes(x = doy, y = chinook_per_hour, color = origin, alpha = origin), linewidth = 1.2)+
  facet_wrap(~year, ncol = 5)+
  scale_color_manual(name = "",labels = c("hatchery", "wild"),values = c("wild" = "salmon","hatchery" = "cadetblue"))+
  scale_alpha_manual(name = "", labels = c("hatchery", "wild"),values = c("wild" = 0.5,"hatchery" = 0.8))+
  theme_classic()+
  labs(x = "Day of Year", y = "Scaled number of Chinook per Hour", color = "Origin")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 14))

ggsave(here("output","dungeness_chinook_wild_hatchery_all_years_linear.png"), width = 7, height = 5, units = "in")


