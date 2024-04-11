#Goal - run the best model from Puyallup with secchi depth

#load libraries

library(MARSS)
library(ggplot2)
library(here)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)
library(GGally)
library(tidyverse)
library(broom)
library(ggpubr)

data <- read.csv(here("puyallup", "data","puyallup_2004-2021_all_days_w_covariates_edited.csv"),
                 na.strings = c("NA",""))

#before interpolating, I need to add 0 to the end of 2018 chinook hatchery column
#if not, there are nans in the 2018 column even after interpolation

data$chinook0_hatchery_perhour_night[dim(data)[1]] <- 0
data$chinook0_hatchery_perhour_day[dim(data)[1]] <- 0

#when secchi_depth_day is NA, replace with secchi_depth values
#when secchi_depth_night is NA, replace with secchi_depth values

data$secchi_depth_day[is.na(data$secchi_depth_day)] <- data$secchi_depth[is.na(data$secchi_depth_day)]

data$secchi_depth_night[is.na(data$secchi_depth_night)] <- data$secchi_depth[is.na(data$secchi_depth_night)]

#interpolating na values for flow_day, flow_night, sechhi_depth_day, secchi_depth_night, hatchery

data <- data %>% 
  mutate(flow_day_inp = na.approx(flow_day, na.rm = FALSE),
         flow_night_inp = na.approx(flow_night, na.rm = FALSE),
         secchi_depth_day_inp = na.approx(secchi_depth_day, na.rm = FALSE),
         secchi_depth_night_inp = na.approx(secchi_depth_night, na.rm = FALSE),
         chinook0_hatchery_perhour_night_inp = na.approx(chinook0_hatchery_perhour_night, na.rm = FALSE),
         chinook0_hatchery_perhour_day_inp = na.approx(chinook0_hatchery_perhour_day, na.rm = FALSE),
         photoperiod_day = photoperiod,
         photoperiod_night = photoperiod,
         lunar_phase_day = lunar_phase,
         lunar_phase_night = lunar_phase)


#making columns for flow diff, photoperiod diff
data <- data %>%
  mutate(flow_diff_day = c(NA,diff(flow_day_inp)),
         flow_diff_night = c(NA,diff(flow_night_inp)),
         photo_diff_day = c(NA,diff(photoperiod_day)),
         photo_diff_night = c(NA,diff(photoperiod_night))
  )

#convert date$Date to as.Date and then get year
#first check where date string is not in the right format
#specify the format and then convert to as.Date


data$Date = as.Date(data$Date, format = "%m/%d/%Y")
data$year = year(data$Date)


# I think the doy limits should be changed from 200 to 218 
covariates_chinook0_puyallup <- arrange(data,doy) %>%
  filter(doy >130 & doy <= 218) %>%
  dplyr::select(year,doy, flow_day_inp, flow_night_inp, 
                photoperiod_day, photoperiod_night, secchi_depth_day_inp, secchi_depth_night_inp,
                lunar_phase_day, lunar_phase_night, flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
                chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_day_inp, flow_night_inp, secchi_depth_day_inp, secchi_depth_night_inp,
    photoperiod_day, photoperiod_night, lunar_phase_day, 
    lunar_phase_night, flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
    chinook0_hatchery_perhour_day_inp, chinook0_hatchery_perhour_night_inp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()



#scaling the variables

num_years = 2021-2004+1
num_rows = num_years*2
num_covariates = 14
total_covariates = dim(covariates_chinook0_puyallup)[1]


for(i in 1:(total_covariates/2)){ # everything except diffs and hatchery
  covariates_chinook0_puyallup[i,] = scale(covariates_chinook0_puyallup[i,])[,1]
}

#just scale

for(i in (total_covariates/2 + 1):(total_covariates)){
  covariates_chinook0_puyallup[i,] = scale(covariates_chinook0_puyallup[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
chinook0_puyallup <- arrange(data,doy) %>%
  filter(doy > 130 & doy <= 218) %>%
  mutate(log.value_day = log(chinook0_wild_perhour_day + 1), 
         log.value_night = log(chinook0_wild_perhour_night + 1)) %>%
  dplyr::select(log.value_day, log.value_night ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_day, log.value_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(chinook0_puyallup)[1]){
  chinook0_puyallup[i,] = scale(chinook0_puyallup[i,])[,1]
}


num_years = 2021-2004+1
num_rows = num_years*2
c <- NULL
name_individual = NULL
for(kk in c(1,5,2,4,7)){
  c = rbind(c,covariates_chinook0_puyallup[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_chinook0_puyallup)[1+(kk-1)*num_rows]
  name_individual = paste(name_individual,substr(name_long,1,nchar(name_long)-9))
  print(name_individual)
  
}

fit.model = c(list(c= c), mod_list(num_rows,5,1, FALSE, FALSE))
fit_chinook0_puyallup <- MARSS(chinook0_puyallup, model=fit.model, silent = TRUE, method = "BFGS",
                               control=list(maxit=2000))

ci_chinook0_puyallup <- tidy(fit_chinook0_puyallup)

chinook0_puyallup_plot <- ggplot(ci_chinook0_puyallup[c(39:44),], 
                                 aes(x = c("Flow", 
                                           "Flow\n difference",
                                           "Secchi depth", 
                                           "Lunar phase",
                                           "Hatchery,\nday", 
                                           "Hatchery,\nnight"),
                                     y = estimate, 
                                     ymin = conf.low, 
                                     ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y =""
       # title = "Dungeness River, Chinook sub yearlings"
  )+
  theme_classic() +
  theme(axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y = element_text(size = 14)) + 
  scale_x_discrete(guide = guide_axis(n.dodge=2),
                   limits = c("Flow\n difference",
                              "Hatchery,\nday", 
                              "Hatchery,\nnight",
                              "Secchi depth", 
                              "Lunar phase",
                              "Flow"))

chinook0_puyallup_plot

#save plot
ggsave(here("output","puyallup_chinook0_secchi_depth.png"),  width = 10, height = 6, units = "in", dpi = 300)
