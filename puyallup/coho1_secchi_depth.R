#load libraries
library(MARSS)
library(ggplot2)
library(here)
library(GGally)
library(tidyverse)
library(zoo)
library(modelr)
library(qpcR)
library(broom)

#load data
data <- read.csv(here("puyallup", "data","puyallup_final.csv"),
                 na.strings = c("NA",""))

data$coho1_hatchery_perhour_night[dim(data)[1]] <- 0
data$coho1_hatchery_perhour_day[dim(data)[1]] <- 0


data <- data %>% 
  mutate(coho1_hatchery_perhour_night = na.approx(coho1_hatchery_perhour_night, na.rm = FALSE),
         coho1_hatchery_perhour_day = na.approx(coho1_hatchery_perhour_day, na.rm = FALSE))



covariates_coho1_puyallup_w_temp <- arrange(data,doy) %>%
  filter(doy >90 & doy <= 160) %>%
  dplyr::select(year,doy, flow_day, flow_night, 
                photoperiod_day, photoperiod_night, secchi_depth_day, secchi_depth_night,
                lunar_phase_day, lunar_phase_night, temp_day, temp_night, atu_day, atu_night,
                flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
                temp_diff_day, temp_diff_night, 
                resid_day, resid_night,
                coho1_hatchery_perhour_day, coho1_hatchery_perhour_night) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_day, flow_night, secchi_depth_day, secchi_depth_night,
    photoperiod_day, photoperiod_night, lunar_phase_day, 
    lunar_phase_night, temp_day, temp_night, atu_day, atu_night,
    flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
    temp_diff_day, temp_diff_night,
    resid_day, resid_night,
    coho1_hatchery_perhour_day, coho1_hatchery_perhour_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2021-2004+1
num_rows = num_years*2
total_covariates = dim(covariates_coho1_puyallup_w_temp)[1]


for(i in 1:(num_rows*6)){ # everything except diffs and hatchery
  covariates_coho1_puyallup_w_temp[i,] = scale(covariates_coho1_puyallup_w_temp[i,])[,1]
}

#just scale

for(i in (num_rows*6 + 1):(total_covariates)){
  covariates_coho1_puyallup_w_temp[i,] = scale(covariates_coho1_puyallup_w_temp[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_coho_summer_perhour <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160) %>%
  mutate(log.value_day = log(coho1_wild_perhour_day + 1), 
         log.value_night = log(coho1_wild_perhour_night + 1)) %>%
  dplyr::select(log.value_day, log.value_night ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_day, log.value_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour)[1]){
  subset_coho_summer_perhour[i,] = scale(subset_coho_summer_perhour[i,])[,1]
}



#trying with hatchery
c<-NULL
for(kk in c(1,2,3,7,11)){
  c = rbind(c,covariates_coho1_puyallup_w_temp[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_coho1_puyallup_w_temp)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  
}
fit.model = c(list(c= c), mod_list(num_rows,5,1, FALSE, TRUE))
fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
ci <- tidy(fit)

coho1_puyallup_plot <- ggplot(ci[c(40:45),], 
                                 aes(x = c("Flow", 
                                           "Secchi depth", 
                                           "Photoperiod",
                                           "Flow\n difference",
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
                              "Photoperiod",
                              "Flow"))

coho1_puyallup_plot

#save plot
ggsave(here("output","puyallup_coho1_secchi_depth.png"),  width = 10, height = 6, units = "in", dpi = 300)
