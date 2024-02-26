#use the output from the code in estimate_figures.R
#use the MARSS output from the best/good models
#to make the prediction figures of the best/good models
#also make additional figure with the hatchery covariate

#load libraries
library(MARSS)
library(ggplot2)
library(tidyverse)
library(here)
library(zoo)

predict_coho1_skagit <- predict(fit_coho1_skagit, type = "ytT", interval = "confidence")
glimpse(predict_coho1_skagit)

head(predict_coho1_skagit$pred)


predict_coho1_skagit$pred$trap <- substring(unique(predict_coho1_skagit$pred$.rownames), 12, 
                                            length(unique(predict_coho1_skagit$pred$.rownames)))

predict_coho1_skagit$pred$year <- substring(unique(predict_coho1_skagit$pred$.rownames), 1,4)

predict_coho1_skagit$pred$daynight_category <- substring(unique(predict_coho1_skagit$pred$.rownames), 6,10)

skagit_covariates_coho1 <- as.data.frame(t(covariates_coho1_skagit_night))
glimpse(skagit_covariates_coho1)

#make skagit_covariates_coho1 long by taking year, daynight_category, trap out of the column names

skagit_covariates_coho1_long <- skagit_covariates_coho1 %>% 
  pivot_longer(cols = -c(doy), names_to = c(".value","year", "daynight_category","trap"),
               names_pattern = "(.*)_(.{4})_(.{5})_(.{5})")
head(skagit_covariates_coho1_long)

#merge skagit_covariates_coho1_long$coho1_hatchery_perhour_inp with predict_coho1_skagit$pred

predict_coho1_skagit$pred <- predict_coho1_skagit$pred %>% 
  left_join(skagit_covariates_coho1_long, by = c("year", "daynight_category", "trap"))


#make column for doy
skagit_covariates_coho1$doy <- as.numeric(rownames(skagit_covariates_coho1))

ggplot(data = predict_coho1_skagit$pred %>% 
         filter(year == 2016 & daynight_category == "night" & trap == 'scoop'))+
  
  geom_line(aes(x = t+100, y = estimate, color = "predicted"))+
  geom_point(aes(x = t+100, y = y, color = "observed"), alpha = 0.5)+
  geom_ribbon(aes(x = t+100, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = skagit_covariates_coho1, aes(x = doy, y = log(coho1_hatchery_perhour_inp_2016_night_scoop+1)), 
            color = "cadetblue", linewidth = 1.2, alpha = 0.7) +
  scale_color_manual(name = "Legend", values = c("predicted" = "salmon", "observed" = "salmon"))+
  labs(x = "Time", y = "log(scaled(coho salmon per hour) + 1)", title = "")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1.2, linetype="solid"
    ),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "none"
  )
#257085


ggplot(data = predict_coho1_skagit$pred)+
  
  geom_line(aes(x = t+100, y = estimate, color = "predicted"))+
  geom_point(aes(x = t+100, y = y, color = "observed"), alpha = 0.5)+
  geom_ribbon(aes(x = t+100, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  facet_wrap(~year+trap, ncol = 4, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(coho salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("predicted" = "salmon", "observed" = "salmon"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA),
                       shape = c(NA, 19))))+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

ggsave(here("output","coho1_skagit_prediction.png"), width = 12, height = 14, units = "in", dpi = 300)

ggplot(data = predict_coho1_skagit$pred)+
  
  geom_line(aes(x = t+100, y = estimate, color = "wild, predicted"))+
  geom_point(aes(x = t+100, y = y, color = "wild, observed"), alpha = 0.8, size = 0.2)+
  geom_ribbon(aes(x = t+100, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_coho1_skagit$pred, aes(x = doy, y = log(coho1_hatchery_perhour_inp+1), 
            color = "hatchery")) +
  facet_wrap(~year+trap, ncol = 4, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log (coho salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", "hatchery" = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA, 1),
                       shape = c(NA, 19, NA),
                       size = c(4,2,4))))+
  scale_y_continuous(breaks = c(-2,0, 2))+
  scale_x_continuous(breaks = c(100, 120,140))+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave(here("output","coho1_skagit_prediction_w_hatchery.png"), width = 6, height = 8, units = "in", dpi = 300)

#make same plots for chinook0_skagit_night

predict_chinook0_skagit_night <- predict(fit_chinook0_skagit, type = "ytT", interval = "confidence")
glimpse(predict_chinook0_skagit_night)

head(predict_chinook0_skagit_night$pred)

predict_chinook0_skagit_night$pred$trap <- substring(predict_chinook0_skagit_night$pred$.rownames, 12, 
                                            length(predict_chinook0_skagit_night$pred$.rownames))

predict_chinook0_skagit_night$pred$year <- substring(predict_chinook0_skagit_night$pred$.rownames, 1,4)

predict_chinook0_skagit_night$pred$daynight_category <- substring(predict_chinook0_skagit_night$pred$.rownames, 6,10)

skagit_covariates_chinook0 <- as.data.frame(t(covariates_chinook0_skagit_night))
glimpse(skagit_covariates_chinook0)

colnames(skagit_covariates_chinook0)

skagit_covariates_chinook0$doy <- as.numeric(rownames(skagit_covariates_chinook0))

skagit_covariates_chinook0_long <- skagit_covariates_chinook0 %>% 
  pivot_longer(cols = -c(doy), names_to = c(".value","year", "daynight_category","trap"),
               names_pattern = "(.*)_(.{4})_(.{5})_(.{5})")
head(skagit_covariates_coho1_long)

predict_chinook0_skagit_night$pred$doy <- predict_chinook0_skagit_night$pred$t+150
predict_chinook0_skagit_night$pred <- predict_chinook0_skagit_night$pred %>% 
  left_join(skagit_covariates_chinook0_long, by = c("doy","year", "daynight_category", "trap"))

head(predict_chinook0_skagit_night$pred)

ggplot(data = predict_chinook0_skagit_night$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "predicted"))+
  geom_point(aes(x = doy, y = y, color = "observed"), alpha = 0.5)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  facet_wrap(~year+trap, ncol = 4, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(Chinook salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("predicted" = "salmon", "observed" = "salmon"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(NA,1),
                       shape = c(19, NA))))+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))

ggsave(here("output","chinook0_skagit_prediction.png"), width = 5, height = 7, units = "in", dpi = 300)


ggplot(data = predict_chinook0_skagit_night$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"))+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_chinook0_skagit_night$pred, aes(x = doy, y = log(chinook0_hatchery_perhour_inp+1), 
                                                  color = "hatchery")) +
  facet_wrap(~year+trap, ncol = 4, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(Chinook salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", hatchery = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave(here("output","chinook0_skagit_prediction_w_hatchery.png"), width = 6, height = 8, units = "in", dpi = 300)

### puyallup

predict_chinook0_puyallup <- predict(fit_chinook0_puyallup, type = "ytT", interval = "confidence")

glimpse(predict_chinook0_puyallup$pred)

head(predict_chinook0_puyallup$pred)

predict_chinook0_puyallup$pred$trap <- 'screw'

#last 4 characters of rowname is year, but length of rowname varies depending on whether it's a day or night trap
predict_chinook0_puyallup$pred$year <-  as.numeric(substr(predict_chinook0_puyallup$pred$.rownames, 
                                               nchar(predict_chinook0_puyallup$pred$.rownames)-3, 
                                               nchar(predict_chinook0_puyallup$pred$.rownames)))

#if rowname has 'day' in it, then it's a day trap, otherwise it's a night trap
predict_chinook0_puyallup$pred$daynight_category <- ifelse(grep("day", 
                                                                predict_chinook0_puyallup$pred$.rownames, 
                                                                value = TRUE) == predict_chinook0_puyallup$pred$.rownames, 
                                                           "day", "night")
predict_chinook0_puyallup$pred$doy <- predict_chinook0_puyallup$pred$t+130
puyallup_covariates_chinook0 <- as.data.frame(t(covariates_chinook0_puyallup))

puyallup_covariates_chinook0$doy <- as.numeric(rownames(puyallup_covariates_chinook0))

puyallup_covariates_chinook0_long <-  puyallup_covariates_chinook0 %>% 
  select(doy, starts_with("chinook")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","daynight_category","year"),
               names_pattern = "(.*)_(.*)_inp_(.{4})") %>%
  mutate(year = as.numeric(year), trap = 'screw')

predict_chinook0_puyallup$pred <- predict_chinook0_puyallup$pred %>%
  left_join(puyallup_covariates_chinook0_long, by = c("doy","year", "daynight_category", "trap"))


ggplot(data = predict_chinook0_puyallup$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"))+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_chinook0_puyallup$pred, aes(x = doy, y = log(chinook0_hatchery_perhour+1), 
                                                           color = "hatchery")) +
  facet_wrap(~year+daynight_category, ncol = 6, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(Chinook salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", hatchery = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  scale_y_continuous(breaks = c(-3,0,3))+
  scale_x_continuous(limit = c(130, 200), breaks = c(140,160,180))+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave(here("output","chinook0_puyallup_prediction_w_hatchery.png"), width = 6, height = 8, units = "in", dpi = 300)


#dungeness

predict_chinook0_dungeness <- predict(fit_chinook0_dungeness, type = "ytT", interval = "confidence")

glimpse(predict_chinook0_dungeness$pred)

head(predict_chinook0_dungeness$pred)

predict_chinook0_dungeness$pred$trap <- 'screw'

predict_chinook0_dungeness$pred$year <-  as.numeric(substr(predict_chinook0_dungeness$pred$.rownames, 
                                               1, 4))

predict_chinook0_dungeness$pred$daynight_category <- ifelse(substr(predict_chinook0_dungeness$pred$.rownames, 
                                               6, nchar(predict_chinook0_dungeness$pred$.rownames)) == 'Day', 'day', 'night')

predict_chinook0_dungeness$pred$doy <- predict_chinook0_dungeness$pred$t+130

dungeness_covariates_chinook0 <- as.data.frame(t(covariates_chinook0_dungeness))

dungeness_covariates_chinook0$doy <- as.numeric(rownames(dungeness_covariates_chinook0))

dungeness_covariates_chinook0_long <-  dungeness_covariates_chinook0 %>% 
  select(doy, starts_with("chinook")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","year","daynight_category"),
               names_pattern = "(.*)_(.{4})_(.*)") %>%
  mutate(year = as.numeric(year), trap = 'screw', daynight_category = ifelse(daynight_category == 'Day', 'day', 'night'))

predict_chinook0_dungeness$pred <- predict_chinook0_dungeness$pred %>%
  left_join(dungeness_covariates_chinook0_long, by = c("doy","year", "daynight_category", "trap"))

ggplot(data = predict_chinook0_dungeness$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"))+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_chinook0_dungeness$pred, aes(x = doy, y = log(chinook0_hatchery_perhour_interpolate+1), 
                                                           color = "hatchery")) +
  facet_wrap(~year+daynight_category, ncol = 5, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(Chinook salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", hatchery = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  scale_y_continuous(breaks = c(-3,0,3))+
  scale_x_continuous(limit = c(130, 200), breaks = c(140,160,180))+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  

ggsave(here("output","chinook0_dungeness_prediction_w_hatchery.png"), width = 6, height = 8, units = "in", dpi = 300)  


#coho

predict_coho1_dungeness <- predict(fit_coho_dungeness, type = "ytT", interval = "confidence")

glimpse(predict_coho1_dungeness$pred)  

head(predict_coho1_dungeness$pred)  

predict_coho1_dungeness$pred$trap <- 'screw'  

predict_coho1_dungeness$pred$year <-  as.numeric(substr(predict_coho1_dungeness$pred$.rownames, 
                                               1, 4))  

predict_coho1_dungeness$pred$daynight_category <- ifelse(substr(predict_coho1_dungeness$pred$.rownames,
                                               6, nchar(predict_coho1_dungeness$pred$.rownames)) == 'Day', 'day', 'night')

predict_coho1_dungeness$pred$doy <- predict_coho1_dungeness$pred$t+120  

dungeness_covariates_coho1 <- as.data.frame(t(covariates_coho1))  

dungeness_covariates_coho1$doy <- as.numeric(rownames(dungeness_covariates_coho1))  

dungeness_covariates_coho1_long <-  dungeness_covariates_coho1 %>%  
  select(doy, starts_with("coho")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","year","daynight_category"),
               names_pattern = "(.*)_(.{4})_(.*)") %>%
  mutate(year = as.numeric(year), trap = 'screw', daynight_category = ifelse(daynight_category == 'Day', 'day', 'night'))

predict_coho1_dungeness$pred <- predict_coho1_dungeness$pred %>%
  left_join(dungeness_covariates_coho1_long, by = c("doy","year", "daynight_category", "trap"))

ggplot(data = predict_coho1_dungeness$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"))+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_coho1_dungeness$pred, aes(x = doy, y = log(coho1_hatchery_perhour_interpolate+1), 
                                                           color = "hatchery")) +
  facet_wrap(~year+daynight_category, ncol = 6, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(Coho salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", hatchery = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  scale_y_continuous(breaks = c(-3,0,3))+
  scale_x_continuous(limit = c(120, 160), breaks = c(120,140,160))+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave(here("output","coho1_dungeness_prediction_w_hatchery.png"), width = 6, height = 8, units = "in", dpi = 300)

  
  
  
  
  #############

# Puyallup new chinook model with temp

### puyallup

predict_chinook0_puyallup_new <- predict(fits.hatchery.all.years[[121]], type = "ytT", interval = "confidence")

glimpse(predict_chinook0_puyallup_new$pred)

head(predict_chinook0_puyallup_new$pred)

predict_chinook0_puyallup_new$pred$trap <- 'screw'

#last 4 characters of rowname is year, but length of rowname varies depending on whether it's a day or night trap
predict_chinook0_puyallup_new$pred$year <-  as.numeric(substr(predict_chinook0_puyallup_new$pred$.rownames, 
                                                          nchar(predict_chinook0_puyallup_new$pred$.rownames)-3, 
                                                          nchar(predict_chinook0_puyallup_new$pred$.rownames)))

#if rowname has 'day' in it, then it's a day trap, otherwise it's a night trap
predict_chinook0_puyallup_new$pred$daynight_category <- ifelse(grep("day", 
                                                                predict_chinook0_puyallup_new$pred$.rownames, 
                                                                value = TRUE) == predict_chinook0_puyallup_new$pred$.rownames, 
                                                           "day", "night")
predict_chinook0_puyallup_new$pred$doy <- predict_chinook0_puyallup_new$pred$t+130
puyallup_covariates_chinook0_new <- as.data.frame(t(covariates_chinook0_puyallup_w_temp))

puyallup_covariates_chinook0_new$doy <- as.numeric(rownames(puyallup_covariates_chinook0_new))

puyallup_covariates_chinook0_long_new <-  puyallup_covariates_chinook0_new %>% 
  select(doy, starts_with("chinook")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","daynight_category","year"),
               names_pattern = "(.*)_(.*)_(.{4})") %>%
  mutate(year = as.numeric(year), trap = 'screw')

predict_chinook0_puyallup_new$pred <- predict_chinook0_puyallup_new$pred %>%
  left_join(puyallup_covariates_chinook0_long_new, by = c("doy","year", "daynight_category", "trap"))


ggplot(data = predict_chinook0_puyallup_new$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"))+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_chinook0_puyallup_new$pred, aes(x = doy, y = log(chinook0_hatchery_perhour+1), 
                                                       color = "hatchery")) +
  facet_wrap(~year+daynight_category, ncol = 6, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(Chinook salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", hatchery = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  scale_y_continuous(breaks = c(-3,0,3))+
  scale_x_continuous(limit = c(130, 200), breaks = c(140,160,180))+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave(here("output","chinook0_puyallup_prediction_w_hatchery_new.png"), width = 6, height = 8, units = "in", dpi = 300)

  
  
  
  
#puyallup coho model


#load data
data <- read.csv(here("puyallup", "data","puyallup_final.csv"),
                 na.strings = c("NA",""))

data$coho1_hatchery_perhour_night[dim(data)[1]] <- 0
data$coho1_hatchery_perhour_day[dim(data)[1]] <- 0


data <- data %>% 
  mutate(coho1_hatchery_perhour_night = na.approx(coho1_hatchery_perhour_night, na.rm = FALSE),
         coho1_hatchery_perhour_day = na.approx(coho1_hatchery_perhour_day, na.rm = FALSE))



covariates_coho1_puyallup_w_temp <- arrange(data,doy) %>%
  filter(doy > 90 & doy <= 160) %>%
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

c<-NULL
for(kk in c(1,3,6,7,11)){
  c = rbind(c,covariates_coho1_puyallup_w_temp[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_coho1_puyallup_w_temp)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  
}
fit.model = c(list(c= c), mod_list(nyears,5,1, FALSE, TRUE))
fit_coho_puyallup_best <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

predict_coho_puyallup <- predict(fit_coho_puyallup_best, type = "ytT", interval = "confidence")

glimpse(predict_coho_puyallup$pred)


predict_coho_puyallup$pred$trap <- 'screw'

#last 4 characters of rowname is year, but length of rowname varies depending on whether it's a day or night trap
predict_coho_puyallup$pred$year <-  as.numeric(substr(predict_coho_puyallup$pred$.rownames, 
                                                              nchar(predict_coho_puyallup$pred$.rownames)-3, 
                                                              nchar(predict_coho_puyallup$pred$.rownames)))

#if rowname has 'day' in it, then it's a day trap, otherwise it's a night trap
predict_coho_puyallup$pred$daynight_category <- ifelse(grep("day", 
                                                                    predict_coho_puyallup$pred$.rownames, 
                                                                    value = TRUE) == predict_coho_puyallup$pred$.rownames, 
                                                               "day", "night")
predict_coho_puyallup$pred$doy <- predict_coho_puyallup$pred$t+90
puyallup_covariates_coho <- as.data.frame(t(covariates_coho1_puyallup_w_temp))

puyallup_covariates_coho$doy <- as.numeric(rownames(puyallup_covariates_coho))

puyallup_covariates_coho_long <-  puyallup_covariates_coho %>% 
  select(doy, starts_with("coho")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","daynight_category","year"),
               names_pattern = "(.*)_(.*)_(.{4})") %>%
  mutate(year = as.numeric(year), trap = 'screw')

predict_coho_puyallup$pred <- predict_coho_puyallup$pred %>%
  left_join(puyallup_covariates_coho_long, by = c("doy","year", "daynight_category", "trap"))



ggplot(data = predict_coho_puyallup$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"))+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_coho_puyallup$pred, aes(x = doy, y = log(coho1_hatchery_perhour+1), 
                                                           color = "hatchery")) +
  facet_wrap(~year+daynight_category, ncol = 6, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(coho salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", hatchery = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  # scale_y_continuous(breaks = c(-3,0,3))+
  scale_x_continuous(limit = c(90, 160), breaks = c(100,120,140))+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave(here("output","coho1_puyallup_prediction_w_hatchery.png"), width = 6, height = 8, units = "in", dpi = 300)





  
  
  