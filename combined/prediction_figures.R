#use the output from the code in estimate_figures.R
#use the MARSS output from the best/good models
#to make the prediction figures of the best/good models
#also make additional figure with the hatchery covariate

#load libraries
library(MARSS)
library(ggplot2)
library(tidyverse)
library(here)

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
  
  geom_line(aes(x = t+100, y = estimate, color = "wild predicted"))+
  geom_point(aes(x = t+100, y = y, color = "wild observed"), alpha = 0.5)+
  geom_ribbon(aes(x = t+100, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_coho1_skagit$pred, aes(x = doy, y = log(coho1_hatchery_perhour_inp+1), 
            color = "hatchery"), linewidth = 1, alpha = 0.5) +
  facet_wrap(~year+trap, ncol = 4, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log (coho salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild predicted" = "salmon", "wild observed" = "salmon", "hatchery" = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA, 1),
                       shape = c(NA, 19, NA))))+
  theme_classic()+
  theme(strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 20)
        )

ggsave(here("output","coho1_skagit_prediction_w_hatchery.png"), width = 12, height = 14, units = "in", dpi = 300)


