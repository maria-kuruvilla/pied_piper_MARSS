

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



ci_chinook0_puyallup_new <- tidy(fits.hatchery.all.years[[121]])


#plotting
#remove grid lines
#dodge the axis labels on the x axis
chinook0_puyallup_plot_new <- ggplot(ci_chinook0_puyallup_new[c(39:45),], 
                                 aes(x = c("Photoperiod\n difference", 
                                           "Flow", 
                                           "Temperature",
                                           "Lunar phase",
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
        axis.text.y = element_text(size = 14))

chinook0_puyallup_plot_new
# chinook0_puyallup_plot2 <- chinook0_puyallup_plot + theme(axis.title.y=element_blank(),
#                                                           axis.text.y = element_text(size = 18),
#                                                           axis.text.x=element_text(size=24)) +
#   scale_x_discrete(guide = guide_axis(n.dodge=3),limits = c("Flow\n difference","Photoperiod\n difference","Hatchery,\nday", 
#                                                             "Hatchery,\nnight", "Flow")) + 
#   scale_y_continuous(breaks = c(-0.25,0,0.25), limits = c(-0.32,0.32))
# 
# chinook0_puyallup_plot2

chinook0_puyallup_plot3 <- chinook0_puyallup_plot_new + theme(axis.title.y=element_text(size=24, family = "Sans", 
                                                                                    margin = margin(t = 10, r = 0, b = 0, l = 10)),
                                                          axis.title.x=element_text(size=24, family = "Sans", 
                                                                                    margin = margin(t = 15, r = 0, b = 0, l = 10)),
                                                          axis.text.y = element_text(size = 24, family = "Sans"),
                                                          axis.text.x=element_text(size=24, family = "Sans")) +
  scale_x_discrete(#guide = guide_axis(n.dodge=3),
    limits = c("Flow","Lunar phase",
               "Hatchery,\nnight","Hatchery,\nday",  "Photoperiod\n difference","Temperature",
               "Flow\n difference"
               
    )) + 
  scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1, 0.2), limits = c(-0.295,0.295))+
  coord_flip()

chinook0_puyallup_plot3

chinook0_dungeness_plot2 <- chinook0_dungeness_plot + theme(axis.title.y=element_text(size=28),
                                                            axis.title.x=element_text(size=28),
                                                            axis.text.y = element_text(size = 18),
                                                            axis.text.x=element_text(size=24)) +
  scale_x_discrete(expand = expansion(add = 1.1), #guide = guide_axis(n.dodge=3),
                   limits = c("Flow\n difference",
                              "Photoperiod\n difference", 
                              "Hatchery,\nday", 
                              "Hatchery,\nnight",
                              "Temperature\n difference")
  )+
  scale_y_continuous(breaks = c(-0.15,0,0.15), limits = c(-0.2,0.2))+
  coord_flip()
chinook0_dungeness_plot3 <- chinook0_dungeness_plot + theme(axis.title.y=element_text(size=24, family = "Sans"),
                                                            axis.title.x=element_text(size=24, family = "Sans", 
                                                                                      margin = margin(t = 15, r = 0, b = 0, l = 10)),
                                                            axis.text.y = element_text(size = 24, family = "Sans"),
                                                            axis.text.x=element_text(size=24, family = "Sans")) +
  scale_x_discrete(#expand = expansion(add = 1.1), #guide = guide_axis(n.dodge=3),
    limits = c("Temperature\n difference",
               "Hatchery,\nnight",
               "Hatchery,\nday", 
               
               "Photoperiod\n difference",
               "Flow\n difference"
    )
  )+
  scale_y_continuous(breaks = c(-0.15,0,0.15), limits = c(-0.2,0.2))+
  coord_flip()


chinook0_dungeness_plot3


chinook0_skagit_plot3 <- chinook0_skagit_plot + theme(axis.text.x=element_text(size=24, family = "Sans"),
                                                      axis.title.y=element_text(size=24, family = "Sans", 
                                                                                margin = margin(t = 10, r = 0, b = 0, l = 10)),
                                                      axis.title.x=element_text(size=24, family = "Sans", 
                                                                                margin = margin(t = 15, r = 20, b = 0, l = 0)),
                                                      axis.text.y = element_text(size = 24, family = "Sans"))+
  scale_y_continuous(breaks = c(-0.1,0,0.1,0.2), limits = c(-0.15,0.15))+
  scale_x_discrete(#expand = expansion(add = 1.15), #guide = guide_axis(n.dodge=3),
    limits = c("Temperature\n difference",
               "Photoperiod", 
               "Residuals",
               
               "Hatchery"
    )) + 
  coord_flip()

#draw rectangle around hatchery covariate
chinook0_skagit_plot4 <- chinook0_skagit_plot3 +
  geom_rect(aes(xmin = 3.5, xmax = 4.2, ymin = -0.15, ymax = 0.15), col = "cadetblue", alpha = 0.0, fill = "cadetblue")+
  geom_text(aes(x = 4.4, y = -0.12, label = "Skagit"), size = 10)

chinook0_puyallup_plot4 <- chinook0_puyallup_plot3+
  geom_rect(aes(xmin = 2.8, xmax = 4.2, ymin = -0.295, ymax = 0.295), col = "cadetblue", alpha = 0.0, fill = "cadetblue")+
  geom_text(aes(x = 7.4, y = -0.225, label = "Puyallup"), size = 10)

chinook0_dungeness_plot4 <- chinook0_dungeness_plot3+
  geom_rect(aes(xmin = 1.8, xmax = 3.2, ymin = -0.2, ymax = 0.2), col = "cadetblue", alpha = 0.0, fill = "cadetblue")+
  geom_text(aes(x = 5.4, y = -0.145, label = "Dungeness"), size = 10)

ggpubr::ggarrange(chinook0_dungeness_plot3,
                  chinook0_puyallup_plot3, 
                  chinook0_skagit_plot3,
                  labels = c("a", "b", "c"), ncol = 1, nrow = 3, font.label = list(size = 28, family = "Sans"),
                  common.legend = TRUE, legend = "right",
                  widths = c(1,1,1), heights = c(1.15,1.4,0.7))


ggsave(here("output","chinook0_covariates_estimates_new.jpeg"), width = 16, height = 15)

ggpubr::ggarrange(chinook0_dungeness_plot4,
                  chinook0_puyallup_plot4, 
                  chinook0_skagit_plot4,
                  labels = c("a", "b", "c"), ncol = 1, nrow = 3, font.label = list(size = 28, family = "Sans"),
                  common.legend = TRUE, legend = "right",
                  widths = c(1,1,1), heights = c(1.15,1.35,1.15))

ggsave(here("output","chinook0_covariates_estimates_new3.jpeg"), width = 16, height = 15)
