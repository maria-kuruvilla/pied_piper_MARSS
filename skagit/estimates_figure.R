#Goal - to combine the figure with chinook estimates
#and the figure with coho estimates into one figure

#load libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)


ci_night
ci_skagit_chinook_w_hatchery
ci_skagit_chinook_wo_hatchery


#combine chinook wo hatchery, chinook w hatchery 
# and coho estimates into one figure

p1 <- ggplot(ci_skagit_chinook_wo_hatchery[c(27:29),], aes(x = c("Residuals",
                                                          "Photoperiod",
                                                          "Temperature\ndifference"),
                                                    y = estimate, 
                                                    ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Estimate of effect") +
  ggtitle("Skagit\nChinook salmon,\nwithout hatchery") + 
  theme_bw() +
  theme(plot.title = element_text(size = 8))+
  theme(axis.text.x=element_text(size=8),axis.title.y=element_text(size=8))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

p2 <- ggplot(ci_skagit_chinook_w_hatchery[c(27:30),], aes(x = c("Residuals",
                                                          "Photo-\nperiod",
                                                          "Temperature\ndifference",
                                                          "Hatchery"),
                                                    y = estimate, 
                                                    ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "") +
  ggtitle("\nChinook salmon,\nwith hatchery")+
  theme_bw() +
  theme(plot.title = element_text(size = 8))+
  theme(axis.text.x=element_text(size=8),axis.title.y=element_text(size=8))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))


p3 <- ggplot(ci_night[c(27:30),], aes(x = c("ATU","Flow","Lunar\nphase",
                                            "Hatchery"),
                                      y = estimate, 
                                      ymin = conf.low, ymax = conf.up)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "") +
  ggtitle("\nCoho salmon\n") + 
  theme_bw() +
  theme(plot.title = element_text(size = 8))+
  theme(axis.text.x=element_text(size=8),axis.title.y=element_text(size=8))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

#combine the three plots into one figure and title it Skagit River
#remove y axis title from the second and third plot
#move title to left
#width should be 7 inches and height should be 3 inches
ggarrange(p1, p2, p3, ncol = 3, nrow = 1) +
  ggtitle("Skagit River")

#save the figure
ggsave(here("skagit", "output","skagit_river_estimates.jpeg"), width = 7, height = 3, units = "in")
