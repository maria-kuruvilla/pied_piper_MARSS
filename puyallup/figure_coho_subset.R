
# Load libraries
library(tidyverse)
library(ggplot2)
library(broom)
library(here)

ci_coho1_puyallup <- tidy(fit_coho_puyallup_best)

coho1_puyallup_plot <- ggplot(ci_coho1_puyallup[c(40:45),], 
                              aes(x = c("Flow", "Photoperiod", 
                                        "ATU", "Flow\ndifference","Hatchery,\nday","Hatchery,\nnight"),
                                  y = estimate, ymin = conf.low, ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")  +
  geom_rect(aes(xmin = 2.5, xmax = 4.5, ymin = -0.28, ymax = 0.28), col = "cadetblue", alpha = 0.0) +
  geom_text(aes(x = 6.35, y = -0.225, label = "Puyallup"), size = 10)+
  labs(x = "", y = "Estimate of effect",
  )+
  theme_classic()+
  theme(axis.title.y=element_text(size=24),
        axis.title.x=element_text(size=24),
        axis.text.y = element_text(size = 24),
        axis.text.x=element_text(size=24)) +
  scale_y_continuous(breaks = seq(-0.2,0.2,0.2), limits = c(-0.28,0.28), n.breaks = 3)+
  scale_x_discrete(limits = c("ATU", "Flow", "Hatchery,\nnight","Hatchery,\nday",
                              "Photoperiod","Flow\ndifference"
  )) +
  coord_flip()
#guide = guide_axis(n.dodge = 2))

coho1_puyallup_plot