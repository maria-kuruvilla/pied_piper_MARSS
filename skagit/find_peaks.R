#Goal - to use a package and find peaks in 
#the chinook0_hatchery_perhour column of the skagit data

#install.packages("pracma")
library(pracma)
library(tidyverse)
library(here)
# library(zoo)
# library(MASS)
# library(modelr)
# library(qpcR)
# library(GGally)
# library(tidyverse)
# library(broom)
library(ggplot2)



data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))

data_day_night$coho1_wild_perhour <- data_day_night$coho1_wild_num/data_day_night$In
data_day_night$coho1_hatchery_perhour <- data_day_night$coho1_hatchery_num/data_day_night$In

data_day_night$chinook0_wild_perhour <- data_day_night$chinook0_wild_num/data_day_night$In
data_day_night$chinook0_hatchery_perhour <- data_day_night$chinook0_hatchery_num/data_day_night$In

#interpolate chinook0_hatchery_perhour column of data_day_night

data_day_night$chinook0_hatchery_perhour_inp <- na.approx(data_day_night$chinook0_hatchery_perhour, 
                                                       na.rm = FALSE)

#first using tidyverse, sum chinook0_hatchery_num for day and night and for 
#scoop and screw

chinook <- data_day_night %>%
  group_by(year, doy) %>%
  filter(doy >= 150 & doy <= 200) %>%
  summarise(chinook0_hatchery_num = sum(chinook0_hatchery_num, na.rm = TRUE),
            chinook0_wild_num = sum(chinook0_wild_num, na.rm = TRUE))


#finding peaks in the chinook0_hatchery_perhour column of data_day_night

peaks <- findpeaks(chinook$chinook0_hatchery_num,
                   threshold = 100, #minimum height of peak
                   minpeakdistance = 50 #minimum distance between peaks
)

#plotting the peaks

plot(chinook$chinook0_hatchery_num, type = "l")
points(peaks[,2], peaks[,1], col = "red", pch = 16)


#using ggplot plot the chinook0_hatchery_num column of chinook and 
#the chinook0_wild_num column of chinook with different colurs
#facet wrap by year
#only plot the doy 25 days before peak and 25 days after peak

ggplot(chinook, aes(x = doy)) +
  xlim(150, 200) +
  geom_line(aes(y = chinook0_hatchery_num), color = "cadetblue", alpha = 0.5,
            linewidth = 1) +
  geom_line(aes(y = chinook0_wild_num), color = "salmon", alpha=0.5,
            linewidth = 1) +
  facet_wrap(~year, scales = "free_y") 

ggsave(here("skagit", "output", "chinook0_hatchery_wild_num.png"), width = 10, height = 10)

#do the same for coho

coho <- data_day_night %>%
  group_by(year, doy) %>%
  filter(doy >= 100 & doy <= 150) %>%
  summarise(coho1_hatchery_num = sum(coho1_hatchery_num, na.rm = TRUE),
            coho1_wild_num = sum(coho1_wild_num, na.rm = TRUE))

#add xlim of doy 100 to 150 
#remove gray background in the facet wrap title
#add legend for hatchery and wild
#show legend
#reduce number of y ticks
#increase size of legend text and x and y axis  labels
ggplot(coho, aes(x = doy)) +
  geom_line(aes(y = coho1_hatchery_num, colour = "Hatchery"), alpha = 0.5,
            linewidth = 1) +
  geom_line(aes(y = coho1_wild_num, colour = "Wild"), alpha=0.5,
            linewidth = 1) +
  facet_wrap(~year, scales = "free_y") +
  xlim(100, 150)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        title = element_text(size = 14, face = "bold")
        )+
  labs(title = "Coho Hatchery and Wild in Skagit", 
       x = "Day of Year", 
       y = "Number of Coho")+
  scale_color_manual(name = "", 
                     values = c("cadetblue", "salmon"),
                     labels = c("Hatchery", "Wild"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(here("skagit", "output", "coho1_hatchery_wild_num.png"), width = 16, height = 10)



#make the same plot for chinook

ggplot(chinook, aes(x = doy)) +
  xlim(150, 200)+
  facet_wrap(~year, scales = "free_y") +
  geom_line(aes(y = chinook0_hatchery_num, colour = "Hatchery"), alpha = 0.5,
            linewidth = 1) +
  geom_line(aes(y = chinook0_wild_num, colour = "Wild"), alpha=0.5,
            linewidth = 1) +
  
  
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        title = element_text(size = 14, face = "bold")
  )+
  labs(title = "Chinook Hatchery and Wild in Skagit", 
       x = "Day of Year", 
       y = "Number of Chinook")+
  scale_color_manual(name = "", 
                     values = c("cadetblue", "salmon"),
                     labels = c("Hatchery", "Wild"))

ggsave(here("skagit", "output", "chinook0_hatchery_wild_num.png"), width = 16, height = 10)

ggplot(data_day_night, aes(x = doy)) +
  geom_point(aes(y=chinook0_wild_num)) +
  facet_wrap(~year, scales = "free_y")

ggplot(data_day_night, aes(x = doy)) +
  geom_point(aes(y=chinook0_hatchery_num)) +
  facet_wrap(~year, scales = "free_y")

