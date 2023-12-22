#find peaks of hatchery chinook subyearlings in the Puyallup data
#plot the hatchery per hour and the wild perhour in separate panels in the same plot
#facet wrap by year

#load libraries
library(pracma)
library(tidyverse)
library(here)
library(ggplot2)
library(cowplot)

data <- read.csv(here("puyallup", "data","puyallup_2004-2021_all_days_w_covariates_edited.csv"),
                 na.strings = c("NA",""))
data$chinook0_hatchery_perhour_inp <- na.approx(data$chinook0_hatchery_perhour, 
                                                          na.rm = FALSE)
#convert date to as.Date

data$date <- as.Date(data$Date, format = "%m/%d/%Y")
data$year <- year(data$date)

chinook <- data %>%
  group_by(year, doy) %>%
  filter(doy >= 150 & doy <= 218) %>%
  summarise(chinook0_hatchery_num = sum(chinook0_hatchery_num, na.rm = TRUE),
            chinook0_wild_num = sum(chinook0_wild_num, na.rm = TRUE))


#finding peaks


peaks <- findpeaks(chinook$chinook0_hatchery_num,
                   threshold = 100, #minimum height of peak
                   minpeakdistance = 50 #minimum distance between peaks
)

#plotting the peaks

plot(chinook$chinook0_hatchery_num, type = "l")
points(peaks[,2], peaks[,1], col = "red", pch = 16)

#make the data dataframe long by adding a origin column which takes 'wild'
#and 'hatchery' as values and a column for cpue which will take values of 
#chinook0_wild_perhour and chinook0_hatchery_perhour

data_long1 <- data %>%
  select(year, doy, chinook0_hatchery_perhour, chinook0_wild_perhour) %>%
  filter(doy >= 130 & doy <= 200 & year <= 2009) %>%
  pivot_longer(cols = c(chinook0_hatchery_perhour, chinook0_wild_perhour),
               names_to = "origin", values_to = "cpue", names_pattern = "chinook0_(.*)_perhour")
data_long2 <- data %>%
  select(year, doy, chinook0_hatchery_perhour, chinook0_wild_perhour) %>%
  filter(doy >= 130 & doy <= 200 & year > 2009 & year <=2015) %>%
  pivot_longer(cols = c(chinook0_hatchery_perhour, chinook0_wild_perhour),
               names_to = "origin", values_to = "cpue", names_pattern = "chinook0_(.*)_perhour")

data_long3 <- data %>%
  select(year, doy, chinook0_hatchery_perhour, chinook0_wild_perhour) %>%
  filter(doy >= 130 & doy <= 200 & year > 2015) %>%
  pivot_longer(cols = c(chinook0_hatchery_perhour, chinook0_wild_perhour),
               names_to = "origin", values_to = "cpue", names_pattern = "chinook0_(.*)_perhour")


#remove legend and combine facet labels
p1 <- ggplot(data_long1, aes(x = doy)) +
  xlim(130, 200) +
  geom_line(aes(y = cpue, color = origin), alpha = 0.5,
            linewidth = 1) +
  facet_wrap(~origin+year, scales = "free_y", nrow = 2, labeller = label_wrap_gen(multi_line=FALSE)) +
  theme(legend.position = "none") +
  labs(x = "", y = "") +
  theme_bw()+
  scale_color_manual(values = c("cadetblue", "salmon")) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        strip.background = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none") 

p2 <- ggplot(data_long2, aes(x = doy)) +
  xlim(130, 200) +
  geom_line(aes(y = cpue, color = origin), alpha = 0.5,
            linewidth = 1) +
  facet_wrap(~origin+year, scales = "free_y", nrow = 2, labeller = label_wrap_gen(multi_line=FALSE)) +
  theme(legend.position = "none") +
  labs(x = "", y = "Number of chinook subyearlings\n caught perhour") +
  theme_bw()+
  scale_color_manual(values = c("cadetblue", "salmon")) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        strip.background = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none")

p3 <- ggplot(data_long3, aes(x = doy)) +
  xlim(130, 200) +
  geom_line(aes(y = cpue, color = origin), alpha = 0.5,
            linewidth = 1) +
  facet_wrap(~origin+year, scales = "free_y", nrow = 2, labeller = label_wrap_gen(multi_line=FALSE)) +
  theme(legend.position = "none") +
  labs(x = "Day of Year", y = "") +
  theme_bw()+
  scale_color_manual(values = c("cadetblue", "salmon")) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        strip.background = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none")

plot_grid(
  p1, p2, p3, ncol = 1
)

#save the plot
ggsave(here("puyallup","output",
            "chinook0_hatchery_wild.png"), width = 12, height = 10, units = "in", dpi = 300)
