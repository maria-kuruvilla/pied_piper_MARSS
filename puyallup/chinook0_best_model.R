#first run chinook0_w_airtemp.R



out.tab.all.years <- NULL
fits.all.years <- NULL
nyears = num_years*2
c <- NULL
for(kk in c(1,4,5,6,7,10)){
  c = rbind(c,covariates_chinook0_puyallup_w_temp[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_chinook0_puyallup_w_temp)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  
}

fit.model = c(list(c= c), mod_list(nyears,6,1, FALSE, FALSE))
fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

ci_chinook0_puyallup_new <- tidy(fit)
chinook0_puyallup_plot_new <- ggplot(ci_chinook0_puyallup_new[c(39:45),], 
                                     aes(x = c("Flow",
                                               "Lunar phase",
                                               "Temperature",
                                               "Flow\n difference",
                                               "Photoperiod\n difference", 
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
