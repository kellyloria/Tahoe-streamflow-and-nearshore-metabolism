##==============================================================================
## Analysis of inter-annual variation WY 22 to WY23
## by Loria et al. 2024
## 09/24/2024
#===============================================================================
library(tidyverse)
library(lubridate)

# plotting packages:
library(ggplot2)
library(ggpubr)
library(ggpattern)

##===========================================
## read data aggregated data for the project:
#============================================
# Temp path setwd("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/Final_Scripts/") 
# NS_analysis_dat.rds"

dat <- readRDS("./NS_analysis_dat.rds")
str(dat)

##===========================================
## average data by shore to account for as many observations as possible:
dat_shore <- dat %>%
  arrange(shore, yday, date, month, WaterYear) %>%
  group_by(shore, yday, date, month, WaterYear) %>%
  summarise(
    GPP_m = mean(middle_GPP, na.rm=T),
    ER_m = mean(c(middle_ER * -1), na.rm=T),
    NEP_m = mean(middle_NEP, na.rm=T),
    temp_m = mean(lake_tempC, na.rm=T),
    LakeDO_m = mean(lake_DO, na.rm=T),
    Kd_m = mean(Kd_fill, na.rm=T),
    par_3m = mean(light_mean, na.rm=T),
    ppt_m = mean(ppt_mm, na.rm=T),
    windsp_m = mean(windsp_mean, na.rm=T)
  ) 

## optional infill ##
dat_fill <- dat_shore %>%
  arrange(shore, WaterYear, yday) %>%
  group_by(shore, WaterYear) %>%
  fill(GPP_m, .direction = "down")%>%
  fill(ER_m, .direction = "down")%>%
  fill(temp_m, .direction = "down")%>%
  fill(LakeDO_m, .direction = "down")%>%
  fill(Kd_m, .direction = "down")%>%
  fill(par_3m, .direction = "down")%>%
  fill(ppt_m, .direction = "down")%>%
  fill(windsp_m, .direction = "down")

##===========================================
# Perform ANOVA for inter annual variation in GPP and ER
##===========================================

dat_23 <- dat_fill%>%
  filter(date>as.Date("2023-02-13"))

# Perform ANOVA for GPP
anova_gpp <- aov(GPP_m ~ shore, data = dat_23)
summary(anova_gpp)

# Post-hoc test for GPP
posthoc_gpp <- TukeyHSD(anova_gpp)
posthoc_gpp

# Perform ANOVA for ER
anova_er <- aov(ER_m ~ shore, data = dat_23)
summary(anova_er)

# Post-hoc test for ER
posthoc_er <- TukeyHSD(anova_er)
posthoc_er

# Perform ANOVA for NEP
anova_NEP <- aov(NEP_m ~ shore, data = dat_23)
summary(anova_NEP)

# Post-hoc test for ER
posthoc_NEP <- TukeyHSD(anova_NEP)
posthoc_NEP

### Metabolism anova 
plot_NEP <- ggplot(dat_23, aes(x = interaction(shore), y = (NEP_m), fill = shore)) +
  geom_jitter(aes(color = shore), width = 0.1, height = 0.01, alpha = 0.35) +  # Adding jittered points with 50% transparency
  geom_boxplot(alpha = 0.9) +  # Pattern angle
  geom_hline(yintercept = 0) +
  #scale_pattern_manual(values = c("2023" = "stripe", "2022" = "none")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  scale_x_discrete(labels = (c("BW", "GB", "SH", "SS"))) +
  labs(y = expression(NEP~(mmol~O[2]~m^-3~d^-1)), x = NULL) 

# ggsave(plot = plot_NEP, filename = paste("./Figures/ANOVA_NEP.png",sep=""),width=3.75,height=5,dpi=300)

###################
## Create WY-days based on data range: 
dat_WYvar22 <- dat_fill%>%
  filter(date>as.Date("2021-09-30") & date<as.Date("2022-09-10")) %>%
  mutate(
    WY_doy = as.numeric(difftime(date, as.Date("2021-09-30"), units = "days")) + 1)

dat_WYvar23 <- dat_fill%>%
  filter(date>as.Date("2022-09-30") & date<as.Date("2023-09-10")) %>%
  mutate(
  WY_doy = as.numeric(difftime(date, as.Date("2022-09-30"), units = "days")) + 1)

dat_WYvar<- rbind(dat_WYvar22, dat_WYvar23)
### 

plot_NEPWY <- ggplot(dat_WYvar, aes(x = interaction(WaterYear, shore), y = (NEP_m), pattern = as.factor(WaterYear), fill = shore)) +
  geom_jitter(aes(color = shore), width = 0.1, height = 0.01, alpha = 0.25) +  # Adding jittered points with 50% transparency
  geom_boxplot_pattern(alpha = 0.9, outlier.alpha = 0.5,
                       pattern_density = 0.25,  
                       pattern_spacing = 0.05,  
                       pattern_fill = "black",  
                       pattern_angle = 45) +  
  geom_hline(yintercept = 0) +
  scale_pattern_manual(values = c("2023" = "stripe", "2022" = "none")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() + #ylim(0,25)+
  theme(legend.position = "bottom") +
  scale_x_discrete(labels = (c("BW", "BW", "GB", "GB", "SH", "SS", "SS"))) +
  labs(y = expression("NEP"~(mmol~O[2]~m^-3~d^-1)), x = NULL)

# ggsave(plot = plot_NEPWY, filename = paste("./Figures/ANOVA_NEP_WYv2.png",sep=""),width=7,height= 5,dpi=300)

###################
## Stat summary: 
dat_WY_S <- dat_WYvar %>%
filter(WaterYear=="2023" & shore=="SS")
mean(na.omit(dat_WY_S$NEP_m))

anova_NEP <- aov(NEP_m ~ as.factor(WaterYear), data = dat_WYvar%>%filter(shore=="SS"))
summary(anova_NEP)


##===========================================
# Perform ANOVA for inter annual variation in environmental conditions 
##===========================================
dat_WYvar$WY <- as.factor(dat_WYvar$WaterYear)

# Perform ANOVA for temp
# sub out sites (BW, GB, SS)

anova_P <- aov(ppt_m ~ as.factor(WaterYear), data =dat_WYvar%>%filter(WY_doy<240, (shore == "SS")))
summary(anova_P)

anova_T <- aov(temp_m ~ as.factor(WaterYear), data = dat_WYvar%>%filter(WY_doy<240, (shore == "GB")))
summary(anova_T)

anova_K <- aov(Kd_m ~ as.factor(WaterYear), data = dat_WYvar%>%filter(WY_doy<240, (shore == "GB")))
summary(anova_K) 

anova_PAR <- aov(par_3m ~ as.factor(WaterYear), data = dat_WYvar%>%filter(WY_doy<240, (shore == "BW")))
summary(anova_PAR) 

anova_W <- aov(windsp_m ~ as.factor(WaterYear), data = dat_WYvar%>%filter(WY_doy<240, (shore == "GB")))
summary(anova_W)

## Stat summary: 
dat_WY_S <- dat_WYvar %>%
  filter(WY_doy<240, WaterYear=="2023" & shore=="BW")
mean(na.omit(dat_WY_S$temp_m))

## Visualization:
plot_Kd <- ggplot(dat_WYvar%>% filter(WY_doy<240,shore!="SH"), aes(x = interaction(WY, shore), y = Kd_m, fill = shore, pattern=WY)) +
  geom_jitter(aes(color = shore), width = 0.1, height = 0.01, alpha = 0.25) +  # Adding jittered points with 50% transparency
  geom_boxplot(alpha = 0.5, outlier.alpha = 0.5) +
  geom_boxplot_pattern(alpha = 0.5, outlier.alpha = 0.5,
                       pattern_density = 0.15,  # Adjust pattern density
                       pattern_spacing = 0.05,  # Adjust pattern spacing
                       pattern_fill = "black",  # Pattern fill color
                       pattern_angle = 45) +  # Pattern angle
  scale_pattern_manual(values = c("2023" = "stripe", "2022" = "none")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_discrete(labels = (c("2022", "2023", "2022", "2023", "2022", "2023"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        plot.subtitle = element_text(size = 14))+
  labs(y = expression(K[d]~(m^-1)), x = NULL) 

plot_Kd <- plot_Kd +
  geom_text(data = data.frame(
    shore = c("BW", "GB"),
    WY = factor(rep("2022", 2)), # Adjust as needed for your specific x-axis values
    y = c(0.33) # Adjust y-values to be above the box plots
  ),
  aes(x = interaction(WY, shore), y = y, label = "*"),
  color = "black", size = 10,
  )

plot_par <- ggplot(dat_WYvar%>%filter(WY_doy<240, shore!="SH"), aes(x = interaction(WY, shore), y = par_3m, fill = shore, pattern=WY)) +
  geom_jitter(aes(color = shore), width = 0.1, height = 0.01, alpha = 0.25) +  # Adding jittered points with 50% transparency
  geom_boxplot(alpha = 0.5, outlier.alpha = 0.25) +
  geom_boxplot_pattern(alpha = 0.5, outlier.alpha = 0.5,
                       pattern_density = 0.15,  # Adjust pattern density
                       pattern_spacing = 0.05,  # Adjust pattern spacing
                       pattern_fill = "black",  # Pattern fill color
                       pattern_angle = 45) +  # Pattern angle
  scale_pattern_manual(values = c("2023" = "stripe", "2022" = "none")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_discrete(labels = (c("2022", "2023", "2022", "2023", "2022", "2023"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        plot.subtitle = element_text(size = 14))+
  labs(y = expression(PAR~(mu~mol~m^-2~s^-1)), x = NULL) 

plot_par <- plot_par +
  geom_text(data = data.frame(
    shore = c("SS"),
    WY = factor(rep("2022", 2)), # Adjust as needed for your specific x-axis values
    y = c(365) # Adjust y-values to be above the box plots
  ),
  aes(x = interaction(WY, shore), y = y, label = "*"),
  color = "black", size = 10,
  )


plot_temp <- ggplot(dat_WYvar%>%filter(WY_doy<240, shore!="SH"), aes(x = interaction(WY, shore), y = temp_m, fill = shore, pattern=WY)) +
  geom_jitter(aes(color = shore), width = 0.1, height = 0.01, alpha = 0.25) +  # Adding jittered points with 50% transparency
  geom_boxplot(alpha = 0.5, outlier.alpha = 0.5) +
  geom_boxplot_pattern(alpha = 0.5, outlier.alpha = 0.5,
                       pattern_density = 0.15,  # Adjust pattern density
                       pattern_spacing = 0.05,  # Adjust pattern spacing
                       pattern_fill = "black",  # Pattern fill color
                       pattern_angle = 45) +  # Pattern angle
  scale_pattern_manual(values = c("2023" = "stripe", "2022" = "none")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_discrete(labels = (c("2022", "2023", "2022", "2023", "2022", "2023"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        plot.subtitle = element_text(size = 14))+
  labs(y= "Water temp. (°C)", x=NULL) 
plot_temp

plot_temp <- plot_temp +
  geom_text(data = data.frame(
    shore = c("GB", "SS"),
    WY = factor(rep("2022", 2)), # Adjust as needed for your specific x-axis values
    y = c(16.5) # Adjust y-values to be above the box plots
  ),
  aes(x = interaction(WY, shore), y = y, label = "*"),
  color = "black", size = 10,
  )


plot_wind <- ggplot(dat_WYvar%>%filter(WY_doy<240, shore!="SH"), aes(x = interaction(WY, shore), y = windsp_m, fill = shore, pattern=WY)) +
  geom_jitter(aes(color = shore), width = 0.1, height = 0.01, alpha = 0.25) +  # Adding jittered points with 50% transparency
  geom_boxplot(alpha = 0.5, outlier.alpha = 0.5) +
  geom_boxplot_pattern(alpha = 0.5, outlier.alpha = 0.5,
                       pattern_density = 0.15,  # Adjust pattern density
                       pattern_spacing = 0.05,  # Adjust pattern spacing
                       pattern_fill = "black",  # Pattern fill color
                       pattern_angle = 45) +  # Pattern angle
  scale_pattern_manual(values = c("2023" = "stripe", "2022" = "none")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_discrete(labels = (c("2022", "2023", "2022", "2023", "2022", "2023"))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        plot.subtitle = element_text(size = 14))+
  labs(y = expression(Wind~(m~s^-1)), x = NULL) 

plot_wind


plot_precip <- ggplot(dat_WYvar%>%filter(WY_doy<240, shore!="SH"), aes(x = interaction(WY, shore), y = log(ppt_m+1), fill = shore, pattern=WY)) +
  geom_jitter(aes(color = shore), width = 0.1, height = 0.01, alpha = 0.25) +  # Adding jittered points with 50% transparency
  geom_boxplot(alpha = 0.5, outlier.alpha = 0.5) +
  geom_boxplot_pattern(alpha = 0.5, outlier.alpha = 0.5,
                       pattern_density = 0.15,  # Adjust pattern density
                       pattern_spacing = 0.05,  # Adjust pattern spacing
                       pattern_fill = "black",  # Pattern fill color
                       pattern_angle = 45) +  # Pattern angle
  scale_pattern_manual(values = c("2023" = "stripe", "2022" = "none")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_discrete(labels = (c("2022", "2023", "2022", "2023", "2022", "2023"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        plot.subtitle = element_text(size = 14))+
  labs(y = expression(log(Precip.+1)~(mm)), x = NULL) 

plot_precip


plot_precip <- plot_precip +
  geom_text(data = data.frame(
      shore = c("BW", "GB", "SS"),
      WY = factor(rep("2022", 3)), # Adjust as needed for your specific x-axis values
      y = c(3.2) # Adjust y-values to be above the box plots
    ),
    aes(x = interaction(WY, shore), y = y, label = "*"),
    color = "black", size = 10,
  )

box_grid <- ggarrange(plot_precip,
                     plot_temp,
                     plot_par,
                     plot_wind,
                     ncol = 4, nrow = 1,
                     common.legend = TRUE, 
                     labels=c("a", "b","c", "d", "e"),
                     legend = "bottom")

# ggsave(plot = box_grid, filename = paste("./Environ_condition.png",sep=""),width=16,height=4.25,dpi=300)


##===========================================
# Create cumulative totals for each site and year based on DOY, keeping yday as an index
##===========================================
cumul_dat <- dat_WYvar %>%
  arrange(shore, WaterYear, WY_doy) %>%
  group_by(shore, WaterYear) %>%  
  mutate(
    GPP_cumulative = cumsum(GPP_m),
    ER_cumulative = cumsum(ER_m),
    NEP_cumulative = cumsum(NEP_m),
    temp_cumulative = cumsum(temp_m),
    LakeDO_cumulative = cumsum(LakeDO_m),
    par_cumulative = cumsum(par_3m),
    kd_cumulative = cumsum(Kd_m),
    ppt_cumulative = cumsum(ppt_m),
    windsp_cumulative = cumsum(windsp_m)
  ) %>%
  ungroup()

temp <- cumul_dat %>%filter(WY_doy<240)

plot_temp_c <- 
  ggplot(cumul_dat%>%filter(WY_doy<240, !shore == "SH" ), aes(x=WY_doy, y=temp_cumulative,  lty=as.factor(WaterYear), colour = shore)) + 
  theme_bw() + 
  geom_line(alpha=0.75, size=1) +
  labs(y="Water temperature (°C)", x=NULL) + 
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        plot.subtitle = element_text(size = 14))+
  theme(legend.position = "bottom") 
plot_temp_c


plot_kd_c <- 
  ggplot(cumul_dat%>%filter(WY_doy<240, !shore == "SH"), aes(x=WY_doy, y=kd_cumulative,  lty=as.factor(WaterYear), colour = shore)) + 
  theme_bw() + 
  geom_line(alpha=0.75, size=1) +
  labs(y=expression(Kd~(m^-1)), x= NULL) + 
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        plot.subtitle = element_text(size = 14))+
  theme(legend.position = "bottom") 
plot_kd_c

plot_par_c <- 
  ggplot(cumul_dat%>%filter(WY_doy<240, !shore == "SH"), aes(x=WY_doy, y=par_cumulative,  lty=as.factor(WaterYear), colour = shore)) + 
  theme_bw() + 
  geom_line(alpha=0.75, size=1) +
  #labs(y=expression(PAR~(m^-1)), x= NULL) + 
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        plot.subtitle = element_text(size = 14))+
  theme(legend.position = "bottom") + ylim(0,30000) +
  labs(y = expression(PAR~(mu~mol~m^-2~s^-1)), x = NULL) 
plot_par_c


plot_ppm_c <- 
  ggplot(cumul_dat%>%filter(WY_doy<240, !(shore == "SH")), aes(x=WY_doy, y=ppt_cumulative,  lty=as.factor(WaterYear), colour = shore)) + 
  theme_bw() + 
  geom_line(alpha=0.75, size=1) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        plot.subtitle = element_text(size = 14))+
  theme(legend.position = "bottom") +
  labs(y = expression(log(Precip.+1)~(mm)), x = NULL) 
plot_ppm_c


plot_wind_c <- 
  ggplot(cumul_dat%>%filter(WY_doy<240, !(shore == "SH")), aes(x=WY_doy, y=windsp_cumulative,  lty=as.factor(WaterYear), colour = shore)) + 
  theme_bw() + 
  geom_line(alpha=0.75, size=1) +
  labs(y=expression(Wind~(m~s^-1)), x= "DOY") + 
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        plot.subtitle = element_text(size = 14))+
  theme(legend.position = "bottom") 
plot_wind_c

covariates24 <- ggarrange(plot_ppm_c,
                          plot_precip,
                          plot_temp_c,
                          plot_temp,
                          plot_kd_c,
                          plot_Kd,
                          plot_par_c,
                          plot_par,
                          plot_wind_c,
                          plot_wind,
                          ncol = 2, nrow = 5,
                          widths = c(0.7, 0.25),
                          common.legend = TRUE, 
                          labels=c("a", "","b", "","c", "","d", "","e",""),
                          legend = "bottom")

# ggsave(plot = covariates24, filename = paste("./Figures/Cumulative_Summer_covarsv3.png",sep=""),width=12,height=14,dpi=300)

##===========================================
# seasonal NEP plots 
##===========================================
### lets look at seasonal differences:

# winter december 1st to march 1st 
dat_WYvar22 <- dat_WYvar%>%
  filter(date>as.Date("2021-12-01") & date<as.Date("2022-03-01")) 
dat_WYvar22$season <- "Winter"
dat_WYvar23 <- dat_WYvar%>%
  filter(date>as.Date("2022-12-01") & date<as.Date("2023-03-01"))
dat_WYvar23$season <- "Winter"
# spring march 1st - to june 1st 
dat_WYvar22sp <- dat_WYvar%>%
  filter(date>as.Date("2022-03-01") & date<as.Date("2022-06-01")) 
dat_WYvar22sp$season <- "Spring"
dat_WYvar23sp <- dat_WYvar%>%
  filter(date>as.Date("2023-03-01") & date<as.Date("2023-06-01"))
dat_WYvar23sp$season <- "Spring"
# summer 1st June- September 10th
dat_WYvar22su <- dat_WYvar%>%
  filter(date>as.Date("2022-06-01") & date<as.Date("2022-09-01")) 
dat_WYvar22su$season <- "Summer"
dat_WYvar23su <- dat_WYvar%>%
  filter(date>as.Date("2023-06-01") & date<as.Date("2023-09-01"))
dat_WYvar23su$season <- "Summer"


season_dat_WY <- rbind(dat_WYvar22, dat_WYvar23,
                         dat_WYvar22sp, dat_WYvar23sp,
                         dat_WYvar22su, dat_WYvar23su)


anova_tep <- aov(NEP_m ~ WY, data = season_dat_WY%>%
                   filter(shore=="SS"& 
                            season=="Summer"))
summary(anova_tep)

SFS_dat_seas_S <- season_dat_WY %>%
  filter(WaterYear=="2022" & 
           season=="Summer" &
           shore=="BW")
mean(na.omit(SFS_dat_seas_S$NEP_m))


### Metabolism anova 
season_dat_WY <- season_dat_WY %>%
  mutate(season = factor(season, levels = c("Winter", "Spring", "Summer")))

plot_NEP <- ggplot(season_dat_WY, aes(x = interaction(WY, shore), y = (NEP_m), pattern = WY, fill = shore)) +
  geom_jitter(aes(color = shore), width = 0.1, height = 0.01, alpha = 0.25) +  # Adding jittered points with 50% transparency
  geom_boxplot_pattern(alpha = 0.7, outlier.alpha = 0.5,
                       pattern_density = 0.20,  # Adjust pattern density
                       pattern_spacing = 0.05,  # Adjust pattern spacing
                       pattern_fill = "black",  # Pattern fill color
                       pattern_angle = 45) +  # Pattern angle
  geom_hline(yintercept = 0) +
  scale_pattern_manual(values = c("2023" = "stripe", "2022" = "none")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  scale_x_discrete(labels = (c("BW", "BW", "GB", "GB", "SH", "SS", "SS"))) +
  labs(y = expression(NEP~(mmol~O[2]~m^-3~d^-1)), x = NULL)  + facet_grid(.~season)

plot_NEP <- plot_NEP +
  geom_text(data = data.frame(
    shore = c("SS"),
    season = factor(c("Winter"), levels = c("Winter", "Spring", "Summer")),  # Ensure correct factor levels
    WY = factor(rep("2022")), # Adjust as needed for your specific x-axis values
    y = c(15) # Adjust y-values to be above the box plots
  ),
  aes(x = interaction(WY, shore), y = y, label = "*"),
  color = "black", size = 10,
  )


plot_NEP1 <- plot_NEP +
  geom_text(data = data.frame(
    shore = c("BW", "GB", "SS"),
    season = factor(c("Spring"), levels = c("Winter", "Spring", "Summer")),  # Ensure correct factor levels
    WY = factor(rep("2022")), # Adjust as needed for your specific x-axis values
    y = c(15) # Adjust y-values to be above the box plots
  ),
  aes(x = interaction(WY, shore), y = y, label = "*"),
  color = "black", size = 10,
  )

plot_NEP2 <- plot_NEP1 +
  geom_text(data = data.frame(
    shore = c("BW", "SS"),
    season = factor(c("Summer"), levels = c("Winter", "Spring", "Summer")),  # Ensure correct factor levels
    WY = factor(rep("2022")), 
    y = c(15) 
  ),
  aes(x = interaction(WY, shore), y = y, label = "*"),
  color = "black", size = 10,
  )

# ggsave(plot = plot_NEP2, filename = paste("./Figures/ANOVA_NEP_season.png",sep=""),width=10,height=4.25,dpi=300)

# end script.