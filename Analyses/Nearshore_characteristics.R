##==============================================================================
## Analysis of spatial variation in nearshore characteristics 
## by Loria et al. 2024
## 09/24/2024
#===============================================================================
library(tidyverse)
library(lubridate)

# plotting packages:
library(ggplot2)
library(reshape2)
library(scales)
# stats packages 
library(PerformanceAnalytics)
library(lme4)
library(lmerTest)
library(MuMIn)
library(car)
library(glmmTMB)

se <- function(dat){
  se <- sd(dat)/sqrt(length(dat))
  return(se)}
##===========================================
## read data aggregated data for the project:
#============================================
# time-series data: 
datQ <- readRDS("./NS_analysis_dat.rds")
summary(datQ)
str(datQ)

# Chemistry data: 
# EDI link: 
chem_dat <- readRDS("./NS_chem_dat.rds")
unique(chem_dat$site)
summary(chem_dat) 

##===========================================
## weekly aggregation  
##===========================================

week_dat <- datQ %>%
  group_by(site, shore, week, WaterYear, position) %>%
  summarise(
    middle_GPP=mean(middle_GPP, na.rm = TRUE),
    lower_GPP=mean(lower_GPP, na.rm = TRUE),
    upper_GPP=mean(upper_GPP, na.rm = TRUE),
    middle_ER=mean(middle_ER, na.rm = TRUE),
    lower_ER=mean(lower_ER, na.rm = TRUE),
    upper_ER=mean(upper_ER, na.rm = TRUE),
    ## Weather 
    tmean_C=mean(tmean_C, na.rm = TRUE),
    light_mean=mean(light_mean, na.rm = TRUE),
    windsp_mean=mean(windsp_mean, na.rm = TRUE),
    ppt_mm=mean(ppt_mm, na.rm = TRUE),
    log_ppmt=mean(log_ppmt, na.rm = TRUE),
    ppmt_sum=sum(precip_bi, na.rm = TRUE),
    ## lake quality 
    lake_tempC=mean(lake_tempC, na.rm = TRUE),
    lake_DO=mean(lake_DO, na.rm = TRUE),
    lake_SPC=mean(lake_SPC, na.rm = TRUE),
    Kd_fill=mean(Kd_fill, na.rm = TRUE),
    par_int_3m=mean(par_int_3m, na.rm = TRUE),
    real_NS_depth=mean(real_NS_depth, na.rm = TRUE),
    ## Stream quality 
    flow_mean_m=mean(flow_mean, na.rm = TRUE),
    log_streamflow=mean(log_streamflow, na.rm = TRUE),
    flow_sum=sum(flow_mean * 86400, na.rm = TRUE),
    stream_SPC=mean(stream_SPC, na.rm = TRUE),
    stream_DO=mean(stream_DO, na.rm = TRUE),
    stream_temp=mean(stream_temp, na.rm = TRUE))
    
summary(week_dat)
##===========================================
## create a new df for complete GPP obs 
GPP_df <- datQ%>%
  dplyr::select( -middle_ER, -upper_ER, -lower_ER)
GPP_df<- GPP_df%>%
  drop_na(middle_GPP)
ER_df <- datQ%>%
  dplyr::select(-middle_GPP, -upper_GPP, -lower_GPP)
ER_df<- ER_df%>%
  drop_na(middle_ER)

ER_df$ER <- c(ER_df$middle_ER *-1)
ER_df$ER_low <- c(ER_df$lower_ER *-1)
ER_df$ER_up <- c(ER_df$upper_ER *-1)
#GPP_df1<- na.omit(GPP_df)
summary(GPP_df)
summary(ER_df)

##===========================================
## Begin analysis 
#============================================
## shore summaries
## BW
BW_dat <- datQ %>%
  filter(shore=="BW" & date> as.Date("2023-02-07") & date < as.Date("2023-09-06"))
round(mean(na.omit(BW_dat$middle_GPP)),2)
round(se(na.omit(BW_dat$middle_GPP)),2)
range(na.omit(BW_dat$middle_GPP))
round(mean(na.omit(BW_dat$middle_ER)),2)
round(se(na.omit(BW_dat$middle_ER)),2)
round(se(na.omit(BW_dat$middle_ER)),2)
range(na.omit(BW_dat$middle_ER))
round(mean(na.omit(BW_dat$middle_NEP)),2)
round(se(na.omit(BW_dat$middle_NEP)),2)
summary(BW_dat)
max_index <- max(BW_dat$middle_GPP)
print(BW_dat[max_index, ])
min_index <- min((BW_dat$middle_GPP))
print(BW_dat[min_index, ])

#============================================
# Figure out auto regressive structure:
#============================================
hist(GPP_df$middle_GPP) # strong skew
hist(log(GPP_df$middle_GPP+1))

gpp_fit2 <- glmmTMB((middle_GPP) ~  lag(middle_GPP, 2) +
                      (1|shore/site), data = GPP_df, family = lognormal)
summary(gpp_fit2)
residuals2 <- residuals(gpp_fit2)
hist(residuals2)
# Plot ACF and PACF of residuals
acf(residuals2, main="ACF of Residuals")
pacf(residuals2, main="PACF of Residuals")


gpp_fit3 <- glmmTMB((middle_GPP) ~  lag(middle_GPP, 1) +
                      (1|shore/site), data = GPP_df, family = lognormal)
summary(gpp_fit3)
residuals2 <- residuals(gpp_fit3)
hist(residuals2)
# Plot ACF and PACF of residuals
acf(residuals2, main="ACF of Residuals")
pacf(residuals2, main="PACF of Residuals")

#
AIC(gpp_fit2, gpp_fit3)
# fit three looks best and is most reasonable 
# yet there is still a high degree of temporal auto correlation not well captured. 

##===========================================
## GPP glms for stream influence
#============================================
# subset for shores with streams BW and GB
GPP_df_stream <- GPP_df %>%
  filter(shore=="BW" |shore=="GB")
# create 1-day lag
GPP_df_stream$middle_GPP_lag <- lag(GPP_df_stream$middle_GPP)

ER_df_stream <- ER_df %>%
  filter(shore=="BW" |shore=="GB")
# create 1-day lag
ER_df_stream$middle_GPP_lag <- lag(ER_df_stream$middle_ER)


##===========================================
## check co-variance 
names(GPP_df_stream)

# scale variables
# GPP_df_stream[, c(7, 19:26)] <- scale(GPP_df_stream[, c(7, 19:26)])
Streamdf_cor <- scale(GPP_df_stream[, c(9, 15, 17, 20, 26:31)])
chart.Correlation(Streamdf_cor, histogram=TRUE, pch=19)

## ER ##
Streamdf_cor <- scale(ER_df_stream[, c(34, 15:17, 19, 20, 22, 23:30 )])
chart.Correlation(Streamdf_cor, histogram=TRUE, pch=19)

##==========================
### MODEL ###
## NULL:
gpp_stream_null <- lmer(log(middle_GPP+1) ~ log(middle_GPP_lag+1) + (1|shore / site) +(1|WaterYear), data = GPP_df_stream)
summary(gpp_stream_null)
residuals3 <- residuals(gpp_stream_null)
hist(residuals3)
r.squaredGLMM(gpp_stream_null)
vif(gpp_stream_null)

##===========================================
## Predictors for stream 
gpp_stream_mod <- lmer(log(middle_GPP+1) ~  log(middle_GPP_lag+1) +
                            scale(stream_temp)  + 
                            scale(flow_mean)+
                            scale(precip_bi)+
                            (1|shore / site) +(1|WaterYear), data = GPP_df_stream)
summary(gpp_stream_mod)
residuals3 <- residuals(gpp_stream_mod)
hist(residuals)
r.squaredGLMM(gpp_stream_mod)
vif(gpp_stream_mod)

gpp_stream_mod1 <- lmer(log(middle_GPP+1) ~  log(middle_GPP_lag+1) +
                         scale(flow_mean)+
                         scale(precip_bi)+
                         (1|shore / site) +(1|WaterYear), data = GPP_df_stream)
summary(gpp_stream_mod1)
residuals3 <- residuals(gpp_stream_mod1)
hist(residuals)
r.squaredGLMM(gpp_stream_mod1)
vif(gpp_stream_mod1)
AIC(gpp_stream_null, gpp_stream_mod, gpp_stream_mod1)

# small improvement with gpp_stream_mod
##===========================================

## Predictors for lake
lakedf_cor <- scale(GPP_df_stream[, c(15, 20, 26, 30, 31)])
chart.Correlation(lakedf_cor, histogram=TRUE, pch=19)
## ALL light correlations are negative -- maybe photo inhibition 
##===========================================
## Predictors for stream + lake
gpp_stream_modl <- lmer(log(middle_GPP+1) ~ log(middle_GPP_lag+1) +
                             scale(lake_tempC)  + 
                             scale(par_int_3m) +
                             scale(flow_mean)+
                             scale(precip_bi)+
                             (1|shore / site) +(1|WaterYear), data = GPP_df_stream)
summary(gpp_stream_modl)
residuals3 <- residuals(gpp_stream_modl)
hist(residuals3)
r.squaredGLMM(gpp_stream_modl)
vif(gpp_stream_modl)

gpp_stream_modl2 <- lmer(log(middle_GPP+1) ~ log(middle_GPP_lag+1) +
                          scale(lake_tempC)  + 
                          scale(par_int_3m) +
                          (1|shore / site) +(1|WaterYear), data = GPP_df_stream)
summary(gpp_stream_modl2)
residuals3 <- residuals(gpp_stream_modl2)
hist(residuals3)
r.squaredGLMM(gpp_stream_modl2)
vif(gpp_stream_modl2)

AIC(gpp_stream_null, gpp_stream_modl, gpp_stream_modl2)
# lake describes better than streams gpp_stream_modl2

##===========================================
## Characterizing biogeochemical stream transport
##===========================================
## Create chem summary table 1
chem_dat_stream_tab <- chem_dat %>%
  group_by(shore, location, substrate) %>%
  summarise(
    NO3_mgL_dl_m = mean(NO3_mgL_dl, na.rm=T),
    NO3_mgL_dl_min = min(NO3_mgL_dl, na.rm=T),
    NO3_mgL_dl_max = max(NO3_mgL_dl, na.rm=T),
    NH4_mgL_dl_m= mean(NH4_mgL_dl, na.rm=T),
    NH4_mgL_dl_min= min(NH4_mgL_dl, na.rm=T),
    NH4_mgL_dl_max= max(NH4_mgL_dl, na.rm=T),
    PO4_ugL_dl_m = mean(PO4_ugL_dl, na.rm=T),
    PO4_ugL_dl_min = min(PO4_ugL_dl, na.rm=T),
    PO4_ugL_dl_max = max(PO4_ugL_dl, na.rm=T),
    DOC_mgL_dl_m = mean(DOC_mgL_dl, na.rm=T),
    DOC_mgL_dl_min = min(DOC_mgL_dl, na.rm=T),
    DOC_mgL_dl_max = max(DOC_mgL_dl, na.rm=T),
    pH_infill_m= mean(pH_infill, na.rm=T),
    pH_infill_min= min(pH_infill, na.rm=T),
    pH_infill_max= max(pH_infill, na.rm=T))

# write.csv(chem_dat_stream_tab, file = "./Chem_table1.csv", row.names = TRUE)

## Aggregate to the shore areas:
chem_dat_shore <- chem_dat %>%
  group_by(shore, location, date, substrate, WaterYear) %>%
  summarise(
    NO3_mgL_dl = mean(NO3_mgL_dl, na.rm=T),
    NH3_mgL_dl= mean(NH3_mgL_dl, na.rm=T),
    NH4_mgL_dl= mean(NH4_mgL_dl, na.rm=T),
    PO4_ugL_dl = mean(PO4_ugL_dl, na.rm=T),
    DOC_mgL_dl = mean(DOC_mgL_dl, na.rm=T))

# Wide to long: 
chem_dat_wide <- chem_dat_shore %>%
  pivot_wider(
    id_cols = c("shore", "date", "substrate", "WaterYear"),
    names_from = c("location"),
    values_from = c("NO3_mgL_dl", "NH4_mgL_dl", "PO4_ugL_dl", "DOC_mgL_dl")
  )

## Three week average grouping: 
chem_dat_filtered <- chem_dat_wide %>%
  mutate(week_period = ceiling(week(date) / 3))%>%
  mutate(grouping_var = paste0(year(date), "-", week_period))

chem_dat_wide_summary <- chem_dat_filtered %>%
  group_by(shore, substrate, grouping_var) %>%
  summarise(
    NO3_mgL_dl_stream = mean(NO3_mgL_dl_stream, na.rm = TRUE),
    NO3_mgL_dl_lake = mean(NO3_mgL_dl_lake, na.rm = TRUE),
    NO3_mgL_dl_interface = mean(NO3_mgL_dl_interface, na.rm = TRUE),
    NH4_mgL_dl_stream = mean(NH4_mgL_dl_stream, na.rm = TRUE),
    NH4_mgL_dl_lake = mean(NH4_mgL_dl_lake, na.rm = TRUE),
    NH4_mgL_dl_interface = mean(NH4_mgL_dl_interface, na.rm = TRUE),
    PO4_ugL_dl_stream = mean(PO4_ugL_dl_stream, na.rm = TRUE),
    PO4_ugL_dl_lake = mean(PO4_ugL_dl_lake, na.rm = TRUE),
    PO4_ugL_dl_interface = mean(PO4_ugL_dl_interface, na.rm = TRUE),
    DOC_mgL_dl_stream = mean(DOC_mgL_dl_stream, na.rm = TRUE),
    DOC_mgL_dl_lake = mean(DOC_mgL_dl_lake, na.rm = TRUE),
    DOC_mgL_dl_interface = mean(DOC_mgL_dl_interface, na.rm = TRUE),
    # Add other variables and summary statistics as needed
    .groups = "drop"
  )

##===========================================
## NO3 ##
##===========================================
# Fit the linear mixed-effects model
stream_NO3 <- lmer(NO3_mgL_dl_lake ~ NO3_mgL_dl_stream + (1 | shore), data = chem_dat_wide_summary%>%filter(shore=="BW" | shore=="GB"))
summary(stream_NO3)
r.squaredGLMM(stream_NO3)
hist(residuals(stream_NO3))
# Extract the fixed effects coefficients
coefficients <- fixef(stream_NO3)
# Extract the slope coefficient
slope <- coefficients["NO3_mgL_dl_stream"]
# Calculate the reciprocal of the slope
inverse_slope <- 1 / slope
# Create a data frame for plotting the trend line
trend_data <- data.frame(NO3_mgL_dl_stream = seq(min(na.omit(chem_dat_wide_summary$NO3_mgL_dl_stream)),
                                                 max(na.omit(chem_dat_wide_summary$NO3_mgL_dl_stream)), length.out = 100),
                         NO3_mgL_dl_lake = inverse_slope * seq(min(na.omit(chem_dat_wide_summary$NO3_mgL_dl_stream)),
                                                               max(na.omit(chem_dat_wide_summary$NO3_mgL_dl_stream)), length.out = 100))
# Plot
NO3_plot <- ggplot(data = chem_dat_wide_summary%>%filter(shore=="BW" | shore=="GB"), 
                   aes(y = NO3_mgL_dl_lake, x = NO3_mgL_dl_stream)) +
  geom_point(size = 2, alpha = 0.75, aes(shape = substrate, color=shore)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  
  xlim(0,0.2) +ylim(0,0.2)+
  geom_line(data = trend_data, aes(y = NO3_mgL_dl_lake, x = NO3_mgL_dl_stream), 
            color = "grey25", linetype = 1) +  # Add trend line
  geom_hline(yintercept = 0.003, color = "grey50") + 
  geom_vline(xintercept = 0.003, color = "grey50") +
  labs(y=expression(Lake~NO[3]^-1~mg~L^-1), x=expression(Stream~NO[3]^-1~mg~L^-1)) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()

##===========================================
## NH4 ##
##===========================================
# Fit the linear mixed-effects model
stream_NH3 <- lmer(NH4_mgL_dl_lake ~ NH4_mgL_dl_stream + (1 | shore), data = chem_dat_wide_summary%>%filter(shore=="BW" | shore=="GB"))
summary(stream_NH3)
r.squaredGLMM(stream_NH3)
hist(residuals(stream_NH3))
# Extract the fixed effects coefficients
coefficients <- fixef(stream_NH3)
# Extract the slope coefficient
slope <- coefficients["NH4_mgL_dl_stream"]
# Calculate the reciprocal of the slope
inverse_slope <- 1 / slope
# Create a data frame for plotting the trend line
trend_data <- data.frame(NH4_mgL_dl_stream = seq(min(na.omit(chem_dat_wide_summary$NH4_mgL_dl_stream)),
                                                 max(na.omit(chem_dat_wide_summary$NH4_mgL_dl_stream)), length.out = 100),
                         NH4_mgL_dl_lake = inverse_slope * seq(min(na.omit(chem_dat_wide_summary$NH4_mgL_dl_stream)),
                                                               max(na.omit(chem_dat_wide_summary$NH4_mgL_dl_stream)), length.out = 100))


# Plot
NH3_plot <- ggplot(data = chem_dat_wide_summary%>%filter(shore=="BW" | shore=="GB") , 
                   aes(y = NH4_mgL_dl_lake, x = NH4_mgL_dl_stream, color = shore)) +
  geom_point(size = 2, alpha = 0.75, aes(shape = substrate)) + 
  geom_line(data = trend_data, aes(y = NH4_mgL_dl_lake, x = NH4_mgL_dl_stream), 
            color = "grey25", linetype = 1) +  # Add trend line
  geom_hline(yintercept = 0.002, color = "grey50") + 
  geom_vline(xintercept = 0.002, color = "grey50") +
  labs(y = expression(Lake~NH[4]^+1~mg~L^-1), x = expression(Stream~NH[4]^+1~mg~L^-1)) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  
  xlim(0,0.1) + ylim(0,0.1) +
  theme_bw()

NH3_plot

##===========================================
## ANOVA NEP + chem 
##===========================================

SFS_c_WYvar22 <- chem_dat%>%
  filter(date>as.Date("2021-10-01") & date<as.Date("2022-09-15"))

SFS_c_WYvar23 <- chem_dat%>%
  filter(date>as.Date("2022-10-01") & date<as.Date("2023-09-15"))

SFS_chem_WYvar<- rbind(SFS_c_WYvar22, SFS_c_WYvar23)
### 

SFS_chem_temp_SH <- SFS_chem_WYvar%>%
  filter(shore!="SH" & site!="BWU" & site!="BWL" &  site!="GBL" & site!="GBU")

str(SFS_chem_temp_SH)
SFS_chem_temp_SH<- water_year(SFS_chem_temp_SH)
SFS_chem_temp_SH$WY <- as.factor(SFS_chem_temp_SH$WaterYear)

SFS_chem_temp_BW <- SFS_chem_temp_SH%>%
  filter(shore=="BW")

SFS_chem_temp_GB <- SFS_chem_temp_SH%>%
  filter(shore=="GB")

SFS_chem_temp_SS <- SFS_chem_temp_SH%>%
  filter(shore=="SS")

anova_NO3 <- aov(NO3_mgL_dl ~ WY, data = SFS_chem_temp_SS%>%filter(NO3_mgL_dl<=0.003))
summary(anova_NO3)

anova_Nh4 <- aov(NH4_mgL_dl ~ WY, data = SFS_chem_temp_SS%>%filter(NH4_mgL_dl<=0.002))
summary(anova_Nh4)

anova_PO4 <- aov(PO4_ugL_dl ~ WY, data = SFS_chem_temp_SS)
summary(anova_PO4)

anova_DOC <- aov(DOC_mgL_dl ~ WY, data = SFS_chem_temp_SS)
summary(anova_DOC)




plot_NO3 <- ggplot(SFS_chem_temp_SH, aes(x = interaction(WY, shore), y = NO3_mgL_dl, fill = shore)) +
  geom_jitter(aes(color = shore, shape=substrate), width = 0.1, height = 0.01, alpha = 0.5) +  # Adding jittered points with 50% transparency
  geom_boxplot(alpha = 0.9, outlier.alpha = 0.5) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_discrete(labels = (c("2022", "2023", "2022", "2023", "2022", "2023"))) +
  theme_bw() +
  geom_hline(yintercept = 0.003, color = "grey50") + 
  labs(y = expression(NO[3]^-1~N~(mg~L^-1)), x = NULL) 

plot_Kd <- plot_Kd +
  geom_text(data = data.frame(
    shore = c("BW", "GB", "SS"),
    WY = factor(rep("2022", 3)), 
    y = c(0.32) 
  ),
  aes(x = interaction(WY, shore), y = y, label = "*"),
  color = "black", size = 10,
  )



plot_NH4 <- ggplot(SFS_chem_temp_SH, aes(x = interaction(WY, shore), y = NH4_mgL_dl, fill = shore)) +
  geom_jitter(aes(color = shore, shape=substrate), width = 0.1, height = 0.01, alpha = 0.5) +  # Adding jittered points with 50% transparency
  geom_boxplot(alpha = 0.9, outlier.alpha = 0.5) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_discrete(labels = (c("2022", "2023", "2022", "2023", "2022", "2023"))) +
  theme_bw() +
  geom_hline(yintercept = 0.002, color = "grey50") + 
  labs(y = expression(NH[4]^+1~-N~(mg~L^-1)), x = NULL) 

plot_Kd <- plot_Kd +
  geom_text(data = data.frame(
    shore = c("BW", "GB", "SS"),
    WY = factor(rep("2022", 3)), # Adjust as needed for your specific x-axis values
    y = c(0.32) # Adjust y-values to be above the box plots
  ),
  aes(x = interaction(WY, shore), y = y, label = "*"),
  color = "black", size = 10,
  )


# Fit the linear mixed-effects model
stream_PO4 <- lmer(PO4_ugL_dl_lake ~ PO4_ugL_dl_stream + (1 | shore), data = chem_dat_wide_summary)
summary(stream_PO4)
# Create a data frame for plotting the trend line
# Plot
PO4_plot <- ggplot(data = chem_dat_wide_summary,
                   aes(x = PO4_ugL_dl_stream, y = PO4_ugL_dl_lake, color = shore)) +
  geom_point(size = 2, alpha = 0.75, aes(shape = substrate)) +
  labs(y=expression(Lake~ophos~ug~L^-1), x=expression(Stream~ophos~ug~L^-1)) +
  geom_hline(yintercept = 0.402, color = "grey50") + 
  geom_vline(xintercept = 0.402, color = "grey50") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 60, by = 12),
                     limits = c(0, 60)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 60, by = 12),
                     limits = c(0, 60)) +  
  theme_bw()
  



library(tidyverse)

# Remove rows with missing values
chem_dat_wide_summary_clean <- chem_dat_wide_summary %>%
  drop_na(DOC_mgL_dl_lake, DOC_mgL_dl_stream)

# Plot
DOC_plot <- ggplot(data = chem_dat_wide_summary_clean,
                   aes(y = DOC_mgL_dl_lake, x = DOC_mgL_dl_stream, color = shore)) +
  geom_point(size = 2, alpha = 0.75, aes(shape = substrate)) +
  geom_hline(yintercept = 0.25, color = "grey50") + 
  geom_vline(xintercept = 0.25, color = "grey50") +
  labs(y = expression(Lake~DOC~mg~L^-1), x = expression(Stream~DOC~mg~L^-1)) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 4.5, by = 0.9),
                     limits = c(0, 4.5)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 4.5, by = 0.9),
                     limits = c(0, 4.5)) +  
  theme_bw()

DOC_plot

# Fit the linear mixed-effects model
stream_DOC <- lmer(DOC_mgL_dl_lake ~ DOC_mgL_dl_stream + (1 | shore), data = chem_dat_wide_summary_clean)
summary(stream_DOC)
r.squaredGLMM(stream_DOC)
hist(residuals(stream_DOC))

chem_grid <- ggarrange(NO3_plot,
                     NH3_plot,
                     DOC_plot,
                     PO4_plot,
                     ncol = 2, nrow = 2,
                     common.legend = TRUE, 
                     labels = c("A", "B", "C", "D"),
                     legend = "bottom")

# ggsave(plot = chem_grid, filename = paste("./Lake_streamChem_plot.png",sep=""),width=5,height=5.05,dpi=300)

# end script.