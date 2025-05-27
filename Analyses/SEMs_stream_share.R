
##==============================================================================
## SEM of streamflow and precipitation dynamics on nearshore metabolism
## by Loria et al. 2024
## 05/27/2025

# inspired by # https://rpubs.com/jebyrnes/brms_bayes_sem
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
library(piecewiseSEM) 
library(DiagrammeR)
library(brms)
library(gridExtra)

se <- function(dat){
  se <- sd(dat)/sqrt(length(dat))
  return(se)}

set.seed(2025)
##===========================================
## read data aggregated data for the project:
#============================================
dat <-readRDS("./NS_analysis_dat.rds")
str(dat)


##===========================================
## metabolism dataframes: 
GPP_df <- dat%>%
  dplyr::select(-middle_ER)

ER_df <- dat%>%
  dplyr::select(-middle_GPP)%>%
  mutate(ER = middle_ER *-1)

NEP_df <- dat%>%
  dplyr::select(-middle_GPP, -middle_ER)

##===========================================
## Shore location dataframes  

#### BW: 

GPP_df_BW <- GPP_df%>%
  filter(shore == "BW")

plot_data_DO <- GPP_df_BW %>%
  filter(!is.na(lake_DO)) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_GPP = !is.na(middle_GPP),
         plot_GPP = ifelse(is.na(middle_GPP), 0, middle_GPP),
         yday=yday(date),
         year=year(date))

# ER
ER_df_BW <- ER_df%>%
  filter(shore == "BW")
plot_data_DO_BW_ER <- ER_df_BW %>%
  filter(!is.na(lake_DO)) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_ER = !is.na(ER),
         plot_ER = ifelse(is.na(ER), 0, ER),
         yday=yday(date),
         year=year(date))

#### GB: 

GPP_df_GB <- GPP_df%>%
  filter(shore == "GB")
plot_data_DO_gb <- GPP_df_GB %>%
  filter(!is.na(lake_DO)) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_GPP = !is.na(middle_GPP),
         plot_GPP = ifelse(is.na(middle_GPP), 0, middle_GPP),
         yday=yday(date),
         year=year(date))

# ER
ER_df_GB <- ER_df%>%
  filter(shore == "GB")
plot_data_DO_GB_ER <- ER_df_GB %>%
  filter(!is.na(lake_DO)) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_ER = !is.na(ER),
         plot_ER = ifelse(is.na(ER), 0, ER),
         yday=yday(date),
         year=year(date))

#### SS: 
GPP_df_SS <- GPP_df%>%
  filter(shore == "SS")
plot_data_DO_SS <- GPP_df_SS %>%
  filter(!is.na(lake_DO)) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_GPP = !is.na(middle_GPP),
         plot_GPP = ifelse(is.na(middle_GPP), 0, middle_GPP),
         yday=yday(date),
         year=year(date))


# ER
ER_df_SS <- ER_df%>%
  filter(shore == "SS")
plot_data_DO_SS_ER <- ER_df_SS %>%
  filter(!is.na(lake_DO)) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_ER = !is.na(ER),
         plot_ER = ifelse(is.na(ER), 0, ER),
         yday=yday(date),
         year=year(date))

##============================================
## Stream SEM dataframe organization:
##============================================
## BW POE 20 % = 0.04334224 - log = 0.04242925
## GB POE 20 % = 0.005187275,  - log = 0.005173867

columns_to_scale <- c("lake_tempC", "light_mean", "windsp_mean","flow_mean", "ppt_mm")

plot_data_DO1<-plot_data_DO %>%
  mutate(
    yday=lubridate::yday(date))

## quick check of flow with POE threshold:
BW_flow_threshold <- ggplot(plot_data_DO1,
                       aes(x = yday, y = flow_mean, color = as.factor(year))) +
  geom_point() +  
  scale_color_viridis_d(name = "Year") +  
  theme_classic() + 
  ylab(expression(Streamflow~(m^3~s^-1))) +
  xlab("Day of year") +
  geom_hline(aes(yintercept =0.42),linetype = "dashed", size = 0.5, alpha = 0.9)


## =====
## BW:
#
GPP_df_BW <- plot_data_DO%>%
  tidyr::drop_na(plot_GPP, lake_tempC, light_mean, flow_mean, ppt_mm) #%>% filter(flow_mean<0.42) 

hist(GPP_df_BW$plot_GPP)
hist(log(GPP_df_BW$plot_GPP +1))
hist(GPP_df_BW$lake_tempC)
hist(log(GPP_df_BW$lake_tempC))
hist(GPP_df_BW$light_mean)
hist(log(GPP_df_BW$lake_tempC))
hist(GPP_df_BW$flow_mean)
hist(log(GPP_df_BW$flow_mean))
hist(scale(GPP_df_BW$flow_mean))
hist(scale(log(GPP_df_BW$flow_mean)))
hist(GPP_df_BW$ppt_mm)
hist(log(GPP_df_BW$ppt_mm+1))


##
##Look at response data structure-- 
### Log-transform the data, because the Shapiro-Wilk test only checks for normality, not log-normality directly
summary(GPP_df_BW)
## Subset only the positive values
# GPP_positive <- GPP_df_BW$plot_GPP[GPP_df_BW$plot_GPP > 0]
## Log-transform
# log_GPP_positive <- log(GPP_positive)
## Shapiro-Wilk test
# shapiro.test(log_GPP_positive)

GPP_df_BW<- as.data.frame(GPP_df_BW)
str(GPP_df_BW)

GPP_df_BW1 <-GPP_df_BW %>%
  mutate(across(all_of(columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_GPP_lag = lag(plot_GPP),
         flow_mean = log(flow_mean)) %>%
  dplyr::select(site, shore, date, plot_GPP, plot_GPP_lag, flow_mean, lake_tempC, light_mean, ppt_mm, ppt_mm) %>%
  tidyr::drop_na(plot_GPP, plot_GPP_lag, lake_tempC, light_mean, flow_mean, ppt_mm, ppt_mm) 
  
# double check flow shapes post log then post scale
summary(GPP_df_BW1)
hist(GPP_df_BW1$flow_mean)

ER_df_BW <- plot_data_DO_BW_ER%>%
  mutate(across(all_of(columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_ER_lag = lag(plot_ER),
         flow_mean = log(flow_mean)) %>%
  dplyr::select(site, shore, date, plot_ER, plot_ER_lag, flow_mean, lake_tempC,light_mean, ppt_mm, ppt_mm) %>%
  drop_na(plot_ER, plot_ER_lag, lake_tempC, light_mean, flow_mean, ppt_mm, ppt_mm) 

hist(ER_df_BW$flow_mean)

## =====
## GB:
# 

## quick check of flow with POE threshold:
GB_flow_threshold <- ggplot(plot_data_DO_gb,
                            aes(x = yday, y = flow_mean, color = as.factor(year))) +
  geom_point() +  
  scale_color_viridis_d(name = "Year") +  
  theme_classic() + 
  ylab(expression(Streamflow~(m^3~s^-1))) +
  xlab("Day of year") +
  geom_hline(aes(yintercept =0.005),linetype = "dashed", size = 0.5, alpha = 0.9)


plot_data_DO_gb <- plot_data_DO_gb %>%
  drop_na(plot_GPP, lake_tempC, light_mean, flow_mean, ppt_mm) # %>% filter(flow_mean<0.005)

GPP_df_GB <- plot_data_DO_gb%>%
  mutate(flow_mean = log(flow_mean)) %>%
  mutate(across(all_of(columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_GPP_lag = lag(plot_GPP)) %>%
  dplyr::select(site, shore, date, plot_GPP, plot_GPP_lag, flow_mean, lake_tempC, light_mean,  ppt_mm, ppt_mm) %>%
  drop_na(plot_GPP, plot_GPP_lag, lake_tempC, light_mean, flow_mean,  ppt_mm, ppt_mm) 

hist(GPP_df_GB$plot_GPP)
hist(GPP_df_GB$lake_tempC)
hist(GPP_df_GB$light_mean)
hist(GPP_df_GB$flow_mean)
hist(GPP_df_GB$ppt_mm)

ER_df_GB <- plot_data_DO_GB_ER%>%
  drop_na(plot_ER, lake_tempC, light_mean, flow_mean, ppt_mm) %>% 
 # filter(flow_mean<0.005)%>%
  mutate(flow_mean = log(flow_mean)) %>%
  mutate(across(all_of(columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_ER_lag = lag(plot_ER)) %>%
  dplyr::select(site, shore, date, plot_ER, plot_ER_lag, flow_mean, lake_tempC, light_mean, ppt_mm, ppt_mm) %>%
  drop_na(plot_ER, plot_ER_lag, lake_tempC, light_mean, flow_mean,  ppt_mm, ppt_mm) 

hist(ER_df_GB$flow_mean)

## =====
## SS:
#
GPP_df_SS <- plot_data_DO_SS%>%
  drop_na(plot_GPP, lake_tempC, light_mean, ppt_mm) %>% 
  mutate(across(all_of(columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_GPP_lag = lag(plot_GPP)) %>%
  dplyr::select(site, shore, date, plot_GPP, plot_GPP_lag, lake_tempC, light_mean,  ppt_mm, ppt_mm) %>%
  drop_na(plot_GPP, plot_GPP_lag, lake_tempC, light_mean,  ppt_mm, ppt_mm) 

hist(GPP_df_SS$plot_GPP)
hist(GPP_df_SS$lake_tempC)
hist(GPP_df_SS$light_mean)

ER_df_SS <- plot_data_DO_SS_ER %>%
  drop_na(ER, lake_tempC, light_mean, ppt_mm) %>% 
  mutate(across(all_of(columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_ER_lag = lag(ER)) %>%
  dplyr::select(site, shore, date, plot_ER, plot_ER_lag, lake_tempC, light_mean, ppt_mm, ppt_mm) %>%
  drop_na(plot_ER, plot_ER_lag, lake_tempC, light_mean,  ppt_mm, ppt_mm) 


##============================================
## Start SEM 
##============================================

###########################
#### revised GPP model ###
##########################

## Was:
# Log- transformed GPP input:
#GPP_mod <- bf(plot_GPP ~ plot_GPP_lag + light_mean + lake_tempC + flow_mean + (1|site)) 

# hurdle transformed :
GPP_mod <- bf(
  plot_GPP ~ plot_GPP_lag + light_mean + lake_tempC + flow_mean + (1|site),
  family = hurdle_lognormal()
)
temp_mod <- bf(lake_tempC ~ flow_mean + ppt_mm + (1|site)) 
light_mod <- bf(light_mean ~ flow_mean + ppt_mm + (1|site)) 
Q_mod <- bf(flow_mean ~ ppt_mm + (1|site)) 

### Select WY 22:
BW_GPP_df22 <- GPP_df_BW1 %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

########
####
## BW

BW_GPP_stream_fit <- brm(
  GPP_mod +
    temp_mod +
    light_mod +
    Q_mod +
    set_rescor(FALSE),
  data = BW_GPP_df22,
  iter = 10000, warmup = 5000,
  control = list(adapt_delta = 0.96, max_treedepth = 16), 
  cores = 4, chains = 3
)

# check convergence: 
plot(BW_GPP_stream_fit)

# Gives a measure of explained variance in each main model path
bayes_R2(BW_GPP_stream_fit)

# Evaluate out-of-sample predictive accuracy using:
# p_loo: Effective number of parameters, analogous to the model complexity or effective degrees of freedom.
# too high = over fitting
loo(BW_GPP_stream_fit) 

## Some notes: 
#     - MCSE of elpd_loo: Monte Carlo Standard Error
#     - A very low value (0.1) means the estimate of elpd_loo is very stable.
#     - Pareto k diagnostics: k < 0 , no problematic observations, reliable across all data points.
#     - values  >0.7 or >1, it would signal potential influence points or outliers that affect prediction.


# Extract posterior samples BW:
posterior_samples_BW <- posterior_samples(BW_GPP_stream_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BW <- posterior_samples_BW %>%
 dplyr::select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm, b_flowmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW_22 <- melt(posterior_samples_BW, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

### Optional save for later plotting: 
# saveRDS(posterior_samples_melted_BW_22, "./SEM_output/posterior_samples_melted_BW_22_final_logflow.rds")



########
####
## GB

GB_GPP_df22 <- GPP_df_GB %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

GB_GPP_stream_fit_22 <- brm(GPP_mod +
                              temp_mod +
                              light_mod +
                              Q_mod +
                              set_rescor(FALSE),
                            data=GB_GPP_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.96, max_treedepth = 16), 
                            cores=4, chains = 3)

plot(GB_GPP_stream_fit_22)
summary(GB_GPP_stream_fit_22)
GB_GPP_stream_fit_sum <- summary(GB_GPP_stream_fit_22)

#
bayes_R2(GB_GPP_stream_fit_22)

# Evaluate out-of-sample predictive accuracy using:
loo(GB_GPP_stream_fit_22)

# Optional plot residuals against fitted values or predictors:
# conditional_effects(GB_GPP_stream_fit_22)

# Extract posterior samples GB:
posterior_samples_GB <- posterior_samples(GB_GPP_stream_fit_22)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))
# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  dplyr::select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
                b_laketempC_flow_mean, b_laketempC_ppt_mm,
                b_lightmean_flow_mean, b_lightmean_ppt_mm, b_flowmean_ppt_mm) %>%
  mutate(shore = shore_levels)
         
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB_22 <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_GB_22, "./SEM_output/posterior_samples_melted_GB_22_GPP_logflow.rds")


########
####
## SS

######## GPP model (GPP_mod1) without streamflow: 
# hurdle transformed :
GPP_mod1 <- bf(
  plot_GPP ~ plot_GPP_lag + light_mean + lake_tempC + (1|site),
  family = hurdle_lognormal()
)
temp_mod1 <- bf(lake_tempC ~  ppt_mm + (1|site)) 
light_mod1 <- bf(light_mean ~  ppt_mm + (1|site)) 


SS_GPP_df22 <- GPP_df_SS %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

## SS
SS_GPP_stream_fit_22 <- brm(GPP_mod1 +
                              temp_mod1 +
                              light_mod1 +
                              set_rescor(FALSE),
                            data=SS_GPP_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.96, max_treedepth = 16), 
                            cores=4, chains = 3)

plot(SS_GPP_stream_fit_22)

bayes_R2(SS_GPP_stream_fit_22)

# Evaluate out-of-sample predictive accuracy using:

loo(SS_GPP_stream_fit_22) 
waic(GB_GPP_stream_fit_22) #Is stable and well-behaved (MCSE and Pareto k look good).

summary(SS_GPP_stream_fit_22)
SS_GPP_stream_fit_sum <- summary(SS_GPP_stream_fit_22)
# Extract posterior samples SS:
posterior_samples_SS <- posterior_samples(SS_GPP_stream_fit_22)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  dplyr::select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC,
         b_laketempC_ppt_mm, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SS_22 <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_SS_22, "./SEM_output/posterior_samples_melted_SS_22_finalv2.rds")


########
####
## WY 2023
####
########


########
####
## BW 

BW_GPP_df23 <- GPP_df_BW1 %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-14"))
summary(BW_GPP_df23)

BW_GPP_stream_fit_23 <- brm(GPP_mod +
                              temp_mod +
                              light_mod +
                              Q_mod +
                              set_rescor(FALSE),
                            data = BW_GPP_df23,
                            #  family = hurdle_lognormal(),  # Add this or inside bf() as shown above
                            iter = 10000, warmup = 5000,
                            control = list(adapt_delta = 0.96, max_treedepth = 16), 
                            cores = 4, chains = 3)

plot(BW_GPP_stream_fit_23)
BW_GPP_stream_fit_sum <- summary(BW_GPP_stream_fit_23)

bayes_R2(BW_GPP_stream_fit_23)

# Evaluate out-of-sample predictive accuracy using:
loo(BW_GPP_stream_fit_23) # 
waic(BW_GPP_stream_fit_23) #

# Extract posterior samples BW:
posterior_samples_BW <- posterior_samples(BW_GPP_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BW <- posterior_samples_BW %>%
  dplyr::select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
                b_laketempC_flow_mean, b_laketempC_ppt_mm,
                b_lightmean_flow_mean, b_lightmean_ppt_mm, b_flowmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW_23 <- melt(posterior_samples_BW, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_BW_23, "./SEM_output/posterior_samples_melted_BW_23_final_logflow.rds")

########
####
## GB

########
GB_GPP_df23 <- GPP_df_GB %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-06"))
summary(GB_GPP_df23)

hist(GB_GPP_df23$flow_mean)

## GB
GB_GPP_stream_fit_23 <- brm(GPP_mod +
                              temp_mod +
                              light_mod +
                              Q_mod +
                              set_rescor(FALSE),
                            data=GB_GPP_df23,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.96, max_treedepth = 16), 
                            cores=4, chains = 3)

plot(GB_GPP_stream_fit_23)
summary(GB_GPP_stream_fit_23)
GB_GPP_stream_fit_sum <- summary(GB_GPP_stream_fit_23)

# Gives a measure of explained variance (like R² in frequentist models
bayes_R2(GB_GPP_stream_fit_23)

# Evaluate out-of-sample predictive accuracy using:
loo(GB_GPP_stream_fit_23) 

# Extract posterior samples GB:
posterior_samples_GB <- posterior_samples(GB_GPP_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))
# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  dplyr::select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
                b_laketempC_flow_mean, b_laketempC_ppt_mm,
                b_lightmean_flow_mean, b_lightmean_ppt_mm, b_flowmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB_23 <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_GB_23, "./SEM_output/posterior_samples_melted_GB_23_GPP_logflow.rds")

########
####
## SS

########
SS_GPP_df23 <- GPP_df_SS %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-14"))

## GB
SS_GPP_stream_fit_23 <- brm(GPP_mod1 +
                              temp_mod1 +
                              light_mod1 +
                              set_rescor(FALSE),
                            data=SS_GPP_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.96, max_treedepth = 16), 
                            cores=4, chains = 3)



plot(SS_GPP_stream_fit_23)
summary(SS_GPP_stream_fit_23)
SS_GPP_stream_fit_sum <- summary(SS_GPP_stream_fit_23)

bayes_R2(SS_GPP_stream_fit_23)

loo(SS_GPP_stream_fit_23) #


# Extract posterior samples SS:
posterior_samples_SS <- posterior_samples(SS_GPP_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  dplyr::select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC,
                b_laketempC_ppt_mm, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SS_23 <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_SS_23, "./SEM_output/posterior_samples_melted_SS_23_finalv2.rds")

###===============================
### ER 
###===============================

# hurdle transformed :
ER_mod <- bf(
  plot_ER ~ plot_ER_lag + light_mean + lake_tempC + flow_mean + (1|site),
  family = hurdle_lognormal()
)
temp_mod <- bf(lake_tempC ~ flow_mean + ppt_mm + (1|site)) 
light_mod <- bf(light_mean ~ flow_mean + ppt_mm + (1|site)) 
Q_mod <- bf(flow_mean ~ ppt_mm + (1|site)) 


########
####
## BW 

BW_ER_df22 <- ER_df_BW %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

BW_ER_stream_fit <- brm(
  ER_mod +
    temp_mod +
    light_mod +
    Q_mod +
    set_rescor(FALSE),
  data = BW_ER_df22,
  #  family = hurdle_lognormal(),  # Add this or inside bf() as shown above
  iter = 10000, warmup = 5000,
  control = list(adapt_delta = 0.96, max_treedepth = 16), 
  cores = 4, chains = 3
)

plot(BW_ER_stream_fit)

bayes_R2(BW_ER_stream_fit)

# Evaluate out-of-sample predictive accuracy using:
loo(BW_ER_stream_fit) #
waic(BW_ER_stream_fit) #

# Extract posterior samples BW:
posterior_samples_BW <- posterior_samples(BW_ER_stream_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BW <- posterior_samples_BW %>%
  dplyr::select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
                b_laketempC_flow_mean, b_laketempC_ppt_mm,
                b_lightmean_flow_mean, b_lightmean_ppt_mm, b_flowmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW_22 <- melt(posterior_samples_BW, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_BW_22, "./SEM_output/posterior_samples_melted_BW_22_ER_final_logflow.rds")


########
####
## GB

########
GB_ER_df22 <- ER_df_GB %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

GB_ER_stream_fit_22 <- brm(ER_mod +
                              temp_mod +
                              light_mod +
                              Q_mod +
                              set_rescor(FALSE),
                            data=GB_ER_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.96, max_treedepth = 16), 
                            cores=4, chains = 3)

plot(GB_ER_stream_fit_22)
summary(GB_ER_stream_fit_22)
GB_ER_stream_fit_sum <- summary(GB_ER_stream_fit_22)

bayes_R2(GB_ER_stream_fit_22)

# Evaluate out-of-sample predictive accuracy using:
loo(GB_ER_stream_fit_22)

# Extract posterior samples GB:
posterior_samples_GB <- posterior_samples(GB_ER_stream_fit_22)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))
# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  dplyr::select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
                b_laketempC_flow_mean, b_laketempC_ppt_mm,
                b_lightmean_flow_mean, b_lightmean_ppt_mm, b_flowmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB_22 <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_GB_22, "./SEM_output/posterior_samples_melted_GB_22_ER_final_logflow.rds")


########
####
## SS

########
# hurdle transformed :
ER_mod1 <- bf(
  plot_ER ~ plot_ER_lag + light_mean + lake_tempC + (1|site),
  family = hurdle_lognormal()
)
temp_mod1 <- bf(lake_tempC ~  ppt_mm + (1|site)) 
light_mod1 <- bf(light_mean ~  ppt_mm + (1|site)) 


SS_ER_df22 <- ER_df_SS %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

## SS
SS_ER_stream_fit_22 <- brm(ER_mod1 +
                              temp_mod1 +
                              light_mod1 +
                              set_rescor(FALSE),
                            data=SS_ER_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.96, max_treedepth = 16), 
                            cores=4, chains = 3)

plot(SS_ER_stream_fit_22)

bayes_R2(SS_ER_stream_fit_22)

# Evaluate out-of-sample predictive accuracy using:
loo(SS_ER_stream_fit_22)
waic(GB_ER_stream_fit_22) 

summary(SS_ER_stream_fit_22)
SS_ER_stream_fit_sum <- summary(SS_ER_stream_fit_22)

# Extract posterior samples SS:
posterior_samples_SS <- posterior_samples(SS_ER_stream_fit_22)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  dplyr::select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC,
                b_laketempC_ppt_mm, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SS_22 <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_SS_22, "./SEM_output/posterior_samples_melted_SS_ER_22_final_v2.rds")


########
####
## WY 2023
####
########

########
####
## BW

BW_ER_df23 <- ER_df_BW %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-14"))
summary(BW_ER_df23)

BW_ER_stream_fit_23 <- brm(ER_mod +
                              temp_mod +
                              light_mod +
                              Q_mod +
                              set_rescor(FALSE),
                            data = BW_ER_df23,
                            #  family = hurdle_lognormal(),  # Add this or inside bf() as shown above
                            iter = 10000, warmup = 5000,
                            control = list(adapt_delta = 0.96, max_treedepth = 16), 
                            cores = 4, chains = 3)

plot(BW_ER_stream_fit_23)
BW_ER_stream_fit_sum <- summary(BW_ER_stream_fit_23)

bayes_R2(BW_ER_stream_fit_23)

# Evaluate out-of-sample predictive accuracy using:
loo(BW_ER_stream_fit_23) 
waic(BW_ER_stream_fit_23) 


# Extract posterior samples BW:
posterior_samples_BW <- posterior_samples(BW_ER_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BW <- posterior_samples_BW %>%
  dplyr::select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
                b_laketempC_flow_mean, b_laketempC_ppt_mm,
                b_lightmean_flow_mean, b_lightmean_ppt_mm, b_flowmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW_23 <- melt(posterior_samples_BW, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


# saveRDS(posterior_samples_melted_BW_23, "./SEM_output/posterior_samples_melted_BW_23_ER_final_logflow.rds")

########
####
## GB

########
GB_ER_df23 <- ER_df_GB %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-06"))
summary(GB_ER_df23)

## GB
GB_ER_stream_fit_23 <- brm(ER_mod +
                              temp_mod +
                              light_mod +
                              Q_mod +
                              set_rescor(FALSE),
                            data=GB_ER_df23,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.96, max_treedepth = 16), 
                            cores=4, chains = 3)

plot(GB_ER_stream_fit_23)
summary(GB_ER_stream_fit_23)
GB_ER_stream_fit_sum <- summary(GB_ER_stream_fit_23)

bayes_R2(GB_ER_stream_fit_23)

# Evaluate out-of-sample predictive accuracy using:
loo(GB_ER_stream_fit_23)
waic(GB_ER_stream_fit_23) 

# Extract posterior samples GB:
posterior_samples_GB <- posterior_samples(GB_ER_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))
# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  dplyr::select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
                b_laketempC_flow_mean, b_laketempC_ppt_mm,
                b_lightmean_flow_mean, b_lightmean_ppt_mm, b_flowmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB_23 <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_GB_23, "./SEM_output/posterior_samples_melted_GB_23_ER_final_logflow.rds")


########
####
## SS

########
SS_ER_df23 <- ER_df_SS %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-14"))

## GB
SS_ER_stream_fit_23 <- brm(ER_mod1 +
                              temp_mod1 +
                              light_mod1 +
                              set_rescor(FALSE),
                            data=SS_ER_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.96, max_treedepth = 16), 
                            cores=4, chains = 3)

plot(SS_ER_stream_fit_23)
summary(SS_ER_stream_fit_23)
SS_ER_stream_fit_sum <- summary(SS_ER_stream_fit_23)

# Gives a measure of explained variance (like R² in frequentist models
bayes_R2(SS_ER_stream_fit_23)

loo(GB_ER_stream_fit_23)


# Extract posterior samples SS:
posterior_samples_SS <- posterior_samples(SS_ER_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  dplyr::select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC,
                b_laketempC_ppt_mm, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)


# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SS_23 <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_SS_23, "./SEM_output/posterior_samples_melted_SS_23_ER_final_v2.rds")



##===========================================
## Optional plot code:
## Just GPP
##===========================================

SEM_df_GPP_comb <- rbind(posterior_samples_melted_BW_22, posterior_samples_melted_GB_22, posterior_samples_melted_SS_22, 
                         posterior_samples_melted_BW_23, posterior_samples_melted_GB_23, posterior_samples_melted_SS_23)

ci_summary_overall <- SEM_df_GPP_comb %>%
  dplyr::group_by(Coefficient, shore, water_year) %>%
  dplyr:: summarize(
    median = median(Estimate),
    l_95 = quantile(Estimate, 0.025),
    u_95 = quantile(Estimate, 0.975),
    .groups = "drop"
  )


coef_order <- c(
  "b_plotGPP_plot_GPP_lag",
  "b_plotGPP_lake_tempC",
  "b_plotGPP_light_mean",
  "b_plotGPP_flow_mean",
  "b_laketempC_flow_mean",
  "b_laketempC_ppt_mm",
  "b_lightmean_flow_mean",
  "b_lightmean_ppt_mm",
  "b_flowmean_ppt_mm"
)



SEM_df_GPP_comb1 <- SEM_df_GPP_comb %>%
  mutate(point_shape = case_when(
    shore == "BW" & water_year == 2022 ~ "Large inflow dry (WY 2022)",
    shore == "BW" & water_year == 2023 ~ "Large inflow wet (WY 2023)",
    shore == "GB" & water_year == 2022 ~ "Small inflow wet (WY 2023)",
    shore == "GB" & water_year == 2023 ~ "Small inflow dry (WY 2022)",
    shore == "SS" & water_year == 2022 ~ "No inflow dry (WY 2022)",
    shore == "SS" & water_year == 2023 ~ "No inflow wet (WY 2023)",
    TRUE ~ "Other"
  ))


SEM_df_GPP_comb1 <- SEM_df_GPP_comb1 %>%
  mutate(Coefficient = fct_relevel(Coefficient, rev(coef_order)))

levels(SEM_df_GPP_comb1$Coefficient)


# Violin plot of Estimate distributions per Coefficient
coef_violin_plot  <- ggplot(SEM_df_GPP_comb1, aes(
  y = Coefficient,
  x = Estimate,
  fill = point_shape,
  color = point_shape
)) +
  geom_violin(
    width = 2,
    alpha = 0.85,
    scale = "width",
    trim = FALSE,
    position = position_dodge(width = 0.5)  # Vertically dodges on y-axis
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
    
  ) +
  scale_color_manual(name = "Flow regime and water year", 
                     values = c("Large inflow dry (WY 2022)" = "#216c8f", "Large inflow wet (WY 2023)" = "#216c8f",
                                "Small inflow dry (WY 2022)" = "#a67d17", "Small inflow wet (WY 2023)" = "#a67d17",
                                "No inflow dry (WY 2022)" = "#09573e", "No inflow wet (WY 2023)" = "#09573e")) +
  scale_fill_manual(name = "Flow regime and water year", 
                    values = c("Large inflow dry (WY 2022)" = "#d3e8f2", "Large inflow wet (WY 2023)" = "#3283a8",
                               "Small inflow dry (WY 2022)" = "#f5e5bc", "Small inflow wet (WY 2023)" = "#a67d17",
                               "No inflow dry (WY 2022)" = "#cce3d6", "No inflow wet (WY 2023)" = "#136F63")) +
  
  labs(
    x = "Estimate",
    y = NULL,
    title = "GPP SEM Coefficients"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 5.5, color = "gray40", linetype = "dotted") +
  geom_hline(yintercept = 3.5,  color = "gray40", linetype = "dotted") +
  geom_hline(yintercept = 1.5,  color = "gray40", linetype = "dotted") +
  
  theme(
    axis.text.y = element_text(size = 14),
    legend.position = "right"
  ) +
  # xlim(-0.85,0.85) +
  theme_minimal() +
  scale_y_discrete(labels = c(
    "b_plotGPP_plot_GPP_lag" = "Metab. ~ AR",
    "b_plotGPP_lake_tempC" = "Metab. ~ Temp.",
    "b_plotGPP_light_mean" = "Metab. ~ Light",
    "b_plotGPP_flow_mean" = "Metab. ~ Flow",
    "b_laketempC_flow_mean" = "Temp ~ Flow",
    "b_laketempC_ppt_mm" = "Temp ~ PPT",
    "b_lightmean_flow_mean" = "Light ~ Flow",
    "b_lightmean_ppt_mm" = "Light ~ PPT",
    "b_flowmean_ppt_mm" = "Flow ~ PPT"
  )) 

coef_violin_plot
