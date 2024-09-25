
##==============================================================================
## SEM of streamflow and precipitation dynamics on nearshore metabolism
## by Loria et al. 2024
## 09/24/2024

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

set.seed(2021)
##===========================================
## read data aggregated data for the project:
#============================================
dat <- readRDS("./NS_analysis_dat.rds")
str(dat)

##===========================================
## create a new df for complete GPP obs 
GPP_df <- dat%>%
  dplyr::select( -middle_ER)%>%
  mutate(middle_GPP_lag = lag(middle_GPP),
         logGPPlag = log(middle_GPP_lag +1),
         logGPP = log(middle_GPP+1),
         loglake_tempC = log(lake_tempC),
         log_streamflow = log(flow_mean +1))

ER_df <- dat%>%
  dplyr::select(-middle_GPP)%>%
  mutate(ER = middle_ER *-1,
         ER_lag = lag(ER),
         logER = log(ER +1),
         logERlag = log(ER_lag +1),
         loglake_tempC = log(lake_tempC),
         log_streamflow = log(flow_mean +1))

NEP_df <- dat%>%
  dplyr::select(-middle_GPP, -middle_ER)%>%
  mutate(NEP_lag = lag(middle_NEP),
         loglake_tempC = log(lake_tempC),
         log_streamflow = log(flow_mean +1))

##============================================
## Stream SEM dataframe organization:
##============================================
## BW POE 20 % = 0.04334224 - log = 0.04242925
## GB POE 20 % = 0.005187275,  - log = 0.005173867

## GPP
## Pick out columns to scale
gpp_columns_to_scale <- c("logGPP", "logGPPlag", "lake_tempC", "light_mean", "windsp_mean","log_streamflow", "ppt_mm")
## BW
GPP_df_BW <- GPP_df %>%
  filter(shore=="BW" & log_streamflow < 0.04243) %>%
  mutate(across(all_of(gpp_columns_to_scale), scale)) %>%
  dplyr::select(site, shore, date, logGPP, logGPPlag, log_streamflow, lake_tempC, light_mean)

BW_GPP_df <- as.data.frame(GPP_df_BW)
BW_GPP_df <- na.omit(BW_GPP_df)

## Check normality... 
hist(BW_GPP_df$logGPP)
hist(BW_GPP_df$log_streamflow)
hist(BW_GPP_df$light_mean)
hist(BW_GPP_df$lake_tempC)


## GB
GPP_df_GB <- GPP_df%>%
  filter(shore=="GB"& log_streamflow < 0.00517) %>%
  mutate(across(all_of(gpp_columns_to_scale), scale))%>%
  dplyr::select(site, shore, date, logGPP, logGPPlag, log_streamflow, lake_tempC, light_mean)

GB_GPP_df <- as.data.frame(GPP_df_GB)
GB_GPP_df <- na.omit(GB_GPP_df)
##===========================================
## ER
## Pick out columns to scale
er_columns_to_scale <- c("logER", "logERlag", "lake_tempC", "light_mean", "windsp_mean","log_streamflow", "ppt_mm")

ER_df_BW <- ER_df %>%
  filter(shore=="BW" & log_streamflow < 0.04243) %>%
  mutate(across(all_of(er_columns_to_scale), scale))%>%
  dplyr::select(site, shore, date, logER, logERlag, log_streamflow, lake_tempC, light_mean)
BW_df_ER <- as.data.frame(ER_df_BW)
BW_df_ER <- na.omit(BW_df_ER)

ER_df_GB <- ER_df %>%
  filter(shore=="GB"& log_streamflow < 0.00517) %>%
 mutate(across(all_of(er_columns_to_scale), scale)) %>%
  dplyr::select(site, shore, date, logER, logERlag, log_streamflow, lake_tempC, light_mean)

GB_df_ER <- as.data.frame(ER_df_GB)
GB_df_ER <- na.omit(GB_df_ER)


#===========================================
## GPP glms for stream influence
#============================================
BW_cor <- (BW_GPP_df[, c(6,7,8)])
chart.Correlation(BW_cor, histogram=TRUE, pch=19) # all less than 0.6 

# Quick visualizations of linear relations: 
plot(BW_GPP_df$lake_tempC~BW_GPP_df$log_streamflow)
summary(lmer(lake_tempC ~ log_streamflow + (1|site), data=BW_GPP_df))

plot(BW_GPP_df$light_mean~BW_GPP_df$log_streamflow)
summary(lmer(light_mean ~ log_streamflow + (1|site), data=BW_GPP_df))

plot(GB_GPP_df$light_mean~GB_GPP_df$log_streamflow)
summary(lmer(light_mean ~ log_streamflow + (1|site), data=GB_GPP_df))


plot(BW_GPP_df$logGPP~BW_GPP_df$log_streamflow)
summary(lmer(logGPP ~ log_streamflow + (1|site), data=BW_GPP_df))

plot(BW_GPP_df$logGPP~BW_GPP_df$lake_tempC)
summary(lmer(logGPP ~ lake_tempC + (1|site), data=BW_GPP_df))

plot(BW_GPP_df$logGPP~BW_GPP_df$light_mean)
summary(lmer(logGPP ~ light_mean + (1|site), data=BW_GPP_df))

#===========================================
## SEM for GPP and streamflow
#============================================
GPP_mod <- bf(logGPP ~ logGPPlag + (light_mean) + (lake_tempC) + (log_streamflow) +(1|site)) 
temp_mod <- bf((lake_tempC) ~ (log_streamflow) + (1|site)) #cover_mod
light_mod <- bf((light_mean) ~ (log_streamflow) + (1|site)) #cover_mod

########
## BW
BW_GPP_stream_fit <- brm(GPP_mod +
                    temp_mod +
                    light_mod +
                    set_rescor(FALSE),
                  data=BW_GPP_df,
                  iter=10000, warmup = 5000,
                  control = list(adapt_delta = 0.95, max_treedepth = 15), 
                  cores=4, chains = 3)

plot(BW_GPP_stream_fit)
BW_GPP_stream_fit_sum <- summary(BW_GPP_stream_fit)

# Extract posterior samples BW:
posterior_samples_BW <- posterior_samples(BW_GPP_stream_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BW <- posterior_samples_BW %>%
  select(b_logGPP_logGPPlag, b_logGPP_light_mean, b_logGPP_lake_tempC, b_logGPP_log_streamflow,
         b_XlaketempC_log_streamflow, b_Xlightmean_log_streamflow) %>%
  mutate(shore = shore_levels)
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW <- melt(posterior_samples_BW, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

########
## GB
GB_GPP_stream_fit <- brm(GPP_mod +
                           temp_mod +
                           light_mod +
                           set_rescor(FALSE),
                         data=GB_GPP_df,
                         iter=10000, warmup = 5000,
                         control = list(adapt_delta = 0.95, max_treedepth = 15), 
                         cores=4, chains = 3)

plot(GB_GPP_stream_fit)
summary(GB_GPP_stream_fit)
GB_GPP_stream_fit_sum <- summary(GB_GPP_stream_fit)


# Extract posterior samples GB:
posterior_samples_GB <- posterior_samples(GB_GPP_stream_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))
# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  select(b_logGPP_logGPPlag, b_logGPP_light_mean, b_logGPP_lake_tempC, b_logGPP_log_streamflow,
         b_XlaketempC_log_streamflow, b_Xlightmean_log_streamflow) %>%
  mutate(shore = shore_levels)
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

#============================================
SEM_df_GPP <- rbind(posterior_samples_melted_GB, posterior_samples_melted_BW)

ci_summary_overall <- SEM_df_GPP %>%
  group_by(Coefficient, shore) %>%
  summarize(
    mean = mean(Estimate),
    l_95 = quantile(Estimate, 0.025),
    u_95 = quantile(Estimate, 0.975)
  )

# Plot histograms with credible intervals and mean estimates
site_semplot_GPP <-ggplot(SEM_df_GPP, aes(x = Estimate, fill = shore, color = shore)) +
  geom_density(alpha = 0.75, linewidth = 0.25) +
  geom_vline(xintercept=0, linetype = "solid", size = 0.5, color = "black") +
  geom_vline(data = ci_summary_overall, aes(xintercept = mean, color=shore), linetype = "solid", size = 0.5) +
  scale_fill_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  labs(
    x = 'Estimate',
    y = 'Density') +
  facet_wrap(Coefficient~., ncol=3, scales = 'free') 

# ggsave(plot = site_semplot_GPP, filename = paste("./NS_GPP_stream_semplot.png",sep=""),width=6.75,height=3.75,dpi=300)

#===========================================
## SEM for ER and streamflow
#============================================
ER_mod <- bf(logER ~ logERlag + (light_mean) + (lake_tempC) + (log_streamflow) +(1|site)) 
temp_mod <- bf((lake_tempC) ~ (log_streamflow) + (1|site)) #cover_mod
light_mod <- bf((light_mean) ~ (log_streamflow) + (1|site)) #cover_mod

########
## BW
BW_ER_stream_fit <- brm(ER_mod +
                           temp_mod +
                           light_mod +
                           set_rescor(FALSE),
                         data=BW_df_ER,
                        iter=10000, warmup = 5000,
                        control = list(adapt_delta = 0.95, max_treedepth = 15), 
                        cores=4, chains = 3)


plot(BW_ER_stream_fit)
summary(BW_ER_stream_fit)
BW_ER_stream_fit_sum <- summary(BW_ER_stream_fit)


# Extract posterior samples BW:
posterior_samples_BW <- posterior_samples(BW_ER_stream_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BWq <- posterior_samples_BW %>%
  select(b_logER_logERlag, b_logER_light_mean, b_logER_lake_tempC, b_logER_log_streamflow,
         b_XlaketempC_log_streamflow, b_Xlightmean_log_streamflow) %>%
  mutate(shore = shore_levels)
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW <- melt(posterior_samples_BWq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

########
## ER
GB_ER_stream_fit <- brm(ER_mod +
                         temp_mod +
                         light_mod +
                         set_rescor(FALSE),
                       data=GB_df_ER,
                       iter=10000, warmup = 5000,
                       control = list(adapt_delta = 0.95, max_treedepth = 15), 
                       cores=4, chains = 3)

plot(GB_ER_stream_fit)
summary(GB_ER_stream_fit)
GB_ER_stream_fit_sum <- summary(GB_ER_stream_fit)


# Extract posterior samples GB:
posterior_samples_GB <- posterior_samples(GB_ER_stream_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))
# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  select(b_logER_logERlag, b_logER_light_mean, b_logER_lake_tempC, b_logER_log_streamflow,
         b_XlaketempC_log_streamflow, b_Xlightmean_log_streamflow) %>%
  mutate(shore = shore_levels)
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

##

SEM_df_ER <- rbind(posterior_samples_melted_GB, posterior_samples_melted_BW)

ci_summary_overall <- SEM_df_ER %>%
  group_by(Coefficient, shore) %>%
  summarize(
    mean = mean(Estimate),
    l_95 = quantile(Estimate, 0.025),
    u_95 = quantile(Estimate, 0.975)
  )

# Plot histograms with credible intervals and mean estimates
site_semplot_ER <-ggplot(SEM_df_ER, aes(x = Estimate, fill = shore, color = shore)) +
  geom_density(alpha = 0.75, linewidth = 0.25) +
  geom_vline(xintercept=0, linetype = "solid", size = 0.5, color = "black") +
  geom_vline(data = ci_summary_overall, aes(xintercept = mean, color=shore), linetype = "solid", size = 0.5) +
  scale_fill_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  labs(
    x = 'Estimate',
    y = 'Density') +
  facet_wrap(Coefficient~., ncol=3, scales = 'free') 

# ggsave(plot = site_semplot_ER, filename = paste("./NS_ER_stream_semplot.png",sep=""),width=6.75,height=3.75,dpi=300)

#===========================================
## SEM for NEP and precip.
#============================================
NEP_23 <- NEP_df %>%
  filter(date>as.Date("2023-02-07"))

## model:
ppt_mod <- bf(middle_NEP ~ NEP_lag + scale(light_mean) + scale(lake_tempC) + scale(ppt_mm) + (1|site)) 

#####
## BW
ppt_brms_BW_fit <- brm(ppt_mod, 
                   data=NEP_23%>%filter(shore=="BW"),
                   iter=10000, warmup = 5000,
                   control = list(adapt_delta = 0.95, max_treedepth = 15), 
                   cores=4, chains = 3)
plot(ppt_brms_BW_fit)
summary(ppt_brms_BW_fit)

# Extract posterior samples
posterior_samples_BW <- posterior_samples(ppt_brms_BW_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BWq <- posterior_samples_BW %>%
   select(b_NEP_lag, b_scalelake_tempC, b_scalelight_mean, b_scaleppt_mm) %>%
  mutate(shore = shore_levels)
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW <- melt(posterior_samples_BWq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

#####
## GB
ppt_brms_GB_fit <- brm(ppt_mod, 
                       data=NEP_23%>%filter(shore=="GB"),
                       iter=10000, warmup = 5000,
                       control = list(adapt_delta = 0.95, max_treedepth = 15), 
                       cores=4, chains = 3)
plot(ppt_brms_GB_fit)
summary(ppt_brms_GB_fit)

# Extract posterior samples
posterior_samples_GB <- posterior_samples(ppt_brms_GB_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))
# Prepare the posterior samples data frame with shore levels
posterior_samples_GB <- posterior_samples_GB %>%
  select(b_NEP_lag, b_scalelake_tempC, b_scalelight_mean, b_scaleppt_mm) %>%
  mutate(shore = shore_levels)
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB <- melt(posterior_samples_GB, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

#####
## SS
ppt_brms_SS_fit <- brm(ppt_mod, 
                       data=NEP_23%>%filter(shore=="SS"),
                       iter=10000, warmup = 5000,
                       control = list(adapt_delta = 0.95, max_treedepth = 15), 
                       cores=4, chains = 3)
plot(ppt_brms_SS_fit)
summary(ppt_brms_SS_fit)

# Extract posterior samples
posterior_samples_SS <- posterior_samples(ppt_brms_SS_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SS <- posterior_samples_SS %>%
  select(b_NEP_lag, b_scalelake_tempC, b_scalelight_mean, b_scaleppt_mm) %>%
  mutate(shore = shore_levels)
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SS <- melt(posterior_samples_SS, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

#####
## SH
ppt_brms_SH_fit <- brm(ppt_mod, 
                       data=NEP_23%>%filter(shore=="SH"),
                       iter=10000, warmup = 5000,
                       control = list(adapt_delta = 0.95, max_treedepth = 15), 
                       cores=4, chains = 3)
plot(ppt_brms_SH_fit)
summary(ppt_brms_SH_fit)

# Extract posterior samples
posterior_samples_SH <- posterior_samples(ppt_brms_SH_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SH"), length.out = nrow(ppt_brms_SH_fit))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SH <- posterior_samples_SH %>%
  select(b_NEP_lag, b_scalelake_tempC, b_scalelight_mean, b_scaleppt_mm) %>%
  mutate(shore = shore_levels)
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SH <- melt(posterior_samples_SH, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

## All 
SEM_df <- rbind(posterior_samples_melted_GB, posterior_samples_melted_BW, posterior_samples_melted_SS, posterior_samples_melted_SH)

ci_summary_overall <- SEM_df %>%
  group_by(Coefficient, shore) %>%
  summarize(
    mean = mean(Estimate),
    l_95 = quantile(Estimate, 0.025),
    u_95 = quantile(Estimate, 0.975)
  )

# Plot histograms with credible intervals and mean estimates
ppt_semplot <-ggplot(SEM_df, aes(x = Estimate, fill = shore, color = shore)) +
  geom_density(alpha = 0.75, linewidth = 0.25) +
  geom_vline(xintercept=0, linetype = "solid", size = 0.5, color = "black") +
  geom_vline(data = ci_summary_overall, aes(xintercept = mean, color=shore), linetype = "solid", size = 0.5) +
  scale_fill_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  labs(
       x = NULL,
       y = 'Density') +
  facet_wrap(Coefficient~., ncol=5, scales = 'free') 
  
# ggsave(plot = ppt_semplot, filename = paste("./NS24_ppt_NEP.png",sep=""),width=10,height=2.25,dpi=300)

# end script.