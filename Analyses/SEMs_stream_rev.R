
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

set.seed(2025)
##===========================================
## read data aggregated data for the project:
#============================================
# dat <- readRDS("./NS_analysis_dat.rds")
dat <-readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/Final_Scripts/NS_analysis_dat.rds")
str(dat)


##===========================================
## create a new df for complete GPP obs 
GPP_df <- dat%>%
  dplyr::select(-middle_ER)

ER_df <- dat%>%
  dplyr::select(-middle_GPP)%>%
  mutate(ER = middle_ER *-1)

NEP_df <- dat%>%
  dplyr::select(-middle_GPP, -middle_ER)

#################
### OVER LAP 
# Filter the dataset to only include rows where do.obs_m or GPP_mean are not NA
GPP_df_BW <- GPP_df%>%
  filter(shore == "BW")

summary(GPP_df_BW)
# date range 2021-06-11  to 2023-09-14   


BW_GPP_df21 <- GPP_df_BW %>%
  filter(date > as.Date("2021-06-10") & date < as.Date("2021-09-10"))

summary(BW_GPP_df21)
# date range 2021-06-11  to 2023-09-14   

plot_data_DO <- GPP_df_BW %>%
  filter(!is.na(lake_DO)) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_GPP = !is.na(middle_GPP),
         plot_GPP = ifelse(is.na(middle_GPP), 0, middle_GPP),
         yday=yday(date),
         year=year(date))

# Identify the unique yday values for each year
yday_2021 <- unique(plot_data_DO$yday[plot_data_DO$year == 2021])
yday_2022 <- unique(plot_data_DO$yday[plot_data_DO$year == 2022])
yday_2023 <- unique(plot_data_DO$yday[plot_data_DO$year == 2023])

# Find the common yday values across all three years
common_yday <- Reduce(intersect, list(yday_2021, yday_2022, yday_2023))

# Filter the dataframe for only those common yday values
filtered_data <- plot_data_DO %>%
  filter(yday %in% common_yday & year %in% c(2021, 2022, 2023))

filtered_data_cs <- filtered_data%>%
  arrange(year, yday) %>%
  group_by(year) %>%
  mutate(GPPcumsum = cumsum(plot_GPP))

# Determine the date range
date_range <- range(plot_data_DO$date, na.rm = TRUE)


# Create the number line plot
templt <- ggplot(filtered_data, aes(x = yday)) +
  #geom_segment(aes(x = date_range[1], xend = date_range[2], y = 0, yend = 0), linewidth = 0.2) + # Main line
  geom_point(data = filter(filtered_data, has_do), aes(y = 0.095), alpha=0.8, color = "black", size = 1) + # DO presence tick marks
  geom_point(data = filter(filtered_data, has_GPP), aes(y = 0.090), alpha= 0.8, color = "#579467", size = 1) + # GPP presence tick marks
  scale_y_continuous(name = "", breaks = NULL) +
  theme_classic() + xlab(NULL) +
  labs(title = "GB Upper - GPP") +facet_grid(year~.)
templt

###
# ER
ER_df_BW <- ER_df%>%
  filter(shore == "BW")

plot_data_DO_BW_ER <- ER_df_BW %>%
  filter(!is.na(lake_DO)) %>%
  #dplyr::select(date, DO_mgL_calibrated, middle_GPP) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_ER = !is.na(ER),
         plot_ER = ifelse(is.na(ER), 0, ER),
         yday=yday(date),
         year=year(date))

#############

GPP_df_GB <- GPP_df%>%
  filter(shore == "GB")

summary(GPP_df_GB)

tempcheck <- GPP_df_GB %>%
  filter(date > as.Date("2021-06-10") & date < as.Date("2021-09-10"))
summary(tempcheck)

plot_data_DO_gb <- GPP_df_GB %>%
  filter(!is.na(lake_DO)) %>%
  #dplyr::select(date, DO_mgL_calibrated, middle_GPP) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_GPP = !is.na(middle_GPP),
         plot_GPP = ifelse(is.na(middle_GPP), 0, middle_GPP),
         yday=yday(date),
         year=year(date))

# Identify the unique yday values for each year
yday_2021 <- unique(plot_data_DO_gb$yday[plot_data_DO_gb$year == 2021])
yday_2022 <- unique(plot_data_DO_gb$yday[plot_data_DO_gb$year == 2022])
yday_2023 <- unique(plot_data_DO_gb$yday[plot_data_DO_gb$year == 2023])

# Find the common yday values across all three years
common_yday <- Reduce(intersect, list(yday_2021, yday_2022, yday_2023))

# Filter the dataframe for only those common yday values
filtered_data <- plot_data_DO_gb %>%
  filter(yday %in% common_yday & year %in% c(2021, 2022, 2023))

filtered_data_cs <- filtered_data%>%
  arrange(year, yday) %>%
  group_by(year) %>%
  mutate(GPPcumsum = cumsum(plot_GPP))

# Determine the date range
date_range <- range(plot_data_DO$date, na.rm = TRUE)


# Create the number line plot
templt <- ggplot(filtered_data, aes(x = yday)) +
  #geom_segment(aes(x = date_range[1], xend = date_range[2], y = 0, yend = 0), linewidth = 0.2) + # Main line
  geom_point(data = filter(filtered_data, has_do), aes(y = 0.095), alpha=0.8, color = "black", size = 1) + # DO presence tick marks
  geom_point(data = filter(filtered_data, has_GPP), aes(y = 0.090), alpha= 0.8, color = "#579467", size = 1) + # GPP presence tick marks
  scale_y_continuous(name = "", breaks = NULL) +
  theme_classic() + xlab(NULL) +
  labs(title = "GB Upper - GPP") +facet_grid(year~.)
templt
#############

###
# ER
ER_df_GB <- ER_df%>%
  filter(shore == "GB")

plot_data_DO_GB_ER <- ER_df_GB %>%
  filter(!is.na(lake_DO)) %>%
  #dplyr::select(date, DO_mgL_calibrated, middle_GPP) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_ER = !is.na(ER),
         plot_ER = ifelse(is.na(ER), 0, ER),
         yday=yday(date),
         year=year(date))



############# SS

GPP_df_SS <- GPP_df%>%
  filter(shore == "SS")

summary(GPP_df_SS)

tempcheck <- GPP_df_SS %>%
  filter(date > as.Date("2021-06-10") & date < as.Date("2021-09-10"))
summary(tempcheck)


plot_data_DO_SS <- GPP_df_SS %>%
  filter(!is.na(lake_DO)) %>%
  #dplyr::select(date, DO_mgL_calibrated, middle_GPP) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_GPP = !is.na(middle_GPP),
         plot_GPP = ifelse(is.na(middle_GPP), 0, middle_GPP),
         yday=yday(date),
         year=year(date))

# Identify the unique yday values for each year
yday_2021 <- unique(plot_data_DO_SS$yday[plot_data_DO_SS$year == 2021])
yday_2022 <- unique(plot_data_DO_SS$yday[plot_data_DO_SS$year == 2022])
yday_2023 <- unique(plot_data_DO_SS$yday[plot_data_DO_SS$year == 2023])

# Find the common yday values across all three years
common_yday <- Reduce(intersect, list(yday_2021, yday_2022, yday_2023))


###
# ER
ER_df_SS <- ER_df%>%
  filter(shore == "SS")

plot_data_DO_SS_ER <- ER_df_SS %>%
  filter(!is.na(lake_DO)) %>%
  #dplyr::select(date, DO_mgL_calibrated, middle_GPP) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_ER = !is.na(ER),
         plot_ER = ifelse(is.na(ER), 0, ER),
         yday=yday(date),
         year=year(date))



############# SH

GPP_df_SH <- GPP_df%>%
  filter(shore == "SH")

plot_data_DO_SH <- GPP_df_SH %>%
  filter(!is.na(lake_DO)) %>%
  #dplyr::select(date, DO_mgL_calibrated, middle_GPP) %>%
  mutate(has_do = !is.na(lake_DO), 
         has_GPP = !is.na(middle_GPP),
         plot_GPP = ifelse(is.na(middle_GPP), 0, middle_GPP),
         yday=yday(date),
         year=year(date))

# Identify the unique yday values for each year
yday_2023 <- unique(plot_data_DO_SH$yday[plot_data_DO_SH$year == 2023])

###
# ER
ER_df_SH <- ER_df%>%
  filter(shore == "SH")

plot_data_DO_SH_ER <- ER_df_SH %>%
  filter(!is.na(lake_DO)) %>%
  #dplyr::select(date, DO_mgL_calibrated, middle_GPP) %>%
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

gpp_columns_to_scale <- c("plot_GPP", "lake_tempC", "light_mean", "windsp_mean","flow_mean", "ppt_mm")

er_columns_to_scale <- c("plot_ER", "lake_tempC", "light_mean", "windsp_mean","flow_mean", "ppt_mm")

## =====
## BW:

## Check normality... 
hist(plot_data_DO$middle_GPP)
hist(log(plot_data_DO$middle_GPP +1))
hist(plot_data_DO$lake_tempC)
hist(log(plot_data_DO$lake_tempC +1))
hist(plot_data_DO$light_mean)
hist(log(plot_data_DO$lake_tempC))
hist(plot_data_DO$flow_mean)
hist(log(plot_data_DO$flow_mean+1))
hist(plot_data_DO$ppt_mm)
hist(log(plot_data_DO$ppt_mm+1))


GPP_df_BW_log <- plot_data_DO%>%
  drop_na(plot_GPP, lake_tempC, light_mean, flow_mean, ppt_mm) %>% 
  filter(flow_mean<0.043) %>%
  mutate(plot_GPP= log(plot_GPP+1),
         flow_mean= log(flow_mean+1),
         lake_tempC= log(lake_tempC),
         light_mean= log(light_mean),
         ppt_mm= log(ppt_mm +1)) %>%
  #mutate(across(all_of(gpp_columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_GPP_lag = lag(plot_GPP)) %>%
  dplyr::select(site, shore, date, plot_GPP, plot_GPP_lag, flow_mean, lake_tempC, light_mean, ppt_mm) %>%
  drop_na(plot_GPP, plot_GPP_lag, lake_tempC, light_mean, flow_mean, ppt_mm) 
  

ER_df_BW_log <- plot_data_DO_BW_ER%>%
  drop_na(plot_ER, lake_tempC, light_mean, flow_mean, ppt_mm) %>% 
  filter(flow_mean<0.043) %>%
  mutate(plot_ER= log(plot_ER+1),
         flow_mean= log(flow_mean+1),
         lake_tempC= log(lake_tempC),
         light_mean= log(light_mean),
         ppt_mm= log(ppt_mm +1)) %>%
  #mutate(across(all_of(er_columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_ER_lag = lag(plot_ER)) %>%
  dplyr::select(site, shore, date, plot_ER, plot_ER_lag, flow_mean, lake_tempC, light_mean, ppt_mm) %>%
  drop_na(plot_ER, plot_ER_lag, lake_tempC, light_mean, flow_mean, ppt_mm) 


## =====
## GB:
GPP_df_GB_log <- plot_data_DO_gb%>%
  drop_na(plot_GPP, lake_tempC, light_mean, flow_mean, ppt_mm) %>% 
  filter(flow_mean<0.0052) %>%
  mutate(plot_GPP= log(plot_GPP+1),
         flow_mean= log(flow_mean+1),
         lake_tempC= log(lake_tempC),
         light_mean= log(light_mean),
         ppt_mm= log(ppt_mm +1)) %>%
 # mutate(across(all_of(gpp_columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_GPP_lag = lag(plot_GPP)) %>%
  dplyr::select(site, shore, date, plot_GPP, plot_GPP_lag, flow_mean, lake_tempC, light_mean, ppt_mm) %>%
  drop_na(plot_GPP, plot_GPP_lag, lake_tempC, light_mean, flow_mean, ppt_mm) 

hist(GPP_df_GB_log$plot_GPP)
hist(GPP_df_GB_log$lake_tempC)
hist(GPP_df_GB_log$light_mean)
hist(GPP_df_GB_log$flow_mean)
hist(GPP_df_GB_log$ppt_mm)



ER_df_GB_log <- plot_data_DO_GB_ER%>%
  drop_na(plot_ER, lake_tempC, light_mean, flow_mean, ppt_mm) %>% 
  filter(flow_mean<0.0052) %>%
  mutate(plot_ER= log(plot_ER+1),
         flow_mean= log(flow_mean+1),
         lake_tempC= log(lake_tempC),
         light_mean= log(light_mean),
         ppt_mm= log(ppt_mm +1)) %>%
  #mutate(across(all_of(er_columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_ER_lag = lag(plot_ER)) %>%
  dplyr::select(site, shore, date, plot_ER, plot_ER_lag, flow_mean, lake_tempC, light_mean, ppt_mm) %>%
  drop_na(plot_ER, plot_ER_lag, lake_tempC, light_mean, flow_mean, ppt_mm) 


## =====
## SS:
GPP_df_SS_log <- plot_data_DO_SS%>%
  drop_na(plot_GPP, lake_tempC, light_mean, ppt_mm) %>% 
  mutate(plot_GPP= log(plot_GPP+1),
         lake_tempC= log(lake_tempC),
         light_mean= log(light_mean),
         ppt_mm= log(ppt_mm +1)) %>%
 # mutate(across(all_of(gpp_columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_GPP_lag = lag(plot_GPP)) %>%
  dplyr::select(site, shore, date, plot_GPP, plot_GPP_lag, lake_tempC, light_mean, ppt_mm) %>%
  drop_na(plot_GPP, plot_GPP_lag, lake_tempC, light_mean, ppt_mm) 

hist(GPP_df_SS_log$plot_GPP)
hist(GPP_df_SS_log$lake_tempC)
hist(GPP_df_SS_log$light_mean)
hist(GPP_df_SS_log$ppt_mm)


ER_df_SS_log <- plot_data_DO_SS_ER %>%
  drop_na(plot_ER, lake_tempC, light_mean, ppt_mm) %>% 
  mutate(plot_ER= log(plot_ER+1),
         lake_tempC= log(lake_tempC),
         light_mean= log(light_mean),
         ppt_mm= log(ppt_mm +1)) %>%
  #mutate(across(all_of(er_columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_ER_lag = lag(plot_ER)) %>%
  dplyr::select(site, shore, date, plot_ER, plot_ER_lag, lake_tempC, light_mean, ppt_mm) %>%
  drop_na(plot_ER, plot_ER_lag, lake_tempC, light_mean, ppt_mm) 



## =====
## SH:
GPP_df_SH_log <- plot_data_DO_SH%>%
  drop_na(plot_GPP, lake_tempC, light_mean, ppt_mm) %>% 
  mutate(plot_GPP= log(plot_GPP+1),
         lake_tempC= log(lake_tempC),
         light_mean= log(light_mean),
         ppt_mm= log(ppt_mm +1)) %>%
  #mutate(across(all_of(gpp_columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_GPP_lag = lag(plot_GPP)) %>%
  dplyr::select(site, shore, date, plot_GPP, plot_GPP_lag, lake_tempC, light_mean, ppt_mm) %>%
  drop_na(plot_GPP, plot_GPP_lag, lake_tempC, light_mean, ppt_mm) 

hist(GPP_df_SH_log$plot_GPP)
hist(GPP_df_SH_log$lake_tempC)
hist(GPP_df_SH_log$light_mean)
hist(GPP_df_SH_log$ppt_mm)


ER_df_SH_log <- plot_data_DO_SH_ER %>%
  drop_na(plot_ER, lake_tempC, light_mean, ppt_mm) %>% 
  mutate(plot_ER= log(plot_ER+1),
         lake_tempC= log(lake_tempC),
         light_mean= log(light_mean),
         ppt_mm= log(ppt_mm +1)) %>%
 # mutate(across(all_of(er_columns_to_scale), scale))%>%
  group_by(site) %>%
  mutate(plot_ER_lag = lag(plot_ER)) %>%
  dplyr::select(site, shore, date, plot_ER, plot_ER_lag, lake_tempC, light_mean, ppt_mm) %>%
  drop_na(plot_ER, plot_ER_lag, lake_tempC, light_mean, ppt_mm) 


###########################
#### revised GPP model ###
##########################

GPP_mod <- bf(plot_GPP ~ plot_GPP_lag + light_mean + lake_tempC + flow_mean + (1|site)) 
temp_mod <- bf(lake_tempC ~ flow_mean + ppt_mm + (1|site)) 
light_mod <- bf(light_mean ~ flow_mean + ppt_mm + (1|site)) 

########
BW_GPP_df21 <- GPP_df_BW_log %>%
  filter(date > as.Date("2021-06-10") & date < as.Date("2021-09-30"))


## BW
BW_GPP_stream_fit <- brm(GPP_mod +
                           temp_mod +
                           light_mod +
                           set_rescor(FALSE),
                         data=BW_GPP_df21,
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
  select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW_21 <- melt(posterior_samples_BW, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_BW_21, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/posterior_samples_melted_BW_21.rds")

########
GB_GPP_df21 <- GPP_df_GB_log %>%
  filter(date > as.Date("2021-06-10") & date < as.Date("2021-09-30"))

## GB
GB_GPP_stream_fit <- brm(GPP_mod +
                           temp_mod +
                           light_mod +
                           set_rescor(FALSE),
                         data=GB_GPP_df21,
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
  select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)


# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB_21 <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_GB_21, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/posterior_samples_melted_GB_21.rds")


########
SS_GPP_df21 <- GPP_df_SS_log %>%
  filter(date > as.Date("2021-08-01") & date < as.Date("2021-09-30"))

## GB
SS_GPP_stream_fit <- brm(GPP_mod +
                           temp_mod +
                           light_mod +
                           set_rescor(FALSE),
                         data=SS_GPP_df21,
                         iter=10000, warmup = 5000,
                         control = list(adapt_delta = 0.95, max_treedepth = 15), 
                         cores=4, chains = 3)

plot(SS_GPP_stream_fit)
summary(SS_GPP_stream_fit)
SS_GPP_stream_fit_sum <- summary(SS_GPP_stream_fit)


# Extract posterior samples SS:
posterior_samples_SS <- posterior_samples(SS_GPP_stream_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SS_21 <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")



#============================================
SEM_df_GPP_21 <- rbind(posterior_samples_melted_BW_21, posterior_samples_melted_GB_21)

ci_summary_overall <- SEM_df_GPP_21 %>%
  group_by(Coefficient, shore) %>%
  summarize(
    mean = mean(Estimate),
    l_95 = quantile(Estimate, 0.025),
    u_95 = quantile(Estimate, 0.975)
  )

# Plot histograms with credible intervals and mean estimates
site_semplot_GPP_21 <-ggplot(SEM_df_GPP_21, aes(x = Estimate, fill = shore, color = shore)) +
  geom_density(alpha = 0.75, linewidth = 0.25) +
  geom_vline(xintercept=0, linetype = "solid", size = 0.5, color = "black") +
  geom_vline(data = ci_summary_overall, aes(xintercept = mean, color=shore), linetype = "solid", size = 0.5) +
  scale_fill_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  labs(
    x = 'Estimate',
    y = 'Density') +
  facet_wrap(Coefficient~., ncol=8, scales = 'free') 

# ggsave(plot = site_semplot_GPP_21, filename = paste("./NS_GPP_stream_semplot_21.png",sep=""),width=12.75,height=2.5,dpi=300)


BW_GPP_df22 <- GPP_df_BW_log %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

########
## BW
BW_GPP_stream_fit_22 <- brm(GPP_mod +
                              temp_mod +
                              light_mod +
                              set_rescor(FALSE),
                            data=BW_GPP_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 15), 
                            cores=4, chains = 3)

plot(BW_GPP_stream_fit_22)
BW_GPP_stream_fit_sum <- summary(BW_GPP_stream_fit_22)


# Extract posterior samples BW:
posterior_samples_BW <- posterior_samples(BW_GPP_stream_fit_22)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BW <- posterior_samples_BW %>%
  select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW_22 <- melt(posterior_samples_BW, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


# saveRDS(posterior_samples_melted_BW_22, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/posterior_samples_melted_BW_22.rds")


########
GB_GPP_df22 <- GPP_df_GB_log %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

hist(GB_GPP_df22$flow_mean)
## GB
GB_GPP_stream_fit_22 <- brm(GPP_mod +
                              temp_mod +
                              light_mod +
                              set_rescor(FALSE),
                            data=GB_GPP_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 15), 
                            cores=4, chains = 3)

plot(GB_GPP_stream_fit_22)
summary(GB_GPP_stream_fit_22)
GB_GPP_stream_fit_sum <- summary(GB_GPP_stream_fit_22)


# Extract posterior samples GB:
posterior_samples_GB <- posterior_samples(GB_GPP_stream_fit_22)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))
# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)
         
         
# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB_22 <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


# saveRDS(posterior_samples_melted_GB_22, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/posterior_samples_melted_GB_22.rds")


########
GPP_mod1 <- bf(plot_GPP ~ plot_GPP_lag + light_mean + lake_tempC) 
temp_mod1 <- bf(lake_tempC ~  ppt_mm) 
light_mod1 <- bf(light_mean ~ ppt_mm) 


GPP_mod2 <- bf(plot_GPP ~ plot_GPP_lag + light_mean + lake_tempC + (1|site)) 
temp_mod2 <- bf(lake_tempC ~  ppt_mm + (1|site)) 
light_mod2 <- bf(light_mean ~ ppt_mm + (1|site)) 


SS_GPP_df22 <- GPP_df_SS_log %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

## GB
SS_GPP_stream_fit_22 <- brm(GPP_mod1 +
                              temp_mod1 +
                              light_mod1 +
                              set_rescor(FALSE),
                            data=SS_GPP_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 17), 
                            cores=3, chains = 3)

plot(SS_GPP_stream_fit_22)
summary(SS_GPP_stream_fit_22)
SS_GPP_stream_fit_sum <- summary(SS_GPP_stream_fit_22)


# Extract posterior samples SS:
posterior_samples_SS <- posterior_samples(SS_GPP_stream_fit_22)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, 
         b_laketempC_ppt_mm,
         b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SS_22 <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


# saveRDS(posterior_samples_melted_SS_22, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/posterior_samples_melted_SS_22.rds")


#============================================
SEM_df_GPP_22 <- rbind(posterior_samples_melted_GB_22, posterior_samples_melted_BW_22, posterior_samples_melted_SS_22)

ci_summary_overall <- SEM_df_GPP_22 %>%
  filter(!(shore == "SS" & Coefficient %in% c("b_plotGPP_flow_mean", "b_laketempC_flow_mean", "b_lightmean_flow_mean")))%>%
  group_by(Coefficient, shore) %>%
  summarize(
    mean = mean(Estimate),
    l_95 = quantile(Estimate, 0.025),
    u_95 = quantile(Estimate, 0.975)
  )

# Plot histograms with credible intervals and mean estimates
site_semplot_GPP_22 <-ggplot(SEM_df_GPP_22 %>%
                               filter(!(shore == "SS" & Coefficient %in% c("b_plotGPP_flow_mean", "b_laketempC_flow_mean", "b_lightmean_flow_mean"))),
                             aes(x = Estimate, fill = shore, color = shore)) +
  geom_density(alpha = 0.75, linewidth = 0.25) +
  geom_vline(xintercept=0, linetype = "solid", size = 0.5, color = "black") +
  geom_vline(data = ci_summary_overall, aes(xintercept = mean, color=shore), linetype = "solid", size = 0.5) +
  scale_fill_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  labs(
    x = 'Estimate',
    y = 'Density') +
  facet_wrap(Coefficient~., ncol=8, scales = 'free') 

# ggsave(plot = site_semplot_GPP_22, filename = paste("./NS_GPP_stream_semplot_22.png",sep=""),width=12.75,height=2.5,dpi=300)



####### 2023


BW_GPP_df23 <- GPP_df_BW_log %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-14"))
summary(BW_GPP_df23)


########
## BW
BW_GPP_stream_fit_23 <- brm(GPP_mod +
                              temp_mod +
                              light_mod +
                              set_rescor(FALSE),
                            data=BW_GPP_df23,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 15), 
                            cores=4, chains = 3)

plot(BW_GPP_stream_fit_23)
BW_GPP_stream_fit_sum <- summary(BW_GPP_stream_fit_23)


# Extract posterior samples BW:
posterior_samples_BW <- posterior_samples(BW_GPP_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BW <- posterior_samples_BW %>%
  select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW_23 <- melt(posterior_samples_BW, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


# saveRDS(posterior_samples_melted_BW_23, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/posterior_samples_melted_BW_23.rds")


########
GB_GPP_df23 <- GPP_df_GB_log %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-14"))
summary(GB_GPP_df23)

## GB
GB_GPP_stream_fit_23 <- brm(GPP_mod +
                              temp_mod +
                              light_mod +
                              set_rescor(FALSE),
                            data=GB_GPP_df23,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 16), 
                            cores=4, chains = 3)

plot(GB_GPP_stream_fit_23)
summary(GB_GPP_stream_fit_23)
GB_GPP_stream_fit_sum <- summary(GB_GPP_stream_fit_23)


# Extract posterior samples GB:
posterior_samples_GB <- posterior_samples(GB_GPP_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))
# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB_23 <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_GB_23, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/posterior_samples_melted_GB_23.rds")


########
SS_GPP_df23 <- GPP_df_SS_log %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-14"))

## GB
SS_GPP_stream_fit_23 <- brm(GPP_mod2 +
                              temp_mod2 +
                              light_mod2 +
                              set_rescor(FALSE),
                            data=SS_GPP_df23,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 15), 
                            cores=4, chains = 3)

plot(SS_GPP_stream_fit_23)
summary(SS_GPP_stream_fit_23)
SS_GPP_stream_fit_sum <- summary(SS_GPP_stream_fit_23)


# Extract posterior samples SS:
posterior_samples_SS <- posterior_samples(SS_GPP_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC,
          b_laketempC_ppt_mm,
         , b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)


# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SS_23 <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(posterior_samples_melted_SS_23, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/posterior_samples_melted_SS_23.rds")


########
SH_GPP_df23 <- GPP_df_SH_log %>%
  filter(date > as.Date("2023-01-01") & date < as.Date("2023-09-10"))

## GB
SH_GPP_stream_fit_23 <- brm(GPP_mod +
                              temp_mod +
                              light_mod +
                              set_rescor(FALSE),
                            data=SH_GPP_df23,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 15), 
                            cores=4, chains = 3)

plot(SH_GPP_stream_fit_23)
summary(SH_GPP_stream_fit_23)
SH_GPP_stream_fit_sum <- summary(SH_GPP_stream_fit_23)


# Extract posterior samples SS:
posterior_samples_SH <- posterior_samples(SH_GPP_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SH"), length.out = nrow(posterior_samples_SH))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SHq <- posterior_samples_SH %>%
  select(b_plotGPP_plot_GPP_lag, b_plotGPP_light_mean, b_plotGPP_lake_tempC, b_plotGPP_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)


# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SH_23 <- melt(posterior_samples_SHq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")



#============================================
SEM_df_GPP_23 <- rbind(posterior_samples_melted_GB_23, posterior_samples_melted_BW_23,
                       posterior_samples_melted_SS_23, posterior_samples_melted_SH_23)

ci_summary_overall <- SEM_df_GPP_23 %>%
  group_by(Coefficient, shore) %>%
  filter(!(shore == "SS" & Coefficient %in% c("b_plotGPP_flow_mean", "b_laketempC_flow_mean", "b_lightmean_flow_mean")))%>%
  filter(!(shore == "SH" & Coefficient %in% c("b_plotGPP_flow_mean", "b_laketempC_flow_mean", "b_lightmean_flow_mean")))%>%
  summarize(
    mean = mean(Estimate),
    l_95 = quantile(Estimate, 0.025),
    u_95 = quantile(Estimate, 0.975)
  )

# ggsave(plot = site_semplot_GPP, filename = paste("./NS_GPP_stream_semplot_21.png",sep=""),width=6.75,height=3.75,dpi=300)

site_semplot_GPP_23 <- ggplot(SEM_df_GPP_23, aes(x = Estimate, fill = shore, color = shore)) +
  geom_density(data = SEM_df_GPP_23 %>%
                 filter(!(shore == "SS" & Coefficient %in% c("b_plotGPP_flow_mean", "b_laketempC_flow_mean", "b_lightmean_flow_mean")))%>%
                 filter(!(shore == "SH" & Coefficient %in% c("b_plotGPP_flow_mean", "b_laketempC_flow_mean", "b_lightmean_flow_mean"))), 
               alpha = 0.75, linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "solid", size = 0.5, color = "black") +
  geom_vline(data = ci_summary_overall, aes(xintercept = mean, color = shore), linetype = "solid", size = 0.5) +
  scale_fill_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(name = "shore", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  labs(
    x = 'Estimate',
    y = 'Density') +
  facet_wrap(Coefficient ~ ., ncol = 8, scales = 'free')


library(ggpubr)

ns_grid <- ggarrange(
  site_semplot_GPP_21,
  site_semplot_GPP_22,
  site_semplot_GPP_23,
  common.legend = TRUE, 
  legend = "bottom",
  labels = c("a","b", "c"),
  align = c("v"),
  ncol = 1, nrow = 3,
  label.x = 0.01,  # Move label to the right
  label.y = 1,     # Keep label at the top
  hjust = 1,       # Right-align the label
  vjust = 1        # Top-align the label
)


# ggsave(plot = ns_grid, filename = paste("./NS_GPP_stream_semplot_remake.png",sep=""),width=12.75,height=5.75,dpi=300)





###================

###########################
#### revised ER model ###
##########################

ER_mod <- bf(plot_ER ~ plot_ER_lag + light_mean + lake_tempC + flow_mean + (1|site)) 
temp_mod <- bf(lake_tempC ~ flow_mean + ppt_mm + (1|site)) 
light_mod <- bf(light_mean ~ flow_mean + ppt_mm + (1|site)) 

########
BW_ER_df21 <- ER_df_BW_log %>%
  filter(date > as.Date("2021-06-10") & date < as.Date("2021-09-30"))


## BW
BW_ER_stream_fit <- brm(ER_mod +
                           temp_mod +
                           light_mod +
                           set_rescor(FALSE),
                         data=BW_ER_df21,
                         iter=10000, warmup = 5000,
                         control = list(adapt_delta = 0.95, max_treedepth = 15), 
                         cores=4, chains = 3)

BW_ER_stream_fit_sum <- summary(BW_ER_stream_fit)


# Extract posterior samples BW:
posterior_samples_BW <- posterior_samples(BW_ER_stream_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BW <- posterior_samples_BW %>%
  select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
ER_posterior_samples_melted_BW_21 <- melt(posterior_samples_BW, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(ER_posterior_samples_melted_BW_21, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/ER_posterior_samples_melted_BW_21.rds")

########
GB_ER_df21 <- ER_df_GB_log %>%
  filter(date > as.Date("2021-06-10") & date < as.Date("2021-09-30"))

## GB
GB_ER_stream_fit <- brm(ER_mod +
                           temp_mod +
                           light_mod +
                           set_rescor(FALSE),
                         data=GB_ER_df21,
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
  select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)


# Melt the data for easier plotting with ggplot2
ER_posterior_samples_melted_GB_21 <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(ER_posterior_samples_melted_GB_21, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/ER_posterior_samples_melted_GB_21.rds")


########
SS_ER_df21 <- ER_df_SS_log %>%
  filter(date > as.Date("2021-08-01") & date < as.Date("2021-09-30"))

## GB
SS_ER_stream_fit <- brm(ER_mod +
                           temp_mod +
                           light_mod +
                           set_rescor(FALSE),
                         data=SS_ER_df21,
                         iter=10000, warmup = 5000,
                         control = list(adapt_delta = 0.95, max_treedepth = 15), 
                         cores=4, chains = 3)

plot(SS_ER_stream_fit)
summary(SS_ER_stream_fit)
SS_ER_stream_fit_sum <- summary(SS_ER_stream_fit)


# Extract posterior samples SS:
posterior_samples_SS <- posterior_samples(SS_ER_stream_fit)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SS_21 <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")



BW_ER_df22 <- ER_df_BW_log %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

########
## BW
BW_ER_stream_fit_22 <- brm(ER_mod +
                              temp_mod +
                              light_mod +
                              set_rescor(FALSE),
                            data=BW_ER_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 15), 
                            cores=4, chains = 3)

plot(BW_ER_stream_fit_22)
BW_ER_stream_fit_sum <- summary(BW_ER_stream_fit_22)


# Extract posterior samples BW:
posterior_samples_BW <- posterior_samples(BW_ER_stream_fit_22)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BW <- posterior_samples_BW %>%
  select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
ER_posterior_samples_melted_BW_22 <- melt(posterior_samples_BW, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


# saveRDS(ER_posterior_samples_melted_BW_22, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/ER_posterior_samples_melted_BW_22.rds")


########
GB_ER_df22 <- ER_df_GB_log %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

## GB
GB_ER_stream_fit_22 <- brm(ER_mod +
                              temp_mod +
                              light_mod +
                              set_rescor(FALSE),
                            data=GB_ER_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 15), 
                            cores=4, chains = 3)

plot(GB_ER_stream_fit_22)
summary(GB_ER_stream_fit_22)
GB_ER_stream_fit_sum <- summary(GB_ER_stream_fit_22)


# Extract posterior samples GB:
posterior_samples_GB <- posterior_samples(GB_ER_stream_fit_22)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))
# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)


# Melt the data for easier plotting with ggplot2
ER_posterior_samples_melted_GB_22 <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


# saveRDS(ER_posterior_samples_melted_GB_22, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/ER_posterior_samples_melted_GB_22.rds")


########
ER_mod1 <- bf(plot_ER ~ plot_ER_lag + light_mean + lake_tempC) 
temp_mod1 <- bf(lake_tempC ~  ppt_mm) 
light_mod1 <- bf(light_mean ~ ppt_mm) 


ER_mod2 <- bf(plot_ER ~ plot_ER_lag + light_mean + lake_tempC + (1|site)) 
temp_mod2 <- bf(lake_tempC ~  ppt_mm + (1|site)) 
light_mod2 <- bf(light_mean ~ ppt_mm + (1|site)) 


SS_ER_df22 <- ER_df_SS_log %>%
  filter(date > as.Date("2021-09-30") & date < as.Date("2022-09-14"))

## SS
SS_ER_stream_fit_22 <- brm(ER_mod1 +
                              temp_mod1 +
                              light_mod1 +
                              set_rescor(FALSE),
                            data=SS_ER_df22,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 17), 
                            cores=3, chains = 3)

plot(SS_ER_stream_fit_22)
summary(SS_ER_stream_fit_22)
SS_ER_stream_fit_sum <- summary(SS_ER_stream_fit_22)


# Extract posterior samples SS:
posterior_samples_SS <- posterior_samples(SS_ER_stream_fit_22)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, 
         b_laketempC_ppt_mm,
         b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
ER_posterior_samples_melted_SS_22 <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")
# saveRDS(ER_posterior_samples_melted_SS_22, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/ER_posterior_samples_melted_SS_22.rds")


####### 2023


BW_ER_df23 <- ER_df_BW_log %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-14"))
summary(BW_ER_df23)


########
## BW
BW_ER_stream_fit_23 <- brm(ER_mod +
                              temp_mod +
                              light_mod +
                              set_rescor(FALSE),
                            data=BW_ER_df23,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 16), 
                            cores=4, chains = 3)

plot(BW_ER_stream_fit_23)
BW_ER_stream_fit_sum <- summary(BW_ER_stream_fit_23)


# Extract posterior samples BW:
posterior_samples_BW <- posterior_samples(BW_ER_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))
# Prepare the posterior samples data frame with shore levels
posterior_samples_BW <- posterior_samples_BW %>%
  select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
ER_posterior_samples_melted_BW_23 <- melt(posterior_samples_BW, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


# saveRDS(ER_posterior_samples_melted_BW_23, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/ER_posterior_samples_melted_BW_23.rds")


########
GB_ER_df23 <- ER_df_GB_log %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-14"))
summary(GB_ER_df23)

## GB
GB_ER_stream_fit_23 <- brm(ER_mod +
                              temp_mod +
                              light_mod +
                              set_rescor(FALSE),
                            data=GB_ER_df23,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 16), 
                            cores=4, chains = 3)

plot(GB_ER_stream_fit_23)
summary(GB_ER_stream_fit_23)
GB_ER_stream_fit_sum <- summary(GB_ER_stream_fit_23)


# Extract posterior samples GB:
posterior_samples_GB <- posterior_samples(GB_ER_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))
# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
ER_posterior_samples_melted_GB_23 <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(ER_posterior_samples_melted_GB_23, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/ER_posterior_samples_melted_GB_23.rds")


########
SS_ER_df23 <- ER_df_SS_log %>%
  filter(date > as.Date("2022-09-30") & date < as.Date("2023-09-14"))

## GB
SS_ER_stream_fit_23 <- brm(ER_mod2 +
                              temp_mod2 +
                              light_mod2 +
                              set_rescor(FALSE),
                            data=SS_ER_df23,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 15), 
                            cores=4, chains = 3)

plot(SS_ER_stream_fit_23)
summary(SS_ER_stream_fit_23)
SS_ER_stream_fit_sum <- summary(SS_ER_stream_fit_23)


# Extract posterior samples SS:
posterior_samples_SS <- posterior_samples(SS_ER_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC,
         b_laketempC_ppt_mm,
         , b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)


# Melt the data for easier plotting with ggplot2
ER_posterior_samples_melted_SS_23 <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# saveRDS(ER_posterior_samples_melted_SS_23, "/Users/kellyloria/Documents/Publications/CH2\ Streamflow\ and\ nearshore\ metabolism/SEM_output/ER_posterior_samples_melted_SS_23.rds")


########
SH_ER_df23 <- ER_df_SH_log %>%
  filter(date > as.Date("2023-01-01") & date < as.Date("2023-09-10"))

## GB
SH_ER_stream_fit_23 <- brm(ER_mod +
                              temp_mod +
                              light_mod +
                              set_rescor(FALSE),
                            data=SH_ER_df23,
                            iter=10000, warmup = 5000,
                            control = list(adapt_delta = 0.95, max_treedepth = 15), 
                            cores=4, chains = 3)

plot(SH_ER_stream_fit_23)
summary(SH_ER_stream_fit_23)
SH_ER_stream_fit_sum <- summary(SH_ER_stream_fit_23)


# Extract posterior samples SS:
posterior_samples_SH <- posterior_samples(SH_ER_stream_fit_23)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SH"), length.out = nrow(posterior_samples_SH))
# Prepare the posterior samples data frame with shore levels
posterior_samples_SHq <- posterior_samples_SH %>%
  select(b_plotER_plot_ER_lag, b_plotER_light_mean, b_plotER_lake_tempC, b_plotER_flow_mean,
         b_laketempC_flow_mean, b_laketempC_ppt_mm,
         b_lightmean_flow_mean, b_lightmean_ppt_mm) %>%
  mutate(shore = shore_levels)


# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SH_23 <- melt(posterior_samples_SHq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

##======================

















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