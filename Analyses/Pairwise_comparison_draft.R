##==============================================================================
## Pairwise comparison of environmental conditions near and way from streams 
## 04/10/2025

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
dat <-readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/Final_Scripts/NS_analysis_dat.rds") %>%
  mutate(middle_ER_ab = c(middle_ER*-1))
str(dat)



chem_dat <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/Final_Scripts/NS_chem_dat.rds")
unique(chem_dat$site)
summary(chem_dat) 

OM_dat <- read.csv("/Users/kellyloria/Documents/UNR/MSMmetab/Sediment_OM_dat.csv") %>%
  mutate(date=as.Date(date, format("%m/%d/%y"))) %>%
  filter(Site=="GB_NS1" | Site=="GB_NS2" | Site=="GB_NS3" |
           Site=="BWNS1" |  Site=="BWNS2" | Site=="BWNS3" |
           Site=="SSNS1" |  Site=="SSNS2" | Site=="SSNS3" |
           Site=="SHNS1" |  Site=="SHNS2" | Site=="SHNS3") %>%
  dplyr::group_by(shore, date) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE))

  
  

pdat <-dat %>%
  dplyr::select(-X, -Latitude, -Longitude, elevation_m) %>%  # drop the 'date' column
  left_join(chem_dat, by=c("date", "shore", "site", "WaterYear")) %>%
  left_join(OM_dat, by=c("date", "shore"))


week_dat <- pdat %>%
  dplyr::select(-date, -yday) %>%  # drop the 'date' column
  dplyr::group_by(site, shore, week, WaterYear) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")


#### Supplemental figures:

## So for SH day range is 35 to 250 in 2023. 
## Looking at stream verse non-stream sites for that range 

# Filter the dataset to only include rows where do.obs_m or GPP_mean are not NA

df_BW_p <- pdat %>%
  filter(shore == "BW") %>%
  filter(year==2023) %>%
  filter(week>1 & week<41) %>%
  dplyr::group_by(shore, week, WaterYear) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

  

df_SS_p <- pdat %>%
  filter(shore == "SS") %>%
  filter(year==2023) %>%
  filter(week>1 & week<41)%>%
  dplyr::group_by(shore, week, WaterYear) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

df_west_subset <- rbind(df_BW_p, df_SS_p)%>%
  dplyr::select(shore, week, middle_GPP, middle_ER_ab, lake_tempC, lake_SPC, ppt_mm, 
                light_mean, windsp_mean, OM,
                NO3_mgL_dl, NH4_mgL_dl, PO4_ugL_dl, DOC_mgL_dl)


df_wide_w <- df_west_subset %>%
  pivot_wider(
    id_cols = week,
    names_from = shore,
    values_from = c(middle_GPP:PO4_ugL_dl, DOC_mgL_dl),  # use your actual variable names here
    names_glue = "{.value}_{shore}"
  )


df_diff_w <- df_wide_w %>%
  mutate(across(ends_with("_BW"), 
                ~ . - get(sub("_BW", "_SS", cur_column())),
                .names = "diff_{.col}"))


df_only_diffs <- df_diff_w %>%
  select(week, starts_with("diff_"))
  
df_long <- df_only_diffs %>%
  pivot_longer(
    cols = -week,
    names_to = "variable",
    values_to = "difference"
  )


ggplot(df_long, aes(x = week, y = difference)) +
  geom_point(color = "steelblue", shape=19) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  theme_minimal(base_size = 12) +
  labs(title = "Time Series of Differences (BW - SS)",
       x = "Date",
       y = "Difference") +
  theme(strip.text = element_text(size = 10, face = "bold"))



West_TS_plots <- ggplot(week_dat, aes(y = (middle_GPP), x = week, color = as.factor(shore), shape=site)) +
  geom_point(data = df_BW_p, aes(x = week, y = middle_GPP), color = "blue", size = 2) +
  geom_point(data = df_SS_p, aes(x = week, y = middle_GPP), color = "red", size = 2) +
  #scale_y_continuous(trans = "log", breaks = scales::log_breaks(base = 10)) +
  geom_smooth(data = df_BW_p, 
              se = TRUE, linewidth = 1, method = "gam", formula = y ~ s(x, k = 10), alpha = 0.1, color="blue") +
  geom_smooth(data = df_SS_p, 
              se = TRUE, linewidth = 1, method = "gam", formula = y ~ s(x, k = 10), alpha = 0.1, color="red") +
  theme_classic()


df_GB_p <- week_dat %>%
  filter(shore == "GB") %>%
  filter(year==2023) %>%
  filter(week>4 & week<41) 

df_SH_p <- week_dat %>%
  filter(shore == "SH") %>%
  filter(year==2023) %>%
  filter(week>4 & week<41) 


East_TS_plots <- ggplot(week_dat, aes(y = (middle_GPP), x = week, color = as.factor(shore), shape=site)) +
  geom_point(data = df_GB_p, aes(x = week, y = middle_GPP), color = "blue", size = 2) +
  geom_point(data = df_SH_p, aes(x = week, y = middle_GPP), color = "red", size = 2) +
  #scale_y_continuous(trans = "log", breaks = scales::log_breaks(base = 10)) +
  geom_smooth(data = df_GB_p, 
              se = TRUE, linewidth = 1, method = "gam", formula = y ~ s(x, k = 10), alpha = 0.1, color="blue") +
  geom_smooth(data = df_SH_p, 
              se = TRUE, linewidth = 1, method = "gam", formula = y ~ s(x, k = 10), alpha = 0.1, color="red") +
  theme_classic()


### by just shore

df_BW_p <- week_dat %>%
  filter(shore == "BW") %>%
  # filter(year==2023) %>%
  filter(week>0 & week<54) %>%
  dplyr::group_by(shore, week, WaterYear) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")


 df_SS_p <- week_dat %>%
  filter(shore == "SS") %>%
  # filter(year==2023) %>%
   filter(week>0 & week<54)  %>%
   dplyr::group_by(shore, week, WaterYear) %>%
   dplyr::summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
 


West_TS_plots <- ggplot(week_dat, aes(y = (middle_GPP), x = week, color = as.factor(shore), shape=shore)) +
  geom_point(data = df_BW_p, aes(x = week, y = middle_GPP), color = "#3283a8", size = 2) +
  geom_point(data = df_SS_p, aes(x = week, y = middle_GPP), color = "#136F63", size = 2) +
  #geom_line(data = df_SS_p, aes(x = week, y = middle_GPP),color = "red") +
  #scale_y_continuous(trans = "log", breaks = scales::log_breaks(base = 10)) +
  ylim(-5,35)+
  geom_smooth(data = df_BW_p, 
              se = TRUE, linewidth = 1, method = "gam", formula = y ~ s(x, k = 10), alpha = 0.1, color="#3283a8", fill="#3283a8") +
  geom_smooth(data = df_SS_p, 
              se = TRUE, linewidth = 1, method = "gam", formula = y ~ s(x, k = 10), alpha = 0.1, color="#136F63", fill= "#136F63") +
  theme_classic() + facet_grid(WaterYear~.)


df_GB_p <- week_dat %>%
  filter(shore == "GB") %>%
  #filter(year==2023) %>%
  filter(week>0 & week<54) %>%
  dplyr::group_by(shore, week, WaterYear) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")


df_SH_p <- week_dat %>%
  filter(shore == "SH") %>%
  filter(year==2023) %>%
  filter(week>0 & week<54) %>%
  dplyr::group_by(shore, week, WaterYear) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")



East_TS_plots <- ggplot(week_dat, aes(y = (middle_GPP), x = week, color = as.factor(shore), shape=shore)) +
  geom_point(data = df_GB_p, aes(x = week, y = middle_GPP), color = "#a67d17", size = 2) +
  geom_point(data = df_SH_p, aes(x = week, y = middle_GPP), color = "#c86640", size = 2) +
  #scale_y_continuous(trans = "log", breaks = scales::log_breaks(base = 10)) +
  ylim(-5,35)+
  geom_smooth(data = df_GB_p, 
              se = TRUE, linewidth = 1, method = "gam", formula = y ~ s(x, k = 10), alpha = 0.1, color="#a67d17", fill="#a67d17") +
  geom_smooth(data = df_SH_p, 
              se = TRUE, linewidth = 1, method = "gam", formula = y ~ s(x, k = 10), alpha = 0.1, color="#c76640", fill="#c76640") +
  theme_classic() + facet_grid(WaterYear~.)


df_SS_p <- week_dat %>%
  filter(shore == "SS") %>%
  filter(year==2023) %>%
  filter(yday>34 & yday<251) 

df_GB_p <- week_dat %>%
  filter(shore == "GB") %>%
  filter(year==2023) %>%
  filter(yday>34 & yday<251) 

df_SH_p <- week_dat %>%
  filter(shore == "SH") %>%
  filter(year==2023) %>%
  filter(yday>34 & yday<251) 

  

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
common_yday <- Reduce(intersect, list(#yday_2021, 
  yday_2022, 
  yday_2023))

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
templt_BW <- ggplot(plot_data_DO, aes(x = date)) +
  #geom_segment(aes(x = date_range[1], xend = date_range[2], y = 0, yend = 0), linewidth = 0.2) + # Main line
  geom_point(data = filter(filtered_data, has_do), aes(y = 0.095), alpha=0.8, color = "black", size = 1) + 
  geom_point(data = filter(filtered_data, has_GPP), aes(y = 0.090), alpha= 0.8, color = "#579467", size = 1) + 
  scale_y_continuous(name = "", breaks = NULL) +
  theme_classic() + xlab(NULL) +
  scale_x_date(date_breaks = "8 week",date_labels = "%b-%y", limits = c(as.Date("2021-09-29"), as.Date("2023-10-01"))) +
  labs(title = "BW NS1-3 metab.") #+facet_grid(year~.)





# Filter the dataset to only include rows where do.obs_m or GPP_mean are not NA
GPP_df_GB <- dat%>%
  filter(shore == "GB")

summary(GPP_df_GB)
# date range 2021-06-11  to 2023-09-14   


GB_GPP_df21 <- GPP_df_GB %>%
  filter(date > as.Date("2021-06-10") & date < as.Date("2021-09-10"))

summary(GB_GPP_df21)
# date range 2021-06-11  to 2023-09-14   

plot_data_DO <- GPP_df_GB %>%
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
common_yday <- Reduce(intersect, list(#yday_2021, 
  yday_2022, 
  yday_2023))

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
templt_GB <- ggplot(plot_data_DO, aes(x = date)) +
  #geom_segment(aes(x = date_range[1], xend = date_range[2], y = 0, yend = 0), linewidth = 0.2) + # Main line
  geom_point(data = filter(filtered_data, has_do), aes(y = 0.095), alpha=0.8, color = "black", size = 1) + # DO presence tick marks
  geom_point(data = filter(filtered_data, has_GPP), aes(y = 0.090), alpha= 0.8, color = "#579467", size = 1) + # GPP presence tick marks
  scale_y_continuous(name = "", breaks = NULL) +
  scale_x_date(date_breaks = "8 week",date_labels = "%b-%y", limits = c(as.Date("2021-09-29"), as.Date("2023-10-01"))) +
  theme_classic() + xlab(NULL) +
  labs(title = "GB NS1-3 metab.") #+facet_grid(year~.)




# Filter the dataset to only include rows where do.obs_m or GPP_mean are not NA
GPP_df_SS <- dat%>%
  filter(shore == "SS")

summary(GPP_df_SS)
# date range 2021-06-11  to 2023-09-14   


SS_GPP_df21 <- GPP_df_SS %>%
  filter(date > as.Date("2021-06-10") & date < as.Date("2021-09-10"))

summary(SS_GPP_df21)
# date range 2021-06-11  to 2023-09-14   

plot_data_DO <- GPP_df_SS %>%
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
common_yday <- Reduce(intersect, list(#yday_2021, 
  yday_2022, 
  yday_2023))

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
templt_SS <- ggplot(plot_data_DO, aes(x = date)) +
  #geom_segment(aes(x = date_range[1], xend = date_range[2], y = 0, yend = 0), linewidth = 0.2) + # Main line
  geom_point(data = filter(filtered_data, has_do), aes(y = 0.095), alpha=0.8, color = "black", size = 1) + # DO presence tick marks
  geom_point(data = filter(filtered_data, has_GPP), aes(y = 0.090), alpha= 0.8, color = "#579467", size = 1) + # GPP presence tick marks
  scale_y_continuous(name = "", breaks = NULL) +
  theme_classic() + xlab(NULL) +
  scale_x_date(
    limits = c(as.Date("2021-09-29"), as.Date("2023-10-01")),
    date_breaks = "8 week", date_labels = "%b-%y") +
  labs(title = "SS NS1 & 2 metab.") #+facet_grid(year~.)





# Filter the dataset to only include rows where do.obs_m or GPP_mean are not NA
GPP_df_SH <- dat%>%
  filter(shore == "SH")

summary(GPP_df_SH)
# date range 2021-06-11  to 2023-09-14   


SH_GPP_df21 <- GPP_df_SH %>%
  filter(date > as.Date("2021-06-10") & date < as.Date("2021-09-10"))

summary(SH_GPP_df21)
# date range 2021-06-11  to 2023-09-14   

plot_data_DO <- GPP_df_SH %>%
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

# Find the common yday values acroSH all three years
common_yday <- Reduce(intersect, list(#yday_2021, 
  yday_2022, 
  yday_2023))

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
templt_SH <- ggplot(plot_data_DO, aes(x = date)) +
  #geom_segment(aes(x = date_range[1], xend = date_range[2], y = 0, yend = 0), linewidth = 0.2) + # Main line
  geom_point(data = filter(plot_data_DO, has_do), aes(y = 0.095), alpha=0.8, color = "black", size = 1) + # DO presence tick marks
  geom_point(data = filter(plot_data_DO, has_GPP), aes(y = 0.090), alpha= 0.8, color = "#579467", size = 1) + # GPP presence tick marks
  scale_y_continuous(name = "", breaks = NULL) +
  theme_classic() + xlab(NULL) +
  scale_x_date(
    limits = c(as.Date("2021-09-29"), as.Date("2023-10-01")),
    date_breaks = "8 week", date_labels = "%b-%y") +
  labs(title = "SH NS1-3 metab.") #+facet_grid(year~.)
