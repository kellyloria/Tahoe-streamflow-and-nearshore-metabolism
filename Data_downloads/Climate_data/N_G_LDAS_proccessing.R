#' Processes downloaded NLDAS data
#' @description This function processes downloaded NLDAS incoming shortwave radiation
#' data (w m-2) for a given Latitude and Longitude.
#'
#' @param read_dir The read directory for downloaded files. For example, "C:/myfolder
#' @param Site_IDs Site name(s), for example "BWNS1"
#' @param write_output Logical value indicating whether to write each individual driver
#' file to disk. Default value is FALSE.
#' @param save_dir Optional parameter when write_output = TRUE. The save directory
#' for files to be placed in. For example, "C:/
#'
#' @return Returns a time series of incoming light data
#' @export

#===============================================================================
#Function for processing the downloaded data modified from
# https://rdrr.io/github/psavoy/StreamLightUtils/src/R/NLDAS_proc.R
# updated: 09/18/24 
#===============================================================================

#===============================================================================
# LIGHT
# NLDAS:NLDAS_FORA0125_H.002:DSWRFsfc
#===============================================================================
# *Check file path from .asc files in the _DL scripts

NLDAS_proc <- function(read_dir, Site_IDs, write_output = FALSE, save_dir = NULL) {
  # Get a list of all downloaded NLDAS data
  setwd(read_dir)
  downloaded <- list.files(read_dir)[grep("*_NLDAS.asc", list.files(read_dir))]

  # Get the names of downloaded sites
  downloaded_names <- stringr::str_sub(downloaded, 1, -11)

  # Function for processing each site
  NLDAS_site <- function(file_name, write_output, save_dir) {
    # Extracting the site information from the file name
    site_match <- stringr::str_match(file_name, "([A-Za-z0-9_]+)_NLDAS.asc")

    # Checking if there is a match
    if (!is.null(site_match) && length(site_match) == 2) {
      site <- site_match[2]
    } else {
      site <- "UnknownSite"
    }

    # Reading in the table, skipping the first 40 lines of header information
    # and removing the last row which contains a calculated mean value
    nldas <- read.table(file_name, skip = 40, nrows = length(readLines(file_name, warn = FALSE)) - 41)

    colnames(nldas) <- c("Date", "hour_raw", "light")

    # Adding in date and time information
    # Extracting the hour information
    nldas[, "Time"] <- as.numeric(substr(nldas[, "hour_raw"], 1, 2))

    # Adding a POSIX time column
    nldas[, "pos_time"] <- as.POSIXct(paste(nldas[, "Date"], " ", as.matrix(sprintf("%02d", nldas[, "Time"])), sep = ""),
                                      format = "%Y-%m-%d %H", tz = "UTC")

    # Adding in Year, DOY, and hour information
    nldas[, "Year"] <- as.numeric(format(nldas[, "pos_time"], format = "%Y", tz = "UTC"))
    nldas[, "DOY"] <- as.numeric(format(nldas[, "pos_time"], format = "%j", tz = "UTC"))
    nldas[, "Hour"] <- as.numeric(format(nldas[, "pos_time"], format = "%H", tz = "UTC"))

    # Selecting the final columns
    final <- nldas[, c("Year", "DOY", "Hour", "pos_time", "light")]
    colnames(final)[4] <- "datetime"

    # Adding the site information
    final$site <- site

    # If write_output == TRUE, save the driver to disk
    if (write_output == TRUE) {
      saveRDS(final, paste0(save_dir, "/", site, "_NLDAS_processed.rds"))
    } else {
      return(final)
    }
  }

  # Apply the function to make all driver files
  if (write_output == TRUE) {
    lapply(downloaded, FUN = NLDAS_site, write_output = write_output, save_dir = save_dir)
  } else {
    processed <- lapply(downloaded, FUN = NLDAS_site, write_output = write_output, save_dir = save_dir)
    names(processed) <- downloaded_names

    return(processed)
  }

  # Notify the user with a list of sites that did not have data
  missing <- Site_IDs[!(Site_IDs %in% downloaded_names)]

  if (length(missing) != 0) {
    print(paste("The following sites did not successfully download NLDAS data:",
                paste(missing, sep = "", collapse = ", ")))
  }
}

# Call the function with your parameters


dat1 <- NLDAS_proc("~./NLDAS/NLDAS_light/",
                   "~./NLDAS/processed_light/")

dat2<- do.call(rbind, dat1)
str(dat2)

###
# quick vizualization
library(ggplot2)
library(tidyverse)

plot_dat <- dat2 %>%
  filter(site=="BWNS2")

plot_dats <- dat2 %>%
  filter(datetime>as.POSIXct("2023-01-01 00:00:00")&datetime<as.POSIXct("2023-01-07 00:00:00"))

plot_light <-
  ggplot(dat2, aes(x=datetime, y=light, colour = as.factor(site))) +
  geom_point(alpha=0.1) + geom_line(alpha=0.25)+
  labs(y = 'light dat', x=NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12)) #+ theme(legend.position = "none") + facet_grid(Name~.)
plot_light

# write.csv(x = dat2, file = "./NLDAS/processed_light/nearshore_NLDAS_light.csv", row.names = TRUE)


############
# OFF SHORE
############
# Call the function with your parameters
dat1 <- NLDAS_proc("~./NLDAS/Offshore/light/",
                   "~./NLDAS/processed_light/")

dat2<- do.call(rbind, dat1)
str(dat2)

# write.csv(x = dat2, file = "./NLDAS/processed_light/offshore_NLDAS_light.csv", row.names = TRUE)

###########
# stream 
###########
dat1 <- NLDAS_proc("~/Documents/LittoralMetabModeling/RawData/NLDAS/stream/light/",
                   "~/Documents/LittoralMetabModeling/RawData/stream/")

dat2<- do.call(rbind, dat1)
str(dat2)
# write.csv(x = dat2, file = "./NLDAS/stream/stream_NLDAS_light.csv", row.names = TRUE)


#===============================================================================
# BARO
# GLDAS2:GLDAS_NOAH025_3H_v2.1:Psurf_f_inst
#===============================================================================

GLDAS_proc <- function(read_dir, Site_IDs, write_output = FALSE, save_dir = NULL) {
  # Get a list of all downloaded NLDAS data
  setwd(read_dir)
  downloaded <- list.files(read_dir)[grep("*_GLDAS.asc", list.files(read_dir))]

  # Get the names of downloaded sites
  downloaded_names <- stringr::str_sub(downloaded, 1, -11)

  # Function for processing each site
  GLDAS_site <- function(file_name, write_output, save_dir) {
    # Extracting the site information from the file name
    site_match <- stringr::str_match(file_name, "([A-Za-z0-9_]+)_GLDAS.asc")

    # Checking if there is a match
    if (!is.null(site_match) && length(site_match) == 2) {
      site <- site_match[2]
    } else {
      site <- "UnknownSite"
    }

    # Reading in the table, skipping the first 40 lines of header information
    # and removing the last row which contains a calculated mean value
    gldas <- read.table(file_name, skip = 13, nrows = length(readLines(file_name, warn = FALSE)) - 41)

    colnames(gldas) <- c("pos_time", "baro_Pa")

    # Adding in date and time information
    # Extracting the hour information
    # Adding a POSIX time column
    gldas[, "pos_time"] <- as.POSIXct(gldas[, "pos_time"], format = "%Y-%m-%dT%H:%M:%S")


    # Adding in Year, DOY, and hour information
    gldas[, "Year"] <- as.numeric(format(gldas[, "pos_time"], format = "%Y", tz = "UTC"))
    gldas[, "DOY"] <- as.numeric(format(gldas[, "pos_time"], format = "%j", tz = "UTC"))
    gldas[, "Hour"] <- as.numeric(format(gldas[, "pos_time"], format = "%H", tz = "UTC"))

    # Selecting the final columns
    final <- gldas[, c("Year", "DOY", "Hour", "pos_time", "baro_Pa")]
    colnames(final)[4] <- "datetime"

    # Adding the site information
    final$site <- site

    # If write_output == TRUE, save the driver to disk
    if (write_output == TRUE) {
      saveRDS(final, paste0(save_dir, "/", site, "_NLDAS_processed.rds"))
    } else {
      return(final)
    }
  }

  # Apply the function to make all driver files
  if (write_output == TRUE) {
    lapply(downloaded, FUN = GLDAS_site, write_output = write_output, save_dir = save_dir)
  } else {
    processed <- lapply(downloaded, FUN = GLDAS_site, write_output = write_output, save_dir = save_dir)
    names(processed) <- downloaded_names

    return(processed)
  }

  # Notify the user with a list of sites that did not have data
  missing <- Site_IDs[!(Site_IDs %in% downloaded_names)]

  if (length(missing) != 0) {
    print(paste("The following sites did not successfully download NLDAS data:",
                paste(missing, sep = "", collapse = ", ")))
  }
}

# Call the function with your parameters
gdat <- GLDAS_proc("~/Documents/LittoralMetabModeling/RawData/NLDAS/GLDAS_baro/",
                   "~/Documents/LittoralMetabModeling/RawData/NLDAS/processed_baro/")

baro<- do.call(rbind, gdat)
str(baro)



plot_dat <- baro %>%
  filter(site=="SHNS2"| site=="SSNS2")


plot_dats <- plot_dat %>%
  filter(datetime>as.POSIXct("2023-01-01 00:00:00")&datetime<as.POSIXct("2023-01-07 00:00:00"))

plot_baro <-
  ggplot(baro, aes(x=datetime, y=baro_Pa, colour = as.factor(site))) +
  geom_point(alpha=0.1) + geom_line(alpha=0.25)+
  #labs(y = 'light dat', x=NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12)) #+ theme(legend.position = "none") + facet_grid(Name~.)
plot_baro

# write.csv(x = baro, file = "./NLDAS/processed_baro/nearshore_NLDAS_baro.csv", row.names = TRUE)


############
# OFF SHORE
############
# Call the function with your parameters
dat <- GLDAS_proc("~/Documents/LittoralMetabModeling/RawData/NLDAS/Offshore/baro/",
                   "~/Documents/LittoralMetabModeling/RawData/NLDAS/processed_baro/")

dat2<- do.call(rbind, dat)
str(dat2)

# write.csv(x = dat2, file = "./NLDAS/processed_baro/offshore_NLDAS_baro.csv", row.names = TRUE)


############
# stream
############
# Call the function with your parameters
dat <- GLDAS_proc("~/Documents/LittoralMetabModeling/RawData/NLDAS/stream/baro/",
                  "~/Documents/LittoralMetabModeling/RawData/NLDAS/stream/")

dat2<- do.call(rbind, dat)
str(dat2)

# write.csv(x = dat2, file = "./NLDAS/stream/stream_NLDAS_baro.csv", row.names = TRUE)



#===============================================================================
# WIND SP
# GLDAS2:GLDAS_NOAH025_3H_v2.1:Wind_f_inst
#===============================================================================


GLDAS_proc <- function(read_dir, Site_IDs, write_output = FALSE, save_dir = NULL) {
  # Get a list of all downloaded NLDAS data
  setwd(read_dir)
  downloaded <- list.files(read_dir)[grep("*_GLDAS.asc", list.files(read_dir))]

  # Get the names of downloaded sites
  downloaded_names <- stringr::str_sub(downloaded, 1, -11)

  # Function for processing each site
  GLDAS_site <- function(file_name, write_output, save_dir) {
    # Extracting the site information from the file name
    site_match <- stringr::str_match(file_name, "([A-Za-z0-9_]+)_GLDAS.asc")

    # Checking if there is a match
    if (!is.null(site_match) && length(site_match) == 2) {
      site <- site_match[2]
    } else {
      site <- "UnknownSite"
    }

    # Reading in the table, skipping the first 40 lines of header information
    # and removing the last row which contains a calculated mean value
    gldas <- read.table(file_name, skip = 13, nrows = length(readLines(file_name, warn = FALSE)) - 41)

    colnames(gldas) <- c("pos_time", "windsp_ms")

    # Adding in date and time information
    # Extracting the hour information
    # Adding a POSIX time column
    gldas[, "pos_time"] <- as.POSIXct(gldas[, "pos_time"], format = "%Y-%m-%dT%H:%M:%S")


    # Adding in Year, DOY, and hour information
    gldas[, "Year"] <- as.numeric(format(gldas[, "pos_time"], format = "%Y", tz = "UTC"))
    gldas[, "DOY"] <- as.numeric(format(gldas[, "pos_time"], format = "%j", tz = "UTC"))
    gldas[, "Hour"] <- as.numeric(format(gldas[, "pos_time"], format = "%H", tz = "UTC"))

    # Selecting the final columns
    final <- gldas[, c("Year", "DOY", "Hour", "pos_time", "windsp_ms")]
    colnames(final)[4] <- "datetime"

    # Adding the site information
    final$site <- site

    # If write_output == TRUE, save the driver to disk
    if (write_output == TRUE) {
      saveRDS(final, paste0(save_dir, "/", site, "_NLDAS_processed.rds"))
    } else {
      return(final)
    }
  }

  # Apply the function to make all driver files
  if (write_output == TRUE) {
    lapply(downloaded, FUN = GLDAS_site, write_output = write_output, save_dir = save_dir)
  } else {
    processed <- lapply(downloaded, FUN = GLDAS_site, write_output = write_output, save_dir = save_dir)
    names(processed) <- downloaded_names

    return(processed)
  }

  # Notify the user with a list of sites that did not have data
  missing <- Site_IDs[!(Site_IDs %in% downloaded_names)]

  if (length(missing) != 0) {
    print(paste("The following sites did not successfully download NLDAS data:",
                paste(missing, sep = "", collapse = ", ")))
  }
}

# Call the function with your parameters
gdat1 <- GLDAS_proc("~./NLDAS/GLDAS_windsp/",
                   "~./NLDAS/processed_windsp/")

windsp<- do.call(rbind, gdat1)
str(windsp)



plot_dat <- windsp %>%
  filter(site=="SHNS2"| site=="SSNS2")

plot_dats <- plot_dat %>%
  filter(datetime>as.POSIXct("2023-01-01 00:00:00")&datetime<as.POSIXct("2023-01-07 00:00:00"))

plot_wsp <-
  ggplot(windsp, aes(x=datetime, y=windsp_ms, colour = as.factor(site))) +
  geom_point(alpha=0.1) + geom_line(alpha=0.25)+
  #labs(y = 'light dat', x=NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12)) #+ theme(legend.position = "none") + facet_grid(Name~.)
plot_wsp

# write.csv(x = windsp, file = "./NLDAS/processed_windsp/nearshore_NLDAS_windsp.csv", row.names = TRUE)

############
# OFF SHORE
############
# Call the function with your parameters
# Call the function with your parameters
wind <- GLDAS_proc("~/Documents/LittoralMetabModeling/RawData/NLDAS/Offshore/windsp/",
                  "~/Documents/LittoralMetabModeling/RawData/NLDAS/processed_windsp/")

wind2<- do.call(rbind, wind)
str(wind2)

# write.csv(x = wind2, file = "./NLDAS/processed_windsp/offshore_NLDAS_windsp.csv", row.names = TRUE)

