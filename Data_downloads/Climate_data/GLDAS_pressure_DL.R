#' Downloads GLDAS pressure data
#' @description This function downloads GLDAS pressure for a given Latitude and Longitude
#' @param save_dir The save directory for files to be placed in. For example, "C:/myfolder
#' @param Site The site name, for example "FL_ICHE2700"
#' @param Lat The site Latitude
#' @param Lon The site Longitude
#' Selecting only data from the central sensor in network (NS2).

#' @param startDate The starting date for the download (YYYY-MM-DD)
#'
#' @return Returns a time series of barometric pressure from the start
#' date to the most recent available data
#' @export

#===============================================================================
#Function for downloading GLDAS pressure "Psurf_f_inst" data from 2000 - 2022 via data rods
#https://disc.gsfc.nasa.gov/information/tools?title=Hydrology%20Data%20Rods
#Created 09/24/24
#===============================================================================

DL_GLDAS <- function(save_dir, Site_ID, Lat, Lon, startDate, endDate){
  #The initial string to build the URL
  http_string <- paste("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=GLDAS2:GLDAS_NOAH025_3H_v2.0:Psurf_f_inst")
  # https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=GLDAS2:GLDAS_NOAH025_3H_v2.0:Psurf_f_inst&startDate=2011-01-01T00:00&endDate=2012-12-31T21:00&location=GEOM:POINT(-93.125,%2041.125)&type=asc2
  
    #  http_string <- paste("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:DSWRFsfc")

  #Separating the date information
  start_split <- strsplit(startDate, "-")[[1]]
  end_split <- strsplit(endDate, "-")[[1]]

  #Build individual components of the url
  location_string <- paste0("&location=GEOM:POINT(", Lon, ",%20", Lat, ")")
  # start_string <- paste0("&startDate=", start_split[1], "-", start_split[2], "-",
  #                        start_split[3], "T00")
  # end_string <- paste0("&endDate=", end_split[1], "-", end_split[2], "-",
  #                      end_split[3], "T00")
  
  start_string <- paste0("&startDate=", start_split[1], "-", start_split[2], "-", start_split[3], "T00:00:00")
  end_string <- paste0("&endDate=", end_split[1], "-", end_split[2], "-", end_split[3], "T23:59:59")
  
  #Generating the URL
  url <-paste0(http_string, location_string, start_string, end_string, "&type=asc2")

  #Downloading the data
  destfile <- paste(save_dir,"/", Site_ID, "_GLDAS.asc", sep = "")


  #Error catch in case the page is inaccessible. A little inelegant at present...
  try_result <- try(download.file(url, destfile), silent = FALSE)

  if(class(try_result) == "try-error") {file.remove(destfile)}

} #End DL_NLDAS function

#############################
####### READ IN DATA ########
#############################

## NEARSHORE

# KAL path: ~/Documents/LittoralMetabModeling/RawData/NLDAS/GLDAS_baro
DL_GLDAS( "./NLDAS/GLDAS_baro",
  Site_ID = "BWNS2",
  Lat = "39.113185",
  Lon = "-120.062018",
  startDate = "2021-01-01",
  endDate = "2023-11-01"
)

DL_GLDAS( "./NLDAS/GLDAS_baro",
  Site_ID = "SSNS2",
  Lat = "39.142064",
  Lon = "-120.058892",
  startDate = "2021-01-01",
  endDate = "2023-11-01"
)

DL_GLDAS( "./GLDAS_baro",
  Site_ID = "SHNS2",
  Lat = "39.103443",
  Lon = "-120.035747",
  startDate = "2021-01-01",
  endDate = "2023-11-01"
)


DL_GLDAS( "./NLDAS/GLDAS_baro",
  Site_ID = "GBNS2",
  Lat = "39.088535",
  Lon = "-120.033959",
  startDate = "2021-01-01",
  endDate = "2023-11-01"
)

###############
## OFFSHORE ##

DL_GLDAS(
  save_dir = "./NLDAS/Offshore/baro/",
  Site_ID = "GB10m",
  Lat = "39.08836",
  Lon = "-119.9474",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)


DL_GLDAS(
  save_dir = "./NLDAS/Offshore/baro/",
  Site_ID = "GB15m", 
  Lat = "39.088485",
  Lon = "-119.948979",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)

DL_GLDAS(
  save_dir = "./NLDAS/Offshore/baro/",
  Site_ID = "GB20m", 
  Lat = "39.088271",
  Lon = "-119.950734",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)

DL_GLDAS(
  save_dir = "./NLDAS/Offshore/baro/",
  Site_ID = "BW10m", 
  Lat = "39.10629",
  Lon = "-120.15701",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)

DL_GLDAS(
  save_dir = "./NLDAS/Offshore/baro/",
  Site_ID = "BW20m", 
  Lat = "39.106182", 
  Lon = "-120.156302",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)


###

DL_GLDAS(
  save_dir = "/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NLDAS/stream/baro/",
  Site_ID = "BWL", # 39.1075414	-120.1646811
  Lat = "39.107541", 
  Lon = "-120.164681",
  startDate = "2020-09-20",
  endDate = "2023-08-01"
)

DL_GLDAS(
  save_dir = "/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NLDAS/stream/baro/",
  Site_ID = "BWU", # 39.105291, -120.195904
  Lat = "39.105291", 
  Lon = "-120.195904",
  startDate = "2021-06-01",
  endDate = "2024-09-01"
)

DL_GLDAS(
  save_dir = "/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NLDAS/stream/baro/",
  Site_ID = "GBLv2", # 39.0880435	-119.9389446
  Lat = "39.08804", 
  Lon = "-119.93895",
  startDate = "2021-03-20",
  endDate = "2023-08-01"
)

DL_GLDAS(
  save_dir = "/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NLDAS/stream/baro/",
  Site_ID = "GBU", # 39.086730, -119.931449
  Lat = "39.086730", 
  Lon = "-119.931449",
  startDate = "2021-03-20",
  endDate = "2023-10-10"
)
# end of script.

DL_GLDAS <- function(save_dir, Site_ID, Lat, Lon, startDate, endDate) {
  # Base URL for the API
  http_string <- "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi"
  
  # Full variable path (with prefix if needed)
  variable_param <- "GLDAS2:GLDAS_NOAH025_3H_v2.0:Psurf_f_inst"
#   GLDAS_NOAH025_3H_EP_v2.1
  
  # Update location format to use a comma
  location_string <- paste0("&location=GEOM:POINT(", Lon, ",", Lat, ")")
  
  # Format start and end date strings
  start_string <- paste0("&startDate=", startDate, "T00:00:00")
  end_string <- paste0("&endDate=", endDate, "T23:59:59")
  
  # Construct full URL
  url <- paste0(http_string, "?variable=", variable_param, location_string, start_string, end_string, "&type=asc2")
  
  # Define destination file path
  destfile <- file.path(save_dir, paste0(Site_ID, "_GLDAS.asc"))
  
  # Attempt download
  try_result <- try(download.file(url, destfile, quiet = FALSE), silent = TRUE)
  
  # Handle errors if download fails
  if (class(try_result) == "try-error") {
    file.remove(destfile)
    cat("Download failed: Check URL and parameters.\n")
  } else {
    cat("Download succeeded: File saved to", destfile, "\n")
  }
}





# Load httr library for handling authentication
library(httr)

DL_GLDAS <- function(save_dir, Site_ID, Lat, Lon, startDate, endDate, username, password) {
  # Base URL for the API
  http_string <- "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi"
  
  # Define variable and location
  variable_param <- "GLDAS2:GLDAS_NOAH025_3H_v2.0:Psurf_f_inst"
  location_string <- paste0("&location=GEOM:POINT(", Lon, ",", Lat, ")")
  
  # Define start and end dates
  start_string <- paste0("&startDate=", startDate, "T00:00:00")
  end_string <- paste0("&endDate=", endDate, "T23:59:59")
  
  # Construct the full URL
  url <- paste0(http_string, "?variable=", variable_param, location_string, start_string, end_string, "&type=asc2")
  
  # Define the file path to save
  destfile <- file.path(save_dir, paste0(Site_ID, "_GLDAS.asc"))
  
  # Perform the download with authentication
  try_result <- try(
    GET(url, authenticate(username, password), write_disk(destfile, overwrite = TRUE)),
    silent = TRUE
  )
  
  # Error handling
  if (class(try_result) == "try-error") {
    cat("Download failed: Please check URL, parameters, and authentication.\n")
  } else if (status_code(try_result) != 200) {
    file.remove(destfile)
    cat("Download failed: HTTP error", status_code(try_result), "\n")
  } else {
    cat("Download succeeded: File saved to", destfile, "\n")
  }
}


DL_GLDAS(
  save_dir = "/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NLDAS/stream/baro/",
  Site_ID = "GBU",
  Lat = "39.086730",
  Lon = "-119.931449",
  startDate = "2021-03-20",
  endDate = "2024-10-10",
  username = "kelly.loria",
  password = "TahoePines2022/"
)




DL_GLDAS(
  save_dir = "/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NLDAS/stream/baro/",
  Site_ID = "GBU",
  Lat = "39.086730",
  Lon = "-119.931449",
  startDate = "2021-03-20",
  endDate = "2021-10-10"
)



# Load httr for HTTP requests and authentication
library(httr)

# Function to download NLDAS data for PSurf
download_nldas_psurf <- function(save_dir, Site_ID, Lat, Lon, startDate, endDate, username, password) {
  # Base URL for the API
  base_url <- "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi"
  
  # NLDAS variable for surface pressure (PSurf)
  variable_param <- "NLDAS_FORA0125_H_2.0:PSurf"
  
  # Construct the location in the GEOM:POINT format
  location_string <- paste0("&location=GEOM:POINT(", Lon, "%20", Lat, ")")
  
  # Define the start and end date (format YYYY-MM-DD)
  start_string <- paste0("&startDate=", startDate, "T00:00:00")
  end_string <- paste0("&endDate=", endDate, "T23:59:59")
  
  # Construct the full request URL
  url <- paste0(base_url, "?variable=", variable_param, location_string, start_string, end_string, "&type=asc2")
  
  # Define destination file path
  destfile <- file.path(save_dir, paste0(Site_ID, "_NLDAS_PSurf.asc"))
  
  
  # Check for errors
  if (class(response) == "try-error") {
    cat("Download failed: Please check the URL, parameters, or authentication.\n")
  } else if (status_code(response) != 200) {
    file.remove(destfile)
    cat("Download failed: HTTP error", status_code(response), "\n")
  } else {
    cat("Download succeeded: File saved to", destfile, "\n")
  }
}

# Example of how to call the function
download_nldas_psurf(
  save_dir = "/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NLDAS/stream/baro/",   # Replace with your directory path
  Site_ID = "GBU",                        # Site ID, e.g., "GBU"
  Lat = "39.086730",                      # Latitude
  Lon = "-119.931449",                    # Longitude
  startDate = "2021-03-20",               # Start date
  endDate = "2021-10-10",                 # End date
)
