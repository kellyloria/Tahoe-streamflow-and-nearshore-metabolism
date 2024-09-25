#' Downloads GLDAS wind speed data: Wind_f_inst
#' @description This function downloads GLDAS pressure for a given Latitude and Longitude
#' @param save_dir The save directory for files to be placed in. For example, "C:/myfolder
#' @param Site The site name, for example "FL_ICHE2700"
#' @param Lat The site Latitude
#' @param Lon The site Longitude
#' @param startDate The starting date for the download (YYYY-MM-DD)
#' only for lake array.
#'
#' @return Returns a time series of barometric pressure from the start
#' date to the most recent available data
#' @export

#===============================================================================
#Function for downloading GLDAS "Wind_f_inst" data from 2000 - 2022 via data rods
#https://disc.gsfc.nasa.gov/information/tools?title=Hydrology%20Data%20Rods
#Created 12 December 2022
#===============================================================================

DL_GLDAS <- function(save_dir, Site_ID, Lat, Lon, startDate, endDate){
  #The initial string to build the URL
  http_string <- paste("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=GLDAS2:GLDAS_NOAH025_3H_v2.1:Wind_f_inst")

  #Separating the date information
  start_split <- strsplit(startDate, "-")[[1]]
  end_split <- strsplit(endDate, "-")[[1]]

  #Build individual components of the url
  location_string <- paste0("&location=GEOM:POINT(", Lon, ",%20", Lat, ")")
  start_string <- paste0("&startDate=", start_split[1], "-", start_split[2], "-",
                         start_split[3], "T00")
  end_string <- paste0("&endDate=", end_split[1], "-", end_split[2], "-",
                       end_split[3], "T00")
  #Generating the URL
  url <-paste0(http_string, location_string, start_string, end_string, "&type=asc2")

  #Downloading the data
  destfile <- paste(save_dir,"/", Site_ID, "_GLDAS.asc", sep = "")


  #Error catch in case the page is inaccessible. A little inelegant at present...
  try_result <- try(download.file(url, destfile), silent = FALSE)

  if(class(try_result) == "try-error") {file.remove(destfile)}

} #End DL_NLDAS function


DL_GLDAS( "./RawData/NLDAS/GLDAS_windsp",
          Site_ID = "BWNS2",
          Lat = "39.10697",
          Lon = "-120.15721",
          startDate = "2021-01-01",
          endDate = "2023-11-01"
)

DL_GLDAS( "./RawData/NLDAS/GLDAS_windsp",
          Site_ID = "SSNS2",
          Lat = "39.13904",
          Lon = "-120.15236",
          startDate = "2021-01-01",
          endDate = "2023-11-01"
)

DL_GLDAS( "./RawData/NLDAS/GLDAS_windsp",
          Site_ID = "SHNS2",
          Lat = "39.0944",
          Lon = "-119.9436",
          startDate = "2021-01-01",
          endDate = "2023-11-01"
)


DL_GLDAS( "./RawData/NLDAS/GLDAS_windsp",
          Site_ID = "GBNS2",
          Lat = "39.0880",
          Lon = "-119.9421",
          startDate = "2021-01-01",
          endDate = "2023-11-01"
)

###############
## OFFSHORE ##

DL_GLDAS(
  save_dir = "./RawData/NLDAS/Offshore/windsp/",
  Site_ID = "GB10m",
  Lat = "39.08836",
  Lon = "-119.9474",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)


DL_GLDAS(
  save_dir = "./RawData/NLDAS/Offshore/windsp/",
  Site_ID = "GB15m", 
  Lat = "39.088485",
  Lon = "-119.948979",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)

DL_GLDAS(
  save_dir = "./RawData/NLDAS/Offshore/windsp/",
  Site_ID = "GB20m", 
  Lat = "39.088271",
  Lon = "-119.950734",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)

DL_GLDAS(
  save_dir = "./RawData/NLDAS/Offshore/windsp/",
  Site_ID = "BW10m", 
  Lat = "39.10629",
  Lon = "-120.15701",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)

DL_GLDAS(
  save_dir = "./RawData/NLDAS/Offshore/windsp/",
  Site_ID = "BW20m", 
  Lat = "39.106182", 
  Lon = "-120.156302",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)

# end of script.
