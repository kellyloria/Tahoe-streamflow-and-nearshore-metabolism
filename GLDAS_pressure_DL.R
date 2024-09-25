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
  http_string <- paste("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=GLDAS2:GLDAS_NOAH025_3H_v2.1:Psurf_f_inst")

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
  save_dir = "./NLDAS/stream/baro/",
  Site_ID = "BWL", # 39.1075414	-120.1646811
  Lat = "39.1075414", 
  Lon = "-120.1646811",
  startDate = "2020-01-01",
  endDate = "2024-01-01"
)


DL_GLDAS(
  save_dir = "./NLDAS/stream/baro/",
  Site_ID = "GBL", # 39.0880435	-119.9389446
  Lat = "39.0880435", 
  Lon = "-119.9389446",
  startDate = "2020-01-01",
  endDate = "2024-01-01"
)

# end of script.