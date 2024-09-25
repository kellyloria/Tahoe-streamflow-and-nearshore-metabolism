#' Downloads NLDAS light data
#' @description This function downloads NLDAS incoming shortwave radiation data
#' (w m-2) for a given Latitude and Longitude.
#'
#' @param save_dir The save directory for files to be placed in. For example, "C:/myfolder
#' @param Site_ID The site ID, for example "NC_NHC"
#' @param Lat The site Latitude
#' @param Lon The site Longitude
#' Selecting only data from the central sensor in network (NS2).
#' @param startDate The starting date for the download (YYYY-MM-DD)
#'
#' @return Returns a time series of incoming shortwave solar radiation from the start
#' date to the most recent available data
#' @export

##===============================================================================
## Function for downloading NLDAS light data via data rods
## https://disc.gsfc.nasa.gov/information/tools?title=Hydrology%20Data%20Rods
## Created 01/19/2024
## cite: Xia, Y., et al., NCEP/EMC (2009), NLDAS Primary Forcing Data L4 Hourly 0.125 x 0.125 degree V002, Edited by David Mocko, NASA/GSFC/HSL, Greenbelt, Maryland, USA, Goddard Earth Sciences Data and Information Services Center (GES DISC), Accessed: [Data Access Date], 10.5067/6J5LHHOHZHN4
##===============================================================================
## LIGHT ##
NLDAS_DL <- function(save_dir, Site_ID, Lat, Lon, startDate, endDate){
  #The initial string to build the URL
  http_string <- paste("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:DSWRFsfc")

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
  destfile <- paste(save_dir, "/", Site_ID, "_NLDAS.asc", sep = "")

  #Error catch in case the page is inaccessible. A little inelegant at present...
  try_result <- try(download.file(url, destfile, method = "libcurl"), silent = FALSE)

  if(class(try_result) == "try-error") {file.remove(destfile)}

} #End DL_NLDAS function


##############################
# call dat for 1 site at time
##############################

## NEARSHORE ##
NLDAS_DL(
  save_dir = "./RawData/NLDAS/",
  Site_ID = "BWNS2",
  Lat = "39.10697",
  Lon = "-120.15721",
  startDate = "2020-01-01",
  endDate = "2023-11-01"
)

NLDAS_DL(
  save_dir = "./RawData/NLDAS/",
  Site_ID = "SSNS2",
  Lat = "39.13904",
  Lon = "-120.15236",
  startDate = "2020-01-01",
  endDate = "2023-11-01"
)


NLDAS_DL(
  save_dir = "./RawData/NLDAS/",
  Site_ID = "SHNS2",
  Lat = "39.0944",
  Lon = "-119.9436",
  startDate = "2020-01-01",
  endDate = "2023-11-01"
)


NLDAS_DL(
  save_dir = "./RawData/NLDAS/",
  Site_ID = "GBNS2",
  Lat = "39.0880",
  Lon = "-119.9421",
  startDate = "2020-01-01",
  endDate = "2023-11-01"
)

###################
# offshore sensors 

NLDAS_DL(
  save_dir = "./NLDAS/Offshore/light/",
  Site_ID = "GB10m",
  Lat = "39.08836",
  Lon = "-119.9474",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)


NLDAS_DL(
  save_dir = "./NLDAS/Offshore/light/",
  Site_ID = "GB15m", 
  Lat = "39.088485",
  Lon = "-119.948979",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)

NLDAS_DL(
  save_dir = "./NLDAS/Offshore/light/",
  Site_ID = "GB20m", 
  Lat = "39.088271",
  Lon = "-119.950734",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)

NLDAS_DL(
  save_dir = "./NLDAS/Offshore/light/",
  Site_ID = "BW10m", 
  Lat = "39.10629",
  Lon = "-120.15701",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)

NLDAS_DL(
  save_dir = "./NLDAS/Offshore/light/",
  Site_ID = "BW20m", 
  Lat = "39.106182", 
  Lon = "-120.156302",
  startDate = "2020-01-01",
  endDate = "2023-10-01"
)

###########
## streams

NLDAS_DL(
  save_dir = "./NLDAS/stream/light/",
  Site_ID = "BWL", # 39.1075414	-120.1646811
  Lat = "39.1075414", 
  Lon = "-120.1646811",
  startDate = "2020-01-01",
  endDate = "2024-01-01"
)



NLDAS_DL(
  save_dir = "./NLDAS/stream/light/",
  Site_ID = "GBL", # 39.0880435	-119.9389446
  Lat = "39.0880435", 
  Lon = "-119.9389446",
  startDate = "2020-01-01",
  endDate = "2024-01-01"
)

# end of script.