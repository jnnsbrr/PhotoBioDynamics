# ---------------------------------------------------------------------------- #
# download and process model input of required format and boundaries (time)
# ---------------------------------------------------------------------------- #

# download and process global radiation data ####
.downloadIrradianceNCs = function(ff, url, nodes, path.originaldata) {
  
  page = xml2::read_html(paste0(url, nodes[ff]))
  
  data_nodes = page %>% 
    rvest::html_nodes("td a") %>% 
    .[rvest::html_text(.) == "file"] %>% 
    rvest::html_attr("href")
  
  if (!dir.exists(paste0(path.originaldata, "/orig/irradiance"))) {
    dir.create(paste0(path.originaldata, "/orig/irradiance"))
  }
  
  y_suburl = sub("contents.html","",nodes[ff])
    sapply(data_nodes, function(x) {
      download.file(url = paste0(url, 
                                 y_suburl,
                                 x),
                    destfile = paste0(path.originaldata,
                                      "/orig/irradiance/",
                                      sub(pattern = "srb_rel3.0_lpsa_monthly",
                                          replacement = "irradiance_sw", 
                                          x = x)),
                    mode ="wb")}) %>% 
      invisible()
}


# download and process global FPAR and LAI data ####
.downloadFPARLAINCs = function(yy, years, co2_tab, path.originaldata) {
  
  url_catalog = "https://icdc.cen.uni-hamburg.de/thredds/catalog/ftpthredds/modis_lai_fpar/global/"
  url_download = "https://icdc.cen.uni-hamburg.de/thredds/fileServer/ftpthredds/modis_lai_fpar/global/"
  
  # web scraping
  page = xml2::read_html(paste0(url_catalog,
                                years[yy],
                                "/catalog.html"))
  # find data nodes
  data_nodes = page %>% 
    rvest::html_nodes("td a") %>%
    rvest::html_attr("href")
  
  # estimate and convert required timespan
  date_format = "([2-9][0-9]{7})"
  dates_raw = str_extract(data_nodes,date_format) %>% 
    lubridate::ymd()
  
  
  dates_int = sapply(X = co2_tab$date[lubridate::year(co2_tab$date) == years[yy]], 
                     FUN = function(x){y = dates_raw[which(abs(
                       dates_raw-x) == min(abs(
                         dates_raw - x), na.rm = TRUE))] %>%  
                       ifelse(length(.) > 1,.[1], .)})  
  
  dates_to_extract = dates_int %>% 
    as.Date(origin = "1970-01-01") %>% 
    format(format = "%Y%m%d")
  
  if (!dir.exists(paste0(path.originaldata, "/orig/fpar_lai"))) {
    dir.create(paste0(path.originaldata, "/orig/fpar_lai"))
  }
  
  
  # download required files for time vectors (within one year)
  lapply(X = dates_to_extract,
         FUN = function(x) download.file(url = paste0(
          url_download,
          years[yy],
          "/MODIS-C006_MOD15A2__LAI_FPAR__LPDAAC__GLOBAL_0.5degree__UHAM-ICDC__",
          x,
          "__fv0.02.nc"), 
         destfile = paste0(path.originaldata, 
                           "/orig/fpar_lai/MODIS_LAI_FPAR_",
                           x,
                           ".nc"),
         mode = "wb")
  ) %>% 
    invisible()
  
  # get dates 
  return(dates_to_extract)
  
}


# download and process global temperature data ####
.getTemperatureData = function(years, path.originaldata, 
                               path.inputdata = "./data/input/") {
  
  url_temp = "https://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.CPC/.GHCN_CAMS/.gridded/.deg0p5/data.nc"
  download.file(url = url_temp,
                destfile = paste0(path.originaldata, "/orig/air.mon.mean2.nc"),
                mode="wb") %>% 
    invisible()
  temp = .read.ncdf(paste0(path.originaldata, "/orig/"),"air.mon.mean2.nc")
  
  # read and assign time metadata
  time_dim = .get.ncdf.dimnames(paste0(path.originaldata,"/orig/"),
                               "air.mon.mean2.nc",
                               which_dim = "T")[[1]]
  dimnames(temp)[[3]] = as.character(time_dim)
  
  # extract required time interval
  time_extract = which(!is.na(match(lubridate::year(time_dim), years)))
  
  temp = temp[c(361:720,1:360),360:1,time_extract] %>% 
    aperm(c(2,1,3))
  
  # project data on raster brick
  temp_ras = raster::raster(res=0.5, crs = raster::crs("+init=epsg:4326")) %>% 
    raster::brick(., nl = dim(temp)[3]) %>% 
    raster::setValues(., values = temp)
  
  if (!dir.exists(path.inputdata)) {
    dir.create(path = path.inputdata, recursive = TRUE)
  }
  saveRDS(object=temp_ras, file=paste0(path.inputdata, "temp_ras.Rds"))
}


# download and process CO2 data from Mauna Loa (good representativity) ####
.getCO2Data = function(years, path.originaldata, 
                       path.inputdata = "./data/input/") {
  
  url_co2 = 'ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt'
  co2_dat = read.table(file = url_co2,
                       na.strings = -99.99)
  
  names(co2_dat) =c("year","month", "decimal date", "average", "interpolated",
                    "trend" ,"#days")
  
  co2_tab = data.frame(year = co2_dat$year,
                       date = lubridate::date(
                         lubridate::date_decimal(co2_dat$`decimal date`, 
                                                 tz = "UTC")),
                       average = co2_dat$average,
                       interpolated = co2_dat$interpolated,
                       trend = co2_dat$trend) %>% 
    dplyr::filter(year >= min(years) & year <= max(years))
  
  if (!dir.exists(path.inputdata)) {
    dir.create(path = path.inputdata, recursive = TRUE)
  }

  saveRDS(object=co2_tab, file = paste0(path.inputdata, "CO2_tab.Rds"))
}


# download and process Leaf Area Index (LAI) & Fraction of Absorbed PAR (FPAR) ####
.getLAIFPARData = function(years, cluster, path.originaldata, 
                           path.inputdata = "./data/input/") {
  
  if (file.exists("./data/input/CO2_tab.Rds")) {
    co2_tab = readRDS(file = "./data/input/CO2_tab.Rds")

  } else {
    stop("Please run PhotoBioDynamics::.wranglingCO2Data() first.")
  }

  # parallel foreach with multidim array returning
  dates_to_extract = foreach::foreach(yy = 1:length(years),
                                      .combine='c',
                                      .packages = c("tidyverse"),
                                      .export=c(# "dates_to_extract", "years",
                                        # "co2_tab", "path.originaldata",
                                        ".downloadFPARLAINCs")
  ) %dopar% {
    
    .downloadFPARLAINCs(yy = yy,
                        years = years,
                        co2_tab = co2_tab,
                        path.originaldata = path.originaldata)
  }
  
  
  # parallel foreach with multidim array returning
  lai_ts = foreach::foreach(y = 1:length(dates_to_extract), 
                            .combine = '.acomb', 
                            .multicombine = TRUE,
                            .export=c(".combineNCData")
  ) %dopar% {
    # model call
    .combineNCData(path = paste0(path.originaldata, "/orig/fpar_lai/"),
                   fn_const = "MODIS_LAI_FPAR_",
                   y = y,
                   var = "lai",
                   dates = dates_to_extract)
  }
  
  lai_ts = aperm(lai_ts, c(2,1,3)) %>% 
    .[360:1,,]
  lai_ts[lai_ts == 255] = 0
  
  
  lai_ras = raster::raster(res=0.5, crs =raster:: crs("+init=epsg:4326")) %>% 
    raster::brick(., nl=dim(lai_ts)[3]) %>% 
    raster::setValues(., lai_ts) %>% 
    # mean over 5 years to get rid of outliers and gaps (cells with no data)
    raster::stackApply(., indices = rep(1:12,times = raster::nlayers(.)/12),
               fun = mean,na.rm = TRUE)
  
  saveRDS(object=lai_ras, file=paste0(path.inputdata,"lai_ras.Rds"))
  
  
  # parallel foreach with multidim array returning
  fpar_ts = foreach::foreach(y = 1:length(dates_to_extract), 
                             .combine = '.acomb', 
                             .multicombine = TRUE,
                             .export=c(# "dates_to_extract, path.originaldata",
                               ".combineNCData")
  ) %dopar% {
    
    .combineNCData(path = paste0(path.originaldata, "/orig/fpar_lai/"),
                   fn_const = "MODIS_LAI_FPAR_",
                   y = y,
                   var = "fpar",
                   dates = dates_to_extract)
  }
  
  
  fpar_ts = aperm(fpar_ts, c(2,1,3))
  fpar_ts[360:1,,] = fpar_ts
  fpar_ts[fpar_ts == 255] = 0
  fpar_ras = raster::raster(res=0.5, crs = raster::crs("+init=epsg:4326")) %>% 
    raster::brick(., nl=dim(fpar_ts)[3]) %>% 
    raster::setValues(., fpar_ts) %>% 
    # mean over 5 years to get rid of outliers and gaps (cells with no data)
    raster::stackApply(., indices = rep(1:12,times = raster::nlayers(.)/12), 
                       fun = mean, na.rm = TRUE)
  
  if (!dir.exists(path.inputdata)) {
    dir.create(path = path.inputdata, recursive = TRUE)
  }
  
  saveRDS(object=fpar_ras, file=paste0(path.inputdata, "fpar_ras.Rds"))
  
}


# download and process Photosynthetic Active Radiation (PAR) ####
.getPARData = function(paryears, path.originaldata, cluster,
                       path.inputdata = "./data/input/") {
  
  par_url = "https://opendap.larc.nasa.gov/opendap/SRB/LPSA/SRB_REL3.0_LPSA_MONTHLY_NC/rel3.0/"
  page = xml2::read_html(paste0(par_url, "contents.html"))
  
  data_nodes = page %>% 
    rvest::html_nodes("td a") %>%
    rvest::html_attr("href")
  
  
  date_format = "([1-2][0-9]{3})"
  par_years_nodes = stringr::str_extract(data_nodes,date_format) %>% 
    as.integer() %>% 
    match(paryears,.)
  
  nodes_use = data_nodes[par_years_nodes]
  
  
  # parallel foreach with multidim array returning
  foreach::foreach(ff = 1:length(nodes_use),
                   .packages = "tidyverse",
                   .export=c(# "par_url", "nodes_use", "path.originaldata",
                     ".downloadIrradianceNCs")
  ) %dopar% {
    
    # model call
    .downloadIrradianceNCs(ff = ff,
                           url = par_url,
                           nodes = nodes_use,
                           path.originaldata = path.originaldata)
  }
  
  
  dates_to_extract = sub("/contents.html","",nodes_use) %>%
    outer(X = ., Y = c("01", "02", "03", "04", "05", "06", "07", "08", "09", 
                       10:12), 
          FUN=paste0) %>% 
    t() %>% 
    as.vector()
  
  
  # parallel foreach with multidim array returning
  par_ts = foreach::foreach(y = 1:length(dates_to_extract), 
                            .combine = '.acomb', 
                            .multicombine = TRUE,
                            .export=c(#"path.originaldata, dates_to_extract",
                              ".combineNCData")
  ) %dopar% {
    # model call
    .combineNCData(path = paste0(path.originaldata, "/orig/irradiance/"),
                   fn_const = "irradiance_sw_",
                   y = y,
                   var = "sw_sfc_dn",
                   dates = dates_to_extract)
  }
  
  # calculate photosynthetic active radiation
  par_ts <- (0.5/0.27)*(par_ts*0.0864) %>% # /24/60/60*10^6*0.22 # to par mol/m^2/day from W/m^2/day
    aperm(.,c(2,1,3)) %>% 
    .[180:1,c(181:360,1:180),]
  
  if (file.exists("./data/input/lai_ras.Rds")) {
    lai_ras =   readRDS(file = "./data/input/lai_ras.Rds")
    
  } else {
    stop("Please run photodynamics::.wranglingLAIFPARData() first.")
  }
  
  par_ras = raster::raster(res=c(1,1), crs=raster::crs(lai_ras)) %>%
    raster::brick(., nl=dim(par_ts)[3]) %>% 
    raster::setValues(., values = par_ts) %>% 
    raster::disaggregate(., fact = 2) %>%
    raster::stackApply(.,indices = rep(1:12, times = dim(par_ts)[3]/12), 
               fun = mean, 
               na.rm = TRUE) %>% 
    mask(., mask = lai_ras[[1]])
  
  if (!dir.exists(path.inputdata)) {
    dir.create(path = path.inputdata, recursive = TRUE)
  }
  
  saveRDS(object=par_ras, file=paste0(path.inputdata,"par_ras.Rds"))
}


#' Get input data 
#'
#' Download and process required input data that is open access and currently
#' (year 2020) available. 
#' This may take a while, depending on your connection.
#' Processed model input data is stored at "./data/input
#'
#' @param year integer. Define a reference year for recycled data (fpar/lai)
#'    end for the time series of outputs to be written. 
#' 
#' @param delete.originaldata logical. Delete the original/raw data that could
#'    likely comprises a large amount of memory. Defaults to `TRUE`
#'
#' @return None. Optional raw output is written to 
#'    `rappdirs::user_data_dir("PhotoBioDynamics")`, processed model input into
#'    "./data/input".
#'
#' @examples \dontrun{
#'
#'  getInputData(year = 2000)
#' }
#' 
#' @export
#' 

getInputData = function(year, delete.originaldata = TRUE){


  # calc specific ranges for each dataset 
  time_range = (year-29):year
  fpar_lai_years = (year-4):year
  par_years = 1984:2007

  path_originaldata = rappdirs::user_data_dir("PhotoBioDynamics") %>% 
    gsub("\\\\", "/", .)

  if (!dir.exists(paste0(path_originaldata, "/orig"))) {
    dir.create(paste0(path_originaldata, "/orig"), recursive = TRUE)
  }
  
  # detect cores for parallel backend
  cores = parallel::detectCores()
  
  # start cluster backend for parallelization
  cl = parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(cl)
  
  # read and process data
  
  .getTemperatureData(years = time_range, 
                            path.originaldata = path_originaldata)  
  
  .getCO2Data(years = time_range,
                            path.originaldata = path_originaldata)  
  
  .getLAIFPARData(years = fpar_lai_years, 
                        cluster = cl,
                        path.originaldata = path_originaldata)  
  
  .getPARData(paryears = par_years, 
                    cluster = cl, 
                    path.originaldata = path_originaldata)  
  
  
  # close parallel cluster backend
  parallel::stopCluster(cl)
  
  if(delete.originaldata) unlink(paste0(path_originaldata, "/orig"))
    
}

