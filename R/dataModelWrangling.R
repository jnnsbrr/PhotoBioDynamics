# ---------------------------------------------------------------------------- #
# download and process model input of required format and boundaries (time)
# ---------------------------------------------------------------------------- #

# download and process global temperature data ####
.getTemperatureData = function(years, path.originaldata, 
                               path.inputdata = NULL) {
  cat("Downloading temperature data ...")
  
  url_temp = "https://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.CPC/.GHCN_CAMS/.gridded/.deg0p5/data.nc"
  download.file(url = url_temp,
                destfile = paste0(path.originaldata, "/orig/air.mon.mean2.nc"),
                mode="wb", quiet = TRUE)
  cat("Completed. Now wrangling data ... \n")
  
  temp = .read.ncdf(paste0(paste0(path.originaldata, "/orig/")),"air.mon.mean2.nc")
  
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
  
  if (is.null(path.inputdata)) {
    saveRDS(object = temp_ras, file = paste0(path.originaldata, "/input/", 
                                         "temp_ras.Rds"))
  } else {
    saveRDS(object = temp_ras, file = paste0(path.inputdata, "temp_ras.Rds"))
  }
  cat("Completed. \n")
  
}

# download and process CO2 data from Mauna Loa (good representativity) ####
.getCO2Data = function(years, path.originaldata, 
                       path.inputdata = NULL) {
  cat("Downloading and wrangling CO2 data ... \n")
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
  
  
  if (is.null(path.inputdata)) {
    saveRDS(object=co2_tab, file = paste0(path.originaldata, "/input/", 
                                         "CO2_tab.Rds"))
  } else {
    saveRDS(object=co2_tab, file = paste0(path.inputdata, "CO2_tab.Rds"))
  }
  cat("Completed. \n")
}


# download and process Leaf Area Index (LAI) & Fraction of Absorbed PAR (FPAR) ####
.getLAIFPARData = function(years, path.originaldata, path.inputdata = NULL) {
  
  if (file.exists(paste0(path.originaldata, "/input/CO2_tab.Rds"))) {
    co2_tab = readRDS(file = paste0(path.originaldata, "/input/CO2_tab.Rds"))

  } else {
    stop("Please run PhotoBioDynamics::.wranglingCO2Data() first.")
  }
  
  if (!dir.exists(paste0(path.originaldata, "/orig/fpar_lai"))) {
    dir.create(paste0(path.originaldata, "/orig/fpar_lai"))
  }
  
  url_catalog = "https://icdc.cen.uni-hamburg.de/thredds/catalog/ftpthredds/modis_lai_fpar/global/"
  url_download = "https://icdc.cen.uni-hamburg.de/thredds/fileServer/ftpthredds/modis_lai_fpar/global/"

  # estimate and convert required timespan
  date_format = "([2-9][0-9]{7})"
  
  dates_to_extract = lapply(years, function(x) {
    page = xml2::read_html(paste0(url_catalog, x, "/catalog.html"))

    data_nodes = page %>% 
      rvest::html_nodes("td a") %>%
      rvest::html_attr("href")
    dates_raw = stringr::str_extract(data_nodes,date_format) %>% 
      lubridate::ymd()
    
    dates_int = sapply(co2_tab$date[lubridate::year(co2_tab$date) == x], 
                       FUN = function(xx){y = dates_raw[which(abs(
                         dates_raw-xx) == min(abs(
                           dates_raw - xx), na.rm = TRUE))] %>%  
                         ifelse(length(.) > 1,.[1], .)})
    dates_to_download = dates_int %>% 
      as.Date(origin = "1970-01-01") %>% 
      format(format = "%Y%m%d")
    
    data.frame(year=rep(x,length(dates_to_download)),suburl=dates_to_download)
  }) %>% dplyr::bind_rows()
  

  
  cat("Downloading LAI and FPAR data ... \n")
  
  # detect cores for parallel backend
  cores = parallel::detectCores()-1 # one core to be left
  
  # start cluster backend for parallelization
  cl = parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(cl)
  
  # parallel foreach with multidim array returning
  foreach::foreach(yy = 1:nrow(dates_to_extract)) %dopar% {
    
    download.file(url = paste0(
      url_download,
      dates_to_extract$year[yy],
      "/MODIS-C006_MOD15A2__LAI_FPAR__LPDAAC__GLOBAL_0.5degree__UHAM-ICDC__",
      dates_to_extract$suburl[yy],
      "__fv0.02.nc"), 
      destfile = paste0(path.originaldata, 
                        "/orig/fpar_lai/MODIS_LAI_FPAR_",
                        dates_to_extract$suburl[yy],
                        ".nc"),
      mode = "wb", quiet = TRUE)
  }
  
  cat("Completed. Wrangling LAI data ... \n")
  
  # parallel foreach with multidim array returning
  lai_ts = foreach::foreach(y = 1:nrow(dates_to_extract), 
                            .combine = '.acomb', 
                            .multicombine = TRUE,
                            .export=c(".combineNCData")
  ) %dopar% {
    # model call
    .combineNCData(path = paste0(path.originaldata, "/orig/fpar_lai/"),
                   fn_const = "MODIS_LAI_FPAR_",
                   y = y,
                   var = "lai",
                   dates = dates_to_extract$suburl)
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
  
  if (is.null(path.inputdata)) {
    saveRDS(object=lai_ras, file=paste0(path.originaldata, "/input/", 
                                         "lai_ras.Rds"))
  } else {
    saveRDS(object=lai_ras, file=paste0(path.inputdata,"lai_ras.Rds"))
  }
  
  cat("Completed. Wrangling FPAR data ...\n")
  
  # parallel foreach with multidim array returning
  fpar_ts = foreach::foreach(y = 1:nrow(dates_to_extract), 
                             .combine = '.acomb', 
                             .multicombine = TRUE,
                             .export=c(# "dates_to_extract, path.originaldata",
                               ".combineNCData")
  ) %dopar% {
    
    .combineNCData(path = paste0(path.originaldata, "/orig/fpar_lai/"),
                   fn_const = "MODIS_LAI_FPAR_",
                   y = y,
                   var = "fpar",
                   dates = dates_to_extract$suburl)
  }
  
  # close parallel cluster backend
  parallel::stopCluster(cl)
  
  fpar_ts = aperm(fpar_ts, c(2,1,3))
  fpar_ts[360:1,,] = fpar_ts
  fpar_ts[fpar_ts == 255] = 0
  fpar_ras = raster::raster(res=0.5, crs = raster::crs("+init=epsg:4326")) %>% 
    raster::brick(., nl=dim(fpar_ts)[3]) %>% 
    raster::setValues(., fpar_ts) %>% 
    # mean over 5 years to get rid of outliers and gaps (cells with no data)
    raster::stackApply(., indices = rep(1:12,times = raster::nlayers(.)/12), 
                       fun = mean, na.rm = TRUE)
  
  if (is.null(path.inputdata)) {
    saveRDS(object = fpar_ras, file = paste0(path.originaldata, "/input/", 
                                         "fpar_ras.Rds"))
  } else {
    saveRDS(object = fpar_ras, file = paste0(path.inputdata, "fpar_ras.Rds"))
  }
  cat("Completed.\n")
  
}


# download and process Photosynthetic Active Radiation (PAR) ####
.getPARData = function(paryears, path.originaldata, path.inputdata = NULL) {
  
  cat("Download irradiance data ...\n")
  
  if (!dir.exists(paste0(path.originaldata, "/orig/irradiance"))) {
    dir.create(paste0(path.originaldata, "/orig/irradiance"))
  }
  
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
  
  nodes = lapply(nodes_use, function(x) {
    page = xml2::read_html(paste0(par_url, x))
    suburl = sub("contents.html","",x)
    data_nodes = page %>% 
      rvest::html_nodes("td a") %>% 
      .[rvest::html_text(.) == "file"] %>% 
      rvest::html_attr("href")
    data.frame(year=rep(suburl,length(data_nodes)),suburl=data_nodes)
    }) %>% dplyr::bind_rows()
  
  # detect cores for parallel backend
  cores = parallel::detectCores() - 1# one core to be left
  
  # start cluster backend for parallelization
  cl = parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(cl)
  
  # parallel foreach with multidim array returning
  foreach::foreach(ff = 1:nrow(nodes)) %dopar% {

    download.file(url = paste0(par_url, nodes$year[ff],"/", nodes$suburl[ff]),
                  destfile = paste0(path.originaldata,
                                    "/orig/irradiance/",
                                    sub(pattern = "srb_rel3.0_lpsa_monthly",
                                        replacement = "irradiance_sw", 
                                        x = nodes$suburl[ff])),
                  mode ="wb", quiet = TRUE)
  }

  dates_to_extract = sub("/contents.html","", nodes_use) %>%
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
  
  # close parallel cluster backend
  parallel::stopCluster(cl)
  
  # calculate photosynthetic active radiation
  par_ts <- (0.5/0.27)*(par_ts*0.0864) %>% # /24/60/60*10^6*0.22 # to par mol/m^2/day from W/m^2/day
    aperm(.,c(2,1,3)) %>% 
    .[180:1,c(181:360,1:180),]
  
  if (file.exists(paste0(path.originaldata, "/input/lai_ras.Rds"))) {
    lai_ras = readRDS(file = paste0(path.originaldata, "/input/lai_ras.Rds"))
    
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
  
  if (is.null(path.inputdata)) {
    saveRDS(object = par_ras, file = paste0(path.originaldata, "/input/", 
                                             "par_ras.Rds"))
  } else {
    saveRDS(object=par_ras, file=paste0(path.inputdata, "par_ras.Rds"))
  }
  cat("Completed. \n")
  
}


#' Get input data 
#'
#' Download and process required input data that is open access and currently
#' (year 2020) available. 
#' This may take a while, depending on your connection.
#'
#' @param latestyear integer. Define the latest year for the input data. 
#'    LAI/FPAR and PAR are averaged and recycled since their long-term variations 
#'    play a minor role within this scope. 
#'    
#' @param timerange integer. Define the timerange to be used for Temperature and
#'    CO2 data. Defaults to `30`
#' 
#' @param delete.originaldata logical. Delete the original/raw data that could
#'    likely comprises a large amount of memory. Defaults to `TRUE`
#'
#' @return None. Output is written to 
#'    `rappdirs::user_data_dir("PhotoBioDynamics")`
#'
#' @examples \dontrun{
#'
#'  getInputData(year = 2000)
#' }
#' 
#' @export
#' 

getInputData = function(latestyear, timerange=30, delete.originaldata = TRUE){

  options(warn=2)
  # calc specific ranges for each dataset 
  time_range = (latestyear-timerange+1):latestyear
  latestyear = latestyear
  fpar_lai_years = latestyear
  # fpar_lai_years = (latestyear-4):latestyear # 5 year average - too slow
  par_years = c(2002) # median year of last sunspot cycle 
  # par_years = 1996:2007 # averaging sunspot cycle- too slow
  path_originaldata = rappdirs::user_data_dir("PhotoBioDynamics") %>% 
    gsub("\\\\", "/", .)

  if (!dir.exists(paste0(path_originaldata, "/orig"))) {
    dir.create(paste0(path_originaldata, "/orig"), recursive = TRUE)
  }
  if (!dir.exists(paste0(path_originaldata, "/input"))) {
    dir.create(paste0(path_originaldata, "/input"), recursive = TRUE)
  }

  # read and process data

  .getTemperatureData(years = time_range, path.originaldata = path_originaldata)  

  .getCO2Data(years = time_range, path.originaldata = path_originaldata)  

  .getLAIFPARData(years = fpar_lai_years, path.originaldata = path_originaldata)  

  .getPARData(paryears = par_years, path.originaldata = path_originaldata)  

  if(delete.originaldata) unlink(paste0(path_originaldata, "/orig"))
    
}

