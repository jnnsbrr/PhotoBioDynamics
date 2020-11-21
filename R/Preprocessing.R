
# packages 
# require("raster")
# require("ncdf4")
# require("lubridate")
# require("doParallel")
# require("tidyverse")
# require("rvest")
# # require("stringr")

# FUN ####

.combineNCData = function(path, fn_const, y, var, dates) {
  out = read.ncdf.var(path = path, 
                      fn = paste0(fn_const, 
                                  dates[y],
                                  ".nc"),
                      varname = var)
}

# abind along third dimension
acomb = function(...) abind::abind(..., along = 3)


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
.wranglingTemperatureData = function(years, path.originaldata) {
  
  url_temp = "https://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.CPC/.GHCN_CAMS/.gridded/.deg0p5/data.nc"
  download.file(url = url_temp,
                destfile = paste0(path.originaldata, "/orig/air.mon.mean2.nc"),
                mode="wb")
  temp = read.ncdf(paste0(path.originaldata, "/orig/"),"air.mon.mean2.nc")
  
  # read and assign time metadata
  time_dim = get.ncdf.dimnames(paste0(path.originaldata,"/orig/"),
                               "air.mon.mean2.nc",
                               which_dim = "T")[[1]]
  dimnames(temp)[[3]] = as.character(time_dim)
  
  # extract required time interval
  time_extract = which(!is.na(match(year(time_dim), years)))
  
  temp = temp[c(361:720,1:360),360:1,time_extract] %>% 
    aperm(c(2,1,3))
  
  # project data on raster brick
  temp_ras = raster(res=0.5, crs = crs("+init=epsg:4326")) %>% 
    brick(., nl = dim(temp)[3]) %>% 
    setValues(., values = temp)
  
  saveRDS(object=temp_ras, file="./data/processed_grid/temp3.Rds")
}


# download and process CO2 data from Mauna Loa (good representativity) #####
.wranglingCO2Data = function(years, path.originaldata) {
  
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
    filter(year >= min(years) & year <= max(years))
  
  saveRDS(object=co2_tab, file = "./data/processed_grid/CO2_tab3.Rds")
}


# download and process Leaf Area Index (LAI) & Fraction of Absorbed PAR (FPAR) ####
.wranglingLAIFPARData = function(years, cluster, path.originaldata) {
  
  if (file.exists("./data/processed_grid/CO2_tab3.Rds")) {
    co2_tab =   readRDS(file = "./data/processed_grid/CO2_tab3.Rds")

  } else {
    stop("Please run photodynamics::.wranglingCO2Data() first.")
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
  lai_ts = foreach::foreach(y = 1:length(dates_int), 
                            .combine = 'acomb', 
                            .multicombine = TRUE,
                            .export=c(# "dates_to_extract, path.originaldata",
                              ".combineNCData")
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
  
  
  lai_ras = brick(temp_ras, nl=dim(lai_ts)[3]) %>% 
    setValues(., lai_ts) %>% 
    # mean over 5 years to get rid of outliers and gaps (cells with no data)
    stackApply(., indices = rep(1:12,times = nlayers(.)/12), fun = mean,na.rm = TRUE)
  
  saveRDS(object=lai_1985_ras, file="./data/processed_grid/lai_3.Rds")
  
  
  # parallel foreach with multidim array returning
  fpar_ts = foreach::foreach(y = 1:length(dates_to_extract), 
                             .combine = 'acomb', 
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
  fpar_ras = brick(temp_ras, nl=dim(fpar_ts)[3]) %>% 
    setValues(., fpar_ts) %>% 
    # mean over 5 years to get rid of outliers and gaps (cells with no data)
    stackApply(., indices = rep(1:12,times = nlayers(.)/12), fun = mean,na.rm = TRUE)
  
  saveRDS(object=fpar_1985_ras, file="./data/processed_grid/fpar_3.Rds")
  
}


# download and process Photosynthetic Active Radiation (PAR) ####
.wranglingPARData = function(years, cluster, path.originaldata) {
  
  par_url = "https://opendap.larc.nasa.gov/opendap/SRB/LPSA/SRB_REL3.0_LPSA_MONTHLY_NC/rel3.0/"
  page = xml2::read_html(paste0(par_url, "contents.html"))
  
  data_nodes = page %>% 
    rvest::html_nodes("td a") %>%
    rvest::html_attr("href")
  
  
  date_format = "([1-2][0-9]{3})"
  par_years_nodes = stringr::str_extract(data_nodes,date_format) %>% 
    as.integer() %>% 
    match(par_years,.)
  
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
                            .combine = 'acomb', 
                            .multicombine = TRUE,
                            .export=c(#"path.originaldata, dates_to_extract",
                              ".combineNCData")
  ) %dopar% {
    # model call
    .combineNCData(path = paste0(path.originaldata, "/orig/irradiance/"),
                   fn_const = "irradiance_sw_",
                   y = y,
                   var = "sw_sfc_net",
                   dates = dates_to_extract)
  }
  
  par_ts = aperm(par_ts,c(2,1,3)) %>% 
    .[180:1,c(181:360,1:180),]
  
  par_ras = raster(res=c(1,1), crs=crs(temp_ras)) %>% 
    brick(., nl=dim(par_ts)[3]) %>% 
    setValues(., values = par_ts) %>% 
    disaggregate(., fact = 2) %>%
    stackApply(.,indices = rep(1:12, times = dim(par_ts)[3]/12), 
               fun = mean, 
               na.rm = TRUE) %>% 
    mask(., mask = temp_ras[[1]])
  
  saveRDS(object=par_ras, file="./data/processed_grid/par_3.Rds")
}


#' Climate data wrangler 
#'
#' Download and process required climate data that is open access and currently
#' (year 2020) available.
#' Processed data is written to the package's data folder.
#'
#' @param year integer. Define a reference year for recycled data (fpar/lai)
#'    end for the time series of outputs to be written. 
#' 
#' @param delete.originaldata logical. Delete the original/raw data that could
#'    likely comprises a large amount of memory. Defaults to `TRUE`
#'
#' @return None
#'
#' @examples \dontrun{
#'
#'  wranglingClimateData(year = 2000)
#' }
#' 
#' @export
#' 

wranglingClimateData = function(year, delete.originaldata = TRUE){


  # calc specific ranges for each dataset 
  time_range = (year-29):year
  fpar_lai_years = (year-4):year
  par_years = 1984:2007

  path_originaldata = user_data_dir("SBD") %>% 
    gsub("\\\\", "/", .)

  if (!dir.exists(paste0(path_originaldata, "/orig"))) {
    dir.create(paste0(path_originaldata, "/orig"), recursive = TRUE)
  }
  
  # detect cores for parallel backend
  cores = detectCores()
  
  # start cluster backend for parallelization
  cl = parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(cl)
  
  # read and process data
  
  .wranglingTemperatureData(years = time_range, 
                            path.originaldata = path_originaldata)  
  
  .wranglingCO2Data(years = time_range,
                            path.originaldata = path_originaldata)  
  
  .wranglingLAIFPARData(years = fpar_lai_years, 
                        cluster = cl,
                        path.originaldata = path_originaldata)  
  
  .wranglingPARData(years = time_range, 
                    cluster = cl, 
                    path.originaldata = path_originaldata)  
  
  
  # close parallel cluster backend
  parallel::stopCluster(cl)
  
  if(delete.originaldata) unlink(paste0(path.originaldata, "/orig"))
    
}



# former LAI FPAR ####

# # fraction absorbed photosynthetic radiation (FPAR) ####
# names_phenology_1985 = get.ncdf.varnames("./data/raw/","global_phenology_1985.nc4")
# fpar_1985 = read.ncdf.var("./data/raw/","global_phenology_1985.nc4","fPAR_dom_biome")
# lai_1985 = read.ncdf.var("./data/raw/","global_phenology_1985.nc4","LAI_dom_biome")
# # harmonize dimensions
# fpar_1985 = aperm(fpar_1985,c(2,1,3))
# lai_1985 = aperm(lai_1985,c(2,1,3))
# # project data on raster brick
# fpar_1985_ras = brick(fpar_1985,xmn=-180, xmx=180, ymn=-90, ymx=90)
# lai_1985_ras = brick(lai_1985,xmn=-180, xmx=180, ymn=-90, ymx=90)
# # assign coordinate reference system and save as RDS
# crs(fpar_1985_ras) = crs("+init=epsg:4326")
# crs(lai_1985_ras) = crs("+init=epsg:4326")
# saveRDS(object=fpar_1985_ras, file="./data/processed_grid/fpar_1985.Rds")
# saveRDS(object=lai_1985_ras, file="./data/processed_grid/lai_1985.Rds")


# former PAR ####

# # ray -> photosynthetic active ray (PAR) ####
# r_mean = read.ncdf("./data/raw/","swdown_erainterim_1901-2010.nc")
# # extract required time interval 
# r_mean = r_mean[,280:1,1009:1320]#1020
# r_mean = r_mean[,280:1,1273:1284]#1020
# 
# r_mean = aperm(r_mean,c(2,1,3))
# # compute means for every month across time interval, currently no faster approach avail
# idx = lapply(1:12,FUN = function(x)(seq(0,(dim(r_mean)[3]-12),12)+x))
# for (ii in 1:12){
#   print(ii)
#   r_mean[,,ii] = apply(r_mean[,,idx[[ii]]],c(1,2), mean, na.rm=T)
# }
# # delete redundant years/data slices
# r_mean = r_mean[,,-c(13:312)]
# # project data on raster brick
# r_mean_ras = brick(r_mean,xmn=-180, xmx=180, ymn=-56, ymx=84)
# # harmonize extent
# r_mean_ras = extend(r_mean_ras,extent(c(-180,180,-90,90)))
# # calculate photosynthetic active radiation
# par_mean_ras = (0.5/0.27)*(r_mean_ras*0.0864)# /24/60/60*10^6*0.22 # to par mol/m^2/day to W/m^2/day
# # assign coordinate reference system and save as RDS
# crs(par_mean_ras) = crs("+init=epsg:4326")
# saveRDS(object=par_mean_ras, file="./data/processed_grid/par_mean1985-2010.Rds")

