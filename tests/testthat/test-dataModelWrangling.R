# test download and processing input files functions

ping_url <- function(url_in,t=1){
  con <- url(url_in)
  check <- suppressWarnings(try(open.connection(con,open="rt",timeout=t),silent=T)[1])
  suppressWarnings(try(close.connection(con),silent=T))
  ifelse(is.null(check),TRUE,FALSE)
}


# test download irradiance daa
test_that("Irradiance data available", {

  test_parurl <- "https://opendap.larc.nasa.gov/opendap/SRB/LPSA/SRB_REL3.0_LPSA_MONTHLY_NC/rel3.0/2000/contents.html"
  ping = ping_url(test_parurl)
  testthat::expect_true(ping)
  
  test_pardataurl = "https://opendap.larc.nasa.gov/opendap/SRB/LPSA/SRB_REL3.0_LPSA_MONTHLY_NC/rel3.0/2000/srb_rel3.0_lpsa_monthly_200001.nc"
  ping = ping_url(test_pardataurl)
  testthat::expect_true(ping)
  
})


# test download and process CO2 data
test_that("CO2 data available", {
  
  test_CO2url = 'ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt'
  ping = ping_url(test_CO2url)
  testthat::expect_true(ping)
})


# test downloading of FPAR and LAI data
test_that("FPAR/LAI data available", {
  
  test_FPARLAIurl = 'https://icdc.cen.uni-hamburg.de/thredds/catalog/ftpthredds/modis_lai_fpar/global/catalog.html'
  ping = ping_url(test_FPARLAIurl)
  testthat::expect_true(ping)

  test_dataFPARLAIurl = 'https://icdc.cen.uni-hamburg.de/thredds/fileServer/ftpthredds/modis_lai_fpar/global/2015/MODIS-C006_MOD15A2__LAI_FPAR__LPDAAC__GLOBAL_0.5degree__UHAM-ICDC__20151227__fv0.02.nc'
  ping = ping_url(test_FPARLAIurl)
  testthat::expect_true(ping)
  
})


# test downloading and processing of Temperature data
test_that("Temperature data available", {

  test_tempurl = 'https://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.CPC/.GHCN_CAMS/.gridded/.deg0p5/data.nc'
  ping = ping_url(test_tempurl)
  testthat::expect_true(ping)
  

})

