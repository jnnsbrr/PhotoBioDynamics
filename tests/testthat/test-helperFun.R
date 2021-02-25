# test helper functions

# test .acomb
test_that(".acomb works", {
  test_array <- array(data = 1, dim = c(1,1,2)) %>% 
    .acomb(.,.,.)
  testthat::expect_equal(dim(test_array)[3], 6)
})


# test .combineNCDatac
test_that(".combineNCData works", {

  path_originaldata = rappdirs::user_data_dir("PhotoBioDynamics") %>% 
    gsub("\\\\", "/", .) %>% 
    paste0("/tests")
  

  download.file(url = "https://icdc.cen.uni-hamburg.de/thredds/fileServer/ftpthredds/modis_lai_fpar/global/2015/MODIS-C006_MOD15A2__LAI_FPAR__LPDAAC__GLOBAL_0.5degree__UHAM-ICDC__20151211__fv0.02.nc",
                destfile = paste0(path_originaldata,"MODIS_LAI_FPAR_20151211.nc"),
                mode = "wb", quiet=TRUE)
  
  test = .combineNCData(path = path_originaldata,
                 fn_const = "MODIS_LAI_FPAR_",
                 y = 1,
                 var = "lai",
                 dates = c(20151211))
  unlink(x = paste0(path_originaldata,"MODIS_LAI_FPAR_20151211.nc"))
  testthat::expect_equal(dim(test),c(720,360))
})


