# test .combineNCDatac


test_that(".combineNCData works", {
  
  path_originaldata = rappdirs::user_data_dir("PhotoBioDynamics") %>% 
    gsub("\\\\", "/", .)
  
  .combineNCData(path = paste0(path_originaldata, "/orig/fpar_lai/"),
                 fn_const = "MODIS_LAI_FPAR_",
                 y = y,
                 var = "lai",
                 dates = dates_to_extract)
})



# test .acomb
test_that(".acomb works", {
  test_array <- array(data = 1, dim = c(1,1,2)) %>% 
    .acomb(.,.,.)
  testthat::expect_equal(dim(a1)[3], 6)
})


