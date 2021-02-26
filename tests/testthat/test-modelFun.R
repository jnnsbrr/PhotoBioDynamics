# test modelFun
test_that(".modelFun", {

  temp_ras = readRDS(file="../testdata/temp_ras.Rds")
  temp = raster::as.array(temp_ras)
  
  fpar = readRDS(file="../testdata/fpar_ras.Rds") %>% 
    raster::as.array(.)
  lai = readRDS(file="../testdata/lai_ras.Rds") %>% 
    raster::as.array(.)
  par = readRDS(file="../testdata/par_ras.Rds") %>% 
    raster::as.array(.)
  co2 = readRDS(file="../testdata/CO2_tab.Rds")$interpolated
  
  test_gpp = .modelPrimaryProduction(outvar ="GPP",
                                     i = 1,
                                     temp = temp,
                                     par = par,
                                     fpar = fpar,
                                     lai = lai, 
                                     co2= co2,
                                     paramfile = "../../inst/data/input/model_parameters.YAML")

  testthat::expect_equal(dim(test_gpp),c(360,720))
  
  test_gpp = array(test_gpp,c(dim(test_gpp),1))
  test_value = sum(test_gpp * raster::as.array(area(temp_ras)), 
                   na.rm = TRUE)*1e-9
  
  testthat::expect_gt(test_value, 5)
})