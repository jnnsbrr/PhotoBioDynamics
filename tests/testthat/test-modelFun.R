# test modelFun
test_that(".modelFun", {

  temp = readRDS(file="../../data/tests/temp_ras.Rds")
  fpar = readRDS(file="../../data/tests/fpar_ras.Rds")
  lai = readRDS(file="../../data/tests/lai_ras.Rds")
  par = readRDS(file="../../data/tests/par_ras.Rds")
  co2 = readRDS(file="../../data/tests/CO2_tab.Rds")
  
  test_gpp = .modelPrimaryProduction(outvar ="GPP",
                                     i = 1,
                                     temp = temp,
                                     par = par,
                                     fpar = fpar,
                                     lai = lai, 
                                     co2= co2,
                                     paramfile = "../../data/input/model_parameters.YAML")

  testthat::expect_equal(dim(test_gpp),c(360,720,1))
  
  test_value = cellStats(x = test_gpp * area(test_gpp),
                         stat = sum, 
                         na.rm = TRUE
                         )*1e-9
  
  testthat::expect_gt(test_value, 5)
})