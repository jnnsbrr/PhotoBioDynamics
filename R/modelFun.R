
# parameter initialization ####

CONST = list(
  # average plant albedo
  albedo = 1-0.17,
  # radiation transmission correction factor
  alphaa = 0.5,
  # O2 partial pressure [Pa]
  po2 = 20.9e3,
  # atmospheric pressure [Pa]
  p = 1.0e5,
  # leaf respiration, fraction of Vmax
  bc3 = 0.015,
  # co-limitation (shape) parameter
  theta = 0.7,
  # inter-cellular to ambient CO2 partial pressure (0.6-0.8) [Wong et al.1979; Long and Hutchin 1991]
  lambda = 0.7,
  # molar mass of carbon
  cmass = 12,
  daylength = 24,
  #  intrinsic quantum efficiency of CO2 uptake in C3 plants
  alphac3 = 0.08,
  # q10 for temperature-sensitive parameter ko
  q10ko = 1.2,
  # q10 for temperature-sensitive parameter kc
  q10kc = 2.1,
  # q10 for temperature-sensitive parameter tau
  q10tau = 0.57,
  tau25 = 2600.0,
  # value of tau at 25 deg C
  ko25 = 3.0e4,
  # Michaelis constant for CO2 (Pa) at 25 deg C
  kc25 = 30.0,
  months_length = c(31,(28*3+29)/4,31,30,31,30,31,31,30,31,30,31),
  mean_month_length = mean(c(31,(28*3+29)/4,31,30,31,30,31,31,30,31,30,31)),
  ## autotrophic respiration specific constants
  CONST_ATR = list(
    aa = 1,
    ln = 50/365, # in gC/m2
    # divided by mean_month_length
    kr = 1.67/30.47917, # in gC /kg /day 
    tref = 10, # in degC
    t0 = 46.02, # in degC
    e0 = 308.56, # in K  
    cn = 0.001 # in g C/ m2
  )
)

# model function ####

modelProduction = function(outvar, i, temp, par, fpar, lai, co2, K){

  ## data prep
  
  # possibility to return multi output variables ("GPP", "NPP", "NEP")  
  if (length(outvar) > 1) {
    multivar = array(NA, dim = c(dim(temp)[1:2],1,length(outvar)))
    j = 1
  }
  
  # subset of a month
  temp_i = temp[,,i]
  co2_i = co2[i]
  
  # repeating index vector for input variables par, fpar and lai that only 
  # comprises 1 year/12 month
  idseq = rep(1:12,times=dim(temp)[3]/12)
  par_i = par[,,idseq[i]]
  lai_i = lai[,,idseq[i]]
  fpar_i = fpar[,,idseq[i]]
  
  
  ## GPP PART 1: PAR-limited photosynthesis
  
  # Calculating APAR absorbed photosynthetic active ray (mol/d/m2) after 
  # Schaphoff et al. 2018
  # while for fpar a remote sensing product is used
  apar = par_i * fpar_i * K$alphaa
  
  ko = K$ko25 * exp((log(K$q10ko)) * (temp_i - 25) * 0.1) # O2 (Pa)
  kc = K$kc25 * exp((log(K$q10kc)) * (temp_i - 25) * 0.1) # CO2 (Pa)
  
  # Calculating temperature dependance of CO2 / O2 specific ratio (Pa/Pa)
  tau = K$tau25 * exp((log(K$q10tau)) * (temp_i - 25) * 0.1)
  
  # internal partial pressure of CO2 (Pa) # CO2 partial pressure as average
  pint = K$lambda * co2_i
  
  # CO2 compensation point (Pa) - where Photosynthesis = Respiration
  gammastar = K$po2 / (2 * tau) 
  
  # (Haxeltine and Prentice 1996)
  phitemp = ((1 + exp(0.2*(10-temp_i)))^-1)#+(-1-exp(0.6*(30-temp_i)))^-1 
  
  c1 = K$alphac3 * phitemp * ((pint - gammastar) / (pint + gammastar))
  
  # PAR-limited photosynthesis rate molC/m2/h
  je = c1 * K$cmass * apar
  
  
  ##  GPP PART 2: Rubisco limited rate of photosynthesis - dependence on 
  ## photorespiration
  
  fac = kc * (1 + K$po2/ ko) 
  
  c2 = (pint - gammastar)/(pint + fac) 
  
  # remove not needed objects since vectorized approach uses lot of memory 
  rm(gammastar, pint, fac, tau, par_i, fpar_i, ko, kc)

  # Calculating maximum daily rate of net photosynthis with Rubisco capacity Vm
  s = (24 / K$daylength) * K$bc3 # daily average Rd/Vm
  
  sigma = 1 - (c2 - s) / (c2 - K$theta * s) 
  
  vm = (1.0 / K$bc3) * (c1 / c2) * ((2.0 * K$theta - 1.0) * s - (
    2.0 * K$theta * s - c2) * sigma) * apar * K$cmass
  
  # Calculation of rubisco-activity-limited photosynthesis rate JC, molC/m2/d
  jc = c2 * vm  
  
  # Calculation of daily gross photosynthesis, gpp, gC/d/m2
  gpp = (je+jc-sqrt((je+jc)*(je+jc)-4.0*K$theta*je*jc))/(2.0*K$theta) 
  gpp[gpp<0] = 0
  
  # remove not needed objects since vectorized approach uses lot of memory 
  rm(je, jc, c1, c2, sigma, apar)
  gc()
  
  if ("GPP" %in% outvar && length(outvar) == 1) {
    return(gpp * rep(K$months_length, dim(temp)[3]/12))
    
  } else if ("GPP" %in% outvar && length(outvar) >1) {
    multivar[,,,j] = gpp * rep(K$months_length, dim(temp)[3]/12)
    j = j+1
    
  }
  
  ## NET PRIMARY PRODUCTION (NPP) [gC/d/m2] (Haxeltine and Prentice 1996) 
  
  # Calculation of daily autotrophic respiration gC/d/m2
  cs = lai_i * K$CONST_ATR$cn
  
  rleaf = K$bc3 * vm
  
  rtrans = K$CONST_ATR$kr * cs*exp(K$CONST_ATR$e0*((
    (K$CONST_ATR$tref-K$CONST_ATR$t0)^-1)-(temp_i-K$CONST_ATR$t0)^-1))
  
  rfine = K$CONST_ATR$aa * lai_i * K$CONST_ATR$ln
  
  rgrowth = gpp * 0.2
  
  rauto = rleaf + rtrans + rfine + rgrowth
  
  # Daily net primary production (NPP), And, gC/d/m2 
  npp = gpp - rauto
  npp[npp < 0] = 0
  
  # remove not needed objects since vectorized approach uses lot of memory 
  rm(vm, cs, gpp, rleaf, rtrans, rfine, rgrowth, rauto, lai_i)
  gc()
  
  if ("NPP" %in% outvar && length(outvar) == 1){
    return(npp * rep(K$months_length, dim(temp)[3]/12))
    
  }else if("NPP" %in% outvar && length(outvar) >1) {
    multivar[,,,j] = npp* rep(K$months_length, dim(temp)[3]/12)
    j = j+1
    
  }
    
  ## Net Ecosystem Production (NEP) derived from Arrhenius-equation 
  ## (in mol CO2 /m2/d) (Lloyd and Taylor 1994)
  # 
  # rhet = (0.1*exp(-30/(temp_i-120)*12)-0.99) * 44 *12^-1 * 1e-6 * 60 * 60 *24
  # rhet =  R10 * exp(308.56 *(1 / 56 - 1 / (temp_i + 46))) 
  # rhet = 0.84 # to be parametrized (Cui. Yi-Bin et al. 2020)
  # nep = npp - rhet 
  # nep[nep < 0] = 0
  # 
  # # remove not needed objects since vectorized approach uses lot of memory 
  # rm(rhet, temp_i)
  # gc()
  # 
  # if ("NEP" %in% outvar && length(outvar) == 1){
  #   return(nep * rep(K$months_length, dim(temp)[3]/12))
  #   
  # } else if("NEP" %in% outvar && length(outvar) > 1) {
  #   multivar[,,,j] = nep * rep(K$months_length, dim(temp)[3]/12)
  #   j = j+1
  # }
  
  if(exists("multivar")) return(multivar)
  
  # Returns a multidimensional numerical array including two dimensions
  #    representing latitude and longitude, a time dimension and if set 
  #    a categorical dimensions reflecting possible outputes "GPP",
  #    "NPP" and/or "NEP".
  
}