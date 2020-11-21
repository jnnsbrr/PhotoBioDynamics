# require("doParallel")
# require("parallel")
# require("abind")
# require("raster")

#' A model wrapper to run different models
#'
#' Simulation of C3 photosynthesis in terms of Gross Primary Production (GPP), 
#' as well as Net Primary Production (NPP) and Net Ecosystem Production (NEP) 
#' based on a modified Farquhar model (Schaphoff et al. 2019, Haxeltine and 
#' Prentice 1996, Farquhar et al. 1980). 
#' This is a vectorized global gridd approach that uses directly and currently 
#' (year 2020) available remote sensing products (open access).
#'
#' @param path character. define a path to write output
#' 
#' @param outvar single character or vector. Define the output variables,
#'    choose between gross primary production - "GPP", net primary prodcution 
#'    - "NPP" and/or net ecosytem production "NEP 
#' 
#' @param format single character or vector. Choose between raster time series
#'    for the defined \code{"outvar"} - "GRID", a global sum time series - 
#'    "GLOBAL_TS" and/or a summed time series for each of 14 biomes ("BIOME_TS).
#'
#' @return None
#'
#' @examples \dontrun{
#'
#'  modelRunner(outvar = c("GPP", "NPP"), format = c("GLOBAL_TS", "BIOME_TS"))
#' }
#' 
#' @export
#' 

modelRunner = function(path, 
                      outvar = c("GPP", "NPP", "NEP"),
                      format = c("GLOBAL_TS", "BIOME_TS","GRID")) {
  
  # pb <- txtProgressBar(style = 3)
  # setTxtProgressBar(pb = pb, value = 0.3)
  
  ## load pre-processed RDS (former large netCDF) files as raster
  
  # extract array for model performance and simple dim binding   
  # temperature time series raster/array
  temp_ras = readRDS(file = "./data/input/temp_ras.Rds")-273.15#in °C
  temp = raster::as.array(temp_ras)
  
  # radiation time series raster -> array
  par = readRDS(file = "./data/input/par_ras.Rds") %>% 
    raster::as.array(.)
  
  # fraction of absorbed radiation raster -> array
  fpar = readRDS(file = "./data/input/fpar_ras.Rds") %>% 
    raster::as.array(.)
  
  # leaf area index raster -> array
  lai = readRDS(file = "./data/input/lai_ras.Rds") %>% 
    raster::as.array(.)
  
  # co2 time series
  co2_tab = readRDS(file = "./data/input/CO2_tab.Rds")
  co2 = co2_tab$interpolated
  ## run model in parallel
  
  # providing a parallel backend cluster for parallelisation
  cl = parallel::makeCluster(2,outfile = "")
  doParallel::registerDoParallel(cl)
  
  # parallel foreach with multidim array returning
  pp = foreach::foreach(i = 1:(dim(temp)[3]), 
                        .combine = '.acomb', 
                        .multicombine = TRUE, 
                        .export = c("modelProduction", "CONST")
  ) %dopar% {
    # model call
    modelProduction(outvar = outvar, 
                    i = i, 
                    temp = temp,
                    par = par, 
                    fpar = fpar, 
                    lai = lai, 
                    co2 = co2, 
                    K = CONST)
  }
  
  # close cluster
  parallel::stopCluster(cl)
  
  ## post processing
  
  # write as gridded output
  if ("GRID" %in% format) {
    
    for (cc in 1:dim(pp)[4]) {
      
      # recycle temp_ras brick 
      pp_ras = raster::setValues(x = temp_ras,
                                 values = pp[,,,cc]
      )
      # save as processed raster
      if (missing(path)) {
        saveRDS(object = pp_ras, 
                file = paste0("./data/output/pp_", outvar[cc], "_ras.Rds"))
      } else {
        saveRDS(object = pp_ras, 
                file = paste0(path,"/pp_",outvar[cc],"_ras.Rds"))
      }
    }
  }
  
  # write as global time series
  if ("GLOBAL_TS" %in% format) {
    
    # calculate area of grid cells using raster package functionality
    pp_area = raster::area(x = temp_ras) %>% 
      raster::as.array(.) %>% 
      "*"(1e-9) %>% 
      "*"(30.47917)
    pp_area = sweep(x = pp,
                    MARGIN = 1, 
                    STATS = pp_area, 
                    FUN = "*",
                    check.margin = FALSE)
    
    # calc global sums for timeseries as well as categories    
    global_ts = apply(pp_area, c(3,4), sum, na.rm=T)
    
    dimnames(global_ts) = list(co2_tab$year,
                              outvar)
    
    # write global ts
    if (missing(path)){
      
      saveRDS(object = global_ts,
              file = paste0("./data/output/globalTS_", 
                            paste(outvar, collapse = "-"), ".Rds"))
    } else {
      saveRDS(object = global_ts, 
              file = paste0(path, "/globalTS_",
                            paste(outvar,collapse = "-"),".Rds"))
    }
    
  }
  
  # write as global time series for each biome
  if ("BIOME_TS" %in% format){
    
    # check if pp_area is already defined, then recycle
    if(!exists("pp_area")) {
      
      # calculate area of grid cells using raster package functionality
      pp_area = raster::area(x = temp_ras) %>% 
        raster::as.array(.) %>% 
        "*"(1e-9) %>% 
        "*"(30.47917)
      
      pp_area = sweep(x = pp,
                      MARGIN = 1, 
                      STATS = pp_area, 
                      FUN = "*",
                      check.margin = FALSE)
    }
    
    # read ascii file of 14 global biomes
    biomes_ras <- raster::raster("./data/raw/14biomes.txt")
    
    # harmonize with default spatial settings used in this context
    biomes_ras <- raster::extend(biomes_ras,extent(temp_ras))
    extent(biomes_ras) <- raster::extent(temp_ras)
    
    # for vectorized extraction use matrix of biomes
    biomes = raster::as.matrix(biomes_ras)
    
    # get values of biomes as vector
    biome_vals = as.vector(na.exclude(as.data.frame(
      raster::freq(biomes_ras))$value))
    
    # parameters as grid of same length for mapply function
    gg = expand.grid(time = 1:dim(pp)[3], biome_vals = biome_vals, 
                     cats = dim(pp)[4])
    
    # use mapply and generate array with corresponding dimensions
    biome_ts = array(
      # calculate time series for each biome and category (GPP, NPP, NEP)
      mapply(function(z,y,x) sum(pp[,,y,z][biomes==x],na.rm=T), 
             y=gg$time,
             x=gg$biome_vals,
             z=gg$cats),
      dim=c(dim(pp)[3],length(biome_vals),dim(pp)[4]))
    
    ## to be parallelized - but facing allocation issues
    # cl <- parallel::makeCluster(detectCores())
    # gg = expand.grid(time=1:dim(pp)[3],biome_vals=biome_vals,cats=dim(pp)[4])
    # biome_ts = array(
    #   parallel::clusterMap(cl,function(z,y,x) sum(pp[,,y,z][biomes==x],na.rm=T), 
    #       y=gg$time,
    #       x=gg$biome_vals,
    #       z=gg$cats,pp=pp,biomes=biomes),
    #   dim=c(dim(pp)[3],length(biome_vals),dim(pp)[4]))
    # parallel::stopCluster(cl)
    
    dimnames(biome_ts) = list(co2_tab$year,
                              c("Trop. Subtrop. Moist Broadleaf Forests", 
                                "Trop. Subtrop. Dry Broadleaf Forests", 
                                "Trop. Subtrop. Moist Coniferous Forests", 
                                "Temperate Broadleaf and Mixed Forests",
                                "Temperate Coniferous Forests", 
                                "Boreal Forests/Taiga",
                                "Tropical and Subtropical Grasslands, Savannas", 
                                "Temperate Grasslands, Savannas and Shrublands", 
                                "Flooded Grasslands and Savannas",
                                "Montane Grasslands and Shrublands",
                                "Tundra",
                                "Mediterranean Forests, Woodlands and Shrub",
                                "Deserts and Xeric Shrublands",
                                "Mangroves",
                                "Lakes",
                                "Rock and Ice"),
                              outvar)
    
    # write global ts
    if (missing(path)){
      saveRDS(object=biome_ts, file=paste0("./data/output/biomeTS_",paste(
        outvar,collapse = "-"),".Rds"))
      
    } else {
      saveRDS(object=biome_ts, file=paste0(path,"/biomeTS_",paste(
        outvar,collapse = "-"),".Rds"))
    }
    
  }
  # close(pb)
}
