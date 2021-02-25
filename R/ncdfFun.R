# ---------------------------------------------------------------------------- #
# functions to read ncdf data
# ---------------------------------------------------------------------------- #

.read.ncdf <- function(path,fn){
  nf <- ncdf4::nc_open(paste0(path,fn))
  data <- ncdf4::ncvar_get(nf,varid=names(nf$var)[1])
  ncdf4::nc_close(nf)
  data
}

.read.ncdf.var <- function(path,fn,varname){
  nf <- ncdf4::nc_open(paste0(path,fn))
  data <- ncdf4::ncvar_get(nf,varid=varname)
  ncdf4::nc_close(nf)
  data
}

.get.ncdf.varnames <- function(path,fn){
  nf <- ncdf4::nc_open(paste0(path,fn))
  data <- names(nf$var)
  ncdf4::nc_close(nf)
  data
}

# helper function to extract proper time objects from netcdf file.
# https://stackoverflow.com/questions/46001573/convert-a-netcdf-time-variable-to-an-r-date-object
.get.ncdf.time <- function(path=NULL, fn=NULL, filecon=NULL, timevar=NULL) {
  if(!is.null(path) || !is.null(fn) && is.null(filecon) && is.null(timevar)){
    nf <- ncdf4::nc_open(paste0(path,fn))
  } else if(is.null(path) && is.null(fn) && !is.null(filecon) || !is.null(
    timevar)) {
    nf <- filecon
  } else {
    stop(paste0("No valid get functionalities for requested arguments: path=\"",
                path, "\", ","fn=\"",fn, "\", ","filecon=\"",filecon, "\", ",
                "timevar=\"",timevar))
  }
  if(is.null(timevar)){
    ncdims <- names(nf$dim) #get netcdf dimensions
    timevar <- ncdims[which(ncdims %in% c("T","time", "Time", "datetime", 
                                          "Datetime", "date", "Date"))[1]] 
                                          #find time variable
  }
  times <- ncdf4::ncvar_get(nf, timevar)
  if (length(timevar)==0) stop(
    "ERROR! Could not identify the correct time variable")
  timeatt <- ncdf4::ncatt_get(nf, timevar) #get attributes
  timedef <- strsplit(timeatt$units, " ")[[1]]
  timeunit <- timedef[1]
  tz <- timedef[5]
  timestart <- strsplit(timedef[4], ":")[[1]]
  if (length(timestart) != 3 || timestart[1] > 24 || 
      timestart[2] > 60 || timestart[3] > 60 || any(timestart < 0)) {
    warning(paste("Warning:", timestart, 
                  "not a valid start time. Assuming 00:00:00\n"))
    timedef[4] <- "00:00:00"
  }
  if (! tz %in% OlsonNames()) {
    warning(paste("Warning:", timestart, 
                  "not a valid start time. Assuming 00:00:00\n"))
    tz <- "UTC"
  }

  timestart <- lubridate::ymd_hms(paste(timedef[3], timedef[4]), tz=tz)

  f <- switch(tolower(timeunit), #Find the correct lubridate time function based on the unit
              seconds=lubridate::seconds, second=lubridate::seconds, sec=lubridate::seconds,
              minutes=lubridate::minutes, minute=lubridate::minutes, min=lubridate::minutes,
              hours=lubridate::hours,     hour=lubridate::hours,     h=lubridate::hours,
              days=lubridate::days,       day=lubridate::days,       d=dlubridate::ays,
              months=base::months,        month=base::months,        m=base::months,
              years=lubridate::years,     year=lubridate::years,     yr=lubridate::years,
              NA
  )
  suppressWarnings(if (is.na(f)) stop(
    "Could not understand the time unit format"))

  if (is.integer(times)) {
    timestart + f(times)
  } else {
    timestart + f(as.integer(times))
  }
}

.get.ncdf.dimnames <- function(path,fn, which_dim=NULL){
  # open netcdf file connection
  nf <- ncdf4::nc_open(paste0(path,fn))
  # get dim vals entries
  data <- lapply(nf$dim,function(x) y=x$vals)
  # if which_dim is provided return specific dimnames
  if(!is.null(which_dim)){
    # conditional for empty or non existing entries
    if(which_dim %in% names(nf$dim)){
      if(which_dim %in% c("T","time", "Time", "datetime", "Datetime", "date", 
                          "Date")){
        data = list(.get.ncdf.time(filecon=nf, timevar= which_dim))
      }else{
        # return as list entry to ensure similar returns
        data = list(get(data[which_dim]))
      }
    } else {
      warning(paste0("which_dim \"", which_dim,
                     "\" is either NULL or does not seem to exist."))      
    }
  }
  ncdf4::nc_close(nf)
  data
}
