#' Download static data 
#'
#' Download and process shapefiles for plotting.
#' This may take a while, depending on your connection.
#'
#' @param delete_zip boolean. Defaults to TRUE. Whether to delete zip folder
#'    after download
#' 
#' @return None. Output is written to 
#'    `rappdirs::user_data_dir("PhotoBioDynamics")`.
#'
#' @examples \dontrun{
#'
#'  downloadSpecialNaturalEarthData(check_valid = TRUE, delete_zip = TRUE)
#' }
#' 
#' @export
#' 

downloadSpecialNaturalEarthData = function(delete_zip = TRUE) {
  ne_data = c("ne_50m_wgs84_bounding_box.zip", "ne_50m_graticules_30.zip",
              "ne_50m_ocean.zip")    
  path_originaldata = rappdirs::user_data_dir("PhotoBioDynamics") %>% 
    gsub("\\\\", "/", .)
  if (!dir.exists(paste0(path_originaldata, "/orig/naturalearthdata"))) {
    dir.create(path = paste0(path_originaldata, "/orig/naturalearthdata"), 
               recursive = TRUE)
  }
    sapply(ne_data, function(fn) {
      url_temp = paste0("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/physical/", 
                        fn)
      tryCatch(
        {
          download.file(url = url_temp,
                        destfile = paste0(path_originaldata, 
                                          "/orig/naturalearthdata/", fn),
                        mode="wb")
        },
        error = function(cond) {
          message(paste0(" "))
          message(cond)
          message(paste0(" "))
          message(paste0("! YOU MAY HAVE TO FIND AN ALTERNATIVE SOURCE. !"))
        })
      
      unzip(paste0(path_originaldata, "/orig/naturalearthdata/", fn),
            exdir=paste0(path_originaldata, "/orig/naturalearthdata"))
      if (delete_zip) {
        file.remove(paste0(path_originaldata, "/orig/naturalearthdata/", fn))
      }
    }) %>% 
    invisible()
    # for convenienve
    return(path_originaldata)
}
