# ---------------------------------------------------------------------------- #
# helper functions
# ---------------------------------------------------------------------------- #

# parallel helper function to combine ncdf data
.combineNCData = function(path, fn_const, y, var, dates) {
  out = .read.ncdf.var(path = path, 
                      fn = paste0(fn_const, 
                                  dates[y],
                                  ".nc"),
                      varname = var)
}

# abind along third dimension
.acomb = function(...) abind::abind(..., along = 3)

