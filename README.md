## Tool for modelling primary production on a global scale.
Based on a previous study project, restructured and extended to include modern R programming methods & styles.

#### Installation
To currently use the tool please use
```R
devtools::install_github("jnnsbrr/PhotoBioDynamics")
```

#### Functionality and Usage
To apply the model, recent data ([Irradiance](https://opendap.larc.nasa.gov/opendap/SRB/LPSA/SRB_REL3.0_LPSA_MONTHLY_NC/rel3.0/contents.html), [FPAR/LAI](https://icdc.cen.uni-hamburg.de/thredds/catalog/ftpthredds/modis_lai_fpar/global/catalog.html), [Temperature](http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.CPC/.GHCN_CAMS/.gridded/) and [CO2](ftp://aftp.cmdl.noaa.gov/products/trends/co2)) can be downloaded and processed by using the [getInputData](https://github.com/jnnsbrr/PhotoBioDynamics/blob/54fc2df3c11bfde4691b091d55af1bc20ce1fd0e/R/dataModelWrangling.R#L358) function.  
The [modelRunner](https://github.com/jnnsbrr/PhotoBioDynamics/blob/54fc2df3c11bfde4691b091d55af1bc20ce1fd0e/R/modelRunner.R#L32) uses the core model functions in [modelFun.R](./R/modelFun.R) to simulate gross primary production (GPP) and net primary production (NPP) on the global scale. The outputs come as monthly gridded time series.  

The following visualization was rendered by [PhotoBioDynamics.Rmd](./PhotoBioDynamics.Rmd):  

![./PhotoBioDynamics.png](https://github.com/jnnsbrr/PhotoBioDynamics/blob/main/PhotoBioDynamics.png)
