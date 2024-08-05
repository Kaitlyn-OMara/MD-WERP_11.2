library(stringr)
library(raster)


filelist <- list.files(path="G:/MDB 11.2 data/2023/Master data and scripts", pattern="monthly_water_percent_landsat_mndwi_gt0_cplt_20_slc_no_clip57_yes_\\d{8}\\.tif") 

raster.water.data <- matrix(nrow=length(filelist), ncol=1)  
raster.water.data <- data.frame(raster.water.data)
names(raster.water.data) <- c("water? 0=no")

files <- filelist
files <- str_extract(files, "[0-9]{8,}")

Nrast <- length(filelist)  # Need to have a vector or matrix of monthly raster files
# Nrast <- 360
rownames(raster.water.data) <- files

sink("run check water values.r")
for(i in 1:Nrast) {
  print(quote(source("check water presence in month rasters.r")))
}
sink()
index <- 0

source("run check water values.r")
