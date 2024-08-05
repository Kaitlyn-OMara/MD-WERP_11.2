library(stringr)
library(raster)


filelist <- list.files(path="G:/MDB 11.2 data/2023/Master data and scripts", pattern="monthly_water_percent_landsat_mndwi_gt0_cplt_20_slc_no_clip57_yes_\\d{8}\\.tif") 


files <- filelist
files <- str_extract(files, "[0-9]{8,}")

Nrast <- length(filelist)  # Need to have a vector or matrix of monthly raster files
# Nrast <- 360

raster.value.data <- matrix(nrow=length(filelist), ncol=1)  
raster.value.data <- data.frame(raster.value.data)
names(raster.value.data) <- c("sum raster values")
rownames(raster.value.data) <- files

sink("run check values.r")
for(i in 1:Nrast) {
  print(quote(source("check values of month rasters.r")))
}
sink()
index <- 0

source("run check values.r")
