# Darling River connectivity Master script #
# MD-WERP Project 11.2 2023-2024 #
# Kaitlyn O'Mara and Ben Stewart-Koster #

#Load required packages#
#Install these first if you don't already have them #
library(sf)
library(raster)
library(igraph)
library(smoothr)
library(units)
library(stringr)
library(terra)
library(seegSDM) #use this code to install: devtools::install_github('SEEG-Oxford/seegSDM')


#Set working directory, make sure to change the file path here to where your scripts and raster tif files are stored
setwd("C:/Users/s5042382/OneDrive - Griffith University/MD-WERP 11.2/MD-WERP 11.2 Data files for MDBA/Code")

# First identify the waterbodies (needs to be done once)

source("Persistent waterbodies.r") # A script to identify the waterbodies based on long term persistence
  # "Create waterbodies.r" produces this 
  # waterbodies <- st_read("p-14waterbodies-zoomed-out.shp") #import cropped waterbodies shapefile created using the 'test water persistence' script

## NOTE: OPEN THE PERSISTENT WATERBODIES SCRIPT AND CHECK THE DIST.DAT SOURCE
#  IF YOU DON'T CHANGE IT, THE DISTANCES USED TO CALCULATE CLOSENESS WILL BE BETWEEN THE WATERBODY CENTROIDS
# WE USED OPTIMAL REGIONS CONNECTIONS TOOL IN ARCGIS TO GET THE DISTANCES BETWEEN THE WATERBODY UPSTREAM AND DOWNSTREAM EDGES, AND READ IT IN USING THE LINE OF CODE HASHED OUT (LINE 147)
# IF YOU ARE NOT USING THE CLOSENESS METRIC (IT WAS NOT USED IN THE REPORT), THEN IT DOESN'T MATTER WHICH DISTANCE SOURCE IS USED SO YOU CAN JUST LEAVE AS IS

#Monthly tifs#
filelist <- list.files(path="G:/MDB 11.2 data/2023/Master data and scripts", pattern="monthly_water_percent_landsat_mndwi_gt0_cplt_20_slc_no_clip57_yes_\\d{8}\\.tif") 

#test monthly tifs#
#filelist <- c("Python_API_Project_112_monthly_water_percent_landsat_mndwi_gt0_cloud_lt20_19950101.tif","Python_API_Project_112_monthly_water_percent_landsat_mndwi_gt0_cloud_lt20_20030301.tif", "Python_API_Project_112_monthly_water_percent_landsat_mndwi_gt0_cloud_lt20_19870701.tif" 
 #             )

files <- filelist
files <- str_extract(files, "[0-9]{8,}")
  
Nrast <- length(filelist)  # Need to have a vector or matrix of monthly raster files


sink("run the whole analysis.r")

for(i in 1:Nrast) {
                   print(quote(source("intersecting polygon area.r")))
                   }
sink()


# Create a data frame to store the Graph theoretic indices - REACH LEVEL
reach.graph.data <- matrix(nrow=length(filelist), ncol=4) # ncol is however many reach level indices we decide to run 
reach.graph.data <- data.frame(reach.graph.data)
names(reach.graph.data) <- c("N components", "N links", "Connected area", "Area of connected waterbodies")
rownames(reach.graph.data) <- files

# Create an empty row to append the monthly data
node.month.data <- matrix(ncol=5)
node.month.data <- data.frame(node.month.data)
names(node.month.data) <- c("month", "overlap.area.data", "Node degree centrality", "Node betweenness centrality", "Closeness")
node.month.data <- node.month.data[-1,]

index <- 0

source("run the whole analysis.r")

reach.graph.data
node.graph.data
node.month.data

date<-rep((files),each=437) #change the value here to the number of waterbodies, if definig waterbodies a different way
date.labels<-data.frame(date)
waterbody.num <- rep(1:437,times=381) #times should equal the number of months in the dataset
waterbody.no <- data.frame(waterbody.num)
node.month.data <-cbind(date.labels, waterbody.no, node.month.data) # add columns for date and waterbody number to node level metric dataset

saveRDS(reach.graph.data, file = "reach.graph.data.RDS") 
saveRDS(node.month.data, file = "node.month.data.RDS") 
write.csv(reach.graph.data,"reach.graph.data.csv", row.names = TRUE)
write.csv(node.month.data, file = "node.month.data.csv", row.names = FALSE)
