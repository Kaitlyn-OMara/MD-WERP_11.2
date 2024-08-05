##MD-WERP project 11.2 ##
## STEP 1 - defining persistent waterbodies in the Darling R ##
##Kaitlyn O'Mara and Ben Stewart-Koster##

#Load in the raster file #
Water <- raster("water_persistence_percent_landsat_mndwi_gt0_cplt_20_slc_no_clip57_yes19870525_20230104.tif") #import raster tif

## use these next 4 lines if creating the polygons on the cropped area##
#e <- as(extent(776000, 779000,-3482000, -3478000), 'SpatialPolygons') ##create an area to crop to, to have a smaller area to play with to test analyses
#crs(Water) #check the co-ordinate reference system (crs) of the Water raster file
#crs(e) <- "+proj=utm +zone=54 +datum=WGS84 +units=m +no_defs" #assign the created crop area the same crs
#r <- crop(Water, e) #crop the larger file to the smaller area

W <- (Water >= 80) ##create a new raster layer that only has the values >=80 and they become 1 while the rest become 0

wc <- clump(W, gaps = FALSE, directions = 8) ## clumping water cells

clump.poly <- rasterToPolygons(wc) # convert clumps to polygons

#aggregate initial polygon output into single multi-part features
#this solves some unrealistic feature separation which resulted from
#the raster to polygon conversion
sp.a = aggregate(clump.poly, dissolve = TRUE)

#disaggregate multi-part polygon into separate features
sp.d <- disaggregate(sp.a)

#apply a negligible buffer distance around polygon features to combine features 
#that are only connected by one raster corner (and therefore counted 
#incorrectly as distinct features after the above steps)
#and dissolve
sp.b <- buffer(sp.d, width = 0.001, dissolve = TRUE)

#disaggregate multi-part buffered polygon into separate features
sp.bd <- disaggregate(sp.b)


##drop waterbodies smaller than 5 pixels(0.0044 km2) or 3 pixels (0.0026 km2)
area_thresh <- units::set_units(0.0044, km^2)
waterbodies <- drop_crumbs(sp.bd, threshold = area_thresh)
( p.df <- data.frame( ID=1:length(waterbodies)) ) 
rownames(p.df)


# Extract polygon ID's
( pid <- sapply(slot(waterbodies, "polygons"), function(x) slot(x, "ID")) )

# Create dataframe with correct rownames
( p.df <- data.frame( ID=1:length(waterbodies), row.names = pid) )    

# Try coersion again and check class
p <- SpatialPolygonsDataFrame(waterbodies, p.df)
class(p) 
names(p) #check it has an ID column
crs(p) #check it has a crs
p


#load in a polygon of the reach during a wet month. This was made by saving the cells that had water during a wet month tif
# as a single polygon. Alternatively, a wet month raster can be loaded into the script and just use the water values from that raster

#perc10_poly <- st_read("C:/Users/s5042382/OneDrive - Griffith University/MD-WERP 11.2/MD-WERP 11.2 Data files for MDBA/Code/whole reach as polygon 10perc persist.shp")
perc10_poly <- st_read("C:/Users/s5042382/OneDrive - Griffith University/MD-WERP 11.2/2023/Master data and scripts/whole reach as polygon 10perc persist.shp")
expanded_channel <- st_buffer(perc10_poly, dist = 30) #expand the channel to fit most of the centroids

#convert to raster to make it a cost raster for reordering waterbodies from upstream to downstream
current.raster = terra::ext(expanded_channel)
current.raster = terra::rast(current.raster,
                             nrow=6456, ncol=8497, 
                             crs = crs(expanded_channel), 
                             vals = 1)
current.raster = terra::mask(current.raster, 
                             terra::vect(expanded_channel),
                             updatevalue = NA)

#get a centroids layer of the persistent waterbodies
centroids <- st_centroid(st_as_sf(p))

#find the most northern centroid (most upstream) so that we can calculate the distance from every other waterbody to this one
most_northern_centroid <- centroids[which.max(st_coordinates(centroids)[, 2]), ] #locate the most northern centroid

#check these things looks like they should by plotting a cropped area of upstream section
#e <- as(extent(870000, 900000,-3400000,-3380000), 'SpatialPolygons') ##create an area to crop to, to have a smaller area to play with to test analyses
#e <- as(extent(630000, 700000,-3600000,-3550000), 'SpatialPolygons') ##create an area to crop to, to have a smaller area to play with to test analyses
#crs(e) <- "+proj=utm +zone=54 +datum=WGS84 +units=m +no_defs" #assign the created crop area the same crs
#cost_zoom <- crop(current.raster, e) #crop the larger file to the smaller area
#plot(cost_zoom)
#plot(p, add=T)
#plot(st_geometry(most_northern_centroid), col = "black", bg = "green", pch = 21, cex = 2, add = TRUE)

#if it looks good, all the waterbodies are inside the cost raster and the point is the most upstream, then we can calculate distances
#calculate distances of each of the raster cells to the most upstream waterbody centroid
#first we need to put our centroids into the right cells of the cost raster:
point.locs = terra::cellFromXY(current.raster,
                               sf::st_coordinates(centroids))
waterbody.1 = terra::cellFromXY(current.raster, 
                               sf::st_coordinates(most_northern_centroid))
#mark the most upstream waterbody as a special cell by giving it the value of 2
terra::values(current.raster)[waterbody.1] = 2
#this tells us which cell it's in
which(terra::values(current.raster) == 2) #Confirms we were successful in making a cell somewhere a cost value of 2!
#now we use gridDist to get the distances of our least-cost paths from every cell that isn't NA to all cells holding the "target" value
#this will take a while because it's computing distances for all the cells
all.dists = terra::gridDistance(current.raster, target = 2, scale = 1, maxiter=500)# the terra vignette asks for a target but it sometimes gives an error if that argument is used instead of origin=2 so see which one works for you, can also try reloading the terra package

#now we can get the distances by location the centroids within the correct cells
centroids_matrix <- st_coordinates(centroids)
extracted <- terra::extract(x = all.dists, y = centroids_matrix)
# there's some centroids that have NA values returned, because the centroid wasn't within the cost raster
#so we have to find out which cells they are and give them the value of the nearest non NA cell
#using the nearestLand function from the seegSDM package (this package needs to be installed via the github install)
outside_mask <- is.na(extracted)
outside_pts <- centroids_matrix[outside_mask,]
oustide_pts_df <- as.data.frame(outside_pts)

all_dists_raster <- raster(all.dists) #this converts it to a rasterLayer from a Spatraster

# find the nearest cell
nearest.cells <- seegSDM::nearestLand(outside_pts, all_dists_raster, 10000)
#Extracting the raster values using the new coordinates:
replaced_vals <- extract(x = all_dists_raster, y = nearest.cells)
# Identify indices of NA values in extracted
na_indices <- which(is.na(extracted[,1])) 
# Replace NA values in extracted with the values from the nearest cells
extracted[na_indices, 1] <- replaced_vals


#add as a column to the persistent waterbody layer and call it RANK
ranks <- rank(extracted, ties.method = "first") #reorders the values first, based on dist from the upstream point
p$RANK <- ranks #add this column to our persistent waterbodies layer
rank.wb.dist <- as.data.frame(cbind(extracted, ranks)) #this is just to check the ranks are correct
#Now reorder the waterbodies by using the rank column as the new ID number 
row.names(p@data) <- p@data[,2]
waterbodies.reordered <- st_as_sf(p)
waterbodies.reordered
waterbodies.reordered <- waterbodies.reordered[ order(as.numeric(row.names(waterbodies.reordered))), ]
waterbodies.reordered$clump.ID <- waterbodies.reordered$ID
waterbodies.reordered <- waterbodies.reordered[,2:3]
waterbodies.reordered.df <- as_Spatial(waterbodies.reordered)

#to save the persistent waterbodies layer after reordering and check it's been corrected properly
#library(rgdal)
#writeOGR(obj=waterbodies.reordered.df, dsn="C:/Users/s5042382/OneDrive - Griffith University/MD-WERP 11.2/MD-WERP 11.2 Data files for MDBA/Code", layer="persistent waterbodies test using gridDistance reordering method V2", driver="ESRI Shapefile")

wb <- st_as_sf(waterbodies.reordered)

#dist.dat <- read.csv("Waterbody_distances_2023_80_5.csv") #read in a distance matrix of calculated distanced between the edges of each waterbody (we did this in arcGIS)

#alternatively, use the code below which uses the distances between centroids, for ease of analysis but be 
#aware of these distances if using them to make ecological conclusions
# Sort the dataframe by ranks
rank.wb.dist <- rank.wb.dist[order(rank.wb.dist$ranks), ]
# Calculate the differences between consecutive values in lyr.1
differences <- diff(rank.wb.dist$lyr.1)
## Create the new dataframe dist.dat
dist.dat <- data.frame(
  REGION1 = rank.wb.dist$ranks[1:(nrow(rank.wb.dist) - 1)],
  REGION2 = rank.wb.dist$ranks[2:nrow(rank.wb.dist)],
  DISTANCE_M = differences
)

adj.dist <- matrix(0, nrow=nrow(dist.dat)+1, ncol=nrow(dist.dat)+1) # Set up the adjacency matrix for the graph


for (i in 1:nrow(dist.dat)) {
  adj.dist[min(dist.dat[i, 4:5]), max(dist.dat[i, 4:5])]  <- dist.dat[i,6]  # Fill the adjacency matrix by working through the rows of the distance dataset 
}
