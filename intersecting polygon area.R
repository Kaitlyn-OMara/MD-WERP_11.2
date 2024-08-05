## test intersecting polygons ##

index <- index + 1

month <- raster(filelist[index])   # Need to create a vector of file names of all monthly tifs or think of an alternative for this line

##following for lines for cropped test area only, don't use for full dataset
#e <- as(extent(770000, 779800,-3485000, -3477000), 'SpatialPolygons') 
#crs(e) <- "+proj=utm +zone=54 +datum=WGS84 +units=m +no_defs" #assign the created crop area the same crs
#wb <- as_Spatial(wb)
#wb <- crop(wb, e) #crop the larger waterbodies file to the smaller area of the test raster
#wb <- st_as_sf(wb)
#######################################################################
######raster to polygon for month data where water is present##########
#month <- raster("crop2monthtest_20030301.tif")


W <- (month >= 2)    # Remove non-water values from the raster (i.e. pixels less than 2)
# plot(W)

#create the polygons of water for the month
wc <- clump(W, gaps = FALSE, directions = 8) 
# freq(wc)
clump.poly <- rasterToPolygons(wc)

# plot(clump.poly)

#aggregate initial polygon output into single multi-part features
#this solves some unrealistic feature separation which resulted from
#the raster to polygon conversion
polya = aggregate(clump.poly, dissolve = TRUE)
# polya      #e.g., features: 1

#disaggregate multi-part polygon into separate features
polyd <- disaggregate(polya)
# polyd      #e.g., features: 228

#apply a negligible buffer distance around polygon features to combine features 
#that are only connected by one raster corner (and therefore counted 
#incorrectly as distinct features after the above steps)
#and dissolve
polyb <- buffer(polyd, width = 0.001, dissolve = TRUE)
# polyb      #e.g., features: 1

#disaggregate multi-part buffered polygon into separate features
month.p <- disaggregate(polyb)
# month.p     #e.g., features: 181



##################################################################################
################################ INTERSECT #######################################

m <- st_as_sf(month.p) # Change format for the intersect function

# Identify which persistent waterbodies intersect with the which monthly waterbodies 
inter.dat <- st_intersects(wb, m, sparse = FALSE)

# Create an empty adjancecy matrix for building the graph
adjac <- matrix(0, nrow=nrow(inter.dat), ncol=nrow(inter.dat))

# Convert inter.dat to binary to work with the adjacency matrix step
inter.dat[inter.dat=="TRUE"] <- 1
inter.dat[inter.dat=="FALSE"] <- 0

# Fill out the adjaceny matrix to show which persistent waterbodies are connected for the month
for (i in 1:(nrow(adjac)-1)) {
                              adjac[i, i+1] <- ifelse(max(colSums(inter.dat[i:(i+1),])) == 2, 1, 0) # Fill out the adjacency matrix  
                              }

#create distance between waterbodies weighted adjacency matrix by multiplying
#the binary matrix with the distance matrix
adjac.dist.weighted <- adjac*adj.dist 

# Build the graph from the adjacency matrix
the.graph <- graph.adjacency(adjac.dist.weighted, weighted=T)

######### REACH LEVEL INDICES ##############
# Extract graph indices
n.components <- components(the.graph)
n.edge <- ecount(the.graph)

# Calculate connected area
poly4area <- st_as_sf(month.p[colSums(inter.dat) >=2])
conn.area <- sum(st_area(poly4area))

waterbodies2count <- rowSums(as.data.frame(inter.dat[,colSums(inter.dat) >=2]))
waterbodies2count <- which(waterbodies2count>0)

# Convert to spatial object for analysis
wb.as.sp <- as(st_geometry(wb), "Spatial")

persist.wb.area.conn <- sum(st_area(st_as_sf(wb.as.sp[waterbodies2count])))

#plot(st_as_sf(wb.as.sp[waterbodies2count]), col= "blue")
#plot(poly4area, add=T, col="green", alpha = 0.2)

# Store graph indices in the reach.graph.data object created in the master script
reach.graph.data[index, ] <- c(n.components$no, n.edge, conn.area, persist.wb.area.conn)


########### NODE LEVEL INDICES #############
# Create a data frame to store the Graph theoretic indices - NODE LEVEL
node.graph.data <- matrix(nrow=437, ncol=3) # nrow is number of waterbodies, ncol is however many node level indices we decide to run 
node.graph.data <- data.frame(node.graph.data)
names(node.graph.data) <- c("Node degree centrality", "Node betweenness centrality", "Closeness")

###Node (waterbody) area###
inter.poly <- st_intersection(wb, m) #creates a shapefile of the intersecting area of the waterbody and month shapefile
inter.poly_sp <- as_Spatial(inter.poly) #convert to spatialpolygonsdataframe to aggregate
inter.agg <- aggregate(inter.poly_sp, by = "RANK") #aggregrate based on RANK (waterbody # that it intersects)
inter.agg <- st_as_sf(inter.agg) #convert back to shapefile
overlapping.area <- st_area(inter.agg) #find areas of aggregated polygons, but there is no value of 0 for waterbodies that weren't intersected
node.area.data <- matrix(nrow=437, ncol=1) # nrow is number of waterbodies, ncol is however many node level indices we decide to run 
node.area.data <- data.frame(node.area.data)
node.area.data[inter.agg$RANK, ] <- st_area(inter.agg)
wb.area.data <- matrix(nrow=437, ncol=1) # nrow is number of waterbodies, ncol is however many node level indices we decide to run 
wb.area.data <- data.frame(wb.area.data)
wb.area.data[wb$RANK, ] <- st_area(wb)
overlap.area.data <- matrix(nrow=437, ncol=1) # nrow is number of waterbodies, ncol is however many node level indices we decide to run 
overlap.area.data <- data.frame(overlap.area.data)
overlap.area.data[wb$RANK, ] <- (node.area.data/wb.area.data*100)



###Node degree centrality###
#What is the node degree of the nodes in our graph, 
#which is the sum of the number of both incoming and outgoing links.
deg1 <-  degree(the.graph, v = V(the.graph), mode = c("all")) ## node degree of all nodes in the network 

### Node betweenness centrality ###
#What is the betweenness centrality of the nodes in our graph, 
#which is the number of shortest paths through the network of which a node is a part
bet1 <-  betweenness(the.graph, v = V(the.graph), directed = FALSE, weights=edge.attributes(the.graph)$weight) ## node degree of all nodes in the network 

##Closeness##
#Uses the weight distances to calculate closeness of one node to every other node#
#Closeness centrality measures how many steps is required to access every other vertex from a given vertex#
close1 <- closeness(the.graph, v = V(the.graph), weights=edge.attributes(the.graph)$weight)

#store the node level indices in the node.graph.data object created in the master script
node.graph.data[,] <- c(deg1, bet1, close1)
month <- rep(index, nrow(node.graph.data))

node.graph.data <-cbind(month, overlap.area.data, node.graph.data)

# output the node level data for the month as its own data frame
# assign(paste("node.graph.data", index, sep="."), node.graph.data)

# append this month's data to previous month's data
node.month.data <- rbind(node.month.data, node.graph.data)
