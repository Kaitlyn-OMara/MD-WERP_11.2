# MD-WERP_11.2
R code to run waterbody scale and reach scale predictor models of connectivity of the Darling River, NSW. Output from MD-WERP project 11.2 2024

How to use the code:
1.	Download all files and make sure all the r scripts, shapefile, and rasters (tifs) are in the same folder
2.	Open the ‘Darling River connectivity master script’ and change the working directory so it has the path to your folder that contains all the files
3.	Load the packages at the top of the script, if you don’t have these already they will need to be installed first
4.	If different raster data is being used to the files on the Github page, some lines of code will need to be updated prior to running the master script. 
5.	If water persistence is being defined from a different raster to the one we used, or you are modifying the definition (e.g. 80% or more inundation, clumps of 5 or more pixels), then you will need to edit the ‘Persistent waterbodies’ script and run it individually prior to running the whole analysis via the master script.  After running the ‘Persistent waterbodies’ script, update the number of waterbodies on the master script (lines 71 & 74). Also update lines 104, 114, 117, 120 of the ‘intersecting polygon area’ script with the new number of waterbodies. 
6.	If different month rasters are being used, it is wise to check them first to make sure they contain good data (clear and contain water). Use the ‘check raster values for monthly tifs master script’ and ‘check water presence for monthly tifs master script’ for this. The connectivity analysis takes a long time (around 12-24hrs) to run, so it’s better to check these first if inputting new data.
7.	Now you are ready to run the code in the master file only (‘Darling River connectivity master script’) to complete the connectivity analysis. It will firstly create the persistent waterbodies shapefile layer from the water persistence raster and then reorder the waterbodies so that their ID’s are in ascending order in a downstream direction beginning at the most upstream waterbody. It will then draw on the ‘intersecting polygon area’ script to fill out the adjacency matrix with data on whether or not each waterbody is connected to its immediate neighbours each month and analyse this data with graph network metrics.

Please email me kaitlyn.omara@griffith.edu.au for access to the raster files, they are too large to upload to Github
