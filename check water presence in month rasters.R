index <- index + 1

test <- raster(filelist[index])
newVal <- (test >= 2) #subset to only include water pixels
sum.water.val <- sum(values(newVal))

raster.water.data[index, ] <- c(sum.water.val)
