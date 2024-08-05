index <- index + 1

test <- raster(filelist[index])
sum.val <- sum(values(test), na.rm=TRUE)

raster.value.data[index, ] <- c(sum.val)
