#Merge predictions for all regions for each indicator

#Load libraries
library(raster); library(maptools)
library(gtools); library(sp); library(spdep)
library(rgdal)
library(ggplot2)

#File paths
filePathData <- c("/mainfs/scratch/ceu1c14/Data_quality_analysis/West_Africa/",
      "/mainfs/scratch/ceu1c14/Data_quality_analysis/Southern_Africa/",
      "/mainfs/scratch/ceu1c14/Data_quality_analysis/East_Africa/",
      "/mainfs/scratch/ceu1c14/Data_quality_analysis/Central_Africa/")

filePathData2 <- "/mainfs/scratch/ceu1c14/Data_quality_analysis/All_Africa/"

#Region 
Reg <- c("WA", "SA", "EA", "CA")

indicator <- c("H1", "H4", "W0", "W1", "W2", "K1", "K6")

########
for (jj in indicator){

#mean
xx <- list()
for (i in 1:length(Reg)){
  hh <- raster(paste0(filePathData[i], Reg[i],"_", jj, "_mean.tif"))
  xx[[i]] <- hh
}
xx$tolerance <- 1
xx$filename <- paste0(filePathData2, jj, "_SSA_mean.tif")
xx$overwrite <- TRUE
mm <- do.call(merge, xx)

#sd
xx <- list()
for (i in 1:length(Reg)){
  hh <- raster(paste0(filePathData[i], Reg[i],"_", jj, "_sd.tif"))
  xx[[i]] <- hh
}
xx$tolerance <- 1
xx$filename <- paste0(filePathData2, jj, "_SSA_sd.tif")
xx$overwrite <- TRUE
mm <- do.call(merge, xx)

#low
xx <- list()
for (i in 1:length(Reg)){
  hh <- raster(paste0(filePathData[i], Reg[i], "_", jj, "_low.tif"))
  xx[[i]] <- hh
}
xx$tolerance <- 1
xx$filename <- paste0(filePathData2, jj, "_SSA_low.tif")
xx$overwrite <- TRUE
mm <- do.call(merge, xx)

#up
xx <- list()
for (i in 1:length(Reg)){
  hh <- raster(paste0(filePathData[i], Reg[i], "_", jj, "_up.tif"))
  xx[[i]] <- hh
}
xx$tolerance <- 1
xx$filename <- paste0(filePathData2, jj, "_SSA_up.tif")
xx$overwrite <- TRUE
mm <- do.call(merge, xx)

#median
xx <- list()
for (i in 1:length(Reg)){
  hh <- raster(paste0(filePathData[i], Reg[i], "_", jj, "_median.tif"))
  xx[[i]] <- hh
}
xx$tolerance <- 1
xx$filename <- paste0(filePathData2, jj, "_SSA_median.tif")
xx$overwrite <- TRUE
mm <- do.call(merge, xx)

}

Reg <- c("WA", "SA", "EA", "CA")

for (jj in indicator){

d1 <- read.csv(paste0(filePathData[1], "WA", "_", jj, "_district_estimates.csv"), header = TRUE)
d2 <- read.csv(paste0(filePathData[2], "SA", "_", jj, "_district_estimates.csv"), header = TRUE)
d3 <- read.csv(paste0(filePathData[3], "EA", "_", jj, "_district_estimates.csv"), header = TRUE)
d4 <- read.csv(paste0(filePathData[4], "CA", "_", jj, "_district_estimates.csv"), header = TRUE)

dd <- rbind(d1,d2,d3,d4)

write.csv(dd, paste0(filePathData2, jj, "_district_estimates.csv"))
}












