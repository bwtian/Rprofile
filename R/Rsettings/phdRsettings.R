#' @usage source("~/SparkleShare/Rprofile/R/Rsettings/phdRsettings.R")
#' @author Bingwei Tian <bwtian@gmail.com>
#'
library(sp)
library(gstat)
library(rgdal)
library(gdalUtils)
library(maptools)
library(raster)
library(plotKML)
library(plyr)
library(rgeos)

### Graph
library(ggplot2)
library(ggmap)
library(rasterVis)
library(grid)
library(gridExtra)
library(lattice)
library(latticeExtra)

### packages
source("~/SparkleShare/Rprofile/R/RprofilesAuto/sourceDir.R")
sourceDir("~/SparkleShare/Rprofile/R/")
sourceDir("~/SparkleShare/TIR/R/")
sourceDir("~/SparkleShare/rLandsat8/src/main/R/rLandsat8/R")
sourceDir("~/SparkleShare/geothermaR/R/")
##  Options
if(.Platform$OS.type == "windows"){
        windowsFonts(Times=windowsFont("TT Times New Roman"))
        windowsFonts(times=windowsFont("TT Times New Roman"))
}

driver     <- "~/Share500sda/" # Linux and Windows Symbolink
dir.tmp    <- file.path(driver, "raster_tmp")
rasterOptions(tmpdir = dir.tmp)
raster::removeTmpFiles(h = 24)
gc()
### Colors
rainbow1  <- rainbow(n = 255, start = 2/6)
phd.rainbow <- grDevices::colorRampPalette(c("purple","blue","cyan","green","yellow", "orange","red"))
oceColorsJet  <- function (n)
{
        if (missing(n) || n <= 0)
                colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        else {
                colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(n)
        }
}

# Extent of Study area of Japan in Phd Thesis
xlimJP <- c(128.5, 146.5)
ylimJP <- c(30.2, 45.8)
certerJp <- c(137.5, 38)

# Coordinate Reference Systems of Japan
## Geographic Coordinate Reference Systems
#tokyoGRS <- "+init=epsg:4301"
tokyoGRS  <- "+proj=longlat +ellps=bessel +no_defs"
#wgs84GRS <- "+init=epsg:4326"
wgs84GRS  <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
#jgd2000GRS <- "+init=epsg:4612"
jgd2000GRS  <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"

## Projected Coordinate Reference Systems
lccBessel <- "+proj=lcc +lat_1=32.8 +lat_2=43.2 +lat_0=38 +lon_0=137.5 +x_0=1000000 +y_0#=1000000 +ellps=bessel +towgs84=-146.414,507.337,680.507,0,0,0,0 +units=m +no_defs"
## change the datum to WGS84 20140508
lccWgs84 <- "+proj=lcc +lat_1=32.8 +lat_2=43.2 +lat_0=38 +lon_0=137.5 +x_0=1000000 +y_0=1000000 +datum=WGS84 +units=m +no_defs"
#  +ellps=WGS84 +towgs84=0,0,0
codeDir <- "~/Dropbox/1code"
dataDir <- "~/Dropbox/2data"
dataData <- "~/Dropbox/2data/data"
dataRaw <- "~/Dropbox/2data/dataRaw"
dataPro <- "~/Dropbox/2data/dataProduct"
figsDir <- "~/Dropbox/3figs"
