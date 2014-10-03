### This file is sourced by or symbol linked to ~/.Rprofile
sourceDir <- function(path = ".") {
        for (file in list.files(path, pattern = "\\.[Rr]$")) {
                source(file.path(path,file))
        }
}
sourceDir("~/SparkleShare/Rprofile/R/RprofilesWin")
sourceDir("~/SparkleShare/Rprofile/R/RprofilesAuto")
sourceDir("~/SparkleShare/phd/R")