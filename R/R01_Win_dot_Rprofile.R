### This file is sourced by or symbol linked to ~/.Rprofile
sourceDir <- function(path = ".") {
        for (file in list.files(path, pattern = "\\.[Rr]$")) {
                source(file.path(path,file))
        }
}
sourceDir("C:/Air/Dropbox/config/R/rProfile/RprofilesAuto")
## sourceDir("~/Dropbox/config/R/rProfile/RprofilesWin")
