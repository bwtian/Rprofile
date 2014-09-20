
sourceDir <- function(path = ".") {
         for (file in list.files(path, pattern = "\\.[Rr]$")) {
         source(file.path(path,file))
         }
}
