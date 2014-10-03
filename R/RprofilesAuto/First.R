
## Don't ask me for my CRAN mirror every time
    .First <- function() {
            cat("Current R.version is:", as.character(getRversion()), "\n")
            cat("Current R.repos is:", getOption("repos"), "\n")
            cat("Current R.Library is:", .libPaths(), "\n")
            cat("Current R.Workdir is:", getwd(), "\n")
            cat("\nSuccessfully loaded .Rprofile at:", date(), "\n")
            if (file.exists(".RData")) load(file=".RData")
}
