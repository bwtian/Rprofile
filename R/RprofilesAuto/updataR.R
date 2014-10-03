R_ver  <- as.character(getRversion())
updateR <- function() {
        if(!require(installr)) {
                install.packages("installr")
        } #load / install+load installr
        updateR() # this will only work AFTER R 3.0.0 
        update.packages(checkBuilt=TRUE, ask = FALSE)
}
updatePkgs <- function() {
        update.packages(checkBuilt=TRUE, ask = FALSE, dependencies = c('Suggests'))
}