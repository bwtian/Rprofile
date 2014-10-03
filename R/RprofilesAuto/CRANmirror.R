
## Don't ask me for my CRAN mirror every time
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.ism.ac.jp/"
       options(repos=r)})
