
## Don't ask me for my CRAN mirror every time
## Japan：http://cran.ism.ac.jp/
local({r <- getOption("repos")
       r["CRAN"] <- "https://mirrors.ustc.edu.cn/CRAN/"
       options(repos=r)})
