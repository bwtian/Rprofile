#!/bin/sh
now=$(date +%Y-%m%d-%H%M)
mv ~/.Rprofile ~/.Rprofile.$now
\ln -sfv ~/SparkleShare/Rprofile/R/dofFiles/R00_Linux_dot_Rprofile.R  ~/.Rprofile
