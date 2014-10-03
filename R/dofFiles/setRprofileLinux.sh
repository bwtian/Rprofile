#!/bin/sh
now=$(date +%Y-%m%d-%H%M)
mv ~/.Rprofile ~/.Rprofile.$now
\ln -sfv ~/SparkleShare/Rprofile/R/dotFiles/R00_Linux_dot_Rprofile.R  ~/.Rprofile
