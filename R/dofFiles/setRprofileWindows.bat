echo on
REM Step1.make a environment variable HOME=%USERPROFILE%
REM Step2.make a environment variable R_HOME=C:\Program Files\R\R-3.1.1
REM Step3. mklink "~/.Rprofile" "~/SparkleShare/Rprofile/R/dotFiles/R01_Win_dot_Rprofile.R"
REM Step4.copy the "~/SparkleShare/Rprofile/R/dotFiles/.Rprofile" to the ~/Rprofile
wmic ENVIRONMENT create name="HOME",username="<system>",VariableValue="%USERPROFILE%"
wmic ENVIRONMENT create name="R_HOME",username="<system>",VariableValue="C:\Program Files\R\R-3.2.5"
mklink "~/.Rprofile" "~/SparkleShare/Rprofile/R/dotFiles/R01_Win_dot_Rprofile.R"