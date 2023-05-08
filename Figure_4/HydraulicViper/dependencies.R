# Depedencies ----------------------
library(readxl) #
library(Matrix)
library(data.table)
library(plyr) 
library(dplyr) #
library(stringr) #
library(ggplot2)
library(purrr)
library(xml2) #
`%!in%` <- compose(`!`, `%in%`)
library(deldir) #
library(alphahull) #
library(sp)

# granar mecha short-cut : random forest model
#source("/home/users/a/h/aheymans/GSUA/R/GRANAR/rf/machine_learning.r")
#PkgTest(c("randomForest"))
# next line is not working ...
# how to load .RData
library(caret)
