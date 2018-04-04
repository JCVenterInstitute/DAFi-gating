#!/usr/bin/env Rscript
#####################################################################################################
# Run FCSTrans
# Note:  ipconvert([FCS file directory path without trailing "/"])
#####################################################################################################
args=commandArgs(TRUE)
print(args[1])
library(methods)
library(utils)
source("~/bin/FCSTrans2TXT.R")
ipconvert(args[1], verbose=T)
