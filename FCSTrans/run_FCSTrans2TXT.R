#!/usr/bin/env Rscript
#####################################################################################################
# Run FCSTrans
# Note:  ipconvert([FCS file directory path without trailing "/"])
#####################################################################################################
args=commandArgs(TRUE)
print(args[1])
source("./FCSTransTXT.R")
ipconvert(args[1], verbose=T)
