#!/usr/bin/env Rscript
#####################################################################################################
# Run FCSTrans
# Note:  ipconvert([FCS file directory path without trailing "/"])
#####################################################################################################
args=commandArgs(TRUE)
print(args[1])
library(methods)
library(utils)
source("/var/DAFi-gating/FCSTrans/FCSTrans2TXT.R")
#ipconvert(args[1], verbose=T)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  print("No compensation file specified")
  ipconvert(args[1], verbose=T)
} else if (length(args)==2) {
  print("Compensation file specified")
  print(args[2])
  ipconvert(args[1], args[2], verbose=T)
}
