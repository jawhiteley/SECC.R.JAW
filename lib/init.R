### INITIALIZE
rm(list=ls())       # house-keeping

## Load Functions and Libraries
cat('Loading functions.\n')
source("./lib/fun.R")   # define functions

## Load Data
cat('Loading data.\n')
if (FALSE) {
  ## do not run when source()'d
  source("./lib/load.R")  # (re-)load Data
} else {
  ## OR if output file is produced instead:
  load("./save/SECC_data.R") # load data
}
