### INITIALIZE
rm(list=ls())       # house-keeping
cat('Memory Cleared.\n')

## Working Directory
if (FALSE) {  # do not run automatically
  ## Set Working Directory: path in quotes "".
  setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")  # iMac@McGill
  setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (Mac).
  getwd()  # Check that we're in the right place
}

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

cat('== SECC Project Initialized ==\n')
