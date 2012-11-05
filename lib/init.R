### Working Directory
if (FALSE) {  # do not run automatically
  ## Set Working Directory: path in quotes "".
  setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")  # iMac@McGill
  setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
  setwd("./ SECC/") # relative to my usual default wd in R GUI (Mac).
  setwd("..")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}

### INITIALIZE
rm(list=ls(all = TRUE))       # house-keeping
options(contrasts=c("contr.treatment", "contr.poly")) # reset options - this affects model-fitting and summary() outputs
cat('Memory Cleared.\n')

## Load Functions and Libraries
cat('Loading functions.\n')
source("./lib/fun.R")   # define functions

## Load Data
cat('Loading data.\n')
if (file.exists("./save/SECC_data.R")) 
{
  load("./save/SECC_data.R") # load data
} else {
  source("./lib/load.R")  # (re-)load Data
}

cat('== SECC Project Initialized ==\n')
