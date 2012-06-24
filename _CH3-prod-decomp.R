################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Main control Script - Chapter 3:
### Decomposition & Productivity in Moss layer
### 
### Jonathan Whiteley		R v2.12		2012-06-22
################################################################
## INITIALIZE
if (FALSE) {  # do not run automatically
  ## Set Working Directory: path in quotes "".
  setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")  # iMac@McGill
  setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("./")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place

  source("./lib/load.R")  # (re-)load data
} else {
  source("./lib/init.R")  # Initialize - all analysis scripts should start with this.
}

## DO
##==============================================================
## Separate Univariate analyses
##==============================================================
source("Moss-growth.R")                # Moss Growth
source("Decomposition.R")              # Decomposition: litter bags


##==============================================================
## Balance
##==============================================================
source("Moss-Prod-Decomp.R")           # Prod - Decomposition mass balance



