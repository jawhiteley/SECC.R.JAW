##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Main control Script - Chapter 4:
### Ecosystem-level synthesis of most (all) measured variables
### N-fixation, cyanobacteria, Moss Growth, Fauna, Moisture, Decomposition, Available N, etc.
### Jonathan Whiteley		R v2.12		2012-07-18
##################################################
## INITIALIZE
## Project Directory
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

SECC.prime <- SECC    # save a copy of the original for reference.

## DO
## source("./SECC-Data-Template.R")    # Generate Template Data Frame & files.

source("N-available.R")                # Available Nitrogen (NH4 + NO3): ion resin capsules
source("N-available-time.R")           # Available N, without duration (time) of exposure

source("CH4-model-fitting/0_reg_models.R")    # Regression modelling: Nfixation, Moss Growth

