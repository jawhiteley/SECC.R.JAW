##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Main control Script - Chapter 4:
### Ecosystem-level synthesis of all measured variables
### N-fixation, cyanobacteria, Fauna, Moisture, Moss Growth, Decomposition, Available N, etc.
### Jonathan Whiteley		R v2.12		2012-07-14
##################################################
## INITIALIZE
if (FALSE) {  ## Project Directory
  ## Set Working Directory: path in quotes "".
  setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")  # iMac@McGill
  setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("../")      # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place

  source("./lib/load.R")  # (re-)load data
} else {
  source("./lib/init.R")  # Initialize - all analysis scripts should start with this.
}

SECC.prime <- SECC    # save a copy of the original for reference.
RegDir <- "./CH4-model-fitting/"       # directory for scripts used in these analyses (this folder)
## I could just setwd() to this folder, but other functions or scripts for the larger projects still depend on having the working directory in the root folder for this project, so it's probably safer to avoid.

## DO
## _setup.R scripts also include source("./lib/init.R"); memore is cleared, data re-loaded & processed every time.
cat("N-fixation ~ \n")
source("CH4-model-fitting/Nfix_setup.R")    # N-fixation: settings
cat("N-fixation ~ : Data Exploration\n")
source("CH4-model-fitting/2_reg_Explore.R") # Data Exploration
source("CH4-model-fitting/Nfix-trees.R")    # N-fixation: regression trees
source("CH4-model-fitting/Nfix-models.R")   # N-fixation: regression models

cat("Moss Growth ~ \n")
source("CH4-model-fitting/Growth_setup.R")  # Moss Growth: settings
source("CH4-model-fitting/2_reg_Explore.R") # Data Exploration
source("CH4-model-fitting/Growth-trees.R")  # Moss Growth: regression trees
source("CH4-model-fitting/Growth-Nfix-models.R") # Moss Growth: regression models

if (FALSE)
{
  cat("Decomposition ~ \n")
  source("CH4-model-fitting/Decomp_setup.R")  # Decomposition: settings
  source("CH4-model-fitting/2_reg_Explore.R") # Data Exploration
  source("CH4-model-fitting/Decomp-trees.R")  # Decomposition: regression trees
  source("CH4-model-fitting/Decomp-models.R") # Decomposition: regression models
}


