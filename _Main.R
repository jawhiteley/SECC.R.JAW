##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Main control Script
# Jonathan Whiteley		R v2.12		2011-01-26
##################################################
## INITIALIZE
rm(list=ls())  # clear memory
## Project Directory
if (FALSE) {  # do not run automatically
  ## Set Working Directory: path in quotes "".
  setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")  # iMac@McGill
  setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
  setwd("./ SECC/") # relative to my usual default wd in R GUI (Mac).
  getwd()           # Check that we're in the right place
}


#source("./lib/load.R")  # (re-)load data
#source("./lib/init.R")  # Initialize - all analysis scripts should start with this.

## DO
source("./SECC-Data-Template.R")    # Generate Template Data Frame & files.
## source("do.R")                      # Perform Analyses

source("ARA.R")             # Acetylene Reduction Assay (N-fixation) Analysis
source("Cyanobacteria.R")   # Cyanobcateria density Analysis

source("H2O.R")             # Cyanobcateria density Analysis

## OUTPUT
source("lib/out.R")   # Produce Outputs: graphs, reports, export.
