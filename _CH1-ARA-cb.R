##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Main control Script - Chapter 1:
### N-fixation, cyanobacteria, and Moisture.
### Jonathan Whiteley		R v2.12		2011-03-29
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

## source("./lib/init.R")  # Initialize - all analysis scripts should start with this.
source("./lib/load.R")  # (re-)load data
SECC.prime <- SECC    # save a copy of the original for reference.

## DO
source("./SECC-Data-Template.R")    # Generate Template Data Frame & files.

source("H2O.R")             # Water Contents Analysis

source("ARA.R")             # Acetylene Reduction Assay (N-fixation) Analysis
source("Cyanobacteria.R")   # Cyanobcateria density Analysis

source("ARA~cyanobacteria.R")   # ARA ~ cyanobacteria (GLMM)

## OUTPUT
source("lib/out.R")   # Produce Outputs: graphs, reports, export.