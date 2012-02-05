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
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("./")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}


#source("./lib/load.R")  # (re-)load data
#source("./lib/init.R")  # Initialize - all analysis scripts should start with this.

## DO
source("./SECC-Data-Template.R")       # Generate Template Data Frame & files.
## source("do.R")                      # Perform Analyses

if (TRUE)
{
  source("_CH1-ARA-cb.R")              # Chapter 1 analyses: ARA ~ Cyanobacteria
  source("_CH2-Fauna.R")               # Chapter 2 analyses: Fauna community analysis
  source("_CH3-Synthesis.R")           # Chapter 3 analyses: ARA ~ Everything else

} else {
  source("H2O.R")                      # Water contents density Analysis

  source("ARA.R")                      # Acetylene Reduction Assay (N-fixation) Analysis
  source("Cyanobacteria.R")            # Cyanobcateria density Analysis

  source("N-available.R")              # Available Nitrogen (NH4 + NO3)
  source("N-available-time.R")         # Available Nitrogen, without duration (time) of exposure
}


## OUTPUT
source("lib/out.R")                    # Produce Outputs: graphs, reports, export.
