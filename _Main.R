##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Main control Script
# Jonathan Whiteley		R v2.12		2011-01-26
##################################################
## INITIALIZE
rm(list=ls())  # clear memory
setwd("./ SECC/")	# Project Directory.

#source("./lib/load.R")  # (re-)load data
#source("./lib/init.R")  # Initialize - all analysis scripts should start with this.

## GO
source("./SECC-Data-Template.R")    # Generate Template Data Frame & files.
source("./do.R")    # Perform Analyses

## OUTPUT
source("./lib/out.R")   # Produce Outputs: graphs, reports, export.
