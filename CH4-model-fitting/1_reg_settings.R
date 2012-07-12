################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Initialize, configure default analysis options
### N-fixation synthesis vs. everything else, really
### Jonathan Whiteley     R v2.12     2012-07-12
################################################################
## INITIALISE
################################################################
## Working Directory: see lib/init.R below [\rd in Vim]
if (FALSE) {  # do not run automatically
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("..")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')

## I'm taking a very different approach to the code organization 
##  for this part of the project than I did for ARA~Cyanobacteria
## This time, I'm using a generic settings script (this one)
##  for consistent defaults + a script for labels & data processing (reg_setup)
## Each modelling script (for each response variable) can start by sourcing these, 
##  overriding the defaults, then proceeding individually.
## The only reason to separate the setup for each response variable
##  is to recycle it for different types of analyses (regression trees, mixed models, etc.)
################################################################
## CONFIGURE BASIC ANALYSIS
################################################################
if (!exists('Y.col'))  Y.col  <- 'Nfix'      # Column to analyze as response variable           *****
# explanatory vars for data exploration (and labels)
if (!exists('X.cols')) X.cols <- c("Cells.m", "Hcells.m", "H2O")  

##==============================================================
## SETTINGS 
##==============================================================
## Specify which treatment levels to include (by index is probably easiest)
## I only really have complete data for:
## - t4 patches
## - Ambient & Full Chamber
## - Inner & Outer patch Positions
## - I might also have to drop pseudo-corridor treatments (due to missing fauna data)
Time.use     <- levels(SECC$Time)[3]           # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]  # Chamber treatments to include
Frag.use     <- levels(SECC$Frag)              # Frag treatments to include
Position.use <- levels(SECC$Position)[c(1, 3)] # Patch Positions to include

Save.results  <- TRUE                  # Output Results?

