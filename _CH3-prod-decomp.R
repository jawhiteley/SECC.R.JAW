################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Main control Script - Chapter 3:
### Ecosystem-level synthesis of all measured variables
### N-fixation, cyanobacteria, Fauna, Moisture, Available N, Decomposition, etc.
### Jonathan Whiteley		R v2.12		2012-02-05
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




################################################################
## PRODUCTION - DECOMPOSITION MASS BALANCE
################################################################
## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')

##==============================================================
## Convert Production and Decomposition to common units (t4 only)
##==============================================================
SECC <- within(SECC, 
               {
                 Productivity <- Prod12 + Prod23
                 if (FALSE)
                 {
                 Productivity <- Productivity * 385.8 # shoots / patch (appoximately, on average)
                 } else {
                 Productivity <- Productivity * Patch.dwt / (Cells.dwt / 1000 * 2) # multiply by fraction of patch weight in ~1 shoot (half of sample used for cyanobacteria Cells)
                 }
                 Productivity <- Productivity / 1000 # mg -> g
                 Decompositng <- Decomposition * Patch.dwt/3 # only dead tissue is really decomposing?
                 PD.bal <- Productivity - Decompositng
               }
)

##==============================================================
## Data Exploration & Checking
##==============================================================
if (FALSE)
{
  with(SECC, hist(Productivity) )
  with(SECC, hist(Decompositng) )
  with(SECC, hist(PD.bal) )
}

SECCsub <- subset(SECC, Position %in% c("Inner", "Outer") ) 
SECCsub <- subset(SECCsub, Chamber %in% c("Ambient", "Full Chamber") ) 
SECCsub$Position <- factor(SECCsub$Position)
SECCsub$Chamber  <- factor(SECCsub$Chamber)

library(ggplot2)
PD.plot <- ggplot(SECCsub, 
                  aes(x = Chamber, y = PD.bal, colour = Position)
) +
               stat_summary(fun.data = "mean_cl_boot", geom = "errorbar" ) +
               stat_summary(fun.y = mean, geom = "line" )

print(PD.plot)
