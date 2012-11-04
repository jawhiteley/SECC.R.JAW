##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Main control Script - Chapter 2:
### Fauna: Mesostigmatid mites, predatory prostigs (?), Collembola (Microarthropods)
### Jonathan Whiteley		R v2.12		2012-02-05
##################################################
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

SECC.prime <- SECC    # save a copy of the original for reference.

## DO
## source("./SECC-Data-Template.R")    # Generate Template Data Frame & files.

source("Fauna/Fauna.R")                      # Multivariate analyses of community structure
## source("Fauna-univariate.R")           # Nested ANOVA of biodiversity indices (included in above)

