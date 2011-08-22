##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Aggregate Fauna data (microarthropod morphospecies counts)
### Jonathan Whiteley     R v2.12     2011-08-22
###===============================================
### Species identified to morphospecies (usually family-level)
### Counts / sample converted to # / g dwt of moss substrate (using 'Patch.dwt' column)
##################################################
## INITIALISE
##################################################
## This script is used in a generic way for most univariate analyses
## Working Directory: see lib/init.R below
if (FALSE) {  # do not run automatically
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (Mac).
  getwd()  # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')


##################################################
## DATA EXPLORATION
##################################################

plot(SECC.fauna.sum[, which( sapply(SECC.fauna.sum, is.numeric) )])







##################################################
## CONFIGURE BASIC ANALYSIS
##################################################
## Assumes analysis on 'SECC' data frame, 
## but relevant fauna data is stored in 'SECC.fauna.sum'
SECC.df <- SECC # temporary storage for this script (restored at the end)
SECC <- SECC.fauna.sum

### ANOVA Response Variable *****
Y.col <- ''     # Column to analyze as response variable           *****
Y.use <- 'Y'    # Which transformation is being used (for labels)? *****

##================================================
## CUSTOM CALCULATIONS 
##================================================



### Load default settings (based on response variable) *****
source("./SECCanova/SECC - ANOVA settings.R", echo = FALSE) 

##================================================
## CUSTOM SETTINGS 
##================================================
## delete lines to use the defaults.

## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)[3]      # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include: A, C
Frag.use     <- levels(SECC$Frag)[c(1, 2, 4)]   # Fragmentation treatments to include: 1, 2, 4
Position.use <- levels(SECC$Position)[c(1, 3)]  # Patch Positions to include: Inner, Outer

## Define Labels
Y.units <- bquote( .(Y.units) )     # store as quote(expression)  *****

## Output Results?
Save.results  <- TRUE


### Load default Labels - dependent on above settings. *****
source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) 

##================================================
## CUSTOM LABELS
##================================================



##################################################
### RUN STANDARD nested ANOVA
##################################################

### Run analysis on each Time point in sequence.
for ( Time.i in 1:length(levels(SECC$Time)) ) {
  ## Specify which treatment levels to include (by index is probably easiest)
  Time.use     <- levels(SECC$Time)[Time.i]      # Time (index: 1-3) to include in this run
  cat("\n\n\nProcessing Time:", Time.use, "\n")

  ## Load default Labels - dependent on above settings. *****
  source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) 

  ## RUN STANDARD nested ANOVA
  source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE)

}


###===============================================
### Include Time as a factor in nested ANOVA
###===============================================
## Note that Samples at different times are actually independent
## in this design, due to destructive sampling.

Time.use     <- levels(SECC$Time)      # Include *ALL* Times (as a Treatment)
source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) # Load default Labels. *****
source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE) # RUN STANDARD nested ANOVA



##################################################
### PUBLICATION GRAPHS
##################################################



##################################################
### CLEAN-UP / HOUSEKEEPING
##################################################
SECC <- SECC.df # restore original
