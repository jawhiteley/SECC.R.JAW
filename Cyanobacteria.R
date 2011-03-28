##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Cyanobacteria counts  @ time #s
### Jonathan Whiteley     R v2.12     2011-03-28
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
## CONFIGURE BASIC ANALYSIS
##################################################

### Response Variable *****
Y.col <- 'Cells.m'     # Column to analyze as response variable           *****
Y.use <- 'Y.4rt'       # Which transformation is being used (for labels)? ****

### Load default settings (based on response variable) *****
source("./SECCanova/SECC - ANOVA settings.R", echo = FALSE) 

##================================================
## CUSTOM SETTINGS 
##================================================
## delete lines to use the defaults.

## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)[1]            # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include
Frag.use     <- levels(SECC$Frag)               # Frag treatments to include
Position.use <- levels(SECC$Position)[c(1, 3)]  # Patch Positions to include

## Define Labels
Y.units <- bquote( sqrt(.(Y.units), 4) )  # store as quote(expression) *****

## Output Results?
Save.results  <- TRUE  


##================================================
## CUSTOM CALCULATIONS 
##================================================
SECC.prime <- SECC    # save a copy of the original for reference.

## str(SECC)
sampleA  <- 6   # sample Area, in cm^2:  pi * (2.75/2)^2 ; pi * (2.8 / 2)^2
      #     6 for rough estimate of inner tube diameter (2.8 cm): pi*(2.8/2)^2,
      #  or 6.4 for 20 shoots, based on density survey.
sample.to.m2 <- (100*100)/sampleA   # scale sample area, in cm^2 to m^2
sample_ml    <- 50  # 50 ml sample
ARA.m2   <- sampleA/(100*100)  # ARA sample area,   in (cm^2 to) m^2
patchA   <- pi * (12.5^2)      # patch area
patch.m2 <- patchA/(100*100)   # patch sample area, in (cm^2 to) m^2

SECC <- within( SECC, { 
    ## Calculations
})


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

