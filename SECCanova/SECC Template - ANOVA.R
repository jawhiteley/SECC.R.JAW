##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Template for basic analyses of experimental data
### Response Variable(s)  @ time #s
### Jonathan Whiteley     R v2.12     2011-03-25
##################################################
## INITIALISE
##################################################
## This script is used in a generic way for most univariate analyses

## Set Working Directory: path in quotes "".
## setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")    # iMac@McGill
## setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
## setwd("./ SECC/")  # relative to my usual default wd in R GUI (Mac).
getwd()  # Check that we're in the right place

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')


##################################################
## CONFIGURE BASIC ANALYSIS
##################################################
### Load default settings *****
source("SECC - ANOVA settings.R") 

##================================================
## SETTINGS - edit
##================================================
## delete lines to use the defaults.

### Response Variable *****
Y.col <- 'Nfix'      # Column to analyze as response variable           *****
Y.use <- 'Y.sqrt'    # Which transformation is being used (for labels)? ****

## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)[1]      # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)      # Chamber treatments to include
Frag.use     <- levels(SECC$Frag)         # Frag treatments to include
Position.use <- levels(SECC$Position)     # Patch Positions to include

## Define Labels
Y.units <- bquote( .(Y.units) )  # sqrt(.(Y.units), 4)  # store as quote(expression)

## Output Results?
## Logical switch determines whether output is saved to files, or left in R.  Easier than setting several values to NULL
Out.results  <- TRUE  


### Load default Labels - dependent on above settings. *****
source("SECC - ANOVA labels.R") 


##================================================
## CALCULATIONS - edit
##================================================
## !is.na(SECC$Time) ; NAs in factors are annoying
SECC.prime <- SECC    # save a copy of the original for reference.

# str(SECC.use)
sampleA  <- 6   # sample Area, in cm^2:  pi * (2.75/2)^2 ; pi * (2.8 / 2)^2
      #     6 for rough estimate of inner tube diameter (2.8 cm): pi*(2.8/2)^2,
      #  or 6.4 for 20 shoots, based on density survey.
sample.to.m2 <- (100*100)/sampleA   # scale sample area, in cm^2 to m^2
sample_ml    <- 50  # 50 ml sample
ARA.m2   <- sampleA/(100*100)  # ARA sample area,   in (cm^2 to) m^2
patchA   <- pi * (12.5^2)      # patch area
patch.m2 <- patchA/(100*100)   # patch sample area, in (cm^2 to) m^2
Nfix.ARA.ratio <- 1/3  # ratio of N-fixation : ARA.

SECC <- within( SECC, { 
  ## change negative ARA values to 0 - should I wait until after aggregation?
  ARA.ml[ARA.ml < 0] <- 0
  ARA.m[ ARA.m  < 0] <- 0
  ARA.g[ ARA.g  < 0] <- 0
  Nfix <- ARA.m * Nfix.ARA.ratio
})



##################################################
### RUN STANDARD ANALYSIS
##################################################

source("SECC - ANOVA.R")

