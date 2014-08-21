##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Template for basic analyses of experimental data
### Response Variable(s)  @ time #s
### Jonathan Whiteley     R v2.12     2011-03-25
##################################################
## INITIALISE
##################################################

par(ask = FALSE)    # Stop asking me to hit <Return> to see next plot!
## options(device.ask.default = FALSE) # same as above: does this work?


##================================================
## SETTINGS
##================================================

### Response Variable *****
## This had better be set first in the parent script.  Only set them here if necessary.
if (!exists('Y.col')) Y.col <- 'Nfix' # Column to analyze as response variable           *****
if (!exists('Y.use')) Y.use <- 'Y'    # Which transformation is being used (for labels)? ****


## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)         # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)      # Chamber treatments to include
Frag.use     <- levels(SECC$Frag)         # Frag treatments to include
Position.use <- levels(SECC$Position)     # Patch Positions to include


##================================================
## CALCULATIONS
##================================================
## !is.na(SECC$Time) ; NAs in factors are annoying
SECC.prime <- SECC    # save a copy of the original for reference.

if (FALSE) {

## str(SECC)
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
  ## change negative ARA values to 0
  ARA.ml[ARA.ml < 0] <- 0
  ARA.m[ ARA.m  < 0] <- 0
  ARA.g[ ARA.g  < 0] <- 0
  Nfix <- ARA.m * Nfix.ARA.ratio
})

}

# options(digits=6)                      # level of detail when printing numbers?

##================================================
## Define Labels
##================================================

Y.label <- attr(SECC, "labels")[[Y.col]]  # response variable label for this script.
Y.units <- attr(SECC, "units" )[[Y.col]]  # response variable units: quote(expression).
Y.units <- bquote( .(Y.units) )  # sqrt(.(Y.units), 4)  # store as quote(expression)

## Output Results?
## Logical switch determines whether output is saved to files, or left in R.  Easier than setting several values to NULL
Save.results <- FALSE  


