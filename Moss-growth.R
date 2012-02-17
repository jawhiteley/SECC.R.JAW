##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Moss Growth / Productivity
### Jonathan Whiteley     R v2.12     2012-02-05
##################################################
## INITIALISE
##################################################
## Working Directory: see lib/init.R below
if (FALSE) {  # do not run automatically
  setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")  # iMac@McGill
  setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("./")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}


## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')
library(car)

## Moss growth is a little different than other measure variables:
## - Measured from the same t3 & t4 patches throughout.
##   - True repeated measures, unlike destructive sampling
##   - t3 data not included in SECC data frame
## - Unbalanced data (missing values)
##   - this causes serious problems for model.tables,
##     which is used to calculate the MSD error bars :(
##   - A much older version of the analysis using LSD bars seemed to work, however :/
## All told, it might be easier to try a more sophisticated approach
##  using lme/lmer with a nested error structure
##  and multcomp for multiple comparisons.
## If it works, I might even be able to replace the nested ANOVA script
##  with this approach...


##================================================
## CUSTOM CALCULATIONS 
##================================================
## using full moss growth data (including t3 & t4)
attr(SECC.moss.full, "labels")[["Block"]] <- attr(SECC, "labels")[["Block"]]
attr(SECC.moss.full, "labels")[["Time"]] <- attr(SECC, "labels")[["Time"]]
attr(SECC.moss.full, "labels")[["Chamber"]] <- attr(SECC, "labels")[["Chamber"]]
attr(SECC.moss.full, "labels")[["Frag"]] <- attr(SECC, "labels")[["Frag"]]
attr(SECC.moss.full, "labels")[["Position"]] <- attr(SECC, "labels")[["Position"]]

SECC <- within( SECC.moss.full, 
               {
                 grow02 <- grow01 + grow12
                 grow03 <- grow01 + grow12 + grow23
                 Prod02 <- Prod01 + Prod12
                 Prod03 <- Prod01 + Prod12 + Prod23
               }
)

attr(SECC, "labels")[["grow02"]] <- "Moss growth 18 months"
attr(SECC, "labels")[["grow03"]] <- "Moss growth 24 months"
attr(SECC, "units" )[["grow02"]] <- attr(SECC, "units" )[["grow01"]]
attr(SECC, "units" )[["grow03"]] <- attr(SECC, "units" )[["grow01"]]
attr(SECC, "labels")[["Prod02"]] <- "Moss Productivity 18 months"
attr(SECC, "labels")[["Prod03"]] <- "Moss Productivity 24 months"
attr(SECC, "units" )[["Prod02"]] <- attr(SECC, "units" )[["Prod01"]]
attr(SECC, "units" )[["Prod03"]] <- attr(SECC, "units" )[["Prod01"]]

## SECC <- checkSECCdata(SECC)            # too many non-standard IDs
SECC <- recodeSECC(SECC)               # drops t3
SECC <- within( SECC, 
               {
                 Time <- factor(Time, 
                                levels = c(3, levels(Time)[3]),
                                labels = c("t3", "t4")
                                )
                 Time[is.na(Time)] <- levels(Time)[1]
                 Position <- recode(Pos, 
                                    "1='Inner'; 0='Outer'; 'c'='corridor'; 'c2'='corridor'; else='other'", 
                                    levels=c( "Inner", "other", "Outer", "corridor"),
                                    as.factor.result=TRUE
                                    )
               }
)



##################################################
## CONFIGURE BASIC ANALYSIS
##################################################
## Moss growth was measured in the same t4 plots
##  throughout the experiment, so unlike destructive samples,
##  these are repeated measures.
## The simple approach is to analyze each time point separately,
##  which avoids pseudo-replication, but also does not measure the effect of Time.

### Response Variable *****
Ycols <- c('Prod01', 'Prod12', 'Prod23')
Y.col <- 'grow01'     # Column to analyze as response variable           *****
Y.use <- 'Y'    # Which transformation is being used (for labels)? *****

### Load default settings (based on response variable) *****
source("./SECCanova/SECC - ANOVA settings.R", echo = FALSE) 

##================================================
## CUSTOM SETTINGS 
##================================================
## delete lines to use the defaults.

## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)               # samples collected from t3 AND t4 (t3 not included in SECC)
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include
Position.use <- levels(SECC$Position)[c(1, 3)]   # Chamber treatments to include

## Define Labels
Y.units <- bquote( .(Y.units) )     # store as quote(expression)  *****

## Output Results?
Save.results  <- FALSE


### Load default Labels - dependent on above settings. *****
source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) 

##================================================
## CUSTOM LABELS
##================================================



##################################################
### RUN nested model
##################################################
if (T) {                               # Nested ANOVA won't work as is with unbalanced data 
  source("./SECCanova/SECC - nested lm.R", echo = FALSE)
} else {
  source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE)
}






###===============================================
### Include Time as a factor in nested ANOVA
###===============================================
## Repeated measures of the same samples (patches) across time...

##################################################
### PUBLICATION GRAPHS
##################################################



if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 6, height = 3, scale = 1.5)
}
