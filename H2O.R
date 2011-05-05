##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Moisture Content (H2O % of moss dry wt)
### Jonathan Whiteley     R v2.12     2011-03-29
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
Y.col <- 'H2O'     # Column to analyze as response variable           *****
Y.use <- 'Y'    # Which transformation is being used (for labels)? *****

##================================================
## CUSTOM CALCULATIONS 
##================================================

SECC <- within( SECC, { 
    H2O.asq <- asin(sqrt(H2O.wwt))
    H2O     <- H2O     * 100  # convert to %
    H2O.wwt <- H2O.wwt * 100
})

attr(SECC, "labels")[["H2O.asq"]] <- "Moisture"
attr(SECC, "units" )[["H2O.asq"]] <- quote(asin(sqrt("% "* H[2]*O)))
attr(SECC, "labels")[["H2O"]] <- "Moisture"
attr(SECC, "units" )[["H2O"]] <- quote("("* H[2]*O *" as % of Dry Wt. moss)")


### Load default settings (based on response variable) *****
source("./SECCanova/SECC - ANOVA settings.R", echo = FALSE) 

##================================================
## CUSTOM SETTINGS 
##================================================
## delete lines to use the defaults.

## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)[1]      # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include

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

if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  postscript(file = Save.final, width = 6, height = 2)
}

Plot.Title <- bquote(.(Time.label) * "Patch means " %+-% "95% Comparison Intervals")
Sub.msd <- "95% comparison intervals (MSR)" 

Y.lim <- c(0, 800)
plot.means <- tapply(SECCp$Y.trans, list(SECCp$Chamber, SECCp$Position), mean)
plot.error <- matrix( as.numeric(msd["Chamber:Position"]/2),
                     nrow = length(levels(SECCp$Chamber)),  # rows: x-axis
                     ncol = length(levels(SECCp$Position))  # cols: y-axis
                     )

old.par <- par(mfrow=c(1,3), las = 1, oma = c(3, 2, 3, 1), mar = c(3, 3, 2, 0) +0.1 )

for(Time.for in 1:length(Time.use)) {
with( SECCp[SECCp$Time == Time.use[Time.for], ], {
  plot.means <- tapply(Y.trans, list(Chamber, Position), mean)
  plot.error <- matrix( as.numeric(msd["Chamber:Position"]/2),
                       nrow = length(levels(SECCp$Chamber)),  # rows: x-axis
                       ncol = length(levels(SECCp$Position))  # cols: y-axis
                       )
  ## custom plotMeans function, with custom error bars (LSD)
  plotMeans( Y.trans, Chamber, Position, 
            error.bars = "custom", level = plot.error, ylim = Y.lim,
            cex = 2, lwd = 2, lty = Position.map$lty, pch = Position.map$pch,
            yaxt = ifelse(Time.for > 1, "n", "s"),
            col = as.character(Position.map$col),
            bg  = as.character(Position.map$bg),
            main = Time.use[Time.for],
            sub  = "",
            xlab = "",
            ylab = ""
            )   
  })
}
mtext(Plot.Title, side = 3, outer = TRUE)
mtext(attr(SECC, "labels")[["Chamber"]], side = 1, outer = TRUE)
mtext(Sub.msd, side = 1, padj = 1.5, outer = TRUE)
mtext(Y.plotlab, side = 2, outer = TRUE, las = 0)

par(old.par)
if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()
