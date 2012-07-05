##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Decomposition (mass loss in 1 year)
### Jonathan Whiteley     R v2.12     2012-02-20
##################################################
## INITIALISE
##################################################
## This script is used in a generic way for most univariate analyses
## Working Directory: see lib/init.R below
if (FALSE) {  # do not run automatically
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("./")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')


##################################################
## CONFIGURE BASIC ANALYSIS
##################################################

### Response Variable *****
Y.col <- 'Decomp.asq'     # Column to analyze as response variable           *****
Y.use <- 'Y'    # Which transformation is being used (for labels)? *****

##================================================
## CUSTOM CALCULATIONS 
##================================================

SECC <- within( SECC, { 
    Decomp.asq    <- asin(sqrt(Decomposition))
    Decomposition <- Decomposition * 100  # convert to %
})

attr(SECC, "labels")[["Decomp.asq"]] <- "Decomposition"
attr(SECC, "units" )[["Decomp.asq"]] <- quote(asin(sqrt("% mass loss")))


### Load default settings (based on response variable) *****
source("./SECCanova/SECC - ANOVA settings.R", echo = FALSE) 

##================================================
## CUSTOM SETTINGS 
##================================================
## delete lines to use the defaults.

## Specify which treatment levels to include (by index is probably easiest)
## I only have Decomposition measures from Inner & Outer patches (all chambers).
## There are a few "other" patches, but they're all in Block 5 :/
Time.use     <- levels(SECC$Time)[3]      # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include (all available, but we'll exclude partial for simplicity?)
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
source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE)

if (FALSE) {

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

}

##################################################
### PUBLICATION GRAPHS
##################################################
FragIconList <- list(FragIcon1 = FragIcon1,
                     FragIcon2 = FragIcon2,
                     FragIcon3 = FragIcon3,
                     FragIcon4 = FragIcon4
                     )

decomp.plotcalcs <- function(plot.means)
{
  plot.data <- within(plot.means, 
                      {
                        lower <- x - error
                        upper <- x + error
                        xtr <- sin(x)^2
                        lwr <- sin(lower)^2
                        upr <- sin(upper)^2
                        if (exists("Time")) levels(Time) <- 
                          paste(c("August", "June", "August"), 
                                levels(Time), sep="\n")
                        if (exists("Chamber")) levels(Chamber)[2] <- "Chamber"
                      }
  )
}

Y.lim <- c(0, 0.20)
Y.plotlab <- bquote( .(attr(SECC, "labels")[["Decomposition"]]) * 
                    " (proportion mass loss / year)")
Plot.Title <- bquote("Patch means " %+-% "95% Comparison Intervals")
Sub.msd <- "95% comparison intervals (MSR)" 
Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = levels(SECC$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]

## 3-way interaction plot!!
plot.means <- aggregate(SECCp$Y.trans, list(Chamber=SECCp$Chamber, Position=SECCp$Position, Time=SECCp$Time, Frag=SECCp$Frag), mean)
plot.means$error <- as.numeric(msd["Chamber:Frag:Position"]/2)
plot.data <- decomp.plotcalcs(plot.means)

CxPxF.plot <- qplot(Frag, xtr, data = plot.data, group = Position, 
                    geom = "point", ylim = Y.lim, size = I(3), 
                    colour = Position, shape = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Chamber"]],
                    ylab = Y.plotlab,
                    legend = FALSE,
                    facets = .~ Chamber)
## CxPxF.plot <- CxPxF.plot + geom_point(aes(Chamber, x), size = 2)
CxPxF.plot <- CxPxF.plot + geom_line(aes(group = Position), size = 0.8)
CxPxF.plot <- CxPxF.plot + geom_errorbar(aes(ymin = lwr, ymax = upr), 
                                         width = 0.2, size = 0.5)
CxPxF.plot <- CxPxF.plot + scale_colour_manual(name = Position.label,
                                           values = Position.map$col, 
                                           breaks = Position.map$label)
CxPxF.plot <- CxPxF.plot + scale_fill_manual(name = Position.label,
                                         values = Position.map$bg, 
                                         breaks = Position.map$label)
CxPxF.plot <- CxPxF.plot + scale_shape_manual(name = Position.label,
                                           values = Position.map$pch, 
                                           breaks = Position.map$label)
CxPxF.plot <- CxPxF.plot + jaw.ggplot()
## Add imported graphics as x-axis tick labels :D
## http://stackoverflow.com/questions/2181902/how-to-use-an-image-as-a-point-in-ggplot
CxPxF.plot <- CxPxF.plot + scale_x_discrete(labels = names(FragIconList), # c(1, 2, 3, 4), 
                                            breaks = levels(plot.means$Frag)) +
opts(axis.ticks.margin = unit(0.2, "lines"),
     axis.text.x = picture_axis(FragIconList, icon.size = unit(1.4, "lines")) 
)
print(CxPxF.plot)


## 2-way (significant) interaction ***
plot.means <- aggregate(SECCp$Y.trans, list(Chamber=SECCp$Chamber, 
                                            Position=SECCp$Position), 
                        mean)
plot.means$error <- as.numeric(msd["Chamber:Position"]/2)
plot.data <- decomp.plotcalcs(plot.means)

CxP.plot <- qplot(Chamber, xtr, data = plot.data, group = Position, 
                    geom = "point", ylim = Y.lim, size = I(3), 
                    colour = Position, shape = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Chamber"]],
                    ylab = Y.plotlab,
                    legend = FALSE)
## CxP.plot <- CxP.plot + geom_point(aes(Chamber, x), size = 2)
CxP.plot <- CxP.plot + geom_line(aes(group = Position), size = 0.8)
CxP.plot <- CxP.plot + geom_errorbar(aes(ymin = lwr, ymax = upr), 
                                         width = 0.2, size = 0.5)
CxP.plot <- CxP.plot + scale_colour_manual(name = Position.label,
                                           values = Position.map$col, 
                                           breaks = Position.map$label)
CxP.plot <- CxP.plot + scale_fill_manual(name = Position.label,
                                         values = Position.map$bg, 
                                         breaks = Position.map$label)
CxP.plot <- CxP.plot + scale_shape_manual(name = Position.label,
                                           values = Position.map$pch, 
                                           breaks = Position.map$label)
CxP.plot <- CxP.plot + jaw.ggplot()
print(CxP.plot)



if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 3, height = 3, scale = 1.5)
}

