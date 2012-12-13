##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Moss dry weight of patch
### Jonathan Whiteley     R v2.12     2012-06-26
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
Y.col <- 'Patch.dwt'  # Column to analyze as response variable           *****
Y.use <- 'Y.log'    # Which transformation is being used (for labels)? *****

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
Time.use     <- levels(SECC$Time)              # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1,3)]   # Chamber treatments to include
Frag.use     <- levels(SECC$Frag)              # Frag treatments to include
Position.use <- levels(SECC$Position)[c(1, 3)] # Patch Positions to include

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
if (FALSE) for ( Time.i in 1:length(levels(SECC$Time)) ) {
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
## Stock graphs: only Frag & Position main effects are significant, no interactions.
## Looks like contiguous patches are consistently bigger in terms of biomass, while "other" patches are slightly lower
## The slightly more massive contiguous patches is not terribly surprising, given the slightly different approach to set them up, and the lower level of disturbance needed to create them.
## The difference in "other" patches is surprising, however.
## I expected the outer chamber patches would have slightly lower mass, given the higher decomposition : productivity ratio.
## - The trend is there, but it is weak and not statistically significant.

Y.lim <- c(0, 30)
Plot.Title <- bquote(.(Time.label) * "Patch means " %+-% "95% Comparison Intervals")
Sub.msd <- "95% comparison intervals (MSR)" 

Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = levels(SECC$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]

## data frame of plot values (for ggplot2).
## might be able to accomplish much the same effect with stat_summary using means in ggplot2?
plot.means <- SECCplotDataANOVA(SECCp$Y.trans, 
                                list(Chamber=SECCp$Chamber, 
                                     Position=SECCp$Position, Time=SECCp$Time), 
                                error = msd["Time:Chamber:Position"]
                                )
levels(plot.means$Chamber)[2] <- "Chamber"

CxP.plot <- qplot(Chamber, x, data = plot.means, group = Position, 
                    geom = "line", ylim = Y.lim, size = Position,
                    colour = Position, shape = Position, fill = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Chamber"]],
                    ylab = Y.plotlab,
                    legend = FALSE,
                    facets = .~Time)
CxP.plot <- CxP.plot + geom_errorbar(aes(ymin = lower, ymax = upper), 
                                         width = 0.2, size = 0.5)
CxP.plot <- CxP.plot + geom_point(aes(group = Position), size = 3)
CxP.plot <- CxP.plot + scale_colour_manual(name = Position.label,
                                           values = Position.map$col, 
                                           breaks = Position.map$label)
CxP.plot <- CxP.plot + scale_fill_manual(name = Position.label,
                                         values = Position.map$bg, 
                                         breaks = Position.map$label)
CxP.plot <- CxP.plot + scale_shape_manual(name = Position.label,
                                          values = Position.map$pch, 
                                          breaks = Position.map$label)
CxP.plot <- CxP.plot + scale_size_manual(name = Position.label,
                                         values = Position.map$lwd*0.5, 
                                         breaks = Position.map$label)
CxP.plot <- CxP.plot + jaw.ggplot()
print(CxP.plot)


## Frag x Pos Interaction
plot.means <- SECCplotDataANOVA(SECCp$Y.trans, 
                                list(Chamber=SECCp$Chamber, Frag=SECCp$Frag, 
                                     Position=SECCp$Position, Time=SECCp$Time), 
                                error = msd["Time:Chamber:Frag:Position"]
                                )
levels(plot.means$Chamber)[2] <- "Chamber"

FragIconList <- list(FragIcon1 = FragIcon1,
                     FragIcon2 = FragIcon2,
                     FragIcon3 = FragIcon3,
                     FragIcon4 = FragIcon4
                     )

FxP.plot <- qplot(Frag, x, data = plot.means, group = Position, 
                    geom = "line", ylim = Y.lim, size = Position, 
                    colour = Position, fill = Position, shape = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Frag"]],
                    ylab = Y.plotlab,
                    legend = FALSE,
                    facets = Chamber~Time)
FxP.plot <- FxP.plot + 
geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.5) +
geom_point(aes(shape = Position), size = 3) +
scale_colour_manual(name = Position.label,
                    values = Position.map$col, 
                    breaks = Position.map$label) +
scale_fill_manual(name = Position.label,
                  values = Position.map$bg, 
                  breaks = Position.map$label) +
scale_shape_manual(name = Position.label,
                   values = Position.map$pch, 
                   breaks = Position.map$label) +
scale_size_manual(name = Position.label,
                  values = Position.map$lwd*0.5, 
                  breaks = Position.map$label) +
scale_x_discrete(labels = c(1, 2, 3, 4), 
                 breaks = levels(plot.means$Frag)) + jaw.ggplot() 
## Add imported graphics as x-axis tick labels :D
## http://stackoverflow.com/questions/2181902/how-to-use-an-image-as-a-point-in-ggplot
FxP.plot <- FxP.plot + scale_x_discrete(labels = names(FragIconList), # c(1, 2, 3, 4), 
                                        breaks = levels(plot.means$Frag)) +
opts(axis.ticks.margin = unit(0.2, "lines"),
     axis.text.x = picture_axis(FragIconList, icon.size = unit(1.4, "lines")) 
)
print(FxP.plot)


if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 6, height = 3, scale = 1)
  ggsave(file = paste(Save.final, "- FxP.eps"), plot = FxP.plot, width = 6, height = 4, scale = 1.2)
}
