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
attr(SECC, "units" )[["H2O"]] <- quote("("* H[2]*O *" as % of dry wt. moss)")


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

Y.lim <- c(0, 800)
Plot.Title <- bquote(.(Time.label) * "Patch means " %+-% "95% Comparison Intervals")
Sub.msd <- "95% comparison intervals (MSR)" 

if (FALSE) {
  postscript(file = Save.final, width = 6, height = 2)


Chamber.map <- plotMap( factor = "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Frag.map <- plotMap( factor = "Frag", labels = levels(SECC$Frag) )
Frag.map <- Frag.map[ levels(SECC$Frag) %in% Frag.use, ]
Position.map <- plotMap( factor = "Position", labels = levels(SECC$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]
Position.col <- Position.map$col
names(Position.col) <- Position.map$label

old.par <- par(mfrow=c(1,3), las = 1, oma = c(3, 2, 3, 1), mar = c(3, 3, 2, 0) +0.1 )

for(Time.for in 1:length(Time.use)) {
with( SECCp[SECCp$Time == Time.use[Time.for], ], {
  plot.error <- matrix( as.numeric(msd["Chamber:Position"]/2),
                       nrow = length(levels(Chamber)),  # rows: x-axis
                       ncol = length(levels(Position))  # cols: y-axis
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
dev.off()
}

Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = levels(SECC$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]

## data frame of plot values (for ggplot2).
## might be able to accomplish much the same effect with stat_summary using means in ggplot2?
plot.means <- aggregate(SECCp$Y.trans, list(Chamber=SECCp$Chamber, Position=SECCp$Position, Time=SECCp$Time), mean)
levels(plot.means$Time) <- paste(c("August", "June", "August"), levels(plot.means$Time), sep="\n")
plot.means$error <- as.numeric(msd["Time:Chamber:Position"]/2)
levels(plot.means$Chamber)[2] <- "Chamber"

CxP.plot <- qplot(Chamber, x, data = plot.means, group = Position, 
                    geom = "point", ylim = Y.lim, size = I(3), 
                    colour = Position, shape = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Chamber"]],
                    ylab = Y.plotlab,
                    legend = FALSE,
                    facets = .~Time)
## CxP.plot <- CxP.plot + geom_point(aes(Chamber, x), size = 2)
CxP.plot <- CxP.plot + geom_line(aes(group = Position), size = 0.8)
CxP.plot <- CxP.plot + geom_errorbar(aes(ymin = x - error, ymax = x + error), 
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


## Frag x Pos Interaction

Y.lim <- c(-100, 1000)
plot.means <- aggregate(SECCp$Y.trans, list(Chamber=SECCp$Chamber, Frag=SECCp$Frag, Position=SECCp$Position, Time=SECCp$Time), mean)
levels(plot.means$Time) <- paste(c("August", "June", "August"), levels(plot.means$Time), sep="\n")
plot.means$error <- as.numeric(msd["Time:Chamber:Frag:Position"]/2)
levels(plot.means$Chamber)[2] <- "Chamber"
FragIconList <- list(FragIcon1 = FragIcon1,
                     FragIcon2 = FragIcon2,
                     FragIcon3 = FragIcon3,
                     FragIcon4 = FragIcon4
                     )

FxP.plot <- qplot(Frag, x, data = plot.means, group = Position, 
                    geom = "point", ylim = Y.lim, size = I(3), 
                    colour = Position, shape = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Frag"]],
                    ylab = Y.plotlab,
                    legend = FALSE,
                    facets = Chamber~Time)
FxP.plot <- FxP.plot + 
geom_line(aes(group = Position), size = 0.8) +
geom_errorbar(aes(ymin = x - error, ymax = x + error), width = 0.2, size = 0.5) +
scale_colour_manual(name = Position.label,
                    values = Position.map$col, 
                    breaks = Position.map$label) +
scale_fill_manual(name = Position.label,
                  values = Position.map$bg, 
                  breaks = Position.map$label) +
scale_shape_manual(name = Position.label,
                   values = Position.map$pch, 
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
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 6, height = 3, scale = 1.5)
  ggsave(file = paste(Save.final, "- FxP.eps"), plot = FxP.plot, width = 6, height = 4, scale = 1.5)
}
