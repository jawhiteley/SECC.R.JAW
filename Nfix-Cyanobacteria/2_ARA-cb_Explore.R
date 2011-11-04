################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Data Exploration
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2011-11-02
################################################################
## INITIALISE
################################################################
## Working Directory: see lib/init.R below [\rd in Vim]
if (FALSE) {  # do not run automatically
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (MBP).
  getwd()  # Check that we're in the right place

  ## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
  source('1_ARA-cb_setup.R')
}

## library(lattice)    # ggplot2 with faceting is easier!
library(ggplot2)
theme_set(theme_bw())                  # change global ggplot2 theme
library(rgl)                           # 3D plots
library(car)                           # diagnostic plots & tools

################################################################
## EXPLORE: PLOTS
################################################################
## make some meaningful plots of data to check for predicted (expected) patterns.
DrawExplorationGraphs <- TRUE # Save.results  # Set to FALSE to suppress all this output
if (Save.results == TRUE && is.null(Save.plots) == FALSE) pdf( file = Save.plots )

### Map of point styles for Chamber treatments
Chamber.map <- plotMap( "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Chamber.map$label <- factor(Chamber.map$label)
point <- 21	# 21 for circles with col & bg ; 16 for solid circles
Chamber.map$pch <- c(21, 16)  # use circles for both treatments
Chamber.label <- attr(SECC, "labels")[["Chamber"]]

## Set up some plotting options for each data point
##  only needed for non-ggplot commands
SECCa <- within( SECCa,{
                colr = ifelse(Chamber == Chamber.map$label[1], 
                              Chamber.map$col[1], 
                              Chamber.map$col[2] 
                              )
                fill = ifelse(Chamber == Chamber.map$label[1], 
                              Chamber.map$bg[1], 
                              Chamber.map$bg[2] 
                              )
                pt   = ifelse(Chamber == Chamber.map$label[1], 
                              Chamber.map$pch[1], 
                              Chamber.map$pch[2]
                              )
})

### pairplots of several variables of interest: check for collinearity, patterns, etc.
### see Zuur et al. books
if (DrawExplorationGraphs) {
  ## full dataset: unbalanced with respect to experimental treatments
  pairplot(SECC[, c("ARA.m", "ARA.g", "H2O", "Cells.m", "Cells.g", "Hcells.m", "Hcells.g", "Stigonema", "Nostoc" )])
  ## filtered dataset: balanced, but am I missing useful information about continuous explanatory variables (H2O, cells)?
  pairplot(SECCa[, c("ARA.m", "ARA.g", "H2O", "Cells.m", "Cells.g", "Hcells.m", "Hcells.g", "Stigonema", "Nostoc" )])
  ## look at log transformations & check for colinearity among explanatory variables
  pairplot(SECCa[, c('Y', 'Y.log', 'X', 'X.log', 'H2O', 'Time', 'Chamber', 'Frag', 'Position')],
           labels=c(Y.col, paste("log(", Y.col, ")"), 
                    X.col, paste("log(", X.col, ")"), "H2O"
                    , "Time", "Chamber", "Frag", "Position"
                    )
          )

  ## Cleveland Dotplots (Zuur et al. 2009)
  Dotplot.y <- "Order of observations"
  op <- par(mfcol=c(2,2))
  dotchart(SECCa$X, ylab=Dotplot.y, xlab=X.label)
  dotchart(SECCa$X.log, ylab=Dotplot.y, xlab=paste("log(", X.label, ")"))
  dotchart(SECCa$Y, ylab=Dotplot.y, xlab=Y.label)
  dotchart(SECCa$Y.log, ylab=Dotplot.y, xlab=paste("log(", Y.label, ")"))
  par(op)
}

if (FALSE) {
  ## the old-fashioned way (low-level)
  with( SECCa,{
       ## scatterplot
       plot(X, Y, type="p",
            ylab=Y.plotlab, xlab=X.plotlab, 
            pch=pt, col=colr
            )
  legend("topright", legend=Chamber.map$label, 
         pch=Chamber.map$pch, col=Chamber.map$col 
  )
  ##     identify(X, Y, labels = SampleID)
  })
}

##==============================================================
## The easy way, with ggplot2
##==============================================================
### prepare plot theme for points varying by Chamber
ChamberPts  <- ggPts.SECC(Chamber.map, Chamber.label) 
TopLegend   <- opts(legend.position = "top", legend.direction = "horizontal")
Time.facets <- facet_grid(facets = .~Time)
All.facets  <- facet_grid(facets = Frag~Time*Position)

ARAcb.plot  <- qplot(X, Y, data = SECCa, group = Chamber,
                   geom = "point", size = I(3),
                   colour = Chamber, shape = Chamber,
                   xlab = X.plotlab,
                   ylab = Y.plotlab
                   )
ARAcb.plot  <- ARAcb.plot + jaw.ggplot() + ChamberPts + TopLegend  # order matters!

ARAcb.time   <- ARAcb.plot + Time.facets # Faceted by Time pts
ARAcb.panels <- ARAcb.plot + All.facets # Full Faceting ***
if (DrawExplorationGraphs) {
  print(ARAcb.plot)
  print(ARAcb.time)
  print(ARAcb.panels)
}

## log-log looks encouraging
##  Why can't I do this by adding + coord_trans(x = "log", y = "log") ??
##  (I get an error when I do this to any of the above plots, even without faceting.
ARAcb.log <- qplot(X, Y, data = SECCa, log = "xy", 
                    group = Chamber, geom = "point", size = I(3),
                    colour = Chamber, shape = Chamber,
                    xlab = X.plotlab,
                    ylab = Y.plotlab
                    )
ARAcb.log <- ARAcb.log + ChamberPts + jaw.ggplot() + TopLegend

ARAcb.time.log   <- ARAcb.log + Time.facets 
ARAcb.panels.log <- ARAcb.log + All.facets 
if (DrawExplorationGraphs) {
  print(ARAcb.log)
  print(ARAcb.time.log)
  print(ARAcb.panels.log)
}

## H2O plots
ARA.H2O <- qplot(H2O, Y, data = SECCa, 
                    group = Chamber, geom = "point", size = I(3),
                    colour = Chamber, shape = Chamber,
                    ylab = Y.plotlab
                    )
ARA.H2O <- ARA.H2O + ChamberPts + jaw.ggplot() + TopLegend

ARA.H2O.time   <- ARA.H2O + Time.facets
ARA.H2O.panels <- ARA.H2O + All.facets
if (DrawExplorationGraphs) {
  print(ARA.H2O)
  print(ARA.H2O.time)
  print(ARA.H2O.panels)
}

ARA.H2O.log <- qplot(H2O, Y, data = SECCa, log = "y", 
                    group = Chamber, geom = "point", size = I(3),
                    colour = Chamber, shape = Chamber,
                    ylab = Y.plotlab
                    )
ARA.H2O.log <- ARA.H2O.log + ChamberPts + jaw.ggplot() + TopLegend

ARA.H2O.log.time   <- ARA.H2O.log + Time.facets
ARA.H2O.log.panels <- ARA.H2O.log + All.facets
if (DrawExplorationGraphs) {
  print(ARA.H2O.log)
  print(ARA.H2O.log.time)
  print(ARA.H2O.log.panels)
}


if (FALSE) {    ## Deprecated scatterplots, broken down in various ways.
  ## coplot() deprecated - qplot with faceting is easier & better-looking (but slower)
  coplot( Y ~ X | Frag * Position, data=SECCa, 
         pch=SECCa$pt, col=SECCa$colr	# , bg=Chamber.map$bg
  )  # why does recycling Chamber.map work for bg, but not col?
  ## I think it's because of the way the arguments are passed: coplot has specific arguments for col & pch, but not bg: bg is simply passed on directly to the plotting function (points) within each panel.
  coplot( Y ~ X | Chamber * Frag , data=SECCa, 
         pch=SECCa$pt, col=SECCa$colr	# , bg= Chamber.map$bg
  )
  coplot( Y ~ X | Chamber * Position , data=SECCa, 
         pch=SECCa$pt, col=SECCa$colr	# , bg= Chamber.map$bg
  )
                                        # Moisture?
  coplot( Y ~ H2O | Frag * Position , data=SECCa, 
         pch=SECCa$pt, col=SECCa$colr	# , bg= Chamber.map$bg
  )
}

##==============================================================
## 3D plot with X, Y & H2O axes (requires rgl package)
if (DrawExplorationGraphs) {
  with(SECCa, plot3d(X, H2O, Y, size = 6, pch = pt, col = colr, bg = fill) )
  with(SECCa, plot3d(X.log, H2O, Y.log, size = 6, pch = pt, col = colr, bg = fill) )
}


##==============================================================
## Averaged over time?
##  I know the measurements were not taken directly from the same patch,
##  but I wonder if some of the noise can be reduced by using long-term averages
##  similar to the way DeLuca pooled several sample together for his data (2007)
SECCtavg <- SECC_aggregate(SECCa, drop.trt="Time")

## Higher correlations with Moisture now.
if (DrawExplorationGraphs) {
  pairplot(SECCtavg[, c("ARA.m", "ARA.g", "H2O", "Cells.m", "Cells.g", "Hcells.m", "Hcells.g", "Stigonema", "Nostoc" )])
  ## look at log transformations
  pairplot(SECCtavg[, c('Y', 'Y.log', 'X', 'X.log', 'H2O')],
           labels=c(Y.col, paste("log(", Y.col, ")"), 
                    X.col, paste("log(", X.col, ")"), "H2O"
                    )
          )

  ## Cleveland Dotplots (Zuur et al. 2009)
  Dotplot.y <- "Order of observations"
  op <- par(mfcol=c(2,2))
  dotchart(SECCtavg$X, ylab=Dotplot.y, xlab=X.label)
  dotchart(SECCtavg$X.log, ylab=Dotplot.y, xlab=paste("log(", X.label, ")"))
  dotchart(SECCtavg$Y, ylab=Dotplot.y, xlab=Y.label)
  dotchart(SECCtavg$Y.log, ylab=Dotplot.y, xlab=paste("log(", Y.label, ")"))
  par(op)
}

## Some clear patterns, but still pretty messy
Tavg.plot  <- qplot(X, Y, data = SECCtavg, group = Chamber,
                   geom = "point", size = I(3),
                   colour = Chamber, shape = Chamber,
                   xlab = X.plotlab,
                   ylab = Y.plotlab
                   )
Tavg.plot  <- Tavg.plot + jaw.ggplot() + ChamberPts + TopLegend  # order matters!
Tavg.panels <- Tavg.plot + facet_grid(Frag~Position) # Full Faceting ***
if (DrawExplorationGraphs) {
  print(Tavg.plot)
  print(Tavg.panels)
}

## A lot of flat (or negative) relationships on these log-log plots.
## I wonder if the low ARA values of the June sample are suppressing patterns in the August sample?
Tavg.log.plot  <- qplot(X, Y, data = SECCtavg, log="xy", group = Chamber,
                   geom = "point", size = I(3),
                   colour = Chamber, shape = Chamber,
                   xlab = X.plotlab,
                   ylab = Y.plotlab
                   )
Tavg.log.plot  <- Tavg.log.plot + jaw.ggplot() + ChamberPts + TopLegend  # order matters!
Tavg.log.panels <- Tavg.log.plot + facet_grid(Frag~Position) # Full Faceting ***
if (DrawExplorationGraphs) {
  print(Tavg.log.plot)
  print(Tavg.log.panels)
}


##==============================================================
## Check Variation, Ranges
X.box <- qplot(Time, X, data = SECCa, geom = "boxplot",
                 ylab = X.plotlab)
X.box <- X.box + jaw.ggplot()
if (DrawExplorationGraphs) print(X.box)

Y.logBlock <- qplot(Block, Y.log, data = SECCa, geom = "boxplot",
                 ylab = Y.plotlab) + jaw.ggplot()
Y.Block <- qplot(Block, Y, data = SECCa, geom = "boxplot",
                 ylab = Y.plotlab) + jaw.ggplot()
Y.Time <- qplot(Time, Y.log, data = SECCa, geom = "boxplot",
                 ylab = Y.plotlab) + jaw.ggplot()
Y.Chamber <- qplot(Chamber, Y.log, data = SECCa, geom = "boxplot",
                   ylab = Y.plotlab) + jaw.ggplot() 
Y.Frag <- qplot(Frag, Y.log, data = SECCa, geom = "boxplot",
                ylab = Y.plotlab) + jaw.ggplot() 
Y.Pos  <- qplot(Position, Y.log, data = SECCa, geom = "boxplot",
                ylab = Y.plotlab) + jaw.ggplot() 
if (DrawExplorationGraphs) {
  print(Y.logBlock)
  print(Y.Block)
  print(Y.Time)
  print(Y.Chamber)
  print(Y.Frag)
  print(Y.Pos)
}

##==============================================================
## Check distributions
X.dist <- qplot(X, data = SECCa, geom = "histogram",
                 xlab = X.plotlab
                 )
X.dist <- X.dist + jaw.ggplot() + All.facets
if (DrawExplorationGraphs) print(X.dist)

Y.dist <- qplot(Y, data = SECCa, geom = "histogram",
                 xlab = Y.plotlab
                 )
Y.dist <- Y.dist + jaw.ggplot() + All.facets
if (DrawExplorationGraphs) print(Y.dist)

if (DrawExplorationGraphs) {
  old.par <- par(mfcol=c(2,2))
  for(i in 1:length(vars.ls) ) {
    var <- vars.ls[i]
    label <- labels.ls[i]
    with( SECCa,{
         X.var <- get(var)
    X.max  <- max( X.var )
    ##         cat(var, " ", X.max, "\n")
    freq.max <- length(X.var)/2
    X.maxD <- max( density( X.var )$y )*1.5
    for(Ch.trt in levels(Chamber)) {
      X.trt <- X.var[Chamber==Ch.trt]
      X.density <- density( X.trt )
      hist(X.trt,
           main=Ch.trt, xlab=label,
           xlim=c(0, X.max), 
           ylim=c(0, freq.max),
           breaks=seq( 0, X.max, length.out=16 ),
           col="#CCCCCC"
           )
###       abline( 5, 0, lty=3, col="#666666" ) # reference line
      plot(X.density,
           main=Ch.trt, xlab=label,
           xlim=c(0, X.max),
           ylim=c(0, X.maxD)
           )
      densityplot( X.trt,
                  main=Ch.trt, xlab=label,
                  xlim=c(0, max( X.var ) )
                  )	# ** TRELLIS plot
      }
    })
    ## qqplot(x, y) to compare distributions.
  }
  par(old.par)
}

