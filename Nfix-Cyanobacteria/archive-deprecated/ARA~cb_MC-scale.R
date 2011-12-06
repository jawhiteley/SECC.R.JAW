################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### GLMM (regression)
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2011-11-02
################################################################
## INITIALISE
################################################################
## Working Directory: see lib/init.R [\rd in Vim]
if (FALSE) {  # do not run automatically
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (MBP).
  getwd()  # Check that we're in the right place

  ## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
  ## load & prepare data
  source('1_ARA-cb_setup.R')

  ## library(lattice)    # ggplot2 with faceting is easier!
  library(ggplot2)
  theme_set(theme_bw())                  # change global ggplot2 theme
  library(rgl)                           # 3D plots
  library(car)                           # diagnostic plots & tools
  library(gvlma)                         # Global Validation of Linear Model Assumptions
  library(nlme)                          # GLMMs (older, but still works)
  ## library(lme4)    # I'd rather use lme4 for GLMMs, but all my references use nlme
  ## library\(mgcv)   # Additive Modelling, Cross-Validation

}




################################################################
## META-COMMUNITY SCALE ANALYSIS
################################################################

SECCa <- SECCmc 


SECCa <- within( SECCa, {
                X <- as.numeric( get(X.col) )
                Y <- as.numeric( get(Y.col) )
                Y.use <- Y  # compatibility with older code
                X.log <- log10(X)
                X.log[X <= 0] <- 0
                Y.log <- log10(Y)
                Y.log[Y <= 0] <- 0
})

##================================================
## EXPLORE: PLOTS
##================================================
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
  pairplot(SECCa[, c("ARA.m", "ARA.g", "H2O", "Cells.m", "Cells.g", "Hcells.m", "Hcells.g", "Stigonema", "Nostoc" )])
  ## look at log transformations
  pairplot(SECCa[, c('Y', 'Y.log', 'X', 'X.log', 'H2O')],
           labels=c(Y.col, paste("log(", Y.col, ")"), 
                    X.col, paste("log(", X.col, ")"), "H2O"
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

##==============================================================
### prepare plot theme for points varying by Chamber
ChamberPts  <- ggPts.SECC(Chamber.map, Chamber.label) 
TopLegend   <- opts(legend.position = "top", legend.direction = "horizontal")
Time.facets <- facet_grid(facets = .~Time)
All.facets  <- facet_grid(facets = Frag~Time) # The only real difference from main data exploration code above

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

## log-log: messy, apparently negative relationships? (too much noise from inner/outer patches being averaged together in Chambers?
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

## H2O plots: A few clear differences in log slope among Frag Treatments, at least in August samples. **
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


##==============================================================
## 3D plot with X, Y & H2O axes (requires rgl package)
if (DrawExplorationGraphs) {
  with(SECCa, plot3d(X, H2O, Y, size = 6, pch = pt, col = colr, bg = fill) )
  with(SECCa, plot3d(X.log, H2O, Y.log, size = 6, pch = pt, col = colr, bg = fill) )
}


