################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### GLMM (regression)
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2011-10-05
################################################################
## INITIALISE
################################################################
## Working Directory: see lib/init.R below [\rd in Vim]
if (FALSE) {  # do not run automatically
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (MBP).
  getwd()  # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')

## library(car)
## library(lattice)    # ggplot2 with faceting is easier!
library(ggplot2)
theme_set(theme_bw())                  # change global ggplot2 theme
library(rgl)        # 3D plots
library(nlme)       # GLMMs (older, but still works)
## library(lme4)    # I'd rather use lme4 for GLMMs, but all my references use nlme
## library\(mgcv)   # Cross-Validation

################################################################
## CONFIGURE BASIC ANALYSIS
################################################################

### Response Variable *****
X.col <- 'Cells.m'   # Column to analyze as explanatory variable        *****
Y.col <- 'ARA.m'     # Column to analyze as response variable           *****
## I would prefer to use measurements on a per-gram dwt of moss basis ('.g'),
## but there aren't enough of such measurements for all ARA values (univariate).
## There are certainly enough for this analysis, but I'm worried about inconsistency
## for readers if N-fixation is analyzed on its own on a per-m^2 basis, 
## while everything else is on a (preferred) per-g dwt basis.

# explanatory vars for data exploration?
vars.ls   <- c("ARA.m", "Cells.m", "Hcells.m", "H2O")  

##==============================================================
## SETTINGS 
##==============================================================
## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)      # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]      # Chamber treatments to include
Frag.use     <- levels(SECC$Frag)         # Frag treatments to include
Position.use <- levels(SECC$Position)[c(1, 3)]     # Patch Positions to include


## Output Results?
Save.results  <- FALSE  


##==============================================================
## CALCULATIONS 
##==============================================================
SECC.prime <- SECC    # save a copy of the original for reference.

# str(SECC)
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
  Nfix    <- ARA.m * Nfix.ARA.ratio
  H2O     <- H2O * 100
  H2O.wwt <- H2O.wwt * 100
})


##==============================================================
## LABELS
##==============================================================

X.label <- attr(SECC, "labels")[[X.col]]  # explanatory variable label
X.units <- attr(SECC, "units" )[[X.col]]  # explanatory variable units

Y.label <- attr(SECC, "labels")[[Y.col]]  # response variable label
Y.units <- attr(SECC, "units" )[[Y.col]]  # response variable units

X.plotlab <- bquote( .(X.label) * "  " * .(X.units) *  "" )
Y.plotlab <- bquote( .(Y.label) * "  " * .(Y.units) *  "" )

labels.ls <- c()
for(i in 1:length(vars.ls) ){
  var <- vars.ls[i]
  labels.ls[i] <- attr(SECC, "labels")[[var]]
}


## Save Output to Files - set to NULL to prevent output.
Save.filename <- paste("Results - ", Y.col, "~" , X.col, " - ",
                       paste(which(levels(SECC$Time) == Time.use), collapse=""),
                       sep = ""
                   )
Save.text  <- paste("./output/", Save.filename, ".txt", sep = "")
Save.plots <- paste("./graphs/", Save.filename, ".pdf", sep = "")
Save.final <- Save.plots              # Destination for final plots.

## Output text
Save.divider <-        "================================================================\n"
Save.header  <- paste( "GLM Results for:", Y.label, "(", Y.col, ")",
                     "\n               ~", X.label, "(", X.col, ")",
                     "\nExpt. Time:   ", paste(Time.use,     collapse = ", "),
                     "\nChamber:      ", paste(Chamber.use,  collapse = ", "),
                     "\nFragmentation:", paste(Frag.use,     collapse = ", "),
                     "\nPatches:      ", paste(Position.use, collapse = ", "),
                     paste("\n\n", date(), "\n\n", Save.divider, sep = "")
                     )
Save.end      <- paste("\n", 
					   "<============================= END ============================>",
						sep = "\n"
					  )

###===============================================
### Include Time as a factor 
###===============================================
## Note that Samples at different times are actually independent
## in this design, due to destructive sampling.

## Time.use     <- levels(SECC$Time)      # Include *ALL* Times (as a Treatment)
## source("./R", echo = FALSE) # Load default Labels. *****
## source("./R", echo = FALSE) # RUN STANDARD nested ANOVA



################################################################
## PROCESS DATA: planned
################################################################
## filter & clean SECC data frame.
SECC.use <- SECCclean(SECC, Time.use, Chamber.use, Frag.use, Position.use)

## Summarize data by means across different (or all) positions to prevent unbalanced effects?
## aggregate 'other' patch Positions?
##   unnecessary: only 'inner' & 'outer' included.
SECCp  <- if (FALSE)  SECC_aggregate( SECC.use, trt = 'Position' )  else SECC.use

## Meta-community (regional) -level analyses (ignoring position):
SECCmc <- SECC_aggregate( SECC.use, trt = 'Frag' )

SECC.scale  <- "patch"  # c("patch", "mc")
SECCa <- if(SECC.scale == "patch") SECCp else if (SECC.scale == "mc") SECCmc else SECC


SECCa <- within( SECCa, {
                X <- as.numeric( get(X.col) )
                Y <- as.numeric( get(Y.col) )
                Y.use <- Y  # compatibility with older code
                X.log <- log10(X)
                X.log[X <= 0] <- 0
                Y.log <- log10(Y)
                Y.log[Y <= 0] <- 0
})

################################################################
## CHECK DATA
################################################################
if (FALSE) {  # do not run if source()d
  head(SECCa)     # have a peek at the first 6 rows & columns: is this what you expected?
  str( SECCa)     # check structure: are the appropriate variables factors, numeric, etc.?
  summary(SECCa)  # summary statistics
}


################################################################
## EXPLORE: PLOTS
################################################################
## make some meaningful plots of data to check for predicted (expected) patterns.
DrawExplorationGraphs <- Save.results  # Set to FALSE to suppress all this output
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
                 ylab = X.plotlab
                 )
X.box <- X.box + jaw.ggplot()
if (DrawExplorationGraphs) print(X.box)

Y.box <- qplot(Time, Y, data = SECCa, geom = "boxplot",
                 ylab = Y.plotlab
                 )
Y.box <- Y.box + jaw.ggplot()
if (DrawExplorationGraphs) print(Y.box)

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





################################################################
## ANALYSIS
################################################################
## Goal: Identify which variables have most relative importance,
##       and identify any interactions with experimental treatments.
## Alternative outline (based on Zuur et al. examples):
## progress from simple to "ultimate" (final) model / approach
## Cycle through:
##  1. Model fitting
##  2. Test Assumptions (graphical analysis of residuals)
##  2a) Identify violated assumptions, and possible solutions
##  3. Model Selection on *valid* models
##  3a) Check Model output, ANOVA tables, parameter values, AIC, R^2, etc.
##  3b) Identify possible improvements
##  END:    choose optimal model; produce model output: tables and graphs.
## Start with GLMM, to incorporate known nested treatments.
## Build in extensions as necessary, time permitting.
##==============================================================
## Model Formula
##==============================================================
### Fixed effects
### Include H2O^2 to test for unimodal relationship?
## Including Time as a factor?
if ( length(Time.use) > 1 ) {
  Y.formula <- Y.log ~ X.log * H2O * Time * Chamber * Frag * Position
} else {
  Y.formula <- Y.log ~ X.log * H2O * Chamber * Frag * Position
}

## Random effects for GLMM
Y.Ri  <- ~ 1|Block/Time/Chamber/Frag                 # Random intercept across Blocks
Y.Ris <- ~ 1 + X.log + H2O | Block/Time/Chamber/Frag # Random intercept + slope




if (FALSE) {    # Model II Regression: nice idea in theory
################################################################
## Model II Regression?                                       **
################################################################
## Given that both continuous explanatory variables (cyanobacteria, H2O)
## were measured and not controlled, this is most likely a "Model II regression"
##  Model 1 regression *underestimates the slope*, because both variables contain error
##  Model II regression allows X to vary (no "Fixed X" assumption),
##  but different methods for *prediction* (OLS) or
##  describing the *functional relationship* (structural relation, major axis, etc.)
##  - See Sokal & Rohlf p. 541, 545
##  - See Legendre & Legendre p.504
##  Honestly, I think I'm more interested in the latter (functional relationship) than the former (prediction), although I'm not aware of any methods for model 2 GLMMs with multiple variables and nested data, and I'm not sure I have the time to learn/develop them.

library(lmodel2)

Y.m2 <- lmodel2(Y.log ~ X.log, data=SECCa, nperm=999) # Model II regression (Legendre & Legendre)
print(Y.m2)                            # What does this tell me, exactly?
plot(Y.m2)

Y.m2 <- lmodel2(Y.log ~ X.log * H2O, data=SECCa, nperm=999) # Model II regression (Legendre & Legendre)
print(Y.m2)                            # What does this tell me, exactly?
plot(Y.m2)

## lmodel2 seems to drop all but the first explanatory variable!
## I would have to fit models for each pair of variables across all combinations
## of treatment factors?  Ouch.  How would I compare them?

## MA (Major Axis) is also called "geometric mean regression" (Sokal & Rohlf p.544)
## and can easily be computed as the geometric mean of the slopes of y~x and x~y
## <http://tolstoy.newcastle.edu.au/R/help/05/06/5992.html>
## I'm guessing that the MA also passes through the intersection of the y~x & x~y regression lines (the "centroid" of the data points?)

## Legendre & Legendre (p.515) suggest: 
## "MA may also be used with dimensionally heterogeneous variables when the purpose of the analysis is 
## **(1) to compare the slopes of the relationships between the same two variables measured under different conditions (e.g. at two or more sampling sites), or 
##   (2) to test the hypothesis that the major axis does not significantly differ from a value given by hypothesis
## "
## AND
## "If a straight line is not an appropriate model, polynomial or nonlinear regression should be considered."
## 
## see also the User Guide to the 'lmodel2' R package:
## vignette("mod2user", package="lmodel2")

## Therefore:
## In theory, I could use MA for a simple relationship between ARA & Cells,
## for each combination of experimental treatments (Chamber x2, Frag x4, Position x2).
## and compare them using some sort of non-parametric, 
## bootstapped, or permutation-based test statistic,
## while controlling for experiment-wise multiple comparisons.
## OR
## I could just proceed with non-linear regression (GLMM/GAMM), 
## and accept that the regressions parameters (slopes) are likely underestimated
## (I call this "conservative").
## Furthermore, I'm not convinced the relationship *is* linear,
## or that is dimensionally homogeneous or bivariate normal,
## and I want to include multiple fixed and random variables, with a nested structure,
## which would be rather difficult for me to do in a reasonable amount of time
## within a model 2 approach.

}




if (FALSE) {    # GLM: deprecated. wrapped to allow folding
##################################################
## BASIC GLM - soon to be deprecated by GLMM (below)
##################################################
##================================================
## Model Fitting
##================================================
### "basic" GLM: initial foray into analysis.  Soon to be deprecated by GLMM.
Y.model <- glm( Y.formula, data=SECCa, family="gaussian" )
Y.model.full <- Y.model
Y.model.main <- glm(Y ~ log(X+1) + Time + Chamber + Frag + Position + H2O + I(H2O^2), data = SECCa)  # Main factors only; no interaction terms.

# Should I not be using a mixed-effects model to account for treatments of different sizes? Chamber / Frag / Position
# YES (See below)

##================================================
## MODEL SELECTION
##================================================
drop1(Y.model)
Y.model.selected <- step(Y.model, direction = "backward")
Y.model.main.selected <- step(Y.model.main, direction = "backward")

Y.model  <- Y.model.selected


##################################################
## CHECK ASSUMPTIONS (MODEL VALIDATION)
##################################################
## analyse residuals, standard diagnostic plots
par(mfrow=c(2,2))	 # panel of figures: 3 rows & 2 columns
plot(Y.model)

## Residuals
Model.resid <- resid(Y.model)
## par(mfrow=c(1,1))	 # panel of figures: 1 rows & 1 columns
hist(Model.resid)    # plot residuals
plot(SECCa$X, Model.resid)
plot(SECCa$H2O, Model.resid)
qplot(Frag, Model.resid, data = SECCa, facets = Chamber * Position ~ Time ) + theme_bw()

## Plot Residuals: see Zuur et al. 2007, pg. 61-63
coplot( Model.resid ~ X | Frag * Position, data=SECCa, 
       pch=SECCa$pt, col=SECCa$colr,
       ylab = "Residuals"
)
## much the same breakdown, using TRELLIS xyplot over 3 factors.
print( xyplot( Model.resid ~ X | Chamber * Frag * Position, data=SECCa, 
       pch=point, col=SECCa$colr, bg = SECCa$fill,
       ylab = "Residuals"
) )

## Plot Residuals: see Zuur et al. 2007, pg. 131-133
plot(Y.model$fitted[,2],resid(Y.model,type="p"),xlab="Fitted values",ylab="Residuals")
qqnorm(resid(Y.model,type="p"),cex=1,main="")



##################################################
## ANALYSIS: GET RESULTS
##################################################
anova(Y.model)
summary(Y.model)

}




################################################################
## GLMM - Hierarchical / Multilevel Mixed Model               **
################################################################
##==============================================================
## Model Fitting
##==============================================================
## Using a mixed-effects model to account for treatments of different sizes: 
##  Chamber / Frag / Position
lmd <- lmeControl()                    # save defaults
lmc <- lmeControl(niterEM = 5000, msMaxIter = 1000) # takes a while to converge...
Y.rim  <- lme(Y.formula, data=SECCa, random=Y.Ri , control=lmd)
## Y.rism <- lme(Y.formula, data=SECCa, random=Y.Ris, control=lmc) # R crashes on this :(

Y.model <- Y.rim

# Variance Components Analysis (variance decomposition)?


##==============================================================
## MODEL SELECTION
##==============================================================



################################################################
## CHECK ASSUMPTIONS: MODEL VALIDATION
################################################################
## Regression (ANOVA): Normality; Homogeneity; Independence; Fixed X*
## Fixed X (Model I): *think about how X data was generated*
##                    - fixed values [I] or large mesurement error [II] ?
## Check Normality: histogram, qqnorm of Residuals
## Check Homogeneity: (standardized) Residuals vs. Fitted / vs. X
## Check Independence: (standardized) Residuals vs. X
## Assess Model fit, specification: Look for patterns in graphs of Residuals
##  - Should be no areas of residuals with similar values, gaps in the cloud, etc.
## Residuals should ideally be spread out equally across all graphs (vs. X / Fitted).


## Standard diagnostic plots 
op <- par(mfrow=c(2,2))	 # panel of figures: 3 rows & 2 columns
plot(Y.model$fitted[,1],resid(Y.model,type="p"),xlab="Fitted values",ylab="Residuals")
plot(Y.model)
par(op)

## Residuals
Model.resid <- resid(Y.model, type="p") #  Zuur et al like to use type="p".  What is this?
## Plot Residuals: see Zuur et al. 2007, pg. 131-133
## par(mfrow=c(1,1))
hist(Model.resid)                      # residuals: Normal distribution?
qqnorm(Model.resid)                    # residuals: Normal distribution?
## Homogeneity, Independence: Pattern in residuals indicates violation
plot(SECCa$X.log, Model.resid)             
plot(SECCa$H2O, Model.resid)
qplot(Chamber, Model.resid, data = SECCa, facets = Frag * Position ~ Time, geom="boxplot" ) + jaw.ggplot()

## Compare model to GA(M)M to check that a linear fit is the most appropriate
## see Zuur et al. (2007, 2009: Appendix)




################################################################
## ANALYSIS: GET RESULTS
################################################################
anova(Y.model)
summary(Y.model)










################################################################
## SAVE OUTPUT
################################################################
if (Save.results == TRUE && is.null(Save.text) == FALSE) {
  capture.output(cat(Save.header), 
				 print(Y.formula),              # model
				 anova(Y.model),                # model summary
				 summary(Y.model),              # model summary
				 cat("\n\n"),                   # for output
				 cat(Save.end),                 # END OUTPUT #
				 file = Save.text
				)
}

if (Save.results == TRUE && is.null(Save.plots) == FALSE && Save.plots != Save.final) dev.off()



################################################################
## FINAL GRAPHICS
################################################################
if (Save.results == TRUE && is.null(Save.final) == FALSE && Save.plots != Save.final) pdf( file = Save.final )

## generate grid to add predicted values to (X-values in all combinations of factors).
Y.pred <- expand.grid(Chamber  = levels(SECCa$Chamber) , 
                      Frag     = levels(SECCa$Frag), 
                      Position = levels(SECCa$Position), 
                      X=seq(0, max(SECCa$X), length.out=100 ) 
                      )
Y.pred$predicted <- predict(Y.model, newdata=Y.pred, type="response" )  # newdata must have same explanatory variable name for predict to work.

if (FALSE) {
  pred.Chamber <- expand.grid(Chamber = levels(SECCa$Chamber) , 
                              X=seq(0, max(SECCa$X), length.out=100 ) 
  )
  pred.Chamber$predicted <- predict(Y.model, newdata=pred.Chamber, type="response" )
  pred.Frag <- expand.grid(Frag = levels(SECCa$Frag) , 
                           X=seq(0, max(SECCa$X), length.out=100 )
  )
  pred.Frag$predicted <- predict(Y.model, newdata=pred.Frag, type="response" )
  pred.Position <- expand.grid(Position = levels(SECCa$Position) , 
                               X=seq(0, max(SECCa$X), length.out=100 ) 
  )
  pred.Position$predicted <- predict(Y.model, newdata=pred.Position, type="response" )
  pred.FxP <- expand.grid(Frag     = levels(SECCa$Frag), 
                          Position = levels(SECCa$Position), 
                          X=seq(0, max(SECCa$X), length.out=100 ) 
                          )
  pred.FxP$predicted <- predict(Y.model, newdata=pred.Position, type="response" )
}


Chamber.map <- plotMap( "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Chamber.map$label <- factor(Chamber.map$label)
point <- 21	# 21 for circles with col & bg ; 16 for solid circles
Chamber.map$pch <- c(21, 16)  # use circles for both treatments

SECCa <- within( SECCa,{
	colr = as.character(ifelse(Chamber == Chamber.map$label[1], 
                               Chamber.map$col[1], 
                               Chamber.map$col[2] 
                               )
    )
	fill = as.character(ifelse(Chamber == Chamber.map$label[1], 
                               Chamber.map$bg[1],
                               Chamber.map$bg[2]
                               )
    )
	pt = ifelse(Chamber == Chamber.map$label[1], 
			Chamber.map$pch[1], 
			Chamber.map$pch[2]
		)
})


par(mfrow=c(1,1))
pred.Y <- with( Y.pred, 
               aggregate(cbind(predicted), list(Chamber = Chamber, X = X), mean)
)  # I should be getting direct predictions, not means of predictions. *****
with( SECCa,{
	# pred | augpred | ?
	plot(X, Y, type="p",
		ylab=Y.plotlab, xlab=X.plotlab,
		pch=pt, col=colr, bg=fill
        )
	lines(predicted ~ X, data=subset(pred.Y, Chamber == "Ambient"), 
          col = Chamber.map$col[1], 
          lty = Chamber.map$lty[1]
          )
	lines(predicted ~ X, data=subset(pred.Y, Chamber == "Full Chamber"), 
          col = as.character(Chamber.map$col[2]), 
          lty = Chamber.map$lty[2]
          )
	legend( "topright", legend=Chamber.map$label, pch=point, col=as.character(Chamber.map$col), pt.bg=as.character(Chamber.map$bg) )
})



##==============================================================
## Plot fitted on observed, by factor?
##==============================================================

## lattice panels
print( xyplot( Y ~ X | Frag * Position , data=SECCa, 
              pch = SECCa$pt, col = SECCa$colr, 
              xlab = quote(X.plotlab), ylab = quote(Y.plotlab), 
              panel = function(..., data, subscripts) {
                panel.xyplot(...)  # regular plot of data points
                Frag.lvl <- unique(SECCa$Frag[subscripts]) # get current factor levels
                Pos.lvl  <- unique(SECCa$Position[subscripts])
                preds    <- Y.pred[which(Y.pred$Frag %in% Frag.lvl 
                                       & Y.pred$Position %in% Pos.lvl), ]
##                  browser()
                for( lvl in levels(preds$Chamber) ) {
                  preds.lvl <- subset(preds, preds$Chamber == lvl)
                  panel.xyplot(preds.lvl$X, preds.lvl$predicted, 
                               type = 'l', 
                               col = Chamber.map$col[Chamber.map$label == lvl]
                               )
                }
              },
              subscripts = T
              )
)




if (FALSE) {
  ## Coplots with linear fits (from Zuur et al. 2007 Chapter 22 R code)
  ## individual lm's within each panel.  Not exactly what I want.
  coplot( Y ~ X | Frag * Position, data=SECCa, 
         pch=SECCa$pt, col=SECCa$colr, # , bg=Chamber.map$bg
         panel = panel.lines2
         )

  ## Plotting: Observed and Fitted from GLMM - from Richard & Zofia's GLMM workshop
  df <- coef( lmList(Y ~ X | Chamber * Position, data=SECCa) )
  cc1 <- as.data.frame(coef(Y.model)$Y)
  names(cc1) <- c("A", "B")
  df <- cbind(df, cc1)
  ff <- fixef(Y.model)

  print( xyplot( Y ~ X | Chamber * Position, data = SECCa, 
                aspect = "xy", layout = c(4,3),
                type = c("g", "p", "r"), coef.list = df[,3:4],
                panel = function(..., coef.list) {
                  panel.xyplot(...)
                  panel.abline(as.numeric( coef.list[packet.number(),] ), 
                               col.line = trellis.par.get("superpose.line")$col[2],
                               lty = trellis.par.get("superpose.line")$lty[2]
                               )
                  panel.abline(fixef(Y.model), 
                               col.line = trellis.par.get("superpose.line")$col[4],
                               lty = trellis.par.get("superpose.line")$lty[4]
                               )
                },
                index.cond = function(x,y) coef(lm(y ~ x))[1],
                xlab = X.plotlab,
                ylab = Y.plotlab,
                key = list(space = "top", columns = 3,
                  text = list(c("Within-subject", "Mixed model", "Population")),
                  lines = list(col = trellis.par.get("superpose.line")$col[c(2:1,4)],
                  lty = trellis.par.get("superpose.line")$lty[c(2:1,4)]
                  ))
                )
  )
}


if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()










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


