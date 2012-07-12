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
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("..")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place

  ## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
  source('./Nfix-Cyanobacteria/1_ARA-cb_setup.R')
  Save.results <- FALSE
}

library(lattice)                       # densityplot
library(ggplot2)
theme_set(theme_bw())                  # change global ggplot2 theme
library(rgl)                           # 3D plots
library(car)                           # diagnostic plots & tools


################################################################
## CUSTOM OPTIONS
################################################################
# explanatory vars for data exploration (and labels)
vars.ls   <- c("Nfix", "Cells.m", "Hcells.m", "H2O")  

labels.ls <- c()
for(i in 1:length(vars.ls) ){
  var <- vars.ls[i]
  labels.ls[i] <- attr(SECC, "labels")[[var]]
}

DrawExplorationGraphs <- TRUE # Save.results  # Set to FALSE to suppress all this output
if (Save.results == TRUE && is.null(Save.plots) == FALSE) {
  pdf( file = gsub("Results", "Exploration", Save.plots, fixed=TRUE) )
}




################################################################
## EXPLORE: PLOTS
################################################################
### Make some meaningful plots of data to check for predicted (expected) patterns.
### Map of point styles for Chamber treatments
Chamber.map <- plotMap( "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Chamber.map$label <- factor(Chamber.map$label)
point <- 21	# 21 for circles with col & bg ; 16 for solid circles
Chamber.map$pch <- c(21, 16)  # use circles for both treatments
Chamber.label <- attr(SECC, "labels")[["Chamber"]]

## Set up some plotting options for each data point
##  only needed for non-ggplot commands
SECCa <- within( SECCa,
                {
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

##==============================================================
## Which Variables to use?
##==============================================================
### pairplots of several variables of interest: check for collinearity, patterns, etc.
### see Zuur et al. books
if (DrawExplorationGraphs) {
  ## full dataset: unbalanced with respect to experimental treatments
  pairplot(SECC[, c("ARA.m", "ARA.g", "H2O", "Cells.m", "Cells.g", "Hcells.m", "Hcells.g", "Stigonema", "Nostoc" )])
  ## filtered dataset: balanced, but am I missing useful information about continuous explanatory variables (H2O, cells)?
  pairplot(SECCa[, c("ARA.m", "ARA.g", "H2O", "Cells.m", "Cells.g", "Hcells.m", "Hcells.g", "Stigonema", "Nostoc" )])
  ## look at log transformations & check for colinearity among explanatory variables
  pairplot(SECCa[, c('Y', 'Y.log', 'H2O', 'Block', 'TempC','Chamber', 'Frag', 'Position')],
           labels=c(Y.col, paste("log(", Y.col, ")"), 
                    "H2O", "Block", "Temperature", "Chamber", "Frag", "Position")
          )
  pairplot(SECCa[, c('Y', 'Y.log', 'Block', 'Frag', 'TempC', 'H2O', 'Cells.m', 
                     'grow2', 'Decomposition', 'TAN', 'Patch.dwt')],
           labels=c(Y.col, paste("log(", Y.col, ")"), 
                    "Block", "Frag", "Temperature", "H2O", "Cells.m", 
                    "Moss growth", "Decomposition", "TAN", "Patch dwt.")
          )
  ## Moisture by Block?
  pairplot(SECCa[, c('H2O', 'Block')])

  ## Cleveland Dotplots (Zuur et al. 2009)
  Dotplot.y <- "Order of observations"
  op <- par(mfcol=c(2,2))
  dotchart(SECCa$Y, ylab=Dotplot.y, xlab=Y.label)
  dotchart(SECCa$Y.log, ylab=Dotplot.y, xlab=paste("log(", Y.label, ")"))
  dotchart(SECCa$Cells.m, ylab=Dotplot.y, xlab=attr(SECCa, "labels")[['Cells.m']]) # outliers > 5e+09
  dotchart(SECCa$H2O, ylab=Dotplot.y, xlab=attr(SECCa, "labels")[['H2O']]) # outliers > 800?
  dotchart(SECCa$Richness, ylab=Dotplot.y, xlab=attr(SECCa, "labels")[['Richness']])
  dotchart(SECCa$grow2, ylab=Dotplot.y, xlab="Moss growth - year 2") # outliers >30??
  dotchart(SECCa$Decomposition, ylab=Dotplot.y, xlab="Moss decomposition - year 2") # outliers >0.3??
  dotchart(SECCa$TAN, ylab=Dotplot.y, xlab="Total Available N") # outliers >0.3??
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
## Scatterplots: The easy way, with ggplot2
##==============================================================
### prepare plot theme for points varying by Chamber
ChamberPts  <- ggPts.SECC(Chamber.map, Chamber.label) 
TopLegend   <- opts(legend.position = "top", legend.direction = "horizontal")
Time.facets <- facet_grid(facets = .~Time)
All.facets  <- facet_grid(facets = Frag~Position)

Ycb.plot  <- ggplot(SECCa, aes(Cells.m, Y, group = Chamber) ) + 
                    geom_point(aes(colour = Chamber, shape = Chamber), size = 3) +
                    xlab( make.axis.title(SECCa, "Cells.m") ) +
                    ylab( Y.plotlab )
Ycb.plot  <- Ycb.plot + jaw.ggplot() + ChamberPts + TopLegend  # order matters!

Ycb.panels <- Ycb.plot + All.facets # Full Faceting ***
if (DrawExplorationGraphs) {
  print(Ycb.plot)
  print(Ycb.panels)
}

## log-y may be the best linear model (according to AIC), but essentially implies an exponential relationship (!)
Ycb.logy <- Ycb.plot + scale_y_log10()

Ycb.time.logy   <- Ycb.logy + Time.facets 
Ycb.panels.logy <- Ycb.logy + All.facets 
if (DrawExplorationGraphs) {
  print(Ycb.logy)
  print(Ycb.time.logy)
  print(Ycb.panels.logy)
}

## log-log looks encouraging, and may be theoretically (Biologically) justified.
Ycb.log <- Ycb.logy + scale_x_log10()

Ycb.panels.log <- Ycb.log + All.facets 
if (DrawExplorationGraphs) {
  print(Ycb.log)
  print(Ycb.panels.log)
}

## sqrt-transformation on cells?
Ycb.sqrt <- ggplot(SECCa, aes(sqrt(Cells.m), Y, group = Chamber) ) + 
                    geom_point(aes(colour = Chamber, shape = Chamber), size = 3) +
                    xlab( bquote( sqrt( .(make.axis.title(SECCa, "Cells.m")) ) ) ) +
                    ylab( Y.plotlab )
Ycb.sqrt <- Ycb.sqrt + ChamberPts + jaw.ggplot() + TopLegend

Ycb.panels.sqrt <- Ycb.sqrt + All.facets 
if (DrawExplorationGraphs) {
  print(Ycb.sqrt)
  print(Ycb.panels.sqrt)
}

## H2O plots
Y.H2O <- ggplot(SECCa, aes(H2O, Y, group = Chamber) ) + 
                    geom_point(aes(colour = Chamber, shape = Chamber), size = 3) +
                    xlab( make.axis.title(SECCa, "H2O") ) +
                    ylab( Y.plotlab )
Y.H2O <- Y.H2O + ChamberPts + jaw.ggplot() + TopLegend

Y.H2O.panels <- Y.H2O + All.facets
if (DrawExplorationGraphs) {
  print(Y.H2O)
  print(Y.H2O.panels)
}

Y.H2O.log <- Y.H2O + scale_y_log10()
Y.H2O.log.panels <- Y.H2O.log + All.facets
if (DrawExplorationGraphs) {
  print(Y.H2O.log)
  print(Y.H2O.log.panels)
}





##==============================================================
## Check Variation, Ranges
## with(SECCa, boxplot(Y.log ~ Block, ylab = Y.plotlab ) )    # bquote doesn't work in base graphics?  Need a different format for expression?
Y.logBlock <- qplot(Block, Y.log, data = SECCa, geom = "boxplot",
                 ylab = Y.plotlab) + jaw.ggplot()
Y.Block <- qplot(Block, Y, data = SECCa, geom = "boxplot",
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
  print(Y.Chamber)
  print(Y.Frag)
  print(Y.Pos)
}

##==============================================================
## Check distributions
X.dist <- qplot(Cells.m, data = SECCa, geom = "histogram", 
                xlab = make.axis.title(SECCa, "Cells.m"))
X.dist <- X.dist + jaw.ggplot() + All.facets
if (DrawExplorationGraphs) print(X.dist)

Y.dist <- qplot(Y, data = SECCa, geom = "histogram",
                xlab = Y.plotlab)
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

if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()





################################################################
## MULTIVARIATE EXPLORATION
################################################################
library(vegan)
SECCmv <- SECCa[, c("Y.trans", "Block", "Warming", "Frag", "H2O", "grow2", "Decomposition", "TAN")]
rownames(SECCmv) <- SECCa$SampleID
SECCmv$Block <- as.numeric(as.character(SECCmv$Block))
SECCmv$Frag <- as.numeric(SECCmv$Frag)
SECC.trans <- decostand(SECCmv, method = "standardize")
boxplot(SECC.trans)

SECC.pca <- rda(na.omit(SECC.trans))
SECC.pca
if (F)
  summary(SECC.pca)

# extract eigenvalues, and apply Kaiser-Guttman criterion to select axes
EV <- SECC.pca$CA$eig
EV[EV > mean(EV)]

# Broken-stick model
N <- length(EV)
BSM <- data.frame(j = seq(1:N), p = 0)
BSM$p[1] <- 1/N
for (i in 2:N)
{
  BSM$p[i] = BSM$p[i - 1] + (1/(N + 1 - i))
}
BSM$p <- BSM$p/N * 100
BSM

# Plot eigenvalues
op <- par(mfrow = c(2,1))
barplot(EV, main = "Eigenvalues", las = 2)
abline(h = mean(EV), col = "red")      # average eigenvalue
legend("topright", "Average eigenvalue", lwd = 1, col = "red", bty = "n")
barplot( t(cbind(100*EV/sum(EV), BSM$p[N:1])), beside = TRUE, 
        main = "% variance", col = c("bisque", "darkblue"), las = 2)
legend("topright", c("% eigenvalue", "Broken stick model"),
       pch = 15, col = c("bisque", "darkblue"), bty = "n")
par(op)

# Plot PCA
op <- par(mfrow = c(1,2))
biplot(SECC.pca, scaling = 1, type = c("text", "points"), main = "PCA - scaling 1 (distance)")
biplot(SECC.pca, scaling = 2, type = c("text", "points"), main = "PCA - scaling 2 (correlation)")
par(op)
