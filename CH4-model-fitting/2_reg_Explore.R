################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Data Exploration
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2012-07-15
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
## vars.ls   <- c("Nfix", "Cells.m", "Hcells.m", "H2O")  
vars.ls   <- c(Y.col, X.cols)          # defined in setup

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
  ## look at log transformations & check for colinearity among explanatory variables
  pairplot(SECCa[, c('Y', 'Y.log', 'H2O', 'Block', 'TempC','Chamber', 'Frag', 'Position')],
           labels=c(Y.col, paste("log(", Y.col, ")"), 
                    "H2O", "Block", "Temperature", "Chamber", "Frag", "Position")
          )
  pairplot(SECCa[, c('Y', 'Y.log', 'Cells.m',
                     'Growth', 'Decomposition', 'TAN', 'Patch.dwt',
                     'H2O', 'TempC', 'Frag', 'Block' )],
           labels=c(Y.col, paste("log(", Y.col, ")"), "Cells.m", 
                    "Moss growth", "Decomposition", "TAN", "Patch dwt.",
                    "H2O", "Temperature", "Frag", "Block" )
          )

  ## Cleveland Dotplots (Zuur et al. 2009)
  Dotplot.y <- "Order of observations"
  op <- par(mfrow=c(2,2))
  dotchart(SECCa$Y, ylab=Dotplot.y, xlab=Y.label)
  dotchart(SECCa$Y.log, ylab=Dotplot.y, xlab=paste("log(", Y.label, ")"))
  dotchart(SECCa$Nfix, ylab=Dotplot.y, xlab=make.axis.title(SECCa, "Nfix"))
  dotchart(SECCa$Cells.m, ylab=Dotplot.y, xlab=attr(SECCa, "labels")[['Cells.m']]) # outliers > 5e+09
  dotchart(SECCa$H2O, ylab=Dotplot.y, xlab=attr(SECCa, "labels")[['H2O']]) # outliers > 800?
  dotchart(SECCa$Growth, ylab=Dotplot.y, xlab="Moss growth - year 2") # outliers >30??
  dotchart(SECCa$Decomposition, ylab=Dotplot.y, xlab="Moss decomposition - year 2") # outliers >0.3??
  dotchart(SECCa$TAN, ylab=Dotplot.y, xlab="Total Available N")
  dotchart(SECCa$Richness, ylab=Dotplot.y, xlab=attr(SECCa, "labels")[['Richness']])
  dotchart(SECCa$Evenness, ylab=Dotplot.y, xlab=attr(SECCa, "labels")[['Evenness']])
  par(op)
}

if (FALSE) 
{ ## the old-fashioned way (low-level)
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

for (Xcol in X.cols)
{
  X.plotlab <- SECC.axislab(SECCa, Xcol)
  ## using get() inside a ggplot call seems risky, but it's the best I can think of for now.
  YX.plot  <- ggplot(SECCa, aes(get(Xcol), Y, group = Chamber) ) + 
  geom_point(aes(colour = Chamber, shape = Chamber), size = 3) +
  xlab( X.plotlab ) +
  ylab( Y.plotlab )
  YX.plot  <- YX.plot + jaw.ggplot() + ChamberPts + TopLegend  # order matters!

  YX.panels <- YX.plot + All.facets # Full Faceting ***
  if (DrawExplorationGraphs) {
    print(YX.plot)
    print(YX.panels)
  }

  ## log-y may be the best linear model (according to AIC), but essentially implies an exponential relationship (!)
  YX.logy <- YX.plot + scale_y_log10()
  YX.panels.logy <- YX.logy + All.facets 
  if (DrawExplorationGraphs) {
    print(YX.logy)
    print(YX.panels.logy)
  }

  ## log-log looks encouraging, and may be theoretically (Biologically) justified.
  YX.log <- YX.logy + scale_x_log10()
  YX.panels.log <- YX.log + All.facets 
  if (DrawExplorationGraphs) {
    print(YX.log)
    print(YX.panels.log)
  }

  ## sqrt-transformation?
  ##   X.plotlab <- SECC.axislab(SECCa, Xcol, trans = "sqrt")   # not quite working as expected
  YX.sqrt <- ggplot(SECCa, aes(sqrt(get(Xcol)), Y, group = Chamber) ) + 
  geom_point(aes(colour = Chamber, shape = Chamber), size = 3) +
  ##   xlab( X.plotlab ) +
  xlab( SECC.axislab(SECCa, Xcol, trans="sqrt") ) +
  ylab( Y.plotlab )
  YX.sqrt <- YX.sqrt + ChamberPts + jaw.ggplot() + TopLegend
  YX.panels.sqrt <- YX.sqrt + All.facets 
  if (DrawExplorationGraphs) {
    print(YX.sqrt)
    print(YX.panels.sqrt)
  }

}


if (FALSE)
{
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
}




##==============================================================
## Check Variation, Ranges
## with(SECCa, boxplot(Y.log ~ Block, ylab = Y.plotlab ) )    # bquote doesn't work in base graphics?  Need a different format for expression?
Y.translab <- SECC.axislab(SECC = SECCa, col = Y.col, trans = "log", parens = FALSE) 
Y.Block <- qplot(Block, Y, data = SECCa, geom = "boxplot",
                 ylab = Y.plotlab) + jaw.ggplot()
Y.logBlock <- qplot(Block, Y.trans, data = SECCa, geom = "boxplot",
                 ylab = Y.translab) + jaw.ggplot()
Y.Chamber <- qplot(Chamber, Y.trans, data = SECCa, geom = "boxplot",
                   ylab = Y.translab) + jaw.ggplot() 
Y.Frag <- qplot(Frag, Y.trans, data = SECCa, geom = "boxplot",
                ylab = Y.translab) + jaw.ggplot() 
Y.Pos  <- qplot(Position, Y.trans, data = SECCa, geom = "boxplot",
                ylab = Y.translab) + jaw.ggplot() 
if (DrawExplorationGraphs) {
  print(Y.Block)
  print(Y.logBlock)
  print(Y.Chamber)
  print(Y.Frag)
  print(Y.Pos)
}

##==============================================================
## Check distributions
Y.dist <- qplot(Y, data = SECCa, geom = "histogram",
                xlab = Y.plotlab)
Y.dist <- Y.dist + jaw.ggplot() + All.facets
if (DrawExplorationGraphs) print(Y.dist)

if (DrawExplorationGraphs) {
  old.par <- par(mfcol=c(2,2))
  for(i in 1:length(vars.ls) ) {
    var <- vars.ls[i]
    label <- labels.ls[i]
    with( SECCa,
         {
           X.var <- get(var)
           Xna   <- !is.na(X.var)
           X.var <- na.omit(X.var)
           X.max  <- max( X.var )
           ##         cat(var, " ", X.max, "\n")
           freq.max <- length(X.var)/2
           X.maxD <- max( density( X.var )$y )*1.5
           for(Ch.trt in levels(Chamber)) {
             X.trt <- X.var[Chamber[Xna]==Ch.trt]
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
                         xlim=c(0, max( X.var ) ),
                         na.rm = TRUE
                         )	# ** TRELLIS plot
      }
    })
    ## qqplot(x, y) to compare distributions.
  }
  par(old.par)
}






################################################################
## MULTIVARIATE EXPLORATION
################################################################
library(vegan)
SECCmv <- SECCa[, c("Y.trans", "Block", "Warming", "Frag", "H2O", "Growth", "Decomposition", "TAN")]
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
cleanplot.pca(SECC.pca, point = TRUE)





if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()
