################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Data Exploration
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2016-02-12
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
vars.ls   <- c("ARA.m", "Cells.m", "Hcells.m", "H2O")  

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
  pairplot(SECCa[, c('Y', 'Y.log', 'X', 'X.log', 'H2O', 'Time', 'Chamber', 'Frag', 'Position')],
           labels=c(Y.col, paste("log(", Y.col, ")"), 
                    X.col, paste("log(", X.col, ")"), "H2O"
                    , "Block", "Time", "Chamber", "Frag", "Position")
          )
  ## Moisture by Block?
  pairplot(SECCa[, c('H2O', 'Block')])

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
## Cell Composition: which species, cell types?
##==============================================================
library(mgcv)                          # gam()
Cells.scale <- 1000
Cells.axlab  <- SECC.axislab(SECC, "Cells", unit=Cells.scale, big.mark=",")
CB.axlab     <- parse(text=sub("Cyanobacteria", "Total Cyanobacteria", Cells.axlab))
Genus.axlab  <- parse(text=sub("Cyanobacteria", "Genus", Cells.axlab))
Cells.units  <- bquote("" %*%.(format(Cells.scale, big.mark=",")) * 
                       " " * .(attr(SECC, "units" )[["Cells"]]) ) 
HCells.axlab <- SECC.axislab(SECC, "Hcells", unit=Cells.scale, big.mark=",")
Cells.plotMap <- data.frame(Var  =c("Cells", "Hcells", "Stigonema", "Nostoc", "Other.cells"),
                            Label=c("Total Cells", "Total Heterocysts", "Stigonema", "Nostoc", "Other"),
                            shape=c(1, 16, 15, 21, 4),
                            col  =c("#000000", "#000000", "#666666", "#000000", "#333333"),
                            size =c(0, 2, 3, 3, 3),
                            lty  =c(1, 1, 1, 2, 3),
                            lcol =c("#000000", "#000000", "#999999", "#999999", "#999999"),
                            lsize=c(1, 1, 0.5, 0.5, 0.5),
                            stringsAsFactors=FALSE
                            )
row.names(Cells.plotMap) <- Cells.plotMap$Var
Spp.plotMap    <- Cells.plotMap[Cells.plotMap$Var %in% c("Stigonema", "Nostoc"), ]
Hcells.plotMap <- Cells.plotMap[Cells.plotMap$Var %in% c("Hcells", "Stigonema", "Nostoc"), ]
## Hcells.plotMap <- Hcells.plotMap[c(2, 3, 1), ] # change order to match layer order :(

format_scale <- function(x, num.scale=1000, ...)
{
  format(x/num.scale, trim=TRUE)
}

Square.plot <- opts(aspect.ratio = 1, 
                    axis.title.x = theme_text(vjust = 6, size=12),
                    axis.title.y = theme_text(angle=90, size=12))
Spp.points <- list(geom_point(aes(y=Stigonema, group="Stigonema", colour="Stigonema", 
                               shape="Stigonema", size="Stigonema") ),
                geom_point(aes(y=Nostoc, group="Nostoc", colour="Nostoc", 
                               shape="Nostoc", size="Nostoc"))
)
Spp.lines <- list(stat_smooth(aes(y=Stigonema, linetype="Stigonema"), 
                              colour="#999999", size=Spp.plotMap[1, "lsize"], 
                              method="gam", se=FALSE),
                  stat_smooth(aes(y=Nostoc, linetype="Nostoc"),
                              colour="#999999", size=Spp.plotMap[1, "lsize"], 
                              method="gam", se=FALSE)
                  )

Cyano_scaleName <- "Cyanobacteria\ngenus"
if (F) {                               # using melted df & group variable?
  Cells.df <- melt(SECCa, id=c("SampleID", "Cells"), 
                   measure=c("Hcells", "Stigonema", "Stigonema.H", "Nostoc", "Nostoc.H", 
                             "Other.cells") )
  Cells.df$variable <- as.character(Cells.df$variable)
  ##   Spp.plotMap    <- Spp.plotMap[c(2, 1), ]

  Cells.plot <- ggplot(data=subset(Cells.df, variable %in% c("Stigonema", "Nostoc")), 
                       aes(x=Cells, y=value, group=variable)) +
                   stat_smooth(aes(group=variable, linetype=variable), colour="#999999",
                               method="gam", se=FALSE) +
                      geom_point(aes(group=variable, 
                                     shape=variable, colour=variable, size=variable) ) +
                      scale_x_continuous(expand=c(0.01,0)) + 
                      scale_y_continuous(expand=c(0.01,0)) +
                   scale_colour_manual(name=Cyano_scaleName, 
                                       values=structure(Spp.plotMap[, "col"], 
                                                        names=Spp.plotMap$Var), 
                                       breaks=Spp.plotMap[, "Var"], 
                                       labels=Spp.plotMap[, "Label"]) +
                   scale_size_manual(name=Cyano_scaleName, 
                                       values=structure(Spp.plotMap[, "size"], 
                                                        names=Spp.plotMap$Var), 
                                       breaks=Spp.plotMap[, "Var"], 
                                       labels=Spp.plotMap[, "Label"]) +
                   scale_shape_manual(name=Cyano_scaleName, 
                                       values=structure(Spp.plotMap[, "shape"], 
                                                        names=Spp.plotMap$Var), 
                                       breaks=Spp.plotMap[, "Var"], 
                                       labels=Spp.plotMap[, "Label"]) +
                   scale_linetype_manual(name=Cyano_scaleName, 
                                       values=structure(Spp.plotMap[, "lty"], 
                                                        names=Spp.plotMap$Var), 
                                       breaks=Spp.plotMap[, "Var"], 
                                       labels=Spp.plotMap[, "Label"])

  ## sort of works, but I still can't control the order in which layers are added :(
  ## structure() allows names to be assigned to vector of values,
  ## so that scale values are assigned by group value, not order. :)
}

Cells.plot <- ggplot(data=SECCa, aes(x=Cells, y=Cells)) +
                geom_abline(intercept=0, slope=1, aes(group="Cells"), 
                            colour="#333333", size=1) +
                xlab(CB.axlab) + ylab(Genus.axlab) +
                scale_x_continuous(expand=c(0.01,0), formatter="format_scale", num.scale=Cells.scale) + 
                scale_y_continuous(expand=c(0.01,0), formatter="format_scale", num.scale=Cells.scale) +
                scale_colour_manual(name=Cyano_scaleName,
                                    values=structure(Spp.plotMap[, "col"], 
                                                     names=Spp.plotMap$Var), 
                                    breaks=Spp.plotMap[, "Var"], 
                                    labels=Spp.plotMap[, "Label"]) +
                 scale_shape_manual(name=Cyano_scaleName,
                                    values=structure(Spp.plotMap[, "shape"], 
                                                     names=Spp.plotMap$Var), 
                                    breaks=Spp.plotMap[, "Var"], 
                                    labels=Spp.plotMap[, "Label"]) +
                 scale_size_manual(name=Cyano_scaleName,
                                   values=structure(Spp.plotMap[, "size"], 
                                                    names=Spp.plotMap$Var), 
                                   breaks=Spp.plotMap[, "Var"], 
                                   labels=Spp.plotMap[, "Label"]) +
                 scale_linetype_manual(name=Cyano_scaleName,
                                       values=structure(Spp.plotMap[, "lty"], 
                                                        names=Spp.plotMap$Var), 
                                       breaks=Spp.plotMap[, "Var"], 
                                       labels=Spp.plotMap[, "Label"]) +
                coord_equal(ratio = 1) 
## log-transforming axes doesn't agree with gam - too many 0s?
Cells.logplot <- Cells.plot + xlab(parse(text=sub("Cyanobacteria", "Total Cyanobacteria", SECC.axislab(SECC, "Cells")))) + 
    ylab(parse(text=sub("Cyanobacteria", "Genus", SECC.axislab(SECC, "Cells")))) +
    scale_x_log10() + scale_y_log10() + Spp.points + jaw.ggplot() + Square.plot
Cells.plot <- Cells.plot + Spp.lines + Spp.points + jaw.ggplot() + Square.plot 
## jaw.ggplot() won't work in large call above

## How can I reliably produce different shapes & colours for lines, in the legend?
HCells.plot <- ggplot(data=SECCa, aes(x=Cells, y=Hcells))
if (T) {
  HCells.plot <- HCells.plot +
                    stat_smooth(aes(y=Stigonema.H), linetype=Hcells.plotMap[2, "lty"],
                                colour=Hcells.plotMap[2, "lcol"], 
                                size=Hcells.plotMap[2, "lsize"], 
                                method="gam", se=FALSE) +
                   stat_smooth(aes(y=Nostoc.H), linetype=Hcells.plotMap[3, "lty"],
                               colour=Hcells.plotMap[3, "lcol"], 
                               size=Hcells.plotMap[3, "lsize"], 
                               method="gam", se=FALSE) +
                   stat_smooth(aes(y=Hcells), linetype=Hcells.plotMap[1, "lty"],
                               colour=Hcells.plotMap[1, "lcol"], 
                               size=Hcells.plotMap[1, "lsize"], 
                               method="gam", se=FALSE)
} else {
  HCells.plot <- HCells.plot +
                    stat_smooth(aes(y=Stigonema.H, linetype="Stigonema",
                                colour="Stigonema", 
                                size="Stigonema"), 
                                method="gam", se=FALSE) +
                   stat_smooth(aes(y=Nostoc.H, linetype="Nostoc",
                                colour="Nostoc", 
                                size="Nostoc"), 
                               method="gam", se=FALSE) +
                   stat_smooth(aes(y=Hcells, linetype="Hcells",
                                colour="Hcells", 
                                size="Hcells"), 
                               method="gam", se=FALSE) +
                   scale_linetype_manual(name="Fitted vs.\nTotal Cells",
                                         values=structure(Hcells.plotMap[, "lty"], 
                                                          names=Hcells.plotMap$Var), 
                                         breaks=Hcells.plotMap[, "Var"], 
                                         labels=Hcells.plotMap[, "Label"]) +
                   scale_colour_manual(name="Fitted vs.\nTotal Cells",
                                       values=structure(Hcells.plotMap[, "lcol"], 
                                                        names=Hcells.plotMap$Var), 
                                       breaks=Hcells.plotMap[, "Var"], 
                                       labels=Hcells.plotMap[, "Label"]) +
                   scale_size_manual(name="Fitted vs.\nTotal Cells",
                                     values=structure(Hcells.plotMap[, "lsize"], 
                                                      names=Hcells.plotMap$Var), 
                                     breaks=Hcells.plotMap[, "Var"], 
                                     labels=Hcells.plotMap[, "Label"])
}
HCells.plot <- HCells.plot +
geom_point(aes(y=Stigonema.H, colour="Stigonema", 
               shape="Stigonema", size="Stigonema") ) +
                   geom_point(aes(y=Nostoc.H, colour="Nostoc", 
                                  shape="Nostoc", size="Nostoc") ) +
                   geom_point(aes(y=Hcells, colour="Hcells", 
                                  shape="Hcells", size="Hcells") ) +
                   xlab(CB.axlab) + ylab(HCells.axlab) +
                   scale_x_continuous(expand=c(0.01,0), formatter="format_scale", num.scale=Cells.scale) + 
                   scale_y_continuous(expand=c(0.01,0), formatter="format_scale", num.scale=Cells.scale) +
                   scale_shape_manual(name=Cyano_scaleName,
                                      values=structure(Hcells.plotMap[, "shape"], 
                                                       names=Hcells.plotMap$Var), 
                                      breaks=Hcells.plotMap[, "Var"], 
                                      labels=Hcells.plotMap[, "Label"]) +
                   scale_colour_manual(name=Cyano_scaleName,
                                       values=structure(Hcells.plotMap[, "col"], 
                                                        names=Hcells.plotMap$Var), 
                                       breaks=Hcells.plotMap[, "Var"], 
                                       labels=Hcells.plotMap[, "Label"]) +
                   scale_size_manual(name=Cyano_scaleName,
                                     values=structure(Hcells.plotMap[, "size"], 
                                                      names=Hcells.plotMap$Var), 
                                     breaks=Hcells.plotMap[, "Var"], 
                                     labels=Hcells.plotMap[, "Label"])
HCells.plot <- HCells.plot + jaw.ggplot() + Square.plot

if (DrawExplorationGraphs) {
  print(Cells.logplot)
  print(Cells.plot)
  print(HCells.plot)
}

## manual scale only seems to work because elements are added in the same *order*: not matched by *values*?
if (F) {
  ## from ggplot2 documentation http://had.co.nz/ggplot2/book/scales.r
  huron <- data.frame(year = 1875:1972, level = LakeHuron)
  ## legend values matched by value, not order
  ## lines added in arbitrary order, total control over legend
  ggplot(huron, aes(year)) +
  geom_line(aes(y = level - 5, colour = "below")) + 
  geom_line(aes(y = level + 5, colour = "above")) + 
  scale_colour_manual("Direction", values=c("above" = "red", "below" = "blue"), 
                      breaks=c("above", "below"), labels=c("Above", "Below"))

  ## legend values matched by order, not value.
  ggplot(huron, aes(year)) +
  geom_line(aes(y = level - 5, colour = "below")) + 
  geom_line(aes(y = level + 5, colour = "above")) + 
  scale_colour_manual("Direction", values=c("red", "blue"), 
                      breaks=c("above", "below"), labels=c("Above", "Below"))
  ## legend values are matched by value only if the 'values' argument to the scale() function is *named* (c(name="value"))
  ## how can I quickly extract a df column, with names == row.names?
  structure(Cells.plotMap$col, names=Cells.plotMap$Var)
}




##==============================================================
## Scatterplots: The easy way, with ggplot2
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

## log-y may be the best linear model (according to AIC), but essentially implies an exponential relationship (!)
ARAcb.logy <- qplot(X, Y, data = SECCa, log = "y", 
                    group = Chamber, geom = "point", size = I(3),
                    colour = Chamber, shape = Chamber,
                    xlab = X.plotlab,
                    ylab = Y.plotlab
                    )
ARAcb.logy <- ARAcb.logy + ChamberPts + jaw.ggplot() + TopLegend

ARAcb.time.logy   <- ARAcb.logy + Time.facets 
ARAcb.panels.logy <- ARAcb.logy + All.facets 
if (DrawExplorationGraphs) {
  print(ARAcb.logy)
  print(ARAcb.time.logy)
  print(ARAcb.panels.logy)
}

## log-log looks encouraging, and may be theoretically ($ Biologically) justified.
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

## sqrt-transformation on cells?
ARAcb.sqrt <- qplot(sqrt(X), Y, data = SECCa, log = "y", 
                    group = Chamber, geom = "point", size = I(3),
                    colour = Chamber, shape = Chamber,
                    xlab = bquote(sqrt(.(X.plotlab))),
                    ylab = Y.plotlab
                    )
ARAcb.sqrt <- ARAcb.sqrt + ChamberPts + jaw.ggplot() + TopLegend

ARAcb.time.sqrt   <- ARAcb.sqrt + Time.facets 
ARAcb.panels.sqrt <- ARAcb.sqrt + All.facets 
if (DrawExplorationGraphs) {
  print(ARAcb.sqrt)
  print(ARAcb.time.sqrt)
  print(ARAcb.panels.sqrt)
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

if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()
