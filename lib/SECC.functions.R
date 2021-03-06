##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# standard functions used in most analyses
# Jonathan Whiteley      R v2.12      2012-11-04
##################################################
## access these functions in another file by using: 
## 	 source("path/to/this/file.R")

##================================================
## Calculations
##================================================
ara2nfix <- function (ara, 
                      d2yr = 4*30,     # growing season: # days / year
                      A2N = 3          # default 3 moles of Acetylene reduced for every mole of N2
                      )
{
  ## convert ARA in umol / m^2 / day to N-fixation in kg / ha / year
  Nmwt   <- 14.07                      # g / mol / atom of N
  Natoms <- 2                          # N atoms / molecule of N2
  u2mol  <- 1000 * 1000                # umol / mol
  g2kg   <- 1000                       # g / kg
  m2ha   <- 100 * 100                  # m2 / ha
  nfix <- (ara/A2N) * ( (Nmwt * Natoms * m2ha * d2yr) / (u2mol * g2kg) )
  nfix
}

nfix2ara <- function (nfix, 
                      d2yr = 4*30,     # growing season: # days / year
                      A2N = 3          # default 3 moles of Acetylene reduced for every mole of N2
                      )
{
  ## convert N-fixation in kg / ha / year to ARA in umol / m^2 / day 
  Nmwt   <- 14.07                      # g / mol / atom of N
  Natoms <- 2                          # N atoms / molecule of N2
  u2mol  <- 1000 * 1000                # umol / mol
  g2kg   <- 1000                       # g / kg
  m2ha   <- 100 * 100                  # m2 / ha
  ara <- (nfix * A2N) / ( (Nmwt * Natoms * m2ha * d2yr) / (u2mol * g2kg) )
  ara
}


## Conversion factors (constants)
sampleA  <- 6   # sample Area, in cm^2:  pi * (2.75/2)^2 ; pi * (2.8 / 2)^2
      #     6 for rough estimate of inner tube diameter (2.8 cm): pi*(2.8/2)^2,
      #  or 6.4 for 20 shoots, based on density survey.
sample.to.m2 <- (100*100)/sampleA   # scale sample area, in cm^2 to m^2
sample_ml    <- 50  # 50 ml sample
ARA.m2   <- sampleA/(100*100)  # ARA sample area,   in (cm^2 to) m^2
patchA   <- pi * ((12.5 / 2)^2)      # patch area, in cm^2 (12.5 cm diameter patch)
patch.m2 <- patchA/(100*100)   # patch sample area, in (cm^2 to) m^2
Nfix.ARA.ratio <- 1/3  # ratio of N-fixation : ARA.



##################################################
## ANOVA: standardized analysis of individual response variables
##################################################
# Standard univariate Analysis of Variance on NESTED treatments: 
# No fancy mixed effects or error distribution families (GLMM) yet.
# The same experimental design applies for each response variable, 
# therefore a standard analytical procedure applies.
# Major differences include details such as:
#   - which Transformation is appropriate
#   - Graphs of significant Main Effects

##================================================
## DATA PROCESSING
##================================================
# See `/lib/fun.R` for data processing functions
# used to produce a data frame with standard columns
# for standardized analysis steps.

SECC_sampleID <- function( factors=NULL ) {
  ## Auto-concatenate factor codes into standard SampleIDs.
  ## the supplied data frame had BETTER contain all necessary columns.
  factors.names <- names(factors)
  factors.present <- c("Chamber", "Time", "Frag", "Pos", "Position") %in% names(factors)
  factors.present[4] <- any( factors.present[c(4,5)] )
  factors.present <- factors.present[1:4]
  if ( all( factors.present ) == FALSE ) {
    stop( paste("ID columns missing:", 
                factors.names[which(factors.present==FALSE)]
    ))
  }
  
  pos.sep <- "." # default
  if( "Pos" %in% names(factors) ) {
    pos.col <- factors$Pos
  } else if( "Position" %in% names(factors) ) {
    pos.col <- factors$Position  # would I ever really do this?
  } else {
    pos.sep <- ""
    pos.col <- ""
  }

  SampleIDs <-  with(factors,
                     paste(Block, Time, Chamber,
                           "-", Frag, pos.sep, pos.col,
                           sep = ""
                           )
                     )

  return(SampleIDs)
}

SECC_aggregate <- function( dataf = NULL , trt = "Frag", drop.trt = FALSE, FUN = mean, ... ) {
  ## aggregate full SECC data to means for each level of specified treatment
  ## - default: Frag (for Meta-Community scale data)
  ## optional: specify a drop argument to aggregate across all other levels, 
  ##		   except for the dropped level(s)
  
  if (is.null(drop.trt)) {
	dropt.trt <- FALSE
  } else {
	if (drop.trt==TRUE) {
	  drop.trt <- trt                     # 'trt' will be dropped
	} else if (!is.logical(drop.trt)) {
	  drop.trt -> trt                     # use 'drop.trt' for 'trt'
	  drop.trt <- TRUE
	}
  }
  data.cols <- names(dataf)
  col.types <- lapply(dataf, class)
  agg.cols  <- names(dataf)[(col.types=="numeric" | col.types=="integer")]  # not sure if integer is necessary, but ...
  agg.by    <- with(dataf, list(Block    = Block,
								Time     = Time,
								Chamber  = Chamber,
								Frag     = Frag,
								Position = Position
								)
                    )
  if (trt %in% names(agg.by) == FALSE) {
    stop( paste("Treatment \'", trt, "\' is not a valid column.  Choose one of:\n",
                names(agg.by),
                sep = "" )
         )
  }
  agg.lvl <- match(TRUE, names(agg.by) %in% trt ) # == or %in%? %in% allows vector of trts
  if (drop.trt==FALSE) {
	agg.by <- agg.by[1:agg.lvl[1]]        # use only the first value, if more than one (avoids Warning messages).
  } else {
	agg.by <- agg.by[-agg.lvl]
  }
  
  data.agg <- with(dataf,
                  aggregate(dataf[, agg.cols], by=agg.by, FUN, na.rm = TRUE, ... )
                  )
  ## note that the mean of only NAs is NaN.
  return(data.agg)
}


SECCclean <- function(data=NULL,
                      Time.lvls     = 1:3,
                      Chamber.lvls  = c("A", "B", "C"),
                      Frag.lvls     = 1:4,
                      Position.lvls = c("Inner", "intermediate", "Outer")
                      )
{
  ## clean empty and unused data from SECC data frame

  ## strip empty rows (rows with only NAs)
  SECC.lvls <- strip_empty_dims( data, dim = 1, col.class = "numeric" )  

  ## Filter data by treatment levels.
  SECC.lvls <- SECC.lvls[SECC.lvls$Time     %in% Time.use     &
                         SECC.lvls$Chamber  %in% Chamber.use  &
                         SECC.lvls$Frag     %in% Frag.use     &
                         SECC.lvls$Position %in% Position.use 
                         , ]
  ## drop unused factor levels (for plotting)
  SECC.lvls <- within( SECC.lvls, {
                      Time     <- factor(Time,     levels = Time.lvls)
                      Chamber  <- factor(Chamber,  levels = Chamber.lvls)
                      Frag     <- factor(Frag,     levels = Frag.lvls)
                      Position <- factor(Position, levels = Position.lvls)
                      })
  return(SECC.lvls)
}


##================================================
## DATA EXPLORATION
##================================================



##================================================
## DEFINE MODEL FORMULA
##================================================




##================================================
## CHECK ASSUMPTIONS: residuals, standard diagnostic plots
##================================================
## see Zuur et al. 2007, 2009 books, 
##  and [Quick-R](http://www.statmethods.net/stats/rdiagnostics.html)
## Regression (ANOVA): Normality; Homogeneity; Independence; Fixed X*
## Fixed X (Model I): *think about how X data was generated*
##                    - fixed values [I] or large mesurement error [II] ?
## Check Normality: histogram, qqnorm of Residuals
## Check Homogeneity: (standardized) Residuals vs. Fitted / vs. X
## Check Independence: (standardized) Residuals vs. X
## Assess Model fit, specification: Look for patterns in graphs of Residuals
## - Should be no areas of residuals with similar values, gaps in the cloud, etc.
## Residuals should ideally be spread out equally across all graphs (vs. X / Fitted).

diagnostics <- function(Y.model=NULL, resType="pearson", label=Y.model$call, X.cols = c("X.trans", "H2O"), more=FALSE) 
{
  ## the defaults for X.cols are really only there for the original use of this function: in ARA~Cyanobacteria (More general defaults would be better for re-use).
  require(car)                         # diagnostic plots
  if (more) require(ggplot2)           # facetted plots
  ## Standard diagnostic plots 
  if (FALSE) {                         # buggy
    op <- par(mfrow=c(2,2))	 # panel of figures
    print(plot(Y.model))
    mtext(label, 3, adj=0.5, line=-2, outer=TRUE)
    par(op)
  }

  if (FALSE) {                         # type="p" or type="normalized"?
    ## type = "normalized" residuals 
    ##  Zuur et al. 2009 use this for 'standardized' residuals, but it actually does something more complicated.  see ?residuals.lme
    ##  - in the examples in Zuur et al. 2009, there is no difference between normalized and standardised ("pearson"), but it might matter with more complicated correlation structures.
    ## type = "pearson" for 'standardised' residuals (see ?residuals.lme)
    ##  - this is what Zuur et al. 2009 are actually referring to, and is what Zuur et al. 2007 use in their R-code (type="p")
    ## see also rstudent (studentized) and rstandard (standardized) residuals.
    ## Zuur et al. use stdres() & studres() from the MASS library - what's the difference?
  }

  ## Residuals ##
  if (resType=="r" | resType=="response")   RE.lab <- "Raw Residuals"
  if (resType=="p" | resType=="pearson")    RE.lab <- "Standardized Residuals"
  if (resType=="n" | resType=="normalized") RE.lab <- "Normalized Residuals"
  RE <- resid(Y.model, type=resType)   # c("response", "pearson", "normalized")

  ## Plot REsiduals: see Zuur et al. 2007, pg. 131-133
  ## REsiduals: Normal distribution?
  op <- par(mfrow=c(2,2))
  if (more) {
    hist(RE)
    hist(RE, freq=FALSE, xlab=RE.lab, main="")
    Xnorm <- seq(min(RE), max(RE), length=40)
    lines(Xnorm, dnorm(Xnorm), col="grey70", lwd=2) # reference Normal distribution
    qqnorm(RE)
  }
  qqPlot(RE, ylab=RE.lab)

  ## Homogeneity, Independence: Pattern in residuals indicates violation
  Fit  <- fitted(Y.model)
  Fit0 <- fitted(Y.model, level=0)
  Fit1 <- fitted(Y.model, level=1)
  plot(Fit, RE, xlab="Fitted values", ylab=RE.lab)
  ## Plot residuals vs. predictor (explanatory) variables
  for (Xcol in X.cols) plot(SECCa[[Xcol]], RE, ylab=RE.lab, xlab = Xcol)             
  ## spreadLevelPlot(Y.model)                  # library(car)
  mtext(label, 3, adj=0.5, line=-2, outer=TRUE)
  par(op)

  if (more) {
    print( qplot(Block, RE, data = SECCa, ylab=RE.lab, facets = . ~ Time, geom="boxplot" ) + jaw.ggplot() )
    print( qplot(Chamber, RE, data = SECCa, ylab=RE.lab, facets = Frag * Position ~ Time, geom="boxplot" ) + jaw.ggplot() )

    for (Xcol in X.cols) 
      print( qplot(get(Xcol), RE, data = SECCa, ylab=RE.lab, xlab = Xcol, facets = Block ~ Time) + jaw.ggplot() )
  }

  ## Global Validation of Linear Model Assumptions (gvlma package)
  if (FALSE & ( "lm" %in% class(Y.model) )) {
    require(gvlma)                         # diagnostic plots
    validation <- gvlma(Y.model)
    plot(validation)
    summary(validation)
  }
  return(invisible(RE))
}



##================================================
## ANALYSIS: GET RESULTS
##================================================

## functions for processing effects() output for graphing
effect2df <- function(eff, column = "effect", 
                         response.var = "Response", fac1.label="fac1", fac2.label = "fac2")
{   # Deprecated: use effect.to.df?
  fac1.label <- names(eff$variables)[1]
  fac2.label <- names(eff$variables)[2]
  effw <- as.data.frame( summary(eff)[[column]] )
  fac2 <- colnames(effw)
  EffCols   <- paste(response.var, colnames(effw), sep=".")
  colnames(effw) <- EffCols
  effw[[fac1.label]] <- rownames(effw)
  effw$Var  <- response.var
  effl <- reshape(effw, varying = list(EffCols), direction = "long")
  effl[[fac2.label]] <- factor(effl$time, labels = fac2)
  effl$Response <- effl[, EffCols[1]]
  effl <- effl[, c("Var", fac1.label, fac2.label, "Response")]
  effl
}

effect.to.df <- function(eff, fun.trans = NULL)
{   # extract plotting data from an `eff` object (& back-transform if necessary)
  vars <- names(eff$variables) 
  eff.df <- cbind(eff$x, effect = eff$fit, lower = eff$lower, upper = eff$upper)
  ## You will have to back-transform any variables in eff.df$x yourself (they may be different)
  cols.trans <- c("effect", "lower", "upper")
  if (!is.null(fun.trans)) eff.df[, cols.trans] <- do.call(fun.trans, list(eff.df[, cols.trans])) 
  eff.df
}

intermean <- function (vec) 
{   # for getting the midpoint between breaks (used for converting continuous variables into discrete groups)
  vec1 <- rep(NA, length(vec) -1)
  for (i in 1:length(vec1) ) {
    vec1[i] <- mean(vec[c(i, i+1)])
  }
  vec1
}

eff.layer <- function(eff.df = NULL, conf.int = TRUE, ...)
{   # assumes certain columns in eff.df (using effect.to.df())
  require(ggplot2)
  ## the ... doesn't really work as planned, because these function calls aren't evaluated until plotting time :(
  result <- list( geom_line(data=eff.df, aes(y=effect, ... ) ) )
  if (conf.int == TRUE) {
    result <- c(result, 
                geom_line(data=eff.df, aes(y=lower, lty=2) ), 
                geom_line(data=eff.df, aes(y=upper, lty=2) ) 
                )
  }
  result
}


##================================================
## ANALYSIS: OUTPUT
##================================================
SaveDir.obj   <- function () "./save/"
SaveDir.plots <- function () "./graphs/"
SaveDir.text  <- function () "./output/"
Save.div  <- function () {
              "================================================================\n" 
}
Save.header  <- function (head.txt="") {
  paste(head.txt,
        paste("", version$version.string, date(), 
              Save.div(), "", sep = "\n"),
        sep="\n"
  )
}
Save.end <-  function () {
  paste(      "\n", 
              "<============================= END ============================>",
        sep = "\n"
        )
}



##================================================
## PLOTS & GRAPHS
##================================================

plotMap <- function (factor = c("Chamber", "Frag", "Position"), 
                     labels = c() ) 
{
  ## allow labels to be a subset of factor levels, and auto-drop labels here?
  ## strings should *NOT* be stored as factors, or they will not be recognized as strings, and cause problems when trying to use them in graphing functions.  Alternatively, I could probably store these as matrices or tables, rather than data frames, but data frames make more sense to me for this purpose.
  
  factor <- match.arg(factor)

  if (factor == "Chamber") {
    if (is.null(labels)) labels <- c("A", "B", "C")
    PlotMap <- 
      data.frame(label = labels, 
                 col = c("#000000","#000099","#990000"), 
                 bg  = c("#FFFFFF","#FFFFFF","#FFFFFF"), 
                 pch = c(21, 23, 18), lty = c(3, 2, 1), lwd = c(2, 1.5, 1),
                 stringsAsFactors = FALSE
                 )
    ## A) Ambient = black, open circles with dotted line ; 
    ## B) Partial = blue, open diamonds with dashed line ; 
    ## C) Full    = red, solid diamond with solid line.
  }

  if (factor == "Frag") {
    if (is.null(labels)) labels <- c("1", "2", "3", "4")
    PlotMap <- 
      data.frame(label = labels, 
                 col = c("#000000", "#666666", "#000099", "#990000"), 
                 bg  = c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF"), 
                 pch = c(19, 15, 22, 21),
                 lty = c(1, 2, 3, 3),
                 stringsAsFactors = FALSE
                 )
    ## 1) Contiguous         = black, filled circles with solid line ; 
    ## 2) Full Corridors     =  grey, filled squares with dashed line ; 
    ## 3) Pseudo-Corridors   =  blue, open squares with dotted line.
    ## 4) Isolated           =   red, open circles with dotted line.
  }

  if (factor == "Position") {
    if (is.null(labels)) labels <- c("I", "*", "O")
    PlotMap <- 
      data.frame(label = labels, 
                 col = c("#000000", "#000099", "#990000"), 
                 bg  = c("#FFFFFF", "#FFFFFF", "#FFFFFF"), 
                 pch = c(19, 17, 21), lty = c(2, 3, 1), lwd = c(2, 1.5, 1),
                 stringsAsFactors = FALSE
                 )
    ## 1) Inner         = black filled circles with dotted line ; 
    ## 2) intermediate  = blue           stars with dashed line ; 
    ## 3) Outer         = red   open   circles with solid  line.
  }  

  return(PlotMap)
}

SECC.axislab <- function(SECC.df=SECC, col.name="ARA.m", parens=TRUE, multiline=FALSE,
                         trans.func=NULL, unit.mult = 1, ...) 
{
  ax.label <- attr(SECC.df, "labels")[[col.name]] 
  ax.units <- attr(SECC.df, "units" )[[col.name]] 
  if (FALSE)
  { # wrap it as an expression() instead - I should have done this from the very beginning.  Let's hope this works anyway ...
    if (!is.null(trans.func)) {          # Not working :(
      ## dynamically wrap the units in a plotmath 'function' (passed as text)?
      ##     ax.units <- paste(trans.func, "(", deparse(ax.units), ")", sep="")
      if (length(grep("\\(", trans.func)) <1 ) 
        trans.func <- paste(trans.func, "(%s)", sep="")
  ax.units <- sprintf(trans.func, deparse(ax.units))
  ax.units <- parse(text = ax.units)
    }
    if (unit.mult!=1) ax.units <- bquote("" %*% .(format(unit.mult, ...)) ~ .(ax.units))
    sep <- if (multiline==TRUE) "\n" else " "
    axis.lab <- if (parens==TRUE) {
      bquote( .(ax.label) ~ .(paste(sep, "(", sep="")) * .(ax.units) * ")" )
    } else {
      bquote( .(ax.label) ~ .(paste(sep, " ", sep="")) * .(ax.units) )
    }
  } else {
    ax.units <- deparse(ax.units)    # -> string
    if (!is.null(trans.func)) {
      if (length(grep("\\(", trans.func)) <1 ) 
        trans.func <- paste(trans.func, "(%s)", sep="")
      ax.units <- sprintf(trans.func, ax.units)
    }
    unit.pre <- ""                     # default
    if (unit.mult!=1) unit.pre <- paste("%*%", deparse( paste( format(unit.mult, ...), "") ) )
    sepc <- if (multiline==TRUE) "\n" else " "
    axis.text <- if (parens==TRUE) {
      paste( deparse(paste(ax.label, sepc, "(", sep="")), 
            sprintf(" %s * %s * ", unit.pre, ax.units), 
            deparse(")"), sep = "" )
    } else {
      paste(deparse(paste(ax.label, sepc, sep="")), 
            sprintf(" %s * %s", unit.pre, ax.units), 
            sep = "" )
    }
    axis.lab <- parse(text = axis.text)
  }
  axis.lab
}

ggPts.SECC <- function (ptMap = plotMap("Chamber"), name = "Chamber Treatment") {
  ## Add Chamber shapes & colours to a ggplot scatterplot
  require(ggplot2)
  plotTemplate <- list(scale_colour_manual(name = name,
                                           values = ptMap$col, 
                                           breaks = levels(ptMap$label)
                                           ),
                       scale_shape_manual(name = name,
                                          values = ptMap$pch, 
                                          breaks = levels(ptMap$label)
                                          ),
                       scale_size_manual(name = name,
                                         values = ptMap$lwd*0.5, 
                                         breaks = levels(ptMap$label)
                                         )
                       )
  ## should just be breaks = Chamber.map$label, but that produces right-aligned text :(
}

SECCicons <- function() 
{ ## Load eps graphics for plot labels
  require(grImport)
  setwd("./save/")     # The loading functions create temp files and I want them to go here.
  FragIcons <- PostScriptTrace("Frag-Black-icons.eps")
  FragIcons <- readPicture("Frag-Black-icons.eps.xml")
  Hex <- PostScriptTrace("hexagon.eps")
  Hex <- readPicture("hexagon.eps.xml")
  setwd("..")
  FragIcon1 <- FragIcons[49:50]        # 1. Contiguous 
  FragIcon2 <- FragIcons[29:48]        # 2. Corridors
  FragIcon3 <- FragIcons[ 9:28]        # 3. Pseudo-corridors
  FragIcon4 <- FragIcons[ 1:8 ]        # 4. Isolated
  IconList <- list(FragIcon1 = FragIcon1,
                   FragIcon2 = FragIcon2,
                   FragIcon3 = FragIcon3,
                   FragIcon4 = FragIcon4,
                   Hex = Hex
                   )
}


SECCplotDataANOVA <- function(SECCdata, 
                              by.agg = list(Chamber = SECCdata$Chamber), 
                              FUN.agg = mean, 
                              error.msd, backtransform = TRUE, ...)
{
  if (class(SECCdata) == "data.frame") SECCdata <- SECCdata$Y.trans
  plot.means <- aggregate(SECCdata, by = by.agg, FUN = FUN.agg, ... )
  if ("Time" %in% colnames(plot.means))
    levels(plot.means$Time) <- paste(c("August", "June", "August"), levels(plot.means$Time), sep="\n")
  plot.means <- within(plot.means, 
                       {
                         error <- as.numeric(error.msd/2)
                         upper <- x + error
                         lower <- x - error
                         if (backtransform)
                         {
                           if (Y.use == "Y.log")
                           {           # Y.log = log(x + 1, 10); subtract 1 from result after?
                             x <- (10^x) -1
                             lower <- (10^lower) -1
                             upper <- (10^upper) -1
                           }
                           if (Y.use == "Y.sqrt")
                           {
                             x <- x^2
                             lower <- lower^2
                             upper <- upper^2
                           }
                           if (Y.use == "Y.4rt")
                           {
                             x <- x^4
                             lower <- lower^4
                             upper <- upper^4
                           }
                           if (Y.use %in% c("Y.sqrt", "Y.4rt"))
                           {  # correct for negative values in limits :(
                             ## lower[which(x - error < 0)] <- 0
                             lower[which(x - error < 0)] <- lower[which(x - error < 0)] * -1
                           }
                         }
                       })
  plot.means
}


## (back)transformation functions
log0 <- function(x, base = exp(1))
{   # log-transform, but 0 stays as 0
  y <- log(x, base)
  y[x <= 0] <- 0
  ##   y[x <  1] <- 0
  y
}
alog0 <- function(y)
{   # back-transform log(X), where 0s stay as 0s
  x <- 10^y
  ##   x[y == 0] <- NA                      # MISSING (these cause glitches in transformed axes in ggplot)
  ##   x[y == 0] <- 0
  x
}
alog1 <- function(y, base = 10)
{   # back-transform log(X+1): this is problematic for log-axes if it results in -ve values :(
  x <- (base^y) -1
  x
}


## Plotting Regression Tree results
RegTreePlot.SECC <- function(Reg.tree, minsplit) 
{
  require(ggplot2)
  require(ggdendro)
  Y.cp <- min(Reg.tree$cptable[, "CP"])
  Suppl.text  <- paste("cp = ", Y.cp, ", split min. n = ", minsplit, sep="")

  Yvals       <- as.numeric(Reg.tree$frame[, "yval"])
  Ylabels     <- paste(" ", formatC(Yvals, digits=1, flag="-", format="f"), sep="")
  Nvals       <- as.numeric(Reg.tree$frame[, "n"])
  Nlabels     <- paste("", formatC(Nvals), sep="") # "n="
  Droot       <- as.numeric(Reg.tree$frame[1, "dev"]) 
  Dvals       <- as.numeric(Reg.tree$frame[, "dev"]) / Droot # proportion of total deviance
  Dlabels     <- paste("", formatC(Dvals, digits=2, big.mark=",", flag="-", format="f"), sep="")
  Leaf.labels <- paste(Ylabels, " (", Dlabels, ")", sep="") #
  Info.labels <- paste(" n=", Nlabels, sep="") # Y.tree.labels
  Var.labels  <- labels(Reg.tree, pretty=F) # includes full factor levels, decision criteria, etc.
  ## clean up labels
  Var.labels <- gsub(" months,", ",", Var.labels) #  duplicate month labels
  Var.labels <- gsub("H2O(..[0-9.]+)", "Moisture\\1%", Var.labels)
  Var.labels <- gsub("X.trans", "Cells", Var.labels)
  Var.labels <- gsub("Full Corridors", "Corridors", Var.labels)
  Var.labels <- gsub("Pseudo-Corridors", "Pseudo-Cs", Var.labels)
  if (FALSE) { # ggplot2 doesn't support math expressions (yet)?
    Var.labels <- gsub("H2O(..[0-9.]+)", "H[2]*O\\1%", Var.labels)
    Var.labels <- gsub("\\de\\+(\\d+)", " %*% 10^{\\1}", Var.labels)
    Var.labels <- expression(Var.labels)      # expressions
  } else {
  }
  Var.labels  <- paste(" ", Var.labels, sep="")
  Node.labels <- paste(Var.labels, "\n", Dlabels, sep="") # space in front of each line for padding (hack) ; blank line for leaf labels in between
  Full.labels <- paste(" ", Var.labels, "\n ", Ylabels, " (n=", Nlabels, ") ", Dlabels, sep="")
  Y.tree.labels <- data.frame(label=Full.labels, 
                                node=Node.labels, leaf=Leaf.labels, info=Info.labels,
                                var=Var.labels, dev=Dlabels, n=Nlabels, yval=Ylabels)
  ## x=Y.plot$x, y=Y.plot$y, 

  class(Reg.tree) <- c("rpart", "tree") # largely the same, but no methods for "rpart"
  dend_data <- dendro_data(Reg.tree)
  ## dend.labels <- Y.tree.labels$label[as.numeric(row.names(dend_data$labels)) +1]
  numLabels <- nrow(Y.tree.labels) -1  # drop "root"
  branch_labels <- segment(dend_data)[1:numLabels, c("x", "y", "xend", "yend")]     # coordinates for segments / branches
  ## names(branch_labels) <- c("x", "y")    # rename columns ;)
  branch_labels <- cbind(branch_labels, Y.tree.labels[-1, 3:ncol(Y.tree.labels)])
  dend_data$branch_label <- branch_labels # package it all up
  dend.leaf   <- Y.tree.labels$leaf[ as.numeric(row.names(dend_data$leaf_label)) ]
  leaf_labels <- dend_data$leaf_label
  leaf_labels$label <- dend.leaf
  dend_data$leaf_label <- leaf_labels
  ## move leaves from edge to leave space for labels?
  ## ylim doesn't work too well, but I did discover that the plotting area will expand
  ## to include all coordinates (on all layers), so I can just add some (blank)
  ## text at the limits to ensure the axes go at least that far.
  Y.lim <- c( min(segment(dend_data)$y), max(segment(dend_data)$y) )
  Y.lim[1] <- Y.lim[1] *0
  Suppl.label <- data.frame(label=Suppl.text, x=(max(segment(dend_data)$x) +0.5), 
                            y=max(segment(dend_data)$y), stringsAsFactors=FALSE)
  Suppl.label <- rbind(Suppl.label, 
                       data.frame(label=" ", x=1, y=Y.lim[1], stringsAsFactors=FALSE)
  ) # add some space 'below' bottom leaves for labels.

  RegTree.plot <- ggplot(segment(dend_data)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=1, color="grey50") +
  scale_size_identity() + scale_color_identity()
  RegTree.plot <- RegTree.plot + geom_text(data=branch_labels, 
                                             aes(label=var, x=xend, y=yend),
                                             hjust=0, vjust=-0.7, size=3.1, colour="grey20") + 
                                   scale_size_identity() + scale_colour_identity()
  RegTree.plot <- RegTree.plot + geom_text(data=branch_labels, 
                                             aes(label=info, x=xend, y=yend),
                                             hjust=0, vjust=1.75, size=3.1, colour="grey20") + 
                                   scale_size_identity() + scale_colour_identity()
  RegTree.plot <- RegTree.plot + geom_text(data=branch_labels, 
                                             aes(label=leaf, x=x, y=y),
                                             hjust=0, vjust=0.5, size=3) + 
                                   scale_size_identity()
  ## RegTree.plot <- RegTree.plot + geom_text(data=leaf_labels, aes(label=label, x=x, y=y), 
  ##                                            hjust=-0.1, size=4) + scale_size_identity()
  RegTree.plot <- RegTree.plot + geom_text(data=Suppl.label, aes(label=label, x=x, y=y), 
                                             hjust=0, size=4) + scale_size_identity()
  RegTree.plot <- RegTree.plot + coord_flip() + scale_y_reverse() + 
  theme_dendro() + opts(panel.border = theme_blank())
  ## expand=c(0.2, 0) within the sacle_y_reverse() adds space to *both* sides, to make room for long labels.  Not ideal, but it works
  ## expand=c(0.2, 0)

  ## these ggplot objects are not entirely stand-alone if they depend on the dendro objects!!
  ## dendro_data()
  RegTree.plot                          # return
}




##################################################
## MULTIVARIATE ORDINATION ANALYSES
##################################################

## SIMPER
# Similarity Percentages: contribution of individual variables to pariwise dissimilarity between all pairs of sites between to groups defined a priori.
## see: Clarke, K.R.  1993.  Non-parametric multivarite analysis of changes in community structure.  Australian journal of ecology, 18, 117-143.

simper <- function( data, groups, method="bray") {
	require(vegan)	# uses the vegan package for ecologically-relevant distance metrics.
	# data is the data matrix: n samples x p variables / species
	# groups is a vector with length = # of samples, specifying a priori groups for each observation / sample
	# method is the distance method used.
	if (length(groups)!=dim(data)[1] || is.null(dim(groups)) != TRUE)
		stop("invalid dimensions of groups: must be a vector with length equal to the number of rows in data")
	if (method != "bray" )
      stop("simper() currently only works with bray-curtis dissimilarity (method= \"bray\")")

	# initialize some variables
	groups.levels <- unique(groups)
	groups.pairs  <- t( combn(groups.levels, 2)	) # all combinations of n elements, taken 2 at a time (all pairs); transposed to 2 columns, 1 with each grouping level in the pair.
    vars <- colnames(data)
    Result <- list(within = list(), between = list())

    ## WITHIN GROUP SIMILARITIES
    for (group in groups.levels)
    {
      gdata <- data[which(groups == group), ]
      group.dist <- mean( vegdist(gdata, method = method), na.rm = TRUE)
      g.pairs <- t( combn( 1:nrow(gdata), 2 ) )
      g.pairs <- as.data.frame(g.pairs)
      colnames(g.pairs) <- c("j", "k")
      bray.j <- gdata[g.pairs$j, ]
      bray.k <- gdata[g.pairs$k, ]
      bray.jk <- g.pairs
      bray.jk$Yjk <- apply( bray.j + bray.k, 1, sum )
      bray.jk  <- cbind(bray.jk[, 1:3], abs(bray.j - bray.k) / bray.jk$Yjk )
      bray.sim <- pmin(bray.j, bray.k)
      colnames(bray.sim) <- colnames(bray.jk[, -(1:3)]) 
      bray.sim <- 2 * bray.sim / bray.jk$Yjk 
      bray.sim <- cbind(bray.jk[, 1:3], bray.sim)
      group.sim <- mean( 1 - apply(bray.jk[, -(1:3)], 1, sum) )
      group.sim <- mean(    apply(bray.sim[, -(1:3)], 1, sum) )    # should be the same

      Result.group <- data.frame(Variable = "", Mean.Value = 0, 
                                 Avg.Sim = 0, SD.Sim = 0, Avg_SD = 0,
                                 Perc.Sim = 0.0, Cum.Perc = 0.0
      )
      Result.group <- Result.group[0, ] # structure only
      for (vari in vars)
      {
        gvdata <- gdata[vari]       # stays as a data.frame
        ##         var.sim <- mean( vegdist(gvdata, method = method), na.rm = TRUE)
        var.sim <- mean( bray.sim[[vari]] )
        SD.sim  <-   sd( bray.sim[[vari]] )
        Result.group <- rbind(Result.group,
                              data.frame(Variable = vari, 
                                         Mean.Value = mean(gvdata[1], na.rm = TRUE), 
                                         Avg.Sim = var.sim, SD.sim = SD.sim,
                                         Avg_SD = var.sim / SD.sim,
                                         Perc.Sim = var.sim / group.sim, 
                                         Cum.Perc = NA
                                         )
        )
      }
      attr(Result.group, "Average Bray-Curtis Similarity") <- group.sim
      colnames(Result.group)[2] <- paste("Mean", group, sep = "_")
      ## sort result
      Result.group <- Result.group[ order(Result.group$Avg.Sim, decreasing = TRUE), ]
      Result.group$Cum.Perc  <- cumsum(Result.group$Perc.Sim)
      Result$within[[group]] <- Result.group
    }
	
    ## BETWEEN GROUP DISTANCES
    for (g in 1:nrow(groups.pairs))
    {
      pair <- groups.pairs[g, ]
      pair.lab <- paste(pair, collapse = ":")
      gdata <- data[which(groups %in% pair), ]
      group.dist <- mean( vegdist(gdata, method = method), na.rm = TRUE)
      ## the hardest part is probably getting the list of sample pairs; the rest is the same as above
      g.pairs <- t( combn( 1:nrow(gdata), 2 ) )
      g.pairs <- as.data.frame(g.pairs)
      colnames(g.pairs) <- c("j", "k")
      g.pairs$group.j <- groups[g.pairs$j]
      g.pairs$group.k <- groups[g.pairs$k]
      g.pairs <- g.pairs[g.pairs$group.j != g.pairs$group.k, 1:2]
      ## Compute Distance between all pairs of samples, between all group pairs, using one variable/species at a time.
      bray.j <- gdata[g.pairs$j, ]
      bray.k <- gdata[g.pairs$k, ]
      bray.jk <- g.pairs
      bray.jk$Yjk <- apply( bray.j + bray.k, 1, sum )
      bray.jk  <- cbind(bray.jk[, 1:3], abs(bray.j - bray.k) / bray.jk$Yjk )
      pair.dist <- mean( apply(bray.jk[, -(1:3)], 1, sum) )

      Result.pair <- data.frame(Variable = "", Mean.Value.1 = 0, Mean.Value.2 = 0,  
                                 Avg.Dist = 0, SD.Dist = 0, Avg_SD = 0,
                                 Perc.Dist = 0.0, Cum.Perc = 0.0
      )
      Result.pair <- Result.pair[0, ] # structure only
      for (vari in vars)
      {
        gvdata <- gdata[vari]       # stays as a data.frame
        ##         var.sim <- mean( vegdist(gvdata, method = method), na.rm = TRUE)
        var.dist <- mean( bray.jk[[vari]] )
        SD.dist  <-   sd( bray.jk[[vari]] )
        Result.pair <- rbind(Result.pair,
                              data.frame(Variable = vari, 
                                         Mean.Value.1 = mean(gdata[which(groups == pair[1]), vari], 
                                                             na.rm = TRUE), 
                                         Mean.Value.2 = mean(gdata[which(groups == pair[2]), vari],
                                                             na.rm = TRUE), 
                                         Avg.Dist = var.dist, SD.dist = SD.dist,
                                         Avg_SD = var.dist / SD.dist,
                                         Perc.Dist = var.dist / pair.dist, 
                                         Cum.Perc = NA
                                         )
        )
      }
      attr(Result.pair, "Average Bray-Curtis Distance") <- pair.dist
      colnames(Result.pair)[2:3] <- paste("Mean", pair, sep = "_")
      ## sort result
      Result.pair <- Result.pair[ order(Result.pair$Avg.Dist, decreasing = TRUE), ]
      Result.pair$Cum.Perc  <- cumsum(Result.pair$Perc.Dist)
      Result$between[[pair.lab]] <- Result.pair
    }
	
	## Similarity percent = mean distance along variable i / mean distance using all variables = % of mean distance contributed by variable i.
	
	# Result:
	# for each group pair:
	# a table (data frame) with each species / variable as a row with the following columns:
	# - mean value of variable in group 1
	# - mean value of variable in group 2
	# - mean distance using only this variable
	# - SD(of distances between all pairs of samples between groups)
	# - mean distance of variable / SD(of distances between all pairs of samples)
	# - mean distance of variable / mean overall distance (%)
	# - cumulative percent
	#
	# Same columns for each single group: to account for variable contributing to similarity within a particular group.
    class(Result) <- c(class(Result), "jaw.simper")
	
    return(Result)

    if (FALSE)
    {                                  # Testing code
      library(vegan)
      data(mite)
      data(mite.env)
      simper(mite, groups = mite.env$Topo)
    }
}

print.jaw.simper <- function(simp, ...)
{
  print.default(simp, ...)
}
