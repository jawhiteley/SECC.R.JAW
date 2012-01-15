##################################################
# functions for graphing with custom error bars
# Jonathan Whiteley     R v2.12     2011-10-05
##################################################
## access these functions in another file by using: 
## 	source("path/to/this/file.R")

##================================================
## MULTIPLE COMPARISONS
##================================================
MSD <- function( model, alpha=0.05, method=c("unplanned", "planned", "MSR", "LSD"), mode.df=c("pairwise", "MSE", "manual"), n=NULL ){
  ## Planned Multiple Comparisons using Least Significant Differences (LSD)
  ## -> comparison intervals for graphical display.
  ## Expects a model of class 'aov' and/or 'lm'? (aovlist)
  ## Design must be balanced (easier for now).

  method  <- match.arg(method)      # default "planned"; synonymous with "MSR"
  mode.df <- match.arg(mode.df)     # default "pairwise"
  conf.level <- 1 - (alpha/2)	# for a 2-tailed test

  ##   vars <- attr( attr(model$terms, "factors"), "dimnames")[[1]]  # dig in for the variable names ...
  ##   data <- data[, vars]  # keep only relevant columns
  ## re-factor columns to drop unused levels.
  ##   for( col in names(data) ) {
  ##     if (class(data[[col]]) == "factor") {
  ##       data[[col]] <- factor(data[[col]])
  ##     }
  ##   }

  if(mode.df == "manual"){
    n <- as.numeric(n)	# relevant group ?
  }
  if(mode.df == "pairwise"){
    ##  n <- replications( model$terms , data )	# sample sizes, according to model structure (formula).  
    ## I couldn't find an easy way to derive this directly from an aov object passed in, so this is the only reason that data is required as an argument (& formula?).  replications() returns a nasty list if the data are unbalanced :(
    tables <- model.tables(model, "means")
    means <- tables$tables
    n <- tables$n	# sample sizes, according to model structure (formula).  
    ### model.tables(model)$n  # does what I want!  Only with *full* model
    ## if the item name is a problem, use as.numeric() to convert to a pure number.  If no group specified, a (named) vector is produced with values for all treatment combinations.
  }

  model.ls <- model  # in case of aovlist

  msd.width <- c()  # container variable
  for (lvl in 1:length(n)) {  # loop through each treatment combination
    
    if ("aovlist" %in% class(model.ls)) {
      msd.lvl <- names(n)[lvl]
      model.lvl <- grep(msd.lvl, names(model.ls))  # find matching name
      ## use last level if no match
      if( length(model.lvl) == 0 ) model.lvl <- length(names(model.ls))
      model.lvl <- model.lvl[1]  # keep only first match
      model <- model.ls[[model.lvl]]
    }

    ## Calculate MSE & se for each level of nesting
    MSE <- sum(resid(model)^2)/model$df.residual	# MSE from model ( SS / df )
    msd.se <- sqrt(2*MSE/n)	# se of a difference (Crawley 2007, pg.465; Sokal & Rohlf 2003/1995, pg. 243)
    ## unbalanced: sqrt( (var[1]/n[1]) + (var[2]/n[2]) )
    
    ## A Minimum Significant Difference (MSD) = (critical value) * SE
    ## - Sokal & Rohlf, pg.243-244
    ## Student's t = a difference / se of the difference
    ## - Crawley, pg. 464
    if( method == "unplanned" || method == "MSR") {
      ngroups <-  length( means[[lvl +1]] )  # skip "Grand Mean"
      msd.df <- ngroups * (n -1)
      width <- qtukey(conf.level, ngroups, msd.df) * sqrt(MSE/n) 
      ## for unplanned comparisons, it should really be Minimum Significant Range
      ## MSR = Q[alpha,k, v] * sqrt(MSE/n)    (see Sokal & Rohlf, pg.244).
    } else {
      if(mode.df=="MSE") {
        msd.df <- model$df	# for MSE from fitted model?  This is probably completely incorrect.
      } else {
        msd.df <- (2*n)-2   # for differences between any 2 levels. 
      }

      width <- qt(conf.level, msd.df)*msd.se  # LSD based on error rate (alpha).
    }

    msd.width[lvl] <- width[lvl]
  }
  names(msd.width) <- names(n)
  return(msd.width)
}



##================================================
## PLOTS & GRAPHS
##================================================
##________________________________________________
## text_to_plot: mtext() wrapper.  Outputs text to a plotting device.
## useful for getting some text output into a pdf
## file with graphs.

text_to_plot <- function( txt = "", auto.size = TRUE,
                         cex = 0.6, family = "mono" , ... ) {
  ## can fit about 60 (56) letters @ cex = 1
  ## auto-word-wrap?  Insert "\n" every X characters, depending on cex?
  plot.new()
  mtext(txt,
        adj = 0, side = 3, line = 0,
        cex = cex, family = family, col = "black",
        ...
        )
}



##________________________________________________
## plotMeans: modified version from 'R commander' package (Rcmdr).  I just hate loading the entire package & GUI, just for this one convenient function.  
# I added an option for custom error bar widths: set 'error.bars' to "custom", and pass the custom width via the 'level' argument (sort of hack).  TO DO: allow a vector of upper & lower error bars for each point.
# I also added the ability to accept more graphing parameters: bg, cex, lwd.  Could access to these & other graphical parameters be provided more easily with the '...' (or '…') argument?
plotMeans <- function (response, factor1, factor2, 
  error.bars = c("se", "sd", "conf.int", "none", "custom"), 
  level = 0.95, xlab = deparse(substitute(factor1)), 
  ylab = paste("mean of", deparse(substitute(response))), 
  legend.lab = deparse(substitute(factor2)),
  main = "Plot of Means", pch = NULL, lty = NULL,
  col = palette(), bg=NULL, cex=2, lwd=NULL, ylim=NULL, ... ) 
{
    if (!is.numeric(response)) 
        stop(gettextRcmdr("Argument response must be numeric."))
    xlab
    ylab
    legend.lab
    error.bars <- match.arg(error.bars)
    if (missing(factor2)) {
        if (!is.factor(factor1)) 
            stop(gettextRcmdr("Argument factor1 must be a factor."))
        valid <- complete.cases(factor1, response)
        factor1 <- factor1[valid]
        response <- response[valid]
        means <- tapply(response, factor1, mean)
        sds <- tapply(response, factor1, sd)
        ns <- tapply(response, factor1, length)
        if (error.bars == "se") 
            sds <- sds/sqrt(ns)
        if (error.bars == "conf.int") 
            sds <- qt((1 - level)/2, df = ns - 1, lower.tail = FALSE) * 
                sds/sqrt(ns)
        if (error.bars == "custom") 
            sds <- level ## custom error bar widths.  Added by JAW.  
            # custom widths passed as an appropriate matrix, or a single value replicated: matrix( level, nrow = dim(sds)[1], ncol = dim(sds)[2] )
        sds[is.na(sds)] <- 0
        yrange <- if (error.bars != "none") 
            c(min(means - sds, na.rm = TRUE), max(means + sds, 
                na.rm = TRUE))
        else range(means, na.rm = TRUE)
        yrange <- if ( is.null(ylim)==FALSE ) ylim else yrange # added by JAW
        levs <- levels(factor1)
        n.levs <- length(levs)
        if (is.null(pch)) pch <- 1
        if (is.null(lty)) lty <- 1
        plot(c(1, n.levs), yrange, type = "n", xlab = xlab, ylab = ylab, 
            axes = FALSE, main = main, ... )	# `...` argument added by JAW
        points(1:n.levs, means, type = "b", pch = pch, cex = cex, bg = bg, lty=lty, lwd=lwd)	# originally: pch=16, cex=2.  Modified by JAW.
        box()
        axis(2)
        axis(1, at = 1:n.levs, labels = levs)
        if (error.bars != "none") 
            arrows(1:n.levs, means - sds, 1:n.levs, means + sds, 
                angle = 90, lty = 2, code = 3, length = 0.125)
    }
    else {
        if (!(is.factor(factor1) | is.factor(factor2))) 
            stop(gettextRcmdr("Arguments factor1 and factor2 must be factors."))
        valid <- complete.cases(factor1, factor2, response)
        factor1 <- factor1[valid]
        factor2 <- factor2[valid]
        response <- response[valid]
        means <- tapply(response, list(factor1, factor2), mean)
        sds <- tapply(response, list(factor1, factor2), sd)
        ns <- tapply(response, list(factor1, factor2), length)
        if (error.bars == "se") 
            sds <- sds/sqrt(ns)
        if (error.bars == "conf.int") 
            sds <- qt((1 - level)/2, df = ns - 1, lower.tail = FALSE) * 
                sds/sqrt(ns)
        if (error.bars == "custom") 
            sds <- level ## custom error bar widths.  Added by JAW.
        sds[is.na(sds)] <- 0
        yrange <- if (error.bars != "none") 
            c(min(means - sds, na.rm = TRUE), max(means + sds, 
                na.rm = TRUE))
        else range(means, na.rm = TRUE)
        yrange <- if ( is.null(ylim)==FALSE ) ylim else yrange # added by JAW
        levs.1 <- levels(factor1)
        levs.2 <- levels(factor2)
        n.levs.1 <- length(levs.1)
        n.levs.2 <- length(levs.2)
        if (is.null(pch)) pch <- 1:n.levs.2
        if (is.null(lty)) lty <- 1:n.levs.2
        if (length(pch) == 1) 
            pch <- rep(pch, n.levs.2)
        if (length(col) == 1) 
            col <- rep(col, n.levs.2)
        if (length(lty) == 1) 
            lty <- rep(lty, n.levs.2)
        ## More graphical parameters added by JAW.
        if (length(bg) == 1) 
            bg <- rep(bg, n.levs.2)
        if (length(lwd) == 1) 
            lwd <- rep(lwd, n.levs.2)
        if (length(cex) == 1) 
            cex <- rep(cex, n.levs.2)
        ## ##
        if (n.levs.2 > length(col)) 
            stop(sprintf(gettextRcmdr("Number of groups for factor2, %d, exceeds number of distinct colours, %d."), 
                n.levs.2, length(col)))
        plot(c(1, n.levs.1 * 1.4), yrange, type = "n", xlab = xlab, 
            ylab = ylab, axes = FALSE, main = main, ... )	# `...` argument added by JAW
        box()
        axis(2)
        axis(1, at = 1:n.levs.1, labels = levs.1)
        for (i in 1:n.levs.2) {
            points(1:n.levs.1, means[, i], type = "b", pch = pch[i], 
                cex = cex[i], col = col[i], lty = lty[i], lwd = lwd[i], bg = bg[i] )
            if (error.bars != "none") 
                arrows(1:n.levs.1, means[, i] - sds[, i], 1:n.levs.1, 
                  means[, i] + sds[, i], angle = 90, code = 3, 
                  col = col[i], lty = lty[i], length = 0.125)
        }
        x.posn <- n.levs.1 * 1.1
        y.posn <- sum(c(0.1, 0.9) * par("usr")[c(3, 4)])
        text(x.posn, y.posn, legend.lab, adj = c(0, -0.5))
        legend(x.posn, y.posn, levs.2, pch = pch, col = col, 
            lty = lty, pt.bg=bg, pt.lwd=lwd, lwd=lwd, pt.cex=cex*0.75 )
    }
    invisible(NULL)
}


##================================================
pairplot <- function(x, panel=points, upper.plots=TRUE, mirror.panels=TRUE,
                     correlations=TRUE, histograms=TRUE, add.smooth=TRUE, ...) 
{
  ## pairplot: wrapper for pairs().
  ## bivariate plots of variable pairs with a few extra useful bits of info
  ## * Mixed Effects Models with Extensions in Ecology with R Figure A.2 p.534 (Appendix)
  ## * Analyzing Ecological Data Figure 4.9 p.33 (CH4)
  ## Requires:
  ## + panel.cor      copied here from books' R code
  ## + panel.hist     copied here from books' R code
  ## + panel.smooth   the books defined a crippled version of the same built-in function
  ##   - see graphics::panel.smooth()
  ## + panel.ablines  There is a version in the lattice package (panel.lines).  
  ##   - I copied the one from the books' code (panel.lines2) and renamed it.
  ## ARGUMENTS ##
  ## x          object with data for pairplots (usually a data frame).  
  ##            Passed to pairs()
  ## panel      panel function, same default as pairs().
  ## upper.plots    TRUE:  correlation values in lower panels, plots in upper panels
  ##                FALSE: vice versa.
  ## mirror.panels  TRUE:  same panel function in upper & lower panels
  ##                FALSE: allow smoothers on one side and default panel on the other
  ##                Only relevant if correlations=FALSE
  ## correlations   TRUE:  Show Correlation values in font size proportional 
  ##                       to value in panels of one side.  
  ##                FALSE: do not show correlation values 
  ##                       (panel function will be used instead).
  ## histograms     TRUE:  Add histograms in diagonal panels.  FALSE: do not ...
  ## add.smooth     TRUE:  Scatterplots with smoother curves in panels.
  ##                FALSE: Use default (or specified) panel function instead.
  panel.plot <- if (add.smooth && mirror.panels) panel.smooth else panel
  panel1     <- if (add.smooth)   panel.smooth else panel.plot
  panel2     <- if (correlations) panel.cor    else panel.plot
  up.panel   <- if (upper.plots)  panel1       else panel2
  low.panel  <- if (upper.plots)  panel2       else panel1
  d.panel    <- if (histograms)   panel.hist   else NULL
  pairs(x, panel=panel,
        upper.panel =  up.panel,
         diag.panel =   d.panel,
        lower.panel = low.panel,
        ...
        )

  if (FALSE) {  ## test code
    pairplot(CO2)
    pairplot(CO2, add.smooth=FALSE)
    pairplot(CO2, hist=FALSE)
    pairplot(CO2, cor=FALSE)
    pairplot(CO2, upper=FALSE)
    pairplot(CO2, panel=panel.smooth, cor=FALSE)
    pairplot(CO2, cor=FALSE, mirror=FALSE)
  }
}

##================================================
## PANEL FUNCTIONS
##================================================
panel.cor <- function(x, y, digits=1, prefix="", cex.cor)
{
  ## "panel.cor" in "MyLibrary.R" from:
  ##   Analysing Ecological Data. (2007). Zuur, Ieno and Smith. Springer, 680p.
  ##   This function was produced by Alain Zuur (highstat@highstat.com)
  ##   www.highstat.com
  ## put correlations on the panels,
  ## with size proportional to the correlations.
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r <- abs(cor(x, y,use="pairwise.complete.obs"))

  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}

panel.hist <- function(x, ...)
{
  ## "panel.hist" in "MyLibrary.R" from:
  ##   Analysing Ecological Data. (2007). Zuur, Ieno and Smith. Springer, 680p.
  ##   This function was produced by Alain Zuur (highstat@highstat.com)
  ##   www.highstat.com
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
}

panel.ablines <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                           cex = 1, ...)
{
  ## "panel.lines2" in "MyLibrary.R" from:
  ##   Analysing Ecological Data. (2007). Zuur, Ieno and Smith. Springer, 680p.
  ##   This function was produced by Alain Zuur (highstat@highstat.com)
  ##   www.highstat.com
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)){
    tmp=lm(y[ok]~x[ok])
    abline(tmp)}

}

## Coplots with linear fits (from Zuur et al. 2007 Chapter 22 R code)
## Deprecated?  Do I still use Coplots rather than faceted ggplots anymore?
panel.lines2 <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
    cex = 1, ...)
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)){
        tmp=lm(y[ok]~x[ok])
        abline(tmp)}

}


##================================================
## GGPLOT2 custom options
##================================================

jaw.ggplot <- function(x) {
  require(ggplot2)
  if (FALSE) {
    if (class(x)!="ggplot") return(x)
    ## add common custom options to ggplot2 objects.
    ## I got tired of copying the same code over & over again.
    x <- x + theme_bw()                  # change to white background colour scheme.
    x <- x + opts(legend.key = theme_rect(colour = NA)) # remove rectangle from legend keys.
    return(x)
  } else {
    ## add common custom options to ggplot2 objects as a template.
    ## useage: MYggPlot + jaw.ggplot()
    ## I got tired of copying the same code over & over again.
    plotTemplate <- list(theme_bw(),   # change to white background colour scheme.
                         opts(legend.key = theme_rect(colour = NA), # remove rectangle from legend keys.
                              ## move axis labels away from axis text :/
                              axis.title.x = theme_text(size = 12, vjust=0, hjust=0.55),
                              axis.title.y = theme_text(size = 12, vjust=0.3, hjust=0.5, angle = 90)
                              )
    )
    return(plotTemplate)
  }
}
